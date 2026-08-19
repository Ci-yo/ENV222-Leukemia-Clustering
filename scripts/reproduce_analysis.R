args <- commandArgs(trailingOnly = FALSE)
file_arg <- gsub("\\\\", "/", sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
script_dir <- dirname(file_arg)
repo_root <- file.path(script_dir, "..")
check_only <- "--check-only" %in% commandArgs(trailingOnly = TRUE)

data_path <- file.path(repo_root, "data", "leukemia_train_cleaned.csv")
figures_dir <- file.path(repo_root, "figures")
results_dir <- file.path(repo_root, "results")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

expr_clean <- read.csv(data_path, check.names = FALSE)
stopifnot(nrow(expr_clean) == 38L, ncol(expr_clean) == 7131L)
stopifnot(identical(sort(unique(expr_clean$cancer)), c("ALL", "AML")))

expr_mat <- expr_clean[, !(names(expr_clean) %in% c("patient", "cancer"))]
expr_mat <- data.frame(lapply(expr_mat, as.numeric), check.names = FALSE)
rownames(expr_mat) <- paste0("P", expr_clean$patient, "_", expr_clean$cancer)

gene_variance <- apply(expr_mat, 2, var, na.rm = TRUE)
top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:100]
expr_scaled <- scale(expr_mat[, top_genes])
rownames(expr_scaled) <- rownames(expr_mat)

distance_matrix <- dist(expr_scaled, method = "euclidean")
hc_samples <- hclust(distance_matrix, method = "ward.D2")
cluster_assign <- cutree(hc_samples, k = 2)
composition <- as.data.frame.matrix(table(Cluster = cluster_assign, Cancer = expr_clean$cancer))
composition$cluster <- as.integer(rownames(composition))
composition <- composition[, c("cluster", "ALL", "AML")]
composition$majority_label <- ifelse(composition$ALL >= composition$AML, "ALL", "AML")

cluster_majority <- setNames(composition$majority_label, composition$cluster)
agreement <- sum(cluster_majority[as.character(cluster_assign)] == expr_clean$cancer)

write.csv(composition, file.path(results_dir, "cluster_composition.csv"), row.names = FALSE)
write.csv(
  data.frame(
    samples = nrow(expr_clean),
    expression_features = ncol(expr_mat),
    selected_genes = length(top_genes),
    linkage = "ward.D2",
    distance = "euclidean",
    post_hoc_label_agreement = agreement,
    agreement_denominator = nrow(expr_clean)
  ),
  file.path(results_dir, "analysis_summary.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(gene = top_genes, variance = unname(gene_variance[top_genes])),
  file.path(results_dir, "top_variable_genes.csv"),
  row.names = FALSE
)

if (!check_only) {
  png(file.path(figures_dir, "leukemia-dendrogram.png"), width = 1600, height = 960, res = 150)
  label_order <- hc_samples$order
  label_text <- rownames(expr_scaled)[label_order]
  label_type <- expr_clean$cancer[label_order]
  label_col <- ifelse(label_type == "ALL", "#4C78A8", "#E15759")
  old_par <- par(mar = c(6.8, 4.8, 4.5, 2))
  plot(
    hc_samples,
    labels = FALSE,
    hang = -1,
    main = "Hierarchical Clustering of Leukemia Samples",
    xlab = "",
    sub = "",
    ylab = "Distance",
    ylim = c(-3.5, max(hc_samples$height) * 1.08),
    cex.main = 1.2,
    cex.lab = 1.05,
    font.main = 2,
    font.lab = 2
  )
  points(seq_along(label_text), rep(-1.0, length(label_text)), pch = 15, col = label_col, cex = 1, xpd = NA)
  text(
    seq_along(label_text), -1.9,
    labels = gsub("_.*", "", label_text),
    srt = 90, adj = 1, cex = 0.55, col = label_col, xpd = NA, font = 2
  )
  legend(
    "topright", legend = c("ALL", "AML"), col = c("#4C78A8", "#E15759"),
    pch = 15, bty = "n", cex = 0.9, text.font = 2
  )
  par(old_par)
  dev.off()

  annotation_row <- data.frame(Cancer = expr_clean$cancer, row.names = rownames(expr_scaled))
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    png(file.path(figures_dir, "leukemia-heatmap.png"), width = 1600, height = 1280, res = 150)
    pheatmap::pheatmap(
      expr_scaled,
      cluster_rows = hc_samples,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      annotation_row = annotation_row,
      main = "Heatmap of Top 100 Variable Genes"
    )
    dev.off()
  } else {
    warning("Package 'pheatmap' is not installed; result tables and dendrogram were produced, but the heatmap was left unchanged.")
  }
}

cat("Analysis complete:", agreement, "of", nrow(expr_clean), "samples align with cluster majority labels.\n")

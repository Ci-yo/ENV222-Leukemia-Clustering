args <- commandArgs(trailingOnly = FALSE)
script_path <- gsub("\\\\", "/", sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
repo_root <- file.path(dirname(script_path), "..")

required <- c(
  "README.md",
  "METHODS.md",
  "data/README.md",
  "data/leukemia_train_cleaned.csv",
  "figures/project-banner.svg",
  "figures/method-pipeline.svg",
  "figures/leukemia-dendrogram.png",
  "figures/leukemia-heatmap.png",
  "presentation/project-walkthrough.qmd",
  "presentation/audio-walkthrough.mp3",
  "results/cluster_composition.csv",
  "results/analysis_summary.csv",
  "results/top_variable_genes.csv",
  "scripts/reproduce_analysis.R"
)
missing <- required[!file.exists(file.path(repo_root, required))]
if (length(missing) > 0) stop("Missing required files: ", paste(missing, collapse = ", "))

d <- read.csv(file.path(repo_root, "data", "leukemia_train_cleaned.csv"), check.names = FALSE)
stopifnot(nrow(d) == 38L, ncol(d) == 7131L)
stopifnot(sum(d$cancer == "ALL") == 27L, sum(d$cancer == "AML") == 11L)
stopifnot(sum(is.na(d)) == 0L, anyDuplicated(d$patient) == 0L)

x <- data.frame(
  lapply(d[, !(names(d) %in% c("patient", "cancer"))], as.numeric),
  check.names = FALSE
)
gene_variance <- apply(x, 2, var, na.rm = TRUE)
top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:100]
clusters <- cutree(hclust(dist(scale(x[, top_genes])), method = "ward.D2"), k = 2)
observed <- table(Cluster = clusters, Cancer = d$cancer)
stopifnot(observed[1, "ALL"] == 27L, observed[1, "AML"] == 1L)
stopifnot(observed[2, "ALL"] == 0L, observed[2, "AML"] == 10L)

majority <- tapply(d$cancer, clusters, function(values) names(sort(table(values), decreasing = TRUE))[1])
agreement <- sum(majority[as.character(clusters)] == d$cancer)
stopifnot(agreement == 37L)

summary_table <- read.csv(file.path(repo_root, "results", "analysis_summary.csv"))
stopifnot(summary_table$post_hoc_label_agreement[1] == 37L)
stopifnot(summary_table$agreement_denominator[1] == 38L)

forbidden_extensions <- c("docx", "pdf", "pptx", "zip", "m4a", "RData", "Rhistory")
repo_files <- list.files(repo_root, recursive = TRUE, all.files = TRUE, full.names = TRUE)
repo_files <- repo_files[!grepl("[\\/]\\.git([\\/]|$)", repo_files)]
bad_files <- repo_files[tolower(tools::file_ext(repo_files)) %in% tolower(forbidden_extensions)]
if (length(bad_files) > 0) stop("Forbidden draft/archive files found: ", paste(basename(bad_files), collapse = ", "))

text_files <- repo_files[grepl("\\.(md|R|qmd|yml|yaml|txt|cff)$", repo_files, ignore.case = TRUE)]
text_files <- gsub("\\\\", "/", text_files)
text_files <- text_files[basename(text_files) != basename(script_path)]
private_tokens <- c(
  paste0("225", "5737"),
  paste0("Xiyao", " Chen"),
  paste0("BEGIN OPENSSH ", "PRIVATE KEY")
)
for (path in text_files) {
  content <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (any(vapply(private_tokens, grepl, logical(1), x = content, fixed = TRUE))) {
    stop("Private identifier or key material found in: ", path)
  }
}

audio_size <- file.info(file.path(repo_root, "presentation", "audio-walkthrough.mp3"))$size
stopifnot(!is.na(audio_size), audio_size < 50 * 1024^2)

cat("Repository checks passed: data, clustering, results, assets and privacy guard are valid.\n")

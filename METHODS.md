# Methods and analytical guardrails

## Question

Can major structure in leukemia gene-expression profiles be recovered without using the known ALL/AML labels during clustering?

## Dataset

The analysis uses 38 training samples from a cleaned distribution of the Golub leukemia dataset:

- 27 acute lymphoblastic leukemia (ALL) samples
- 11 acute myeloid leukemia (AML) samples
- 7,129 numerical gene-expression features
- no missing values in the supplied cleaned table

The `patient` column is a sample index. `cancer` is retained only for the post-clustering interpretation audit.

## Analysis pipeline

1. Remove `patient` and `cancer` from the expression matrix.
2. Calculate variance for every expression feature.
3. Retain the 100 most variable genes.
4. Standardise each retained feature to place genes on a comparable scale.
5. Calculate pairwise Euclidean distances between samples.
6. Apply agglomerative hierarchical clustering with Ward's `ward.D2` linkage.
7. Cut the dendrogram into two clusters for a transparent label-composition audit.
8. Compare cluster membership with ALL/AML labels only after the unsupervised analysis is complete.

## Result audit

The two-cluster cut produces one branch containing 27 ALL samples and 1 AML sample, and a second branch containing 10 AML samples. Assigning each branch its majority label yields 37/38 label agreement.

That number is deliberately described as **post-hoc concordance**, not accuracy. The labels did not train a classifier, and the same small cohort is used for exploration and interpretation.

## Why top-variance filtering?

With thousands of features but only 38 samples, low-variation genes can dilute visible between-sample structure. The top-100 rule is a simple, auditable filter suited to this exploratory demonstration. It is not claimed to be the uniquely correct biological feature set.

## Limitations

- The sample size is small and no external cohort is evaluated.
- Gene selection, scaling, distance and linkage choices influence the dendrogram.
- A two-cluster cut is chosen to compare with two known labels; other cuts reveal different levels of the hierarchy.
- Cluster concordance does not establish diagnostic sensitivity, specificity or generalisability.
- The analysis does not control for batch effects or other technical and clinical covariates.

## Appropriate use

This repository supports teaching, portfolio review and reproducible exploration. It is not intended for patient diagnosis, prognosis, treatment selection or any other clinical decision.

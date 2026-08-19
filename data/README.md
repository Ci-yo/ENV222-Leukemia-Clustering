# Data note

## File

`leukemia_train_cleaned.csv` contains 38 sample rows and 7,131 columns:

- `patient`: sample index from 1 to 38
- `cancer`: known label (`ALL` or `AML`), used only after clustering
- 7,129 numerical gene-expression features

The table has 27 ALL samples, 11 AML samples and no missing values.

## Provenance

The cleaned training table was prepared from a public distribution of the Golub leukemia gene-expression dataset:

- Dataset page: https://www.kaggle.com/datasets/crawford/gene-expression/data
- Original study: https://doi.org/10.1126/science.286.5439.531

The repository includes the cleaned analysis table so that the exact figures and cluster audit can be reproduced. Consult the original dataset source for its terms and full documentation.

## Privacy and use boundary

The `patient` values are dataset sample indices rather than names or hospital identifiers. This educational repository does not contain direct personal identifiers. The data and results are not suitable for clinical decision-making.

## Transformations used here

The original preparation removed non-expression metadata and detection-call columns before arranging samples in rows and genes in columns. The repository analysis then selects the 100 highest-variance genes and standardises them before clustering; the checked-in CSV itself is not modified by the analysis script.

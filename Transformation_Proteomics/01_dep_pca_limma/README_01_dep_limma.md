# Proteomics Analysis Pipeline: Script 1 (DEP & Limma)

This script performs the foundational proteomics analysis for the study. It utilizes the `DEP` package to filter, normalize, and impute missing values from raw mass spectrometry TPM data. It then calculates differential expression and extracts Bayesian statistics using `limma`.

### Prerequisites
* R (>= 4.0)
* Bioconductor packages: `DEP`, `limma`, `SummarizedExperiment`
* CRAN packages: `dplyr`

### Inputs
1. **Raw TPM Data:** A CSV file containing raw protein intensities (e.g., `Cal_Transformation_ProteomicsTPM.csv`). Missing values should be represented as `NA` or `0`.
2. **Experimental Design:** A CSV file mapping sample names to experimental conditions (e.g., `CalProteomicsExpSummary.csv`).

### Outputs
* `_DEP_TestvControl.csv`: Standard differential expression results.
* `_limma_topTable_[Contrast].csv`: Bayesian statistics and posterior probabilities for each defined contrast.
* `_PCA_[Group].csv`: Raw coordinates for Principal Component Analysis plotting.
* `_Pearson_Correlation_Matrix.csv`: Sample-to-sample correlation matrix.

### Usage
Update the `FILE_RAW_TPM` and `FILE_EXP_DESIGN` variables in the User Configuration block, then run the script.
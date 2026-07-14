# Proteomics Analysis Pipeline: Script 1 (DEP & Limma)

This script performs the foundational proteomics analysis for the study. It utilizes the `DEP` package to filter, normalize, and impute missing values from raw mass spectrometry TPM data. It then calculates differential expression and extracts Bayesian statistics using `limma`.

### Prerequisites
* R (>= 4.0)
* Bioconductor packages: `DEP`, `limma`, `SummarizedExperiment`
* CRAN packages: `dplyr`

### Inputs
1. **Raw TPM Data:** A CSV file containing raw protein intensities (e.g., `Cal_Transformation_ProteomicsTPM.csv`). Missing values should be represented as `NA`.
2. **Experimental Design:** A CSV file mapping sample names to experimental conditions (e.g., `CalProteomicsExpSummary.csv`).

### Outputs
* `_DEP_TestvControl.csv`: DEP differential expression results comparing transformed samples to respective untreated preculture samples. VSN-normalized protein intensities labeled as "centered" values are used to create DEP_Centered.csv which is used as the input for 03_correlation_heatmap.R.
* `_limma_topTable_[Contrast].csv`: Bayesian statistics and posterior probabilities for each defined contrast. Up- and downregulated proteins are defined as Log2FC >= 1 or <= -1 with p.adj < 0.05. Up- and downregulated proteins from these files are used as inputs for 02_venn_diagrams.R and YNB vs YPD DAVID analysis for 04_go_dotplots.R.
* `_PCA_[Group].csv`: Raw coordinates for Principal Component Analysis plotting (Figure 3B, Figure S9B-C).
* `_Pearson_Correlation_Matrix.csv`: Sample-to-sample correlation matrix (Figure S9A).

### Usage
Update the `FILE_RAW_TPM` and `FILE_EXP_DESIGN` variables in the User Configuration block, then run the script.

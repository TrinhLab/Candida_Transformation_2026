# Proteomics Analysis Pipeline: Script 3 (Correlation Heatmap)

This script correlates individual protein expression profiles against phenotypic transformation efficiency data. It bins proteins based on their Pearson correlation coefficient and generates a clustered, annotated heatmap.

### Prerequisites
* R (>= 4.0)
* Bioconductor packages: `ComplexHeatmap`
* CRAN packages: `circlize`, `RColorBrewer`, `svglite`

### Inputs
* **Centered DEP Data:** A CSV file containing the centered protein intensities for the treated conditions (e.g., `DEP_Centered_Final.csv`).

### Outputs
* `_Correlation_Heatmap.svg`: A high-resolution, publication-ready heatmap clustered by correlation bins.
* `_Protein_Cluster_Assignments.csv`: A data frame listing every protein, its Pearson correlation value, and its assigned bin (e.g., "Strong Positive", "Moderate Negative").

### Usage
Update the `FILE_DEP_CENTERED` variable in the User Configuration block. If applying this to a new dataset, update the `trans_eff` vector with the corresponding phenotypic transformation efficiencies. Run the script to generate the heatmap and cluster assignments.
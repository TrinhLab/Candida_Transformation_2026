# Proteomics Analysis Pipeline: Script 4 (GO Term Dotplots)

This script visualizes Gene Ontology (GO) enrichment data. It takes raw CSV outputs from DAVID or Panther and generates publication-ready bubble/dotplots, scaling the points by Fold Enrichment and coloring them by statistical significance (-log10 P-value).

### Prerequisites
* R (>= 4.0)
* CRAN packages: `ggplot2`, `stringr`, `egg`, `svglite`

### Inputs
* **GO Enrichment CSVs:** Output files from DAVID Bioinformatics Resources (e.g., `StrongPositive_DAVID_DotPlot.csv`). The script expects columns for `Term`, `GeneRatio`, `FoldEnrichment`, and `Pvalue`.

### Outputs
* `_[Condition]_Dotplot.svg`: High-resolution, standardized SVG plots for each input file.

### Usage
Update the input file variables (e.g., `FILE_GO_POS`, `FILE_GO_NEG`) in the User Configuration block. The script contains two functions (`generate_dotplot` and `generate_dotplot_YNBvYPD`) with slightly different scaling limits to accommodate different enrichment magnitudes. Run the script to export the plots.
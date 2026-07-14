# Proteomics Analysis Pipeline: Correlation Heatmap
# Author: Christopher Cotter
# Manuscript: Cellular Stress Tolerance Governs Genetic Transformability in Recalcitrant Candida Species

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(svglite)
})

# User Configuration
FILE_DEP_CENTERED <- 'DEP_Centered.csv'
OUT_PREFIX <- "Cal_Proteomics_"

print("Generating Correlation Heatmap...")
centered_mat <- read.table(FILE_DEP_CENTERED, header = TRUE, sep = ',', na.strings = "NA")
protein_ids_all <- centered_mat$Protein.ID
centered_mat_numeric <- centered_mat[, setdiff(colnames(centered_mat), "Protein.ID")]

# Transformation efficiencies for 4 treated conditions
trans_eff <- c(YNB_pH5 = 968.7, YNB_TE = 289.3, YPD_pH5 = 57.3, YPD_TE = 19.3)
trans_eff <- log10(trans_eff)
treated_mat <- as.matrix(centered_mat_numeric[, names(trans_eff)])

# Calculate Pearson correlation
prot_vs_eff <- apply(treated_mat, 1, function(profile) cor(as.numeric(profile), trans_eff, method = "pearson"))

bin_breaks <- c(-1, -0.5, -0.3, 0.3, 0.5, 1)
bin_labels <- c("Strong Negative", "Moderate Negative", "Weak", "Moderate Positive", "Strong Positive")
protein_bins <- cut(prot_vs_eff, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)

valid_proteins <- !is.na(protein_bins)
treated_mat <- treated_mat[valid_proteins, ]
protein_bins <- protein_bins[valid_proteins]
prot_vs_eff <- prot_vs_eff[valid_proteins]
protein_ids <- protein_ids_all[valid_proteins]

heatmap_mat_scaled <- t(scale(t(treated_mat)))
column_order <- names(sort(trans_eff, decreasing = FALSE))
heatmap_mat_scaled_ordered <- heatmap_mat_scaled[, column_order]

annotation_col <- data.frame(TransformationEfficiency = trans_eff)
rownames(annotation_col) <- names(trans_eff)
annotation_col_ordered <- annotation_col[column_order, , drop = FALSE]

row_order <- unlist(lapply(levels(protein_bins), function(bin) {
  inds <- which(protein_bins == bin)
  inds[order(prot_vs_eff[inds], decreasing = FALSE)]
}))

heatmap_mat_clustered <- heatmap_mat_scaled_ordered[row_order, ]
protein_bins_ordered <- protein_bins[row_order]
protein_ids_ordered <- protein_ids[row_order]

col_fun <- colorRamp2(breaks = c(min(heatmap_mat_clustered), 0, max(heatmap_mat_clustered)), colors = c("#d5bb59", "#d3d3d1", "#e43628"))
bin_colors <- brewer.pal(length(levels(protein_bins_ordered)), "Set3")
row_anno <- rowAnnotation(Cluster = protein_bins_ordered, col = list(Cluster = setNames(bin_colors, levels(protein_bins_ordered))))

blue_col_fun <- colorRamp2(breaks = range(annotation_col_ordered$TransformationEfficiency, na.rm = TRUE), colors = c("white", "#3182bd"))
col_anno <- HeatmapAnnotation(TransformationEfficiency = annotation_col_ordered$TransformationEfficiency, col = list(TransformationEfficiency = blue_col_fun))
column_split <- factor(colnames(heatmap_mat_clustered), levels = column_order)

svglite::svglite(paste0(OUT_PREFIX, "_Correlation_Heatmap.svg"), width = 7, height = 5)
ht <- Heatmap(
  heatmap_mat_clustered, name = "z-score", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = protein_bins_ordered, column_split = column_split, column_gap = unit(1, "mm"),
  show_row_names = FALSE, show_column_names = TRUE, top_annotation = col_anno, left_annotation = row_anno,
  heatmap_legend_param = list(title = "Protein Intensity"), row_title = NULL,
  column_title = "Proteins Binned by Correlation with Transformation Efficiency",
  column_names_rot = 0, column_names_side = "bottom", column_names_gp = gpar(fontsize = 12, just = "center")
)
draw(ht)
dev.off()

cluster_df <- data.frame(ProteinID = protein_ids_ordered, CorrelationBin = as.character(protein_bins_ordered), CorrelationValue = prot_vs_eff[row_order])
write.csv(cluster_df, file = paste0(OUT_PREFIX, "_Protein_Cluster_Assignments.csv"), row.names = FALSE)
print("Script 3 Complete.")

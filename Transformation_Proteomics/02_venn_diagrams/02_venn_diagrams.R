# Proteomics Analysis Pipeline: Venn Diagrams
# Author: Christopher Cotter
# Manuscript: Cellular Stress Tolerance Governs Genetic Transformability in Recalcitrant Candida Species

suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
})

# User Configuration
FILE_UP_DEPS <- 'UpregulatedDEPs_limma.csv'
FILE_DOWN_DEPS <- 'DownregulatedDEPs_limma.csv'
OUT_PREFIX <- "Cal_Proteomics_"

print("Generating Venn Diagrams...")

generate_venn_and_sets <- function(input_file, out_name, palette) {
  prots <- read.csv(input_file, header = TRUE, sep = ',')
  
  A <- prots$YPD_TE_vs_YPD
  B <- prots$YNB_pH5_vs_YNB_A
  C <- prots$YPD_pH5_vs_YPD
  D <- prots$YNB_TE_vs_YNB_A
  
  A <- A[!is.na(A) & A != ""]
  B <- B[!is.na(B) & B != ""]
  C <- C[!is.na(C) & C != ""]
  D <- D[!is.na(D) & D != ""]
  
  # Draw Venn
  vennplot <- draw.quad.venn(
    area1 = length(A), area2 = length(B), area3 = length(C), area4 = length(D),
    n12 = length(intersect(A, B)), n13 = length(intersect(A, C)), n14 = length(intersect(A, D)),
    n23 = length(intersect(B, C)), n24 = length(intersect(B, D)), n34 = length(intersect(C, D)),
    n123 = length(Reduce(intersect, list(A, B, C))), n124 = length(Reduce(intersect, list(A, B, D))),
    n134 = length(Reduce(intersect, list(A, C, D))), n234 = length(Reduce(intersect, list(B, C, D))),
    n1234 = length(Reduce(intersect, list(A, B, C, D))),
    category = c("YPD TE pH 8","YNB AA pH 5", "YPD AA pH 5", "YNB TE pH 8"),
    fill = palette, cex = 3, cex.axis = 3, cex.lab = 3, cat.cex = 2.5,
    cat.fontface = "plain", cat.fontfamily = "Arial", fontface = "plain", fontfamily = "Arial",
    cat.just = list(c(0.4, 0.6), c(0.65, 0.6), c(0.5, 0.8), c(0.5, 0.9))
  )
  
  png(paste0(OUT_PREFIX, "_", out_name, "_Venn.png"), width = 14, height = 10, units = "in", res = 300)
  grid.draw(vennplot)
  dev.off()
  
  # Extract Unique Sets
  all_sets <- list(
    A_only = setdiff(A, union(B, union(C, D))),
    B_only = setdiff(B, union(A, union(C, D))),
    C_only = setdiff(C, union(A, union(B, D))),
    D_only = setdiff(D, union(A, union(B, C))),
    AB_only = setdiff(intersect(A, B), union(C, D)),
    AC_only = setdiff(intersect(A, C), union(B, D)),
    AD_only = setdiff(intersect(A, D), union(B, C)),
    BC_only = setdiff(intersect(B, C), union(A, D)),
    BD_only = setdiff(intersect(B, D), union(A, C)),
    CD_only = setdiff(intersect(C, D), union(A, B)),
    ABC_only = setdiff(intersect(intersect(A, B), C), D),
    ABD_only = setdiff(intersect(intersect(A, B), D), C),
    ACD_only = setdiff(intersect(intersect(A, C), D), B),
    BCD_only = setdiff(intersect(intersect(B, C), D), A),
    ABCD = intersect(intersect(intersect(A, B), C), D)
  )
  
  max_len <- max(sapply(all_sets, length))
  all_sets <- lapply(all_sets, function(x) c(x, rep("", max_len - length(x))))
  write.csv(as.data.frame(all_sets), paste0(OUT_PREFIX, "_", out_name, "_UniqueSets.csv"), row.names = FALSE)
}

# Run for Up and Down
generate_venn_and_sets(FILE_UP_DEPS, "Upregulated", c("#BB5566", "#B69D71", "#E3DED4", "#F49E8B"))
generate_venn_and_sets(FILE_DOWN_DEPS, "Downregulated", c("#2E58A4", "#B69D71", "#E3DED4", "#70AFC7"))
print("Script 2 Complete.")

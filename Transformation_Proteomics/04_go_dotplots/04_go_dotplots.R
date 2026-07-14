# Proteomics Analysis Pipeline: GO Term Dotplots
# Author: Christopher Cotter
# Manuscript: Cellular Stress Tolerance Governs Genetic Transformability in Recalcitrant Candida Species

suppressPackageStartupMessages({
  library(ggplot2)
  library(stringr)
  library(egg)
  library(svglite)
})

# User Configuration
FILE_GO_POS <- 'StrongPositive_DAVID_DotPlot.csv'
FILE_GO_NEG <- 'StrongNegative_DAVID_DotPlot.csv'
FILE_YNBvYPD_Up <- 'YNBvYPD_DAVID_Up_DotPlot.csv'
FILE_YNBvYPD_Down <- 'YNBvYPD_DAVID_Down_DotPlot.csv'
OUT_PREFIX <- "Cal_Proteomics_"

print("Generating GO Term Dotplots...")

generate_dotplot <- function(input_file, out_name, high_color) {
  panther <- read.csv(input_file, header = TRUE)
  panther_ordered <- panther[order(-panther$GeneRatio), ]
  panther_top <- panther_ordered[panther_ordered$Pvalue < 0.05, ]
  panther_top$Term <- stringr::str_wrap(panther_top$Term, width = 40)
  panther_top$Term <- factor(panther_top$Term, levels = panther_top$Term[order(panther_top$GeneRatio)])
  
  p <- ggplot(panther_top, aes(x = GeneRatio, y = Term)) +
    geom_point(aes(size = FoldEnrichment, color = -log10(Pvalue))) +
    scale_color_gradient(low = "#d3d3d1", high = high_color, limits = c(1, 2.5), oob = scales::squish) +
    scale_size_area(breaks = c(1, 2, 3, 4), limits = c(1, 4), max_size = 8) +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    labs(title = "Biological Process", x = "Protein Ratio", y = "", size = "Fold Enrichment") +
    guides(color = guide_colorbar(title = "-log10(P-value)")) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(t = 50, r = 30, b = 60, l = 100),
      legend.text = element_text(size = 12), legend.title = element_text(size = 13)
    ) +
    coord_cartesian(clip = "off")
  
  fixed_plot <- set_panel_size(p, width = unit(4, "cm"), height = unit(7, "cm"))
  ggsave(paste0(OUT_PREFIX, "_", out_name, "_Dotplot.svg"), fixed_plot, width = 8, height = 8, dpi = 300, bg = "transparent")
}

generate_dotplot_YNBvYPD <- function(input_file, out_name, high_color) {
  panther <- read.csv(input_file, header = TRUE)
  panther_ordered <- panther[order(-panther$GeneRatio), ]
  panther_top <- panther_ordered[panther_ordered$Pvalue < 0.05, ]
  panther_top$Term <- stringr::str_wrap(panther_top$Term, width = 40)
  panther_top$Term <- factor(panther_top$Term, levels = panther_top$Term[order(panther_top$GeneRatio)])
  
  p <- ggplot(panther_top, aes(x = GeneRatio, y = Term)) +
    geom_point(aes(size = FoldEnrichment, color = -log10(Pvalue))) +
    scale_color_gradient(low = "#d3d3d1", high = high_color, limits = c(1, 4), oob = scales::squish) +
    scale_size_area(breaks = c(1, 3, 6, 9), limits = c(1, 9), max_size = 8) +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    labs(title = "Biological Process", x = "Protein Ratio", y = "", size = "Fold Enrichment") +
    guides(color = guide_colorbar(title = "-log10(P-value)")) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(t = 50, r = 30, b = 60, l = 100),
      legend.text = element_text(size = 12), legend.title = element_text(size = 13)
    ) +
    coord_cartesian(clip = "off")
  
  fixed_plot <- set_panel_size(p, width = unit(4, "cm"), height = unit(7, "cm"))
  ggsave(paste0(OUT_PREFIX, "_", out_name, "_Dotplot.svg"), fixed_plot, width = 8, height = 8, dpi = 300, bg = "transparent")
}

generate_dotplot(FILE_GO_POS, "Positive_Correlation", "#e43628")
generate_dotplot(FILE_GO_NEG, "Negative_Correlation", "#d5bb59")
generate_dotplot_YNBvYPD(FILE_YNBvYPD_Up, "Upregulated", "#e43628")
generate_dotplot_YNBvYPD(FILE_YNBvYPD_Down, "Downregulated", "#d5bb59")
print("Script 4 Complete.")

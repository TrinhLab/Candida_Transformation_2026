# Proteomics Analysis Pipeline: DEP, PCA, and Limma
# Author: Christopher Cotter
# Manuscript: Cellular Stress Tolerance Governs Genetic Transformability in Recalcitrant Candida Species

suppressPackageStartupMessages({
  library(DEP)
  library(limma)
  library(SummarizedExperiment)
  library(dplyr)
})

# User Configuration
FILE_RAW_TPM <- 'Cal_Transformation_ProteomicsTPM.csv'
FILE_EXP_DESIGN <- 'CalProteomicsExpSummary.csv'
OUT_PREFIX <- "Cal_Proteomics_"

print("Starting DEP Analysis...")
Cal_prot <- read.table(FILE_RAW_TPM, header = TRUE, sep = ',', na.strings = c("NA"))
expDesign <- read.csv2(FILE_EXP_DESIGN, header = TRUE, sep = ',')

data_unique <- make_unique(Cal_prot, 'Accession', 'Accession', delim = ';')
data_se <- make_se(data_unique, 2:19, expDesign)

data_filt <- filter_missval(data_se, thr = 2)
data_norm <- normalize_vsn(data_filt)
set.seed(123)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Test for differential expression
data_diff_all_contrasts <- test_diff(data_imp, type = "manual", test = c('YNB_pH5_vs_YNB_A', 'YNB_TE_vs_YNB_A', 'YPD_pH5_vs_YPD', 'YPD_TE_vs_YPD'))
data_diff_all_contrasts_media <- test_diff(data_imp, type = "manual", test = c('YNB_A_vs_YPD'))
data_diff_ImprovedvStandard <- test_diff(data_imp, type = "manual", test = c('YNB_pH5_vs_YPD_TE'))

# Denote significant proteins
dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(1))
dep_media <- add_rejections(data_diff_all_contrasts_media, alpha = 0.05, lfc = log2(1))
dep_ImprovedvStandard <- add_rejections(data_diff_ImprovedvStandard, alpha = 0.05, lfc = log2(1))

# Extract PCA Data for external plotting
pca_all <- plot_pca(dep, x = 1, y = 2, n = nrow(dep), plot = FALSE)
write.csv(pca_all, paste0(OUT_PREFIX, "_PCA_All_Samples.csv"), row.names = FALSE)

dep_ynb <- dep[, grep("YNB", dep$condition)]
pca_ynb <- plot_pca(dep_ynb, x = 1, y = 2, n = nrow(dep_ynb), plot = FALSE)
write.csv(pca_ynb, paste0(OUT_PREFIX, "_PCA_YNB_Only.csv"), row.names = FALSE)

dep_ypd <- dep[, grep("YPD", dep$condition)]
pca_ypd <- plot_pca(dep_ypd, x = 1, y = 2, n = nrow(dep_ypd), plot = FALSE)
write.csv(pca_ypd, paste0(OUT_PREFIX, "_PCA_YPD_Only.csv"), row.names = FALSE)

# Calculate and print PCA Percent Variance for Axis Labels
get_pca_variance <- function(se_object) {
  exprs <- as.matrix(assay(se_object))
  prot_var <- apply(exprs, 1, var)
  select <- order(prot_var, decreasing = TRUE)[seq_len(nrow(exprs))]
  pca <- prcomp(t(exprs[select, ]))
  percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  return(c(PC1 = percentVar[1], PC2 = percentVar[2]))
}

var_all <- get_pca_variance(dep)
cat("--- ALL SAMPLES ---\nPC1:", var_all["PC1"], "%\nPC2:", var_all["PC2"], "%\n\n")

var_ynb <- get_pca_variance(dep_ynb)
cat("--- YNB ONLY ---\nPC1:", var_ynb["PC1"], "%\nPC2:", var_ynb["PC2"], "%\n\n")

var_ypd <- get_pca_variance(dep_ypd)
cat("--- YPD ONLY ---\nPC1:", var_ypd["PC1"], "%\nPC2:", var_ypd["PC2"], "%\n")

# Extract Pearson Correlation Matrix for external plotting
cor_matrix <- plot_cor(dep, significant = TRUE, lower = 0, upper = 1, plot = FALSE)
cor_df <- as.data.frame(cor_matrix)
cor_df$Sample <- rownames(cor_df)
cor_df <- cor_df[, c("Sample", setdiff(colnames(cor_df), "Sample"))]
write.csv(cor_df, paste0(OUT_PREFIX, "_Pearson_Correlation_Matrix.csv"), row.names = FALSE)

# Summarize and write results
data_results <- get_results(dep)
data_results_media <- get_results(dep_media)
data_results_ImprovedvStandard <- get_results(dep_ImprovedvStandard)

write.csv(data_results, paste0(OUT_PREFIX, "_DEP_TestvControl.csv"), row.names = FALSE)
write.csv(data_results_media, paste0(OUT_PREFIX, "_DEP_YNBvYPD.csv"), row.names = FALSE)
write.csv(data_results_ImprovedvStandard, paste0(OUT_PREFIX, "_DEP_YNBpH5vYPDTE.csv"), row.names = FALSE)

# Limma Bayesian Stats
print("Running Limma Bayesian Stats...")
exprs_imp <- as.matrix(assay(data_imp))
stopifnot(!any(dim(exprs_imp) == 0))

cd <- as.data.frame(colData(data_imp))
cond_col_candidates <- c("condition","Condition","group","replicate","sample_group")
cond_name <- intersect(names(cd), cond_col_candidates)
condition <- factor(cd[[cond_name[1]]])
design <- model.matrix(~0 + condition)
colnames(design) <- sub("^condition", "", colnames(design))

contrasts_to_make <- c(
  "YNB_A_vs_YPD"     = "YNB_A - YPD",
  "YNB_pH5_vs_YNB_A" = "YNB_pH5 - YNB_A",
  "YNB_TE_vs_YNB_A"  = "YNB_TE - YNB_A",
  "YPD_pH5_vs_YPD"   = "YPD_pH5 - YPD",
  "YPD_TE_vs_YPD"    = "YPD_TE - YPD", 
  "YNB_pH5_vs_YPD_TE" = "YNB_pH5 - YPD_TE"
)
present_groups <- colnames(design)
contrasts_to_make <- contrasts_to_make[sapply(contrasts_to_make, function(x) {
  groups_in_formula <- unique(gsub(" *[\\+-] *", " ", gsub("\\s+"," ", x)))
  all(sapply(strsplit(groups_in_formula, " ")[[1]], function(g) g %in% present_groups))
})]

fit <- lmFit(exprs_imp, design)
contr.mat <- makeContrasts(contrasts = contrasts_to_make, levels = design)
fit2 <- contrasts.fit(fit, contr.mat)
fit2 <- eBayes(fit2)

for(ct in colnames(contr.mat)){
  tt <- topTable(fit2, coef = ct, number = Inf, sort.by = "none")
  tt$protein <- rownames(tt)
  if("B" %in% colnames(tt)){
    tt$postProb_limmaB <- plogis(tt$B)
    tt$lfdr_limmaB <- 1 - tt$postProb_limmaB
  } else {
    tt$postProb_limmaB <- NA
    tt$lfdr_limmaB <- NA
  }
  write.csv(tt, paste0(OUT_PREFIX, "_limma_topTable_", ct, ".csv"), row.names = FALSE)
}
print("Script 1 Complete.")

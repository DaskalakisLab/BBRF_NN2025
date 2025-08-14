########################################
# Region specific analysis for MDD
# Author: Artemis Iatrou
# Date: 2025-08-08
# Description:
#   This script performs differential expression, pathway enrichment,
#   and spatial registration for the region-specific analysis.
########################################

library(dplyr)
library(readr)
library(purrr)
library(variancePartition)

# ===== USER SETTINGS =====
regions <- c("DLPFC", "mPFC", "dACC")  
base_dir <- "/data"
spatial_code <- file.path(base_dir, "/SciencePaper/DLPFCSpatialEnrichment/spatialDLPFC/code/analysis/10_clinical_gene_set_enrichment")
source(file.path(spatial_code, "gene_set_enrichment_plot_complex.R"))
source(file.path(base_dir, "/BBRF/Scripts/Combined_analysis_functions_sex.R"))

# Parameters
k_list <- c(2, 7, 9, 16, 28)
names(k_list) <- paste0("k", sprintf("%02d", k_list))
spatial_reg_enrich_threshold <- 0.05
gsea_categories <- list(
  c("C2", "REACTOME"),
  c("C5", "BP"),
  c("C5", "MF"),
  c("C5", "CC")
)

# ===== FUNCTION TO PROCESS EACH REGION =====
process_region <- function(region) {
  message("Processing region: ", region)
  
  # === Load demo data ===
  counts_demo_file <- file.path(base_dir, "/BBRF/Bulk", region, "counts_demo.rdata")
  load(counts_demo_file)  # Should load e.g. Demo_dlpfc
  demo_var <- ls(pattern = "^Demo_")[1]
  demo_df <- get(demo_var)
  
  demo_df$Race <- demo_df$Race...34
  
  # === Load cell type proportions ===
  ctp_file <- file.path(base_dir, "/Science/Cell_Proportion_Y1234/results", 
                        paste0("deconv_prop_full_", region, ".csv"))
  ctp <- read_csv(ctp_file)
  colnames(ctp)[1] <- "SampleID"
  
  common_SampleID <- intersect(demo_df$SampleID, ctp$SampleID)
  merged_demo <- merge(
    demo_df %>% filter(SampleID %in% common_SampleID),
    ctp %>% filter(SampleID %in% common_SampleID),
    by = "SampleID", all.x = TRUE
  )
  
  # === Load SVs ===
  sv_file <- file.path(base_dir, "/BBRF/Bulk", region, paste0("SVs_15.RDS")) #here make sure to load the respective SV file.
  sv <- readRDS(sv_file)
  merged_demo <- cbind(merged_demo, sv)
  
  # Adjust diagnosis factor
  merged_demo$gDx <- as.factor(merged_demo$PrimaryDx)
  levels(merged_demo$gDx) <- c("Control", "MDD", "MDD")
  
  # === Load expression data ===
  expr_file <- file.path(base_dir, "/BBRF/Bulk", region, paste0("Voom", region, ".rds"))
  expr <- readRDS(expr_file)$E
  expr_ss <- expr[, colnames(expr) %in% merged_demo$SampleID]
  
  demo <- merged_demo[merged_demo$SampleID %in% colnames(expr_ss), ]
  
  # === BayesSpace files ===
  bayesSpace_registration_fn <- map(k_list, ~ 
                                      file.path(base_dir, "/SciencePaper/DLPFCSpatialEnrichment/spatialDLPFC/processed-data/rdata/spe/07_layer_differential_expression",
                                                paste0("modeling_results_BayesSpace_k", sprintf("%02d", .x), ".Rdata"))
  )
  
  bayes_anno_file <- file.path(base_dir, "/SciencePaper/DLPFCSpatialEnrichment/spatialDLPFC/processed-data/rdata/spe/08_spatial_registration/bayesSpace_layer_annotations.csv")
  
  bayesSpace_registration <- lapply(bayesSpace_registration_fn[-5], function(x) {
    get(load(x))
  })
  
  bayes_anno <- read.csv(bayes_anno_file) %>%
    dplyr::select(layer_combo, test = cluster, Annotation = bayesSpace)
  
  # === Run pipeline ===
  out_dir <- file.path(base_dir, "/BBRF/Bulk", region, "CombinedResults_svs")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  run_pipeline(
    demo_data = demo,
    expr_data = expr_ss,
    covariates = c("gDx", "Age", "Sex", "sv1", "sv2", "sv3"),
    covariates2 = c("gDx", "Age", "sv1", "sv2", "sv3"),
    out_dir = out_dir,
    bayesSpace_registration = bayesSpace_registration,
    bayes_anno = bayes_anno,
    spatial_reg_enrich_threshold = spatial_reg_enrich_threshold,
    gsea_categories = gsea_categories
  )
  
  # === Subanalysis ===
  sub_dir <- file.path(base_dir, "/BBRF/Bulk", region, "Subanalysis")
  dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)
  
  limma_analysis_split(
    anno = demo,
    expr = expr_ss,
    covariates = c("gDx", "Age", "Sex", "sv1", "sv2", "sv3"),
    out_dir = sub_dir
  )
}

# ===== RUN FOR ALL REGIONS =====
walk(regions, process_region)
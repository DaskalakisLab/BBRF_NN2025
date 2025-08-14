########################################
# Spatial Enrichment Analysis for MDD
# Author: Artemis Iatrou
# Date: 2025-08-08
# Description:
#   This script takes differential expression results for MDD_RE,
#   identifies significant genes, and performs spatial domain
#   enrichment analysis using BayesSpace registrations.
########################################

# =====================
# 0. Load Required Packages
# =====================
required_packages <- c(
  "dplyr", "limma", "fgsea", "msigdbr", "spatialLIBD",
  "ggplot2", "jaffelab", "here", "purrr", "bacon", "ggrepel"
)

lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
})

# =====================
# 1. User Parameters
# =====================
script_dir <- here("code", "analysis", "10_clinical_gene_set_enrichment")
data_dir   <- here("data", "BBRF", "Bulk", "AcrossRegions")
spatial_dir <- here("data", "SciencePaper", "DLPFCSpatialEnrichment", "spatialDLPFC")

# Script that contains gene_set_enrichment_plot_complex() and related helpers
source(file.path(script_dir, "gene_set_enrichment_plot_complex.R"))

# Differential expression results file
de_results_file <- file.path(data_dir, "results_gDxMDD_All.RDS")

# BayesSpace registration files
k_values <- c(2, 7, 9, 16, 28)
names(k_values) <- sprintf("k%02d", k_values)

bayesSpace_registration_files <- purrr::map(
  k_values,
  ~ file.path(spatial_dir, "processed-data", "rdata", "spe",
              "07_layer_differential_expression",
              paste0("modeling_results_BayesSpace_k", sprintf("%02d", .x), ".Rdata"))
)

# BayesSpace annotation file
bayes_anno_file <- file.path(spatial_dir, "processed-data", "rdata", "spe",
                             "08_spatial_registration", "bayesSpace_layer_annotations.csv")

# Output directories
output_dir <- file.path(data_dir, "SpEnrich")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis threshold
spatial_reg_enrich_threshold <- 0.05

# =====================
# 2. Load Data
# =====================
results <- readRDS(de_results_file)
results[["ensID"]] <- rownames(results)

# Load BayesSpace registration objects (skipping k=28 as per original code)
bayesSpace_registration <- lapply(bayesSpace_registration_files[-5], function(file) {
  get(load(file))
})

# Read annotation file
bayes_anno <- read.csv(bayes_anno_file) %>%
  select(layer_combo, test = cluster, Annotation = bayesSpace)

# =====================
# 3. Identify Significant Genes
# =====================
sig_indices <- which(results$adj.P.Val <= spatial_reg_enrich_threshold)

if (length(sig_indices) <= 2) {
  message("No significant genes found for spatial enrichment.")
  write("No significant genes found for spatial enrichment.",
        file.path(output_dir, "spatial_enrichment_log.txt"),
        append = TRUE)
  quit(save = "no")
}

sigGeneSign <- data.frame(
  gene = sub("\\..*", "", as.character(results$ensID[sig_indices])),
  change = sign(results$logFC[sig_indices]),
  stringsAsFactors = FALSE
)

mdd_gene_list <- list(MDD = sigGeneSign)

# =====================
# 4. Run Enrichment Analysis
# =====================
enriched <- map(mdd_gene_list, function(gene_list) {
  map(bayesSpace_registration, function(modeling_res) {
    gene_set_enrichment(
      gene_list = gene_list,
      modeling_results = modeling_res,
      model_type = "enrichment"
    ) %>%
      left_join(bayes_anno, by = "test") %>%
      mutate(test = factor(layer_combo,
                           levels = bayes_anno$layer_combo[
                             bayes_anno$layer_combo %in% layer_combo
                           ])) %>%
      select(-c(layer_combo, Annotation, fdr_cut, model_type))
  })
})

# Adjust p-values for all enrichment results
enriched <- map_depth(enriched, 2, ~ {
  .x$padj <- p.adjust(.x$Pval, method = "BH")
  .x
})

# Save enrichment object
saveRDS(enriched, file.path(output_dir, "enriched_MDD_RE.rds"))

# =====================
# 5. Summarize Gene Counts
# =====================
gene_enrichment_count <- map(bayesSpace_registration, function(r) {
  en_count <- get_gene_enrichment_count(r)
  rownames(en_count) <- bayes_anno$layer_combo[match(rownames(en_count), bayes_anno$test)]
  layer_order <- bayes_anno$layer_combo[bayes_anno$layer_combo %in% rownames(en_count)]
  en_count[layer_order, , drop = FALSE]
})

gene_list_count <- map(mdd_gene_list, get_gene_list_count)

# =====================
# 6. Generate Plots
# =====================
# Example: k09 plot
pdf(file.path(output_dir, "Enrich_MDD_k09_adj.pdf"), height = 8, width = 10)
gene_set_enrichment_plot_complex_mod(
  enriched$MDD$k09,
  gene_count_col = gene_list_count[["MDD"]],
  gene_count_row = gene_enrichment_count[["k09"]],
  anno_title_col = "n DE Genes",
  anno_title_row = "n Domain\nGenes"
)
dev.off()

# All datasets/ks
walk2(enriched, names(enriched), function(enrich_obj, ds_name) {
  pdf(file.path(output_dir, paste0("Enrich_MDD_", ds_name, ".pdf")),
      height = 8, width = 9)
  map2(enrich_obj, names(enrich_obj), function(data_k, k_name) {
    message(ds_name, " - ", k_name)
    print(
      gene_set_enrichment_plot_complex(
        data_k,
        gene_count_col = gene_list_count[[ds_name]],
        gene_count_row = gene_enrichment_count[[k_name]],
        anno_title_col = "n DE Genes",
        anno_title_row = "n Domain\nGenes"
      )
    )
  })
  dev.off()
})
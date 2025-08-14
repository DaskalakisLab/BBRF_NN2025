########################################
# LogFC Correlation Analysis: Iatrou vs Girgenti
# Author: Artemis Iatrou
# Date: 2025-08-08
# Description:
#   This script computes correlations between log fold changes
#   from Iatrou et al. and Girgenti et al. datasets across 
#   multiple brain regions, thresholds, and correlation methods.
#   Results are summarized in CSV and heatmap formats.
########################################

# =====================
# 0. Load Required Packages
# =====================
library(readr)
library(dplyr)
library(ggplot2)

# =====================
# 1. Define Parameters
# =====================
region_mcl <- c("RE", "DLPFC", "mPFC", "dACC")
girg_region <- c("OFC", "dACC", "sgPFC", "dlPFC")
method <- c("pearson", "spearman")
threshold <- c("genome", "nominal_Iatrou", "fdr_Iatrou", "fdr_both")

out_dir <- "output/Correlations" 
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =====================
# 2. Define Correlation Function
# =====================
corr_regions <- function(region_mcl, girg_region, method = "pearson", threshold = "genome") {
  
  # Load main dataset
  if(region_mcl == "RE") {
    df1 <- readRDS("data/AcrossRegions/results_gDxMDD_All.RDS")
  } else {
    df1 <- readRDS(paste0("data/Bulk/", region_mcl, "/CombinedResults_svs/all/results_gDxMDD_All.RDS"))
  }
  
  # Load Girgenti dataset
  df2 <- read_csv("data/Girgenti/MDD_Girgenti.csv", show_col_types = FALSE)
  df2 <- df2[, c(2, grep(pattern = girg_region, colnames(df2)))]
  colnames(df2) <- c("symbol", "logFC", "adj.P.Val")
  df2 <- df2[complete.cases(df2$logFC), ]
  
  # Apply thresholds
  if(threshold == "fdr_both") {
    df1 <- df1 %>% filter(adj.P.Val < 0.05)
    df2 <- df2 %>% filter(adj.P.Val < 0.05)
  } else if(threshold == "nominal_Iatrou") {
    df1 <- df1 %>% filter(P.Value < 0.05)
  } else if(threshold == "fdr_Iatrou") {
    df1 <- df1 %>% filter(adj.P.Val < 0.05)
  }
  
  # Keep overlapping genes
  keep_genes <- intersect(df1$symbol, df2$symbol)
  df1_ss <- df1[match(keep_genes, df1$symbol), ]
  df2_ss <- df2[match(keep_genes, df2$symbol), ]
  
  # Compute correlation
  logFC_corr <- cor(df1_ss$logFC, df2_ss$logFC, method = method)
  
  # Return results as list
  list(
    our_brain = region_mcl,
    Girg_brain = girg_region,
    threshold = threshold,
    method = method,
    logFC_corr = logFC_corr,
    nr_genes = length(keep_genes)
  )
}

# =====================
# 3. Loop Over All Combinations
# =====================
all_results <- list()
counter <- 1

for (mcl in region_mcl) {
  for (girg in girg_region) {
    for (meth in method) {
      for (thresh in threshold) {
        all_results[[counter]] <- corr_regions(mcl, girg, meth, thresh)
        counter <- counter + 1
      }
    }
  }
}

# =====================
# 4. Save Results as CSV
# =====================
results_df <- do.call(rbind, lapply(all_results, as.data.frame))
write_csv(results_df, file.path(out_dir, "Girgenti_Artemis_Main.csv"))

# =====================
# 5. Prepare Heatmap Data
# =====================
correlation_results <- results_df %>%
  filter(method == "spearman", Girg_brain %in% c("OFC", "sgPFC")) %>%
  mutate(comp = paste0(our_brain, "~", Girg_brain),
         threshold = factor(threshold, levels = c("genome", "nominal_Iatrou", "fdr_Iatrou", "fdr_both")))

heatmap_data <- correlation_results %>%
  select(comp, threshold, logFC_corr, nr_genes, our_brain) %>%
  mutate(comp = factor(comp, levels = unique(comp))) %>%
  arrange(threshold, comp)

heatmap_data$group <- factor(heatmap_data$our_brain, levels = c("RE", "mPFC", "DLPFC", "dACC"))

# =====================
# 6. Generate Heatmap
# =====================
heatmap_plot <- ggplot(heatmap_data, aes(x = comp, y = threshold, fill = logFC_corr)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0("r: ", round(logFC_corr, 2), "\nN: ", nr_genes)), size = 4) +
  scale_fill_gradient2(
    low = "white",
    mid = "#ffb3b3",
    high = "#f23a3a",
    midpoint = 0.5,
    limits = c(0, 1),
    breaks = c(0, 0.3, 0.6, 1),
    labels = c("0", "0.3", "0.6", "1"),
    name = "Rs"
  ) +
  labs(x = "", y = "Threshold", title = "LogFC Correlation by Threshold: Iatrou vs Girgenti") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ group, scales = "free", space = "fixed")

# Save heatmap
ggsave(filename = file.path(out_dir, "Girgenti_RE_correlation_heatmap.pdf"),
       plot = heatmap_plot,
       width = 10, height = 6)
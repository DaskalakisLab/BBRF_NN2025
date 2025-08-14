########################################
########## Overlap & Correlation #######
########################################

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# Load your DEG results
df <- readRDS("/data/BBRF/Finalized_output/singlecell/Volcanos/combined_df.rds")
df <- df %>% filter(adj.P.Val < 0.05)

# Load Girgenti et al. DEG table
gf <- read_csv("/data/BBRF/Girgenti_sc/Supp_Data_Table_19_MDD_specific_DEG.csv")

# Inspect cell type distributions
table(gf$Celltype)
table(df$CellType)

# Standardize cell type labels
levels(df$CellType) <- c("EXN", "IN", "AST", "MG", "OLG", "END", "FB", "MURAL", "OPC", "TCELLS", "MACRO")

#########################################################
# Compute gene overlaps per cell type
#########################################################

shared_celltypes <- intersect(unique(df$CellType), unique(gf$Celltype))

overlap_summary <- lapply(shared_celltypes, function(ct) {
  
  df_genes <- df %>% filter(CellType == ct) %>% pull(genes) %>% unique()
  gf_genes <- gf %>% filter(Celltype == ct) %>% pull(Genename) %>% unique()
  
  intersect_genes <- intersect(df_genes, gf_genes)
  
  data.frame(
    CellType = ct,
    DF_Sig = length(df_genes),
    GF_Sig = length(gf_genes),
    Shared = length(intersect_genes),
    Percent_DF = length(intersect_genes) / length(df_genes) * 100,
    Percent_GF = length(intersect_genes) / length(gf_genes) * 100
  )
})

overlap_df <- bind_rows(overlap_summary)
print(overlap_df)

#########################################################
# Merge tables for correlation
#########################################################

# Clean Girgenti table for consistency
gf_clean <- gf %>%
  rename(
    gene = Genename,
    celltype = Celltype,
    gf_logFC = `MAST log2FC`
  )

# Merge with your DEGs
merged_df <- df %>%
  rename(gene = genes, logFC_df = logFC, celltype = CellType) %>%
  inner_join(gf_clean, by = c("gene", "celltype"))

#########################################################
# Correlation analysis
#########################################################

cor_test <- cor.test(merged_df$logFC_df, merged_df$gf_logFC)
print(cor_test)


########################################
# MAGMA Gene-Based Analysis for MDD 2025
# Author: Artemis Iatrou
# Date: 2025-08-14
# Description:
#   Prepares SNP locations and summary statistics for MAGMA,
#   annotates genes, runs gene-based analysis, and performs
#   cell-subclass enrichment analyses.
########################################

# =====================
# 0. Load Packages & Set Working Directory
# =====================
library(R.utils)
library(readr)
library(dplyr)

setwd("/data/BBRF/GWAS/MAGMA_output/")

# =====================
# 1. Load GWAS Summary Statistics
# =====================
mdd2025 <- read_tsv(
  '/data/GWAS/PGC_2025/pgc-mdd2025_no23andMe_eur_v3.49.24.11.tsv.gz',
  comment = '##'
) %>% as.data.frame()

# =====================
# 2. Prepare MAGMA Input Files
# =====================

# SNP locations
snp_loc <- mdd2025 %>%
  dplyr::select(ID, `#CHROM`, POS) %>%
  setNames(c("SNP", "CHR", "BP"))

write.table(
  snp_loc,
  "Input/mdd_snp_loc",
  sep = " ",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# SNP P-values
pval <- mdd2025 %>%
  dplyr::select(ID, PVAL) %>%
  setNames(c("SNP", "P"))

write.table(
  pval,
  "Input/mdd_sumstat_pval",
  sep = " ",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Gene location file
gene_loc <- "~/MAGMA/cc_data/NCBI37.3.gene.loc"

# =====================
# 3. Annotate SNPs to Genes
# =====================
cmd_annotate <- paste0(
  "~/MAGMA/magma --annotate ",
  "--snp-loc Input/mdd_snp_loc ",
  "--gene-loc ", gene_loc, " ",
  "--out Results/snp_loc_annotate"
)
system(cmd_annotate)

# =====================
# 4. Gene-Based Analysis
# =====================
cmd_gene_analysis <- paste0(
  "~/MAGMA/magma --bfile ~/MAGMA/cc_data/g1000_eur ",
  "--pval Input/mdd_sumstat_pval N=412305 ",
  "--gene-annot Results/snp_loc_annotate.genes.annot ",
  "--out Results/mdd_gene_based_analysis"
)
system(cmd_gene_analysis)

# =====================
# 5. Cell-Subclass Enrichment Analyses
# =====================
gene_results <- "Results/mdd_gene_based_analysis.genes.raw"

# Final cell subclasses
cmd_subclasses_final <- paste0(
  "~/MAGMA/magma --gene-results ", gene_results, " ",
  "--set-annot ~/MAGMA/InputFiles/BBRF/Final_cell_subclasses.txt ",
  "--out Results/mdd/mdd_final_subclasses"
)
system(cmd_subclasses_final)

# Finer cell subclasses
cmd_subclasses_finer <- paste0(
  "~/MAGMA/magma --gene-results ", gene_results, " ",
  "--set-annot ~/MAGMA/InputFiles/BBRF/Final_finer_cell_subclasses.txt ",
  "--out Results/mdd/mdd_finer_final_subclasses"
)
system(cmd_subclasses_finer)

# January-specific subclasses
cmd_subclasses_jan <- paste0(
  "~/MAGMA/magma --gene-results ", gene_results, " ",
  "--set-annot ~/MAGMA/InputFiles/BBRF/Final_cell_subclasses_Jan.txt ",
  "--out Results/mdd/mdd_final_subclasses_Jan"
)
system(cmd_subclasses_jan)


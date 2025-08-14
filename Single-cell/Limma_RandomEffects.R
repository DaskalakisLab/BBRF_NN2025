########################################
########## Normalize & DGE #############
########################################

library(calibrate)
library(limma)
library(edgeR)
library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(bacon)
library(sva)
library(variancePartition)

# Load helper functions
source("/data/humgen/daskalakislab/dipietro/SciencePaper/Code/RNA/Limma_Functions.R")
source("/data/BBRF/Scripts/Limma_4_sub.class_Functions.R")
source("/data/Science/Scripts/Functions4Limma_mcGill_2024.R")

# Set up model name and output directory
model_name <- model <- paste("LimmaResults_RE", "sex", "age", "version", "tissue", sep = "_")
output_dir <- paste0("/data/BBRF/Pseudobulk_results_May/", model_name, "/")
dir.create(output_dir, recursive = TRUE)
covariates <- c("sex", "age", "tissue", "version")

# Define cell types for analysis (astrocytes for now)
cts <- c("astro")

#########################################################
# Loop through each cell type
#########################################################

for(celltype in cts) {
  
  # Load pseudobulk metadata
  phenofile <- readRDS(paste0("/data/BBRF/Pheno_Demo/Pseudobulk_Dec/phenofile_subclass_", celltype, ".rds"))
  
  # Identify subclass columns
  start_col <- which(colnames(phenofile) == "classes") + 1
  end_col <- which(colnames(phenofile) == "total") - 1
  ct <- colnames(phenofile)[start_col:end_col]
  
  # Expand metadata for each subclass
  pheno_mod <- data.frame()
  for(i in seq_along(ct)){
    pheno2 <- phenofile
    pheno2$nam <- paste0(ct[i], "_", pheno2$sampleID)
    pheno2$time <- ct[i]
    pheno2$cells_num <- pheno2[, ct[i]]
    setDT(pheno2, keep.rownames = TRUE)
    pheno_mod <- if(i == 1) pheno2 else rbind(pheno_mod, pheno2)
  }
  rownames(pheno_mod) <- pheno_mod$nam
  
  # Load pseudobulk count matrices
  pwd11 <- paste0("/data/BBRF/Pseudobulk_Dec/December/Subclasses/", celltype)
  atble11 <- list.files(pwd11)
  
  for(i in seq_along(atble11)){
    pg <- paste0(pwd11,"/",atble11[i],"/")
    allfiles <- list.files(path = pg)
    dge <- readDGE(allfiles, path = pg, columns = c(1,2))
    
    brid <- sub(".gct", "", allfiles)
    colname <- paste0(atble11[i], "_", brid)
    colnames(dge) <- colname
    
    if(i == 1){
      expfin_all <- dge
    } else {
      expfin_all <- cbind(expfin_all, dge)
    }
  }
  
  # Match metadata with expression columns
  pheno_mod <- pheno_mod[match(colnames(expfin_all), rownames(pheno_mod)),]
  
  # Create DGEList and normalize
  temp <- DGEList(expfin_all)
  temp$samples$SampleID_cl <- rownames(temp$samples)
  temp$samples$group <- paste0(pheno_mod$gDx)
  
  pheno_mod$Dx <- factor(pheno_mod$gDx)
  pheno_mod$sample <- factor(pheno_mod$sampleID)
  pheno_mod$SampleID_cl <- pheno_mod$nam
  
  keep2 <- rowSums(cpm(temp) > 1) >= 70
  dge <- temp[keep2, keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  colnames(dge) <- pheno_mod$SampleID_cl
  
  # Run Limma with random effects
  run_RE_limma_analysis(
    group = "MDD",
    output_dir = output_dir,
    covariates = c("gDx", "time", covariates),
    celltype = celltype
  )
}

#########################################################
# Aggregate breakdown files
#########################################################

bkdown_directory <- paste0("/data/BBRF/Pseudobulk_results_Dec/", model, "/")
bkdown_files <- list.files(path = bkdown_directory, pattern = "breakdown.*\\.csv", full.names = TRUE)
data_list <- list()

for(f in bkdown_files){
  df <- read.csv(f, stringsAsFactors = FALSE)
  celltype <- gsub("DEGs_(.*)_breakdown.csv", "\\1", basename(f))
  df$celltype <- celltype
  data_list[[basename(f)]] <- df
}

combined_df <- do.call(rbind, data_list)
write.csv(combined_df, file = paste0("/data/BBRF/Pseudobulk_results_Dec/", model, "/limma_summary.csv"))
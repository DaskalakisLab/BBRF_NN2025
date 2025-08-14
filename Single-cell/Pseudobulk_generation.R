#########################################################
# Pseudobulk Generation for Classes and Subclasses
#########################################################

library(Seurat)
library(future)

# Increase memory allowance
options(future.globals.maxSize = 300000 * 1024^2)
plan("multicore", workers = 65)

# Define cell types to process
nn_cells <- c("ex0", "inh", "micro", "oligo", "opc", "astro")

#########################################################
# Loop through each cell type
#########################################################

for(cell in nn_cells){
  
  # Create directories for output
  dir.create(paste0("/data/BBRF/Pseudobulk_Dec/December/Subclasses/", cell, "/"), recursive = TRUE)
  dir.create(paste0("/data/BBRF/Pseudobulk_Dec/December/Classes/", cell, "/"), recursive = TRUE)
  
  # Load Seurat object (astro uses _DAS.rds)
  so_file <- if(cell == "astro"){
    paste0("/data/BBRF/Seurat_obj/Subclusters/", cell, "_tissue_final2_DAS.rds")
  } else {
    paste0("/data/BBRF/Seurat_obj/Subclusters/", cell, "_tissue_final2.rds")
  }
  so <- readRDS(so_file)
  
  #########################################################
  # Generate pseudobulk for Classes
  #########################################################
  
  setwd(paste0("/data/BBRF/Pseudobulk_Dec/December/Classes/", cell, "/"))
  
  class_levels <- names(table(so$classes))
  count_table <- table(so$sampleID, so$classes)
  count_mtx <- as.data.frame.matrix(count_table)
  count_mtx$sampleID <- rownames(count_mtx)
  count_mtx[,1:length(class_levels)] <- lapply(count_mtx[,1:length(class_levels)], as.numeric)
  count_mtx$total <- rowSums(count_mtx[,1:length(class_levels)])
  
  # Create metadata for pseudobulk
  metadata <- so@meta.data
  pheno <- metadata[!duplicated(metadata$sampleID),]
  pheno <- merge(pheno, count_mtx, by = "sampleID")
  saveRDS(pheno, paste0("/data/BBRF/Pheno_Demo/Pseudobulk_Dec/phenofile_class_", cell, ".rds"))
  
  # Generate per-sample .gct files for Classes
  mat <- GetAssayData(so, layer = "counts", assay = "RNA")
  for(class_id in class_levels){
    dir.create(paste0(getwd(), "/", class_id))
    d2 <- subset(metadata, metadata$classes == class_id)
    samples <- unique(d2$sampleID)
    
    for(s in samples){
      dtemp <- subset(d2, d2$sampleID == s)
      mat2 <- as.matrix(mat[, row.names(dtemp)])
      summed <- rowSums(mat2)
      dxx <- data.frame(Name = row.names(mat2), Value = summed)
      colnames(dxx)[2] <- paste0("S_", s)
      filename <- paste0(getwd(), "/", class_id, "/", s, ".gct")
      write.table(dxx, filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
    print(class_id)
  }
  
  #########################################################
  # Generate pseudobulk for Subclasses
  #########################################################
  
  setwd(paste0("/data/BBRF/Pseudobulk_Dec/December/Subclasses/", cell, "/"))
  
  subclass_levels <- names(table(so$subclasses))
  count_table <- table(so$sampleID, so$subclasses)
  count_mtx <- as.data.frame.matrix(count_table)
  count_mtx$sampleID <- rownames(count_mtx)
  count_mtx[,1:length(subclass_levels)] <- lapply(count_mtx[,1:length(subclass_levels)], as.numeric)
  count_mtx$total <- rowSums(count_mtx[,1:length(subclass_levels)])
  
  # Create metadata for pseudobulk
  pheno <- metadata[!duplicated(metadata$sampleID),]
  pheno <- merge(pheno, count_mtx, by = "sampleID")
  saveRDS(pheno, paste0("/data/BBRF/Pheno_Demo/Pseudobulk_Dec/phenofile_subclass_", cell, ".rds"))
  
  # Generate per-sample .gct files for Subclasses
  mat <- GetAssayData(so, layer = "counts", assay = "RNA")
  for(subclass_id in subclass_levels){
    dir.create(paste0(getwd(), "/", subclass_id))
    d2 <- subset(metadata, metadata$subclasses == subclass_id)
    samples <- unique(d2$sampleID)
    
    for(s in samples){
      dtemp <- subset(d2, d2$sampleID == s)
      mat2 <- as.matrix(mat[, row.names(dtemp)])
      summed <- rowSums(mat2)
      dxx <- data.frame(Name = row.names(mat2), Value = summed)
      colnames(dxx)[2] <- paste0("S_", s)
      filename <- paste0(getwd(), "/", subclass_id, "/", s, ".gct")
      write.table(dxx, filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
    print(subclass_id)
  }
}
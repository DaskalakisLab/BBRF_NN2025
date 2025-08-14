library(dplyr)
library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
library(tibble)
library(AnnotationDbi)

# Load your combined DE results
df <- readRDS("/data/BBRF/Finalized_output/singlecell/Volcanos/combined_df.rds")

# Prepare MSigDB gene sets — GO:BP or Hallmark as example
msigdb_bp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

# Convert to a list format for fgsea
pathways_bp <- split(msigdb_bp$gene_symbol, msigdb_bp$gs_name)

# Output directory
outdir <- "/data/BBRF/Finalized_output/singlecell/GSEA_MSigDB"
dir.create(outdir, showWarnings = FALSE)

# Loop over cell types
cell_types <- unique(df$CellType)

for (ct in cell_types) {
  message("Processing ", ct, "...")
  
  # Subset for current cell type
  df_ct <- df %>%
    filter(CellType == ct, !is.na(logFC)) 
  
  if (nrow(df_ct) < 100) {
    warning("Skipping ", ct, ": not enough data.")
    next
  }
  
  gsea_input<- df_ct %>%
    arrange(desc(logFC)) %>% 
    dplyr::select(genes, logFC)
  
  ranks<- deframe(gsea_input)
  
  #Run GSEA with 
  fgseaRes<- fgseaMultilevel(pathways_bp, stats = ranks, 
                             minSize=15, maxSize = 500,
                             eps = 0) 
  
  # Save results
  if (!is.null(fgseaRes)) {
    saveRDS(fgseaRes, file = file.path(outdir, paste0("fgsea_", ct, "_GO_BP.rds")))
  }
}

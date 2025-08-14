# Load required libraries
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# Define cell types to iterate over
celltypes <- c("oligo", "astro", "endo", "micro", "inh", "ex0")

# Loop over each cell type
for (celltype in celltypes) {
  df_path <- if(celltype == "astro") {
    paste0("/data/BBRF/Seurat_obj/Subclusters/", celltype, "_tissue_final2_DAS.rds")
  } else {
    paste0("/data/BBRF/Seurat_obj/Subclusters/", celltype, "_tissue_final2.rds")
  }
  df <- readRDS(df_path)
  
  # Print subclass summary
  print(table(df$subclasses))
  
  # Set identity to subclasses
  Idents(df) <- df$subclasses
  
  # Find markers
  all <- FindAllMarkers(df, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
  saveRDS(all, paste0("/data/BBRF/Finalized_outputs/sketch_clusters/Subclasses_markers/", celltype, "_markers_new.rds"))
  
  # Significant genes list by cluster
  significant_genes_list <- all %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC)) %>%
    group_by(cluster) %>%
    summarise(genes = list(gene)) %>%
    pull(genes)
  names(significant_genes_list) <- all %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC)) %>%
    group_by(cluster) %>%
    summarise(cluster = unique(cluster)) %>%
    pull(cluster)
  
  # GO enrichment per cluster
  for (cluster in names(significant_genes_list)) {
    sig_genes <- significant_genes_list[[cluster]]
    bckgrd <- rownames(df)
    ego <- enrichGO(
      gene          = sig_genes,
      universe      = bckgrd,
      keyType       = "SYMBOL",
      OrgDb         = org.Hs.eg.db,
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.01,
      qvalueCutoff  = 0.05,
      readable      = FALSE
    )
    saveRDS(ego@result, paste0("/data/BBRF/GO/ALL_", celltype, "_", cluster, "_new.RDS"))
  }
}

# Special handling: Endothelial and pericytes
celltype <- "endo"
df <- readRDS(paste0("/data/BBRF/Seurat_obj/Subclusters/", celltype, "_tissue_final2.rds"))
print(table(df$subclasses))
Idents(df) <- df$subclasses

Endo_cap <- FindMarkers(df, ident.1 = "Endo_cap", ident.2 = c("Endo_veinous", "Endo_arterial"), min.pct = 0.5, logfc.threshold = 1, only.pos = TRUE)
Endo_cap$cluster <- "Endo_cap"; Endo_cap$gene <- rownames(Endo_cap)
Endo_art <- FindMarkers(df, ident.1 = "Endo_arterial", ident.2 = c("Endo_veinous", "Endo_cap"), min.pct = 0.5, logfc.threshold = 1, only.pos = TRUE)
Endo_art$cluster <- "Endo_arterial"; Endo_art$gene <- rownames(Endo_art)
Endo_v <- FindMarkers(df, ident.1 = "Endo_veinous", ident.2 = c("Endo_arterial", "Endo_cap"), min.pct = 0.5, logfc.threshold = 1, only.pos = TRUE)
Endo_v$cluster <- "Endo_veinous"; Endo_v$gene <- rownames(Endo_v)

Peri_1 <- FindMarkers(df, ident.1 = "Pericytes_1", ident.2 = c("Pericytes_2", "Pericytes_3"), min.pct = 0.5, logfc.threshold = 0.75, only.pos = TRUE)
Peri_1$cluster <- "Pericytes_1"; Peri_1$gene <- rownames(Peri_1)
Peri_2 <- FindMarkers(df, ident.1 = "Pericytes_2", ident.2 = c("Pericytes_1", "Pericytes_3"), min.pct = 0.5, logfc.threshold = 1, only.pos = TRUE)
Peri_2$cluster <- "Pericytes_2"; Peri_2$gene <- rownames(Peri_2)
Peri_3 <- FindMarkers(df, ident.1 = "Pericytes_3", ident.2 = c("Pericytes_2", "Pericytes_1"), min.pct = 0.5, logfc.threshold = 1, only.pos = TRUE)
Peri_3$cluster <- "Pericytes_3"; Peri_3$gene <- rownames(Peri_3)

SMC <- FindMarkers(df, ident.1 = "SMC", ident.2 = c("Pericytes_1", "Pericytes_2"), min.pct = 0.5, logfc.threshold = 1, only.pos = TRUE)
SMC$cluster <- "SMC"; SMC$gene <- rownames(SMC)

Fibro_peri <- FindMarkers(df, ident.1 = "Fibro_peri", ident.2 = c("Fibro_men"), min.pct = 0.5, logfc.threshold = 1, only.pos = TRUE)
Fibro_peri$cluster <- "Fibro_peri"; Fibro_peri$gene <- rownames(Fibro_peri)
Fibro_men <- FindMarkers(df, ident.1 = "Fibro_men", ident.2 = c("Fibro_peri"), min.pct = 0.5, logfc.threshold = 1, only.pos = TRUE)
Fibro_men$cluster <- "Fibro_men"; Fibro_men$gene <- rownames(Fibro_men)

all <- rbind(Endo_art, Endo_cap, Endo_v, Peri_1, Peri_2, Peri_3, SMC, Fibro_men, Fibro_peri)
saveRDS(all, paste0("/data/BBRF/Finalized_outputs/sketch_clusters/Subclasses_markers/", celltype, "2_markers.rds"))

significant_genes_list <- all %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  summarise(genes = list(gene)) %>%
  pull(genes)
names(significant_genes_list) <- all %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  summarise(cluster = unique(cluster)) %>%
  pull(cluster)

for (cluster in names(significant_genes_list)) {
  sig_genes <- significant_genes_list[[cluster]]
  bckgrd <- rownames(df)
  ego <- enrichGO(
    gene          = sig_genes,
    universe      = bckgrd,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = FALSE
  )
  saveRDS(ego@result, paste0("/data/BBRF/GO/ALL_", celltype, "_", cluster, "_2.RDS"))
}

# Microglia example
celltype <- "micro"
df <- readRDS(paste0("/data/BBRF/Seurat_obj/Subclusters/", celltype, "_tissue_final2.rds"))
print(table(df$subclasses))
Idents(df) <- df$subclasses

Micro_homeo <- FindMarkers(df, ident.1 = "Micro_homeo", ident.2 = c("Micro_act"), min.pct = 0.5, logfc.threshold = 0.5, only.pos = TRUE)
Micro_homeo$cluster <- "Micro_homeo"; Micro_homeo$gene <- rownames(Micro_homeo)
Micro_act <- FindMarkers(df, ident.1 = "Micro_act", ident.2 = c("Micro_homeo"), min.pct = 0.5, logfc.threshold = 1, only.pos = TRUE)
Micro_act$cluster <- "Micro_act"; Micro_act$gene <- rownames(Micro_act)
Macrophages <- FindMarkers(df, ident.1 = "Macrophages", ident.2 = c("Micro_homeo", "Micro_act"), min.pct = 0.5, logfc.threshold = 1, only.pos = TRUE)
Macrophages$cluster <- "Macrophages"; Macrophages$gene <- rownames(Macrophages)
Tcells <- FindMarkers(df, ident.1 = "Tcells", ident.2 = c("Micro_homeo", "Micro_act", "Macrophages"), min.pct = 0.5, logfc.threshold = 1, only.pos = TRUE)
Tcells$cluster <- "Tcells"; Tcells$gene <- rownames(Tcells)

all <- rbind(Micro_homeo, Micro_act, Macrophages, Tcells)
saveRDS(all, paste0("/data/BBRF/Finalized_outputs/sketch_clusters/Subclasses_markers/", celltype, "2_markers.rds"))

significant_genes_list <- all %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  summarise(genes = list(gene)) %>%
  pull(genes)
names(significant_genes_list) <- all %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  summarise(cluster = unique(cluster)) %>%
  pull(cluster)

for (cluster in names(significant_genes_list)) {
  sig_genes <- significant_genes_list[[cluster]]
  bckgrd <- rownames(df)
  ego <- enrichGO(
    gene          = sig_genes,
    universe      = bckgrd,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = FALSE
  )
  saveRDS(ego@result, paste0("/data/BBRF/GO/ALL_", celltype, "_", cluster, "_2.RDS"))
}
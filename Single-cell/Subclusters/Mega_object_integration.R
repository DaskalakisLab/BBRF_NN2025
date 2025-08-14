########################################
########## Seurat Mega Project #########
########################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(tidyr)

# -----------------------------
# 1. Load metadata and subset
# -----------------------------
cells_keep <- read.csv("/data/BBRF/Pheno_Demo/combined_meta.csv")
mega_prj$CellName <- rownames(mega_prj@meta.data)
mega_prj_ss <- subset(mega_prj, subset = CellName %in% cells_keep$Cell_ID)

# Match metadata to Seurat object
cells_keep <- cells_keep[match(mega_prj_ss$CellName, cells_keep$Cell_ID),]
stopifnot(identical(rownames(mega_prj_ss@meta.data), cells_keep$Cell_ID))

mega_prj_ss$Classes <- cells_keep$Classes
mega_prj_ss$Subclasses <- cells_keep$Subclasses
mega_prj_ss <- DietSeurat(mega_prj_ss)

# -----------------------------
# 2. Normalize, variable features, sketch
# -----------------------------
mega <- NormalizeData(mega_prj_ss)
mega <- FindVariableFeatures(mega)
mega <- SketchData(mega, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(mega) <- "sketch"

mega <- FindVariableFeatures(mega)
mega <- ScaleData(mega, vars.to.regress = c("percent.mt", "nCount_RNA"))
mega <- RunPCA(mega, npcs = 100)

# -----------------------------
# 3. Harmony integration
# -----------------------------
mega <- RunHarmony(
  mega, 
  reduction.save = "Harmony_", 
  group.by.vars = c("tissue", "orBatch", "sampleID", "version"), 
  plot_convergence = FALSE
)

# -----------------------------
# 4. UMAP, neighbors, clusters
# -----------------------------
mega <- RunUMAP(mega, dims = 1:100, reduction = "Harmony_")
mega <- FindNeighbors(mega, reduction = "Harmony_", k.param = 100, nn.eps = 0)
mega <- FindClusters(mega, algorithm = 1, resolution = 0.1)

# -----------------------------
# 5. UMAP plots - categoricals
# -----------------------------
reductions <- grep("^umap", names(mega@reductions), value = TRUE)
categoricals <- c("orBatch", "version", "sex", "primarydx", "tissue", "sketch_snn_res.0.1", "Classes", "Subclasses")

plots_list <- lapply(categoricals, function(cat) {
  DimPlot(mega, reduction = reductions[1], group.by = cat)
})

pdf("/data/UMAPS_categoricals.pdf", width = 12, onefile = TRUE)
plots_list
dev.off()

# -----------------------------
# 6. UMAP plots - continuous
# -----------------------------
continuous_vars <- c("nFeature_RNA", "nCount_RNA", "pmi", "age", "percent.mt")

plots_list <- lapply(continuous_vars, function(var) {
  FeaturePlot(mega, features = var, reduction = reductions[1])
})

pdf("/data/UMAPS_continuous.pdf", width = 12, height = 16, onefile = TRUE)
plots_list
dev.off()

# -----------------------------
# 7. Marker gene plots
# -----------------------------
feat_plot <- c("GFAP","SLC1A2","AQP4","SLC14A1","VWF","CLDN5","FLT1","GAD1","GAD2","SATB2",
               "SLC17A7","NRGN","SNAP25","MOBP","MOG","TF","DOCK8","CSF1R","C3","VCAN","PDGFRA",
               "CSPG4","ACTA2","RGS5","PDGFRB")

pdf("/data/MarkerGenes_canonical.pdf", width = 20, height = 20)
VlnPlot(mega, group.by = "sketch_snn_res.0.1", features = feat_plot, pt.size = 0, ncol = 5)
FeaturePlot(mega, features = feat_plot, ncol = 5)
dev.off()

# -----------------------------
# 8. Save Seurat object
# -----------------------------
saveRDS(mega, file = "/data/combined_subset_sketched.rds")

# -----------------------------
# 9. Projected Seurat object with Azimuth
# -----------------------------
mega_prj <- ProjectData(
  object = mega,
  assay = "RNA",
  full.reduction = "Harmony.full",
  sketched.assay = "sketch",
  sketched.reduction = "Harmony_",
  umap.model = "umap",
  dims = 1:100,
  refdata = list(cluster_01 = "sketch_snn_res.0.1")
)

DefaultAssay(mega_prj) <- "RNA"
Idents(mega_prj) <- mega_prj$cluster_01

mega_prj <- RunAzimuth(mega_prj, reference = "/path/to/Azimuth_annotation")
mega_prj$fin_classes <- Idents(mega_prj)
levels(mega_prj$fin_classes) <- c("Excitatory","Inhibitory","Oligo","Microglia","Astrocytes","Inhibitory","Perivascular","OPC","Excitatory")

# -----------------------------
# 10. UMAP plots - projected annotations
# -----------------------------
colors <- c("#ff2d55","#007aff","#5856d6","#ffcc00","#ff9500","#4cd964","#7D61D6")
pdf("/data/umap_projected.pdf", width = 8, height = 8)
DimPlot(mega_prj, label = TRUE, reduction = "ref.umap", group.by = "fin_classes",
        label.size = 6, label.box = TRUE, repel = TRUE, label.color = "white", cols = colors) + NoLegend() + ggtitle(NULL)
dev.off()

# -----------------------------
# 11. Cell type proportion plots
# -----------------------------
plotCellTypeProps_mod <- function(x, clusters, sample){
  prop.list <- getTransformedProps(clusters, sample)
  Proportions <- as.vector(t(prop.list$Proportions))
  Samples <- rep(colnames(prop.list$Proportions), nrow(prop.list$Proportions))
  Clusters <- rep(rownames(prop.list$Proportions), each = ncol(prop.list$Proportions))
  plotdf <- data.frame(Samples = Samples, Clusters = Clusters, Proportions = Proportions)
  colors <- c("#ff9500","#ff2d55","#007aff","#ffcc00","#5856d6","#7D61D6","#4cd964")
  
  ggplot(plotdf, aes(x = Samples, y = Proportions, fill = Clusters)) + 
    geom_bar(stat = "identity") + 
    theme_bw() +
    scale_fill_manual(values = colors) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 14))
}

pdf("/data/celltype_props_samples.pdf", width = 20)
plotCellTypeProps_mod(mega_prj, clusters = mega_prj$fin_classes, sample = mega_prj$sampleID)
dev.off()

# -----------------------------
# 12. Pathway module scores
# -----------------------------
pathways <- readRDS("/path/to/top20_pathways.rds")
mega_prj <- AddModuleScore(mega_prj, features = pathways, name = "path_")
Idents(mega_prj) <- mega_prj$fin_classes

for(i in seq_along(pathways)){
  pdf(paste0("/data/pathway_", names(pathways)[[i]], ".pdf"))
  VlnPlot(mega_prj, features = paste0("path_", i), pt.size = 0)
  FeaturePlot(mega_prj, features = paste0("path_", i), pt.size = 0, reduction = "ref.umap") & NoLegend()
  dev.off()
}

# -----------------------------
# 13. Save final projected object
# -----------------------------
saveRDS(mega_prj, file = "/data/Dec_mega_prj.rds")
########################################
# Microglia Subcluster Processing - oFC and dlPFC
# Author: Artemis Iatrou
# Date: 2025-08-14
# Description:
#   Subsets microglia by tissue, performs SCTransform, PCA, Harmony integration,
#   UMAP, clustering, Azimuth annotation, marker gene analysis, and saves outputs.
########################################

# =====================
# 0. Load Packages & Set Options
# =====================
library(Seurat)
library(future)
library(speckle)
options(future.globals.maxSize = 1000000 * 1024^2)

# =====================
# 1. Load and Prepare Seurat Object
# =====================
micro_sub <- readRDS("/data/BBRF/Seurat_obj/Subclusters/micro_sub2.RDS")
micro_sub <- DietSeurat(micro_sub)

micro_sub_oFC <- subset(micro_sub, subset = tissue == "oFC")
micro_sub_dlpfc <- subset(micro_sub, subset = tissue == "dlPFC")
mito <- micro_sub_oFC$percent.mt

# Assign dataset and sample IDs
micro_sub_oFC$dataset <- "eb"
micro_sub_dlpfc$dataset <- "ai"
micro_sub_oFC$sampleID <- micro_sub_oFC$Donor
micro_sub_dlpfc$sampleID <- micro_sub_dlpfc$sample

# Assign additional metadata
micro_sub_oFC$tissue <- "oFC"
micro_sub_dlpfc$tissue <- "dlPFC"
micro_sub_oFC$pr_clusters <- micro_sub_oFC$celltypes_final
micro_sub_oFC$pmi <- micro_sub_oFC$PMI
micro_sub_oFC$orBatch <- "Binder_oFC"
micro_sub_oFC$version <- "v3"
micro_sub_oFC$sex <- factor(micro_sub_oFC$Sex, levels = c("F", "M"), labels = c("female", "male"))
micro_sub_oFC$age <- micro_sub_oFC$Age
micro_sub_oFC$gDx <- factor(micro_sub_oFC$Classification, levels = c("Control", "MDD"))

micro_sub_dlpfc$gDx <- factor(micro_sub_dlpfc$primarydx, levels = c("Control", "MDD", "MDD"))

# =====================
# 2. SCTransform, PCA, Harmony
# =====================
micro_sub_oFC <- SCTransform(micro_sub_oFC, vars.to.regress = c("percent.mt", "nCount_RNA"))
micro_sub_dlpfc <- SCTransform(micro_sub_dlpfc, vars.to.regress = c("percent.mt", "nCount_RNA"))

micro_sub_oFC <- RunPCA(micro_sub_oFC, npcs = 50)
micro_sub_dlpfc <- RunPCA(micro_sub_dlpfc, npcs = 50)

set.seed(1234)
micro_sub_oFC <- RunHarmony(micro_sub_oFC, reduction.save = "Harmony_", group.by.vars = c("sampleID", "sex"), plot_convergence = F)
micro_sub_dlpfc <- RunHarmony(micro_sub_dlpfc, reduction.save = "Harmony_", group.by.vars = c("orBatch", "sampleID", "version", "sex"), plot_convergence = F)

# =====================
# 3. UMAP, Neighbors, Clustering
# =====================
micro_sub_oFC <- RunUMAP(micro_sub_oFC, dims = 1:30, reduction = "Harmony_")
micro_sub_dlpfc <- RunUMAP(micro_sub_dlpfc, dims = 1:30, reduction = "Harmony_")

micro_sub_oFC <- FindNeighbors(micro_sub_oFC, reduction = "Harmony_", k.param = 30, nn.eps = 0)
micro_sub_dlpfc <- FindNeighbors(micro_sub_dlpfc, reduction = "Harmony_", k.param = 30, nn.eps = 0)

micro_sub_oFC <- FindClusters(micro_sub_oFC, algorithm = 1, resolution = c(0.1, 0.3))
micro_sub_dlpfc <- FindClusters(micro_sub_dlpfc, algorithm = 1, resolution = c(0.1, 0.3))

# =====================
# 4. Azimuth Annotation
# =====================
micro_sub_oFC <- RunAzimuth(micro_sub_oFC, reference = "/data/BBRF/Azimuth_annotation")
micro_sub_oFC$percent.mt <- mito
micro_sub_dlpfc <- RunAzimuth(micro_sub_dlpfc, reference = "/data/BBRF/Azimuth_annotation")

# =====================
# 5. Subset PVM Microglia & Merge
# =====================
micro_sub_oFC_ss <- subset(micro_sub_oFC, subset = predicted.subclass == "Micro-PVM")
micro_sub_dlpfc_ss <- subset(micro_sub_dlpfc, subset = predicted.subclass == "Micro-PVM")

DefaultAssay(micro_sub_oFC_ss) <- "RNA"
DefaultAssay(micro_sub_dlpfc_ss) <- "RNA"

micro_sub_merge <- merge(x = micro_sub_oFC_ss, y = micro_sub_dlpfc_ss)
micro_sub_merge <- DietSeurat(micro_sub_merge, layers = "counts")
micro_sub_merge@assays$SCT <- NULL

micro_sub_obj <- SplitObject(micro_sub_merge, split.by = "dataset")
features <- SelectIntegrationFeatures(object.list = micro_sub_obj, nfeatures = 3000)

micro_sub_merge <- NormalizeData(micro_sub_merge)
VariableFeatures(micro_sub_merge) <- features
micro_sub_merge <- ScaleData(micro_sub_merge, vars.to.regress = c("percent.mt", "nCount_RNA"))
micro_sub_merge <- RunPCA(micro_sub_merge, npcs = 50)
micro_sub_merge <- RunHarmony(micro_sub_merge, reduction.save = "Harmony_", group.by.vars = c("sampleID", "sex", "tissue"), plot_convergence = F)

micro_sub_merge <- RunUMAP(micro_sub_merge, dims = 1:10, reduction = "Harmony_")
micro_sub_merge <- FindNeighbors(micro_sub_merge, reduction = "Harmony_", k.param = 10, nn.eps = 0)
micro_sub_merge <- FindClusters(micro_sub_merge, algorithm = 1, resolution = c(0.1, 0.3))

micro_sub <- micro_sub_merge
micro_sub <- RunAzimuth(micro_sub, reference = "/data/BBRF/Azimuth_annotation")

# =====================
# 6. UMAP Plots (Categorical & Continuous)
# =====================
categoricals <- c("orBatch", "version", "sex", "gDx", "tissue", "RNA_snn_res.0.1", "RNA_snn_res.0.3", "predicted.subclass", "predicted.cluster")
reductions <- grep(names(micro_sub@reductions), pattern = "^umap", value = TRUE)

plots_list <- list()
for(reduction in reductions){
  for(categorical in categoricals){
    plots_list[[paste(reduction, categorical, sep = "_")]] <- DimPlot(micro_sub, reduction = reduction, group.by = categorical)
  }
}

pdf("/data/BBRF/Finalized_outputs/sketch_clusters/microglia/4_UMAPS_PCS_harmony_categorical_int2.pdf",
    onefile = TRUE, width = 12)
plots_list
dev.off()

# Continuous variables
continuous_vars <- c("nFeature_RNA", "nCount_RNA", "pmi", "age", "percent.mt")
plots_list <- list()
for(reduction in reductions){
  plots_list[[reduction]] <- FeaturePlot(micro_sub, ncol = 2, features = continuous_vars, reduction = reduction)
}

pdf("/data/BBRF/Finalized_outputs/sketch_clusters/microglia/4_UMAPS_PCS_harmony_continuous_int2.pdf",
    onefile = TRUE, width = 12, height = 16)
plots_list
dev.off()

# =====================
# 7. Marker Gene Analysis
# =====================
Idents(micro_sub) <- micro_sub$RNA_snn_res.0.1
mrks <- FindAllMarkers(micro_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dir.create("/data/BBRF/Finalized_outputs/sketch_clusters/microglia/MarkerGenes/", recursive = TRUE)
saveRDS(mrks, file = "/data/BBRF/Finalized_outputs/sketch_clusters/microglia/MarkerGenes/mrks_0.1_2.rds")

# =====================
# 8. UMAPs for Azimuth Predictions
# =====================
pdf("/data/BBRF/Finalized_outputs/sketch_clusters/microglia/UMAPs/azimuth_predictions.pdf",
    width = 15, height = 15)
DimPlot(micro_sub, group.by = "predicted.class", label = TRUE, label.size = 3)
DimPlot(micro_sub, group.by = "predicted.subclass", label = TRUE, label.size = 3)
DimPlot(micro_sub, group.by = "predicted.cluster", label = TRUE, label.size = 3)
dev.off()

# =====================
# 9. Cell Type-Specific Violin Plots
# =====================
homeo.specific <- c("P2RY12", "CX3CR1")
activated.specific <- c("AIF1", "CD68", "CD206", "CD45")
macro.specific <- c("F13A1", "LYVE1", "SIGLEC1")
tcell.specific <- c("BCL11B", "CD247")

pdf("/data/BBRF/Finalized_outputs/sketch_clusters/microglia/MarkerGenes/micro_sub_celltypes_1_2.pdf",
    width = 20)
VlnPlot(micro_sub, features = homeo.specific, pt.size = 0, split.by = "dataset")
VlnPlot(micro_sub, features = activated.specific, pt.size = 0, split.by = "dataset")
VlnPlot(micro_sub, features = macro.specific, pt.size = 0, split.by = "dataset")
VlnPlot(micro_sub, features = tcell.specific, pt.size = 0, split.by = "dataset")
VlnPlot(micro_sub, features = "SNAP25", pt.size = 0, split.by = "dataset")
dev.off()

# =====================
# 10. Save Final Seurat Object & Propeller Tests
# =====================
micro_sub_ss <- subset(micro_sub, idents = c(8), invert = TRUE)

micro_sub_ss <- SCTransform(micro_sub_ss, vars.to.regress = c("percent.mt", "nCount_RNA"))
micro_sub_ss <- RunPCA(micro_sub_ss, npcs = 50)
micro_sub_ss <- RunHarmony(micro_sub_ss, reduction.save = "Harmony_", group.by.vars = c("sampleID", "sex", "tissue"), plot_convergence = F)
micro_sub_ss <- RunUMAP(micro_sub_ss, dims = 1:10, reduction = "Harmony_")
micro_sub_ss <- FindNeighbors(micro_sub_ss, reduction = "Harmony_", k.param = 10, nn.eps = 0)
micro_sub_ss <- FindClusters(micro_sub_ss, algorithm = 1, resolution = c(0.1, 0.3))

# Annotate subclasses/classes
micro_sub_ss$subclasses <- factor(micro_sub_ss$SCT_snn_res.0.1, labels = c("Micro_homeo", "Micro_homeo", "Micro_homeo", "Micro_homeo", "Micro_act", "Macrophages", "Micro_homeo","Tcells"))
micro_sub_ss$classes <- factor(micro_sub_ss$SCT_snn_res.0.1, labels = c("Micro", "Micro", "Micro", "Micro", "Micro", "Macrophages", "Micro","Tcells"))

saveRDS(micro_sub_ss, file = "/data/BBRF/Seurat_obj/Subclusters/micro_tissue_final2.rds")

# Propeller tests (cell type proportions)
dir.create("/data/BBRF/Finalized_outputs/sketch_clusters/microglia/ctp/", recursive = TRUE)

groups <- list("orBatch", "tissue", "version", "sex", "gDx")
clusters <- list(micro_sub_ss$subclasses, micro_sub_ss$classes)
resolutions <- c(0.3, 0.1)
for(i in 1:2){
  clust <- clusters[[i]]
  res <- resolutions[i]
  for(g in groups){
    ctp <- propeller(x = micro_sub_ss, clusters = clust, sample = micro_sub_ss$sampleID, group = micro_sub_ss[[g]], trend = FALSE, robust = TRUE, transform = "logit")
    write.csv(ctp, paste0("/data/BBRF/Finalized_outputs/sketch_clusters/microglia/ctp/ctp_2_umap_", g, "_", res, ".csv"))
  }
}

# Save metadata
meta <- data.frame(
  Cell_ID = rownames(micro_sub_ss@meta.data),
  Classes = as.character(micro_sub_ss$classes),
  Subclasses = as.character(micro_sub_ss$subclasses),
  Cluster = micro_sub_ss$RNA_snn_res.0.1
)
write.csv(meta, "/data/BBRF/Pheno_Demo/combined_meta_clusters_micro2.csv", row.names = FALSE)

########################################
# End of Microglia Subcluster Processing
########################################
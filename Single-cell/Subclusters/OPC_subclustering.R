########################################
# OPC Subcluster Processing - oFC and dlPFC
# Author: Artemis Iatrou
# Date: 2025-08-08
# Description:
#   Subsets OPCs by tissue, performs SCTransform, PCA, Harmony integration,
#   UMAP, clustering, Azimuth annotation, marker gene analysis, and saves outputs.
########################################

# =====================
# 0. Load Packages & Set Options
# =====================
library(Seurat)
library(future)
options(future.globals.maxSize = 1000000 * 1024^2)

# =====================
# 1. Load Seurat Object
# =====================
opc_sub <- readRDS("/data/BBRF/Seurat_obj/Subclusters/opc_sub.RDS")
opc_sub <- DietSeurat(opc_sub)

# =====================
# 2. Subset by Tissue
# =====================
opc_sub_oFC <- subset(opc_sub, subset = tissue == "oFC")
opc_sub_dlpfc <- subset(opc_sub, subset = tissue == "dlPFC")

# Preserve percent.mt for oFC
mito <- opc_sub_oFC$percent.mt

# =====================
# 3. Annotate Metadata
# =====================
opc_sub_oFC$dataset <- "eb"
opc_sub_dlpfc$dataset <- "ai"

opc_sub_oFC$sampleID <- opc_sub_oFC$Donor
opc_sub_dlpfc$sampleID <- opc_sub_dlpfc$sample

opc_sub_oFC$tissue <- "oFC"
opc_sub_dlpfc$tissue <- "dlPFC"

opc_sub_oFC$pr_clusters <- opc_sub_oFC$celltypes_final
opc_sub_oFC$pmi <- opc_sub_oFC$PMI
opc_sub_oFC$orBatch <- "Binder_oFC"
opc_sub_oFC$version <- "v3"

opc_sub_oFC$sex <- factor(opc_sub_oFC$Sex, levels = c("F", "M"), labels = c("female", "male"))
opc_sub_oFC$age <- opc_sub_oFC$Age

opc_sub_oFC$gDx <- factor(opc_sub_oFC$Classification, levels = c("Control", "MDD"))
opc_sub_dlpfc$gDx <- factor(opc_sub_dlpfc$primarydx, levels = c("Control", "MDD", "MDD"))

# =====================
# 4. SCTransform & PCA
# =====================
opc_sub_oFC <- SCTransform(opc_sub_oFC, vars.to.regress = c("percent.mt", "nCount_RNA"))
opc_sub_dlpfc <- SCTransform(opc_sub_dlpfc, vars.to.regress = c("percent.mt", "nCount_RNA"))

opc_sub_oFC <- RunPCA(opc_sub_oFC, npcs = 50)
opc_sub_dlpfc <- RunPCA(opc_sub_dlpfc, npcs = 50)

# =====================
# 5. Azimuth Annotation
# =====================
opc_sub_oFC <- RunAzimuth(opc_sub_oFC, reference = "/data/BBRF/Azimuth_annotation")
opc_sub_dlpfc <- RunAzimuth(opc_sub_dlpfc, reference = "/data/BBRF/Azimuth_annotation")

# Restore percent.mt
opc_sub_oFC$percent.mt <- mito

# =====================
# 6. Subset to OPC predicted subclass
# =====================
opc_sub_oFC_ss <- subset(opc_sub_oFC, subset = predicted.subclass == "OPC")
opc_sub_dlpfc_ss <- subset(opc_sub_dlpfc, subset = predicted.subclass == "OPC")

DefaultAssay(opc_sub_oFC_ss) <- "RNA"
DefaultAssay(opc_sub_dlpfc_ss) <- "RNA"

# =====================
# 7. Merge Datasets & Preprocess
# =====================
opc_sub_comb <- merge(opc_sub_oFC_ss, y = opc_sub_dlpfc_ss)
opc_sub_comb <- DietSeurat(opc_sub_comb, layers = "counts")
opc_sub_comb@assays$SCT <- NULL

opc_sub_list <- SplitObject(opc_sub_comb, split.by = "dataset")
features <- SelectIntegrationFeatures(opc_sub_list, nfeatures = 3000)

opc_sub_comb <- NormalizeData(opc_sub_comb)
VariableFeatures(opc_sub_comb) <- features
opc_sub_comb <- ScaleData(opc_sub_comb, vars.to.regress = c("percent.mt", "nCount_RNA"))
opc_sub_comb <- RunPCA(opc_sub_comb, npcs = 50)
opc_sub_comb <- RunHarmony(opc_sub_comb, reduction.save= "Harmony_", group.by.vars = c("sampleID", "sex", "tissue"))

opc_sub_comb <- RunUMAP(opc_sub_comb, dims = 1:10, reduction = "Harmony_")
opc_sub_comb <- FindNeighbors(opc_sub_comb, reduction = "Harmony_", k.param = 10, nn.eps = 0)
opc_sub_comb <- FindClusters(opc_sub_comb, algorithm = 1, resolution = c(0.1, 0.3))

opc_sub <- opc_sub_comb

# =====================
# 8. Visualizations
# =====================
plots_dir <- "/data/BBRF/Finalized_outputs/sketch_clusters/opc/"
dir.create(plots_dir, recursive = TRUE)

# Categorical variables
categoricals <- c("orBatch", "version", "sex", "gDx", "tissue",
                  "RNA_snn_res.0.1", "RNA_snn_res.0.3",
                  "predicted.subclass", "predicted.cluster")

reductions <- grep(names(opc_sub@reductions), pattern = "^umap", value = TRUE)

plots_list <- list()
for(reduction in reductions){
  for(cat_var in categoricals){
    plots_list[[paste(reduction, cat_var, sep="_")]] <- DimPlot(opc_sub, reduction=reduction, group.by=cat_var)
  }
}
pdf(file = paste0(plots_dir, "UMAPS_categorical.pdf"), onefile = TRUE, width = 12)
plots_list
dev.off()

# Continuous variables
continuous_vars <- c("nFeature_RNA", "nCount_RNA", "pmi", "age", "percent.mt")
plots_list <- list()
for(reduction in reductions){
  plots_list[[reduction]] <- FeaturePlot(opc_sub, features = continuous_vars, reduction = reduction, ncol = 2)
}
pdf(file = paste0(plots_dir, "UMAPS_continuous.pdf"), onefile = TRUE, width = 12, height = 16)
plots_list
dev.off()

# =====================
# 9. Marker Gene Analysis
# =====================
Idents(opc_sub) <- opc_sub$RNA_snn_res.0.1
mrks <- FindAllMarkers(opc_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

marker_dir <- paste0(plots_dir, "MarkerGenes/")
dir.create(marker_dir, recursive = TRUE)
saveRDS(mrks, file = paste0(marker_dir, "mrks_0.1_2.rds"))

# Example violin plots
other_genes <- c("AQP4", "MBP", "PLP1")
pdf(paste0(marker_dir, "opc_sub_celltypes_1_2.pdf"))
VlnPlot(opc_sub, features = other_genes, pt.size = 0, split.by = "dataset")
dev.off()

# =====================
# 10. Subset & Final Processing
# =====================
opc_sub_ss <- subset(opc_sub, idents = c(4:6), invert = TRUE)
opc_sub_ss <- SCTransform(opc_sub_ss, vars.to.regress = c("percent.mt", "nCount_RNA"))
opc_sub_ss <- RunPCA(opc_sub_ss, npcs = 50)
opc_sub_ss <- RunHarmony(opc_sub_ss, reduction.save= "Harmony_", group.by.vars = c("sampleID", "sex", "tissue"))
opc_sub_ss <- RunUMAP(opc_sub_ss, dims = 1:10, reduction = "Harmony_")
opc_sub_ss <- FindNeighbors(opc_sub_ss, reduction = "Harmony_", k.param = 10, nn.eps = 0)
opc_sub_ss <- FindClusters(opc_sub_ss, algorithm = 1, resolution = c(0.1, 0.3))

# Final UMAP plots
pdf(paste0(plots_dir, "UMAPS_final_ss2.pdf"), width = 8, height = 8)
DimPlot(opc_sub_ss, group.by = "SCT_snn_res.0.1", label = TRUE, label.size = 3)
DimPlot(opc_sub_ss, group.by = "sex", label = TRUE, label.size = 3)
DimPlot(opc_sub_ss, group.by = "tissue", label = TRUE, label.size = 3)
dev.off()

# Marker genes for final subset
Idents(opc_sub_ss) <- opc_sub_ss$SCT_snn_res.0.1
mrks_ss <- FindAllMarkers(opc_sub_ss, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(mrks_ss, file = paste0(marker_dir, "mrks_ss_0.1_2.rds"))

# Annotate classes/subclasses
opc_sub_ss$subclasses <- "OPC"
opc_sub_ss$classes <- "OPC"

pdf(paste0(plots_dir, "UMAPS_final_classes_subclasses.pdf"), width = 8, height = 8)
DimPlot(opc_sub_ss, group.by = "classes", label = TRUE, reduction = "umap", label.box = TRUE)
DimPlot(opc_sub_ss, group.by = "subclasses", label = TRUE, reduction = "umap", label.box = TRUE)
dev.off()

# Save final object
saveRDS(opc_sub_ss, file = "/data/BBRF/Seurat_obj/Subclusters/opc_tissue_final2.rds")

# =====================
# 11. Export Metadata
# =====================
meta <- data.frame(
  Cell_ID = rownames(opc_sub_ss@meta.data),
  Classes = as.character(opc_sub_ss$classes),
  Subclasses = as.character(opc_sub_ss$subclasses),
  Cluster = opc_sub_ss$SCT_snn_res.0.1
)
write.csv(meta, "/data/BBRF/Pheno_Demo/combined_meta_clusters_opc2.csv")
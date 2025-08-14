########################################
# Oligodendrocyte (Oligo) Subcluster Processing - oFC and dlPFC
# Author: Artemis Iatrou
# Date: 2025-08-08
# Description:
#   Subsets Oligos by tissue, performs SCTransform, PCA, Harmony integration,
#   UMAP, clustering, Azimuth annotation, marker gene analysis, propeller tests,
#   and saves outputs.
########################################

# =====================
# 0. Load Packages & Set Options
# =====================
library(Seurat)
library(future)
library(speckle)

options(future.globals.maxSize = 1000000 * 1024^2)

# =====================
# 1. Load Subcluster Object & Create Tissue Subsets
# =====================
oligo_sub <- readRDS("/data/BBRF/Seurat_obj/Subclusters/oligo_sub.RDS")
oligo_sub <- DietSeurat(oligo_sub)

oligo_sub_oFC <- subset(oligo_sub, subset = tissue == "oFC")
mito <- oligo_sub_oFC$percent.mt
oligo_sub_dlpfc <- subset(oligo_sub, subset = tissue == "dlPFC")

# =====================
# 2. Assign Metadata
# =====================
# Dataset
oligo_sub_oFC$dataset <- "eb"
oligo_sub_dlpfc$dataset <- "ai"

# Sample ID
oligo_sub_oFC$sampleID <- oligo_sub_oFC$Donor
oligo_sub_dlpfc$sampleID <- oligo_sub_dlpfc$sample

# Tissue
oligo_sub_oFC$tissue <- "oFC"
oligo_sub_dlpfc$tissue <- "dlPFC"

# Original Clusters & PMI
oligo_sub_oFC$pr_clusters <- oligo_sub_oFC$celltypes_final
oligo_sub_oFC$pmi <- oligo_sub_oFC$PMI

# Batch & Version
oligo_sub_oFC$orBatch <- "Binder_oFC"
oligo_sub_oFC$version <- "v3"

# Sex & Age
oligo_sub_oFC$sex <- factor(oligo_sub_oFC$Sex, levels = c("F", "M"), labels = c("female", "male"))
oligo_sub_oFC$age <- oligo_sub_oFC$Age

# Diagnosis
oligo_sub_oFC$gDx <- factor(oligo_sub_oFC$Classification, levels = c("Control", "MDD"))
oligo_sub_dlpfc$gDx <- factor(oligo_sub_dlpfc$primarydx, levels = c("Control", "MDD"))

# =====================
# 3. SCTransform & PCA
# =====================
oligo_sub_oFC <- SCTransform(oligo_sub_oFC, vars.to.regress = c("percent.mt", "nCount_RNA"))
oligo_sub_dlpfc <- SCTransform(oligo_sub_dlpfc, vars.to.regress = c("percent.mt", "nCount_RNA"))

oligo_sub_oFC <- RunPCA(oligo_sub_oFC, npcs = 50)
oligo_sub_dlpfc <- RunPCA(oligo_sub_dlpfc, npcs = 50)

# =====================
# 4. Azimuth Annotation
# =====================
oligo_sub_oFC <- RunAzimuth(oligo_sub_oFC, reference = "/data/BBRF/Azimuth_annotation")
oligo_sub_dlpfc <- RunAzimuth(oligo_sub_dlpfc, reference = "/data/BBRF/Azimuth_annotation")

# Restore original mito percent
oligo_sub_oFC$percent.mt <- mito

# =====================
# 5. Subset Predicted Oligos
# =====================
oligo_sub_oFC_ss <- subset(oligo_sub_oFC, subset = predicted.subclass == "Oligo")
oligo_sub_dlpfc_ss <- subset(oligo_sub_dlpfc, subset = predicted.subclass == "Oligo")

DefaultAssay(oligo_sub_oFC_ss) <- "RNA"
DefaultAssay(oligo_sub_dlpfc_ss) <- "RNA"

# =====================
# 6. Merge & Integration Prep
# =====================
oligo_sub_merged <- merge(x = oligo_sub_oFC_ss, y = oligo_sub_dlpfc_ss)

oligo_sub_list <- SplitObject(oligo_sub_merged, split.by = "dataset", assay = "SCT")
features <- SelectIntegrationFeatures(object.list = oligo_sub_list, nfeatures = 3000, assay = c("SCT", "SCT"))

oligo_sub_merged <- DietSeurat(oligo_sub_merged, layers = "counts")
oligo_sub_merged@assays$SCT <- NULL

# =====================
# 7. Normalize, Scale, PCA, Harmony, UMAP, Clustering
# =====================
oligo_sub_merged <- NormalizeData(oligo_sub_merged)
VariableFeatures(oligo_sub_merged) <- features
oligo_sub_merged <- ScaleData(oligo_sub_merged, vars.to.regress = c("percent.mt", "nCount_RNA"))
oligo_sub_merged <- RunPCA(oligo_sub_merged, npcs = 50)

oligo_sub_merged <- RunHarmony(oligo_sub_merged, reduction.save = "Harmony_", group.by.vars = c("sampleID", "sex", "tissue"), plot_convergence = FALSE)
oligo_sub_merged <- RunUMAP(oligo_sub_merged, dims = 1:50, reduction = "Harmony_")
oligo_sub_merged <- FindNeighbors(oligo_sub_merged, reduction = "Harmony_", k.param = 50, nn.eps = 0)
oligo_sub_merged <- FindClusters(oligo_sub_merged, algorithm = 1, resolution = c(0.1, 0.3))

oligo_sub_merged <- RunAzimuth(oligo_sub_merged, reference = "/data/BBRF/Azimuth_annotation")

oligo_sub <- oligo_sub_merged

# =====================
# 8. Marker Gene Analysis
# =====================
Idents(oligo_sub) <- oligo_sub$RNA_snn_res.0.1
mrks <- FindAllMarkers(oligo_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

dir.create("/data/BBRF/Finalized_outputs/sketch_clusters/oligo/MarkerGenes/", recursive = TRUE)
saveRDS(mrks, file = "/data/BBRF/Finalized_outputs/sketch_clusters/oligo/MarkerGenes/mrks_0.1_2.rds")

# =====================
# 9. Optional Subset for Cleaned Oligos (remove cluster 3)
# =====================
oligo_sub_ss <- subset(oligo_sub, idents = 3, invert = TRUE)

# Re-run integration pipeline on subset
oligo_sub_list_ss <- SplitObject(oligo_sub_ss, split.by = "dataset")
features_ss <- SelectIntegrationFeatures(object.list = oligo_sub_list_ss, nfeatures = 3000)

oligo_sub_ss <- NormalizeData(oligo_sub_ss)
VariableFeatures(oligo_sub_ss) <- features_ss
oligo_sub_ss <- ScaleData(oligo_sub_ss, vars.to.regress = c("percent.mt", "nCount_RNA"))

oligo_sub_ss <- RunHarmony(oligo_sub_ss, reduction.save = "Harmony_", group.by.vars = c("sampleID", "orBatch", "sex", "tissue"), plot_convergence = FALSE)
oligo_sub_ss <- RunUMAP(oligo_sub_ss, dims = 1:10, reduction = "Harmony_")
oligo_sub_ss <- FindNeighbors(oligo_sub_ss, reduction = "Harmony_", k.param = 10, nn.eps = 0)
oligo_sub_ss <- FindClusters(oligo_sub_ss, graph.name = "RNA_snn", algorithm = 1, resolution = c(0.1, 0.3))

oligo_sub <- oligo_sub_ss

# =====================
# 10. Plotting UMAPs (Categorical, Sample, Continuous)
# =====================
reductions <- grep(names(oligo_sub@reductions), pattern = "^umap", value = TRUE)

# Categoricals
categoricals <- c("orBatch", "version", "sex", "gDx", "tissue", "RNA_snn_res.0.1", "RNA_snn_res.0.3", "predicted.subclass", "predicted.cluster")
plots_list <- list()
for (reduction in reductions) {
  for (cat in categoricals) {
    plots_list[[paste(reduction, cat, sep = "_")]] <- DimPlot(oligo_sub, reduction = reduction, group.by = cat)
  }
}
pdf("/data/BBRF/Finalized_outputs/sketch_clusters/oligo/4_UMAPS_PCS_harmony_categorical_int2.pdf", onefile = TRUE, width = 12)
plots_list
dev.off()

# Sample-only plots
plots_list <- list()
categoricals <- c("sampleID")
for (reduction in reductions) {
  for (cat in categoricals) {
    plots_list[[paste(reduction, cat, sep = "_")]] <- DimPlot(oligo_sub, reduction = reduction, group.by = cat) + NoLegend()
  }
}
pdf("/data/BBRF/Finalized_outputs/sketch_clusters/oligo/4_UMAPS_PCS_harmony_samples_int2.pdf", onefile = TRUE, width = 12)
plots_list
dev.off()

# Continuous features
continuous_vars <- c("nFeature_RNA", "nCount_RNA", "pmi", "age", "percent.mt")
plots_list <- list()
for (reduction in reductions) {
  plots_list[[reduction]] <- FeaturePlot(oligo_sub, features = continuous_vars, reduction = reduction, ncol = 2)
}
pdf("/data/BBRF/Finalized_outputs/sketch_clusters/oligo/4_UMAPS_PCS_harmony_continuous_int2.pdf", onefile = TRUE, width = 12, height = 16)
plots_list
dev.off()

# =====================
# 11. Subclass/Classes Assignment & VlnPlot
# =====================
oligo_sub$subclasses <- Idents(oligo_sub)
levels(oligo_sub$subclasses) <- c("Oligo_1", "Oligo_2", "Oligo_3")

oligo_sub$classes <- "Oligo"


# =====================
# 12. Violin Plots for Key Marker Genes
# =====================
marker_genes <- c("MOBP", "MBP", "PLP1", "ST18", "SNAP25")

dir.create("/data/BBRF/Finalized_outputs/sketch_clusters/oligo/MarkerGenes/", recursive = TRUE)
pdf("/data/BBRF/Finalized_outputs/sketch_clusters/oligo/MarkerGenes/oligo_sub_celltypes_1_2.pdf")
VlnPlot(oligo_sub, features = marker_genes, pt.size = 0, split.by = "dataset")
dev.off()

# =====================
# 13. UMAP Plots for Classes/Subclasses
# =====================
dir.create("/data/BBRF/Finalized_outputs/sketch_clusters/oligo/UMAPs/", recursive = TRUE)
pdf("/data/BBRF/Finalized_outputs/sketch_clusters/oligo/UMAPs/final_ss2.pdf", width = 8, height = 8)
DimPlot(oligo_sub, group.by = "classes", label = TRUE, reduction = "umap", label.box = TRUE)
DimPlot(oligo_sub, group.by = "subclasses", label = TRUE, reduction = "umap", label.box = TRUE)
dev.off()

# =====================
# 14. Propeller Tests (Cell Type Proportions)
# =====================
dir.create("/data/BBRF/Finalized_outputs/sketch_clusters/oligo/ctp/", recursive = TRUE)

ctp_1 <- propeller(x = oligo_sub, clusters = oligo_sub$subclasses, sample = oligo_sub$sampleID, group = oligo_sub$orBatch, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_1, "/data/BBRF/Finalized_outputs/sketch_clusters/oligo/ctp/ctp_2_umap_batch_0.3.csv")

ctp_1.1 <- propeller(x = oligo_sub, clusters = oligo_sub$subclasses, sample = oligo_sub$sampleID, group = oligo_sub$tissue, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_1.1, "/data/BBRF/Finalized_outputs/sketch_clusters/oligo/ctp/ctp_2_umap_tissue_0.3.csv")

ctp_1.2 <- propeller(x = oligo_sub, clusters = oligo_sub$subclasses, sample = oligo_sub$sampleID, group = oligo_sub$version, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_1.2, "/data/BBRF/Finalized_outputs/sketch_clusters/oligo/ctp/ctp_2_umap_version_0.3.csv")

ctp_1.3 <- propeller(x = oligo_sub, clusters = oligo_sub$subclasses, sample = oligo_sub$sampleID, group = oligo_sub$sex, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_1.3, "/data/BBRF/Finalized_outputs/sketch_clusters/oligo/ctp/ctp_2_umap_sex_0.3.csv")

ctp_1.4 <- propeller(x = oligo_sub, clusters = oligo_sub$subclasses, sample = oligo_sub$sampleID, group = oligo_sub$gDx, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_1.4, "/data/BBRF/Finalized_outputs/sketch_clusters/oligo/ctp/ctp_2_umap_dx_0.3.csv")

# =====================
# 15. Save Final Seurat Object & Metadata
# =====================
saveRDS(oligo_sub, file = "/data/BBRF/Seurat_obj/Subclusters/oligo_tissue_final2.rds")

meta <- data.frame(
  Cell_ID = rownames(oligo_sub@meta.data),
  Classes = as.character(oligo_sub$classes),
  Subclasses = as.character(oligo_sub$subclasses),
  Cluster = oligo_sub$RNA_snn_res.0.1
)
write.csv(meta, "/data/BBRF/Pheno_Demo/combined_meta_clusters_oligo2.csv", row.names = FALSE)

########################################
# End of Oligo Subcluster Processing
########################################
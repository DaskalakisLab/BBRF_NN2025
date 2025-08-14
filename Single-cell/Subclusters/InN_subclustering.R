########################################
# Inhibitory Neurons Subcluster Processing - oFC and dlPFC
# Author: Artemis Iatrou
# Date: 2025-08-14
# Description:
#   Subsets inhibitory neurons by tissue, performs SCTransform, PCA, Harmony integration,
#   UMAP, clustering, Azimuth annotation, marker gene analysis, categorical & continuous plotting,
#   and propeller cluster proportion testing (CTP), then saves outputs.
########################################

# =====================
# 0. Load Packages & Set Options
# =====================
library(Seurat)
library(future)
library(speckle)
options(future.globals.maxSize = 1000000 * 1024^2)

# =====================
# 1. Load Data & Initial Subsetting
# =====================
inh_sub <- readRDS("/data/BBRF/Seurat_obj/Subclusters/inh_sub2.RDS")
inh_sub <- DietSeurat(inh_sub)

inh_sub_oFC <- subset(inh_sub, subset = tissue == "oFC")
inh_sub_dlpfc <- subset(inh_sub, subset = tissue == "dlPFC")
inh_sub_oFC$dataset <- "eb"
inh_sub_dlpfc$dataset <- "ai"

inh_sub_oFC$sampleID <- inh_sub_oFC$Donor
inh_sub_dlpfc$sampleID <- inh_sub_dlpfc$sample
inh_sub_oFC$tissue <- "oFC"
inh_sub_dlpfc$tissue <- "dlPFC"

inh_sub_oFC$pr_clusters <- inh_sub_oFC$celltypes_final
inh_sub_oFC$pmi <- inh_sub_oFC$PMI
inh_sub_oFC$orBatch <- "Binder_oFC"
inh_sub_oFC$version <- "v3"

inh_sub_oFC$sex <- factor(inh_sub_oFC$Sex, levels = c("F", "M"), labels = c("female", "male"))
inh_sub_oFC$age <- inh_sub_oFC$Age
inh_sub_oFC$gDx <- factor(inh_sub_oFC$Classification, levels = c("Control", "MDD"))

inh_sub_dlpfc$gDx <- factor(inh_sub_dlpfc$primarydx, levels = c("Control", "MDD"))

# =====================
# 2. Azimuth Annotation
# =====================
mito <- inh_sub_oFC$percent.mt
inh_sub_oFC <- RunAzimuth(inh_sub_oFC, reference = "/data/BBRF/Azimuth_annotation")
inh_sub_dlpfc <- RunAzimuth(inh_sub_dlpfc, reference = "/data/BBRF/Azimuth_annotation")
inh_sub_oFC$percent.mt <- mito

# =====================
# 3. Subset Inhibitory Classes & SCTransform
# =====================
inh_sub_oFC_ss <- subset(inh_sub_oFC, subset = predicted.subclass %in% c("Lamp5","Pvalb","Sncg","Sst","Sst Chodl","Vip"))
inh_sub_dlpfc_ss <- subset(inh_sub_dlpfc, subset = predicted.subclass %in% c("Lamp5","Pvalb","Sncg","Sst","Sst Chodl","Vip"))

inh_sub_oFC_ss <- SCTransform(inh_sub_oFC_ss, vars.to.regress = c("percent.mt","nCount_RNA"))
inh_sub_dlpfc_ss <- SCTransform(inh_sub_dlpfc_ss, vars.to.regress = c("percent.mt","nCount_RNA"))

DefaultAssay(inh_sub_oFC_ss) <- "RNA"
DefaultAssay(inh_sub_dlpfc_ss) <- "RNA"

# =====================
# 4. Merge & Integration Prep
# =====================
inh_sub. <- merge(x = inh_sub_oFC_ss, y = inh_sub_dlpfc_ss)
inh_sub.obj <- SplitObject(inh_sub., split.by = "dataset")
features <- SelectIntegrationFeatures(object.list = inh_sub.obj, nfeatures = 3000)

inh_sub. <- DietSeurat(inh_sub., layers = "counts")
inh_sub.@assays$SCT <- NULL

inh_sub. <- NormalizeData(inh_sub.)
VariableFeatures(inh_sub.) <- features
inh_sub. <- ScaleData(inh_sub., vars.to.regress = c("percent.mt","nCount_RNA"))

# =====================
# 5. PCA & Harmony
# =====================
inh_sub. <- RunPCA(inh_sub., npcs = 50)
inh_sub. <- RunHarmony(inh_sub., reduction.save = "Harmony_", group.by.vars = c("sampleID","sex","tissue","orBatch"), plot_convergence = FALSE)

# =====================
# 6. UMAP, Neighbors, Clustering
# =====================
inh_sub. <- RunUMAP(inh_sub., dims = 1:50, reduction = "Harmony_")
inh_sub. <- FindNeighbors(inh_sub., reduction = "Harmony_", k.param = 50, nn.eps = 0)
inh_sub. <- FindClusters(inh_sub., algorithm = 1, resolution = c(0.1, 0.3))
inh_sub. <- RunAzimuth(inh_sub., reference = "/data/BBRF/Azimuth_annotation")

# =====================
# 7. Final GABAergic Subset & Re-Processing
# =====================
inh_sub. <- subset(inh_sub., subset = predicted.class == "GABAergic")
inh_sub. <- SCTransform(inh_sub., vars.to.regress = c("percent.mt","nCount_RNA"))
inh_sub. <- RunPCA(inh_sub., npcs = 50)
inh_sub. <- RunHarmony(inh_sub., reduction.save = "Harmony_", group.by.vars = c("sampleID","sex","tissue","orBatch"))
inh_sub. <- RunUMAP(inh_sub., dims = 1:50, reduction = "Harmony_")
inh_sub. <- FindNeighbors(inh_sub., reduction = "Harmony_", k.param = 50)
inh_sub. <- FindClusters(inh_sub., algorithm = 1, resolution = c(0.1, 0.3))
inh_sub. <- RunAzimuth(inh_sub., reference = "/data/BBRF/Azimuth_annotation")

# =====================
# 8. Assign Subclasses & Classes
# =====================
inh_sub <- inh_sub.
inh_sub$subclasses <- Idents(inh_sub)
levels(inh_sub$subclasses) <- c("Inh_VIP","Inh_PVALB_TAC1","Inh_SST_NPY","Inh_SNCG_RELN","Inh_LAMP5_KIT","Inh_PVALB_Cha","Inh_LAMP5_CHST9")
inh_sub$classes <- "Inh"

# =====================
# 9. Plotting
# =====================
# Categorical
plots_list <- list()
reductions <- grep(names(inh_sub@reductions), pattern = "^umap", value = TRUE)
categoricals <- c("orBatch","version","sex","gDx","tissue","SCT_snn_res.0.1","predicted.subclass","predicted.cluster")
for(reduction in reductions){
  for(cat in categoricals){
    plots_list[[paste(reduction,cat,sep="_")]] <- DimPlot(inh_sub, reduction = reduction, group.by = cat)
  }
}
pdf("/data/BBRF/Finalized_outputs/sketch_clusters/inh/4_UMAPS_PCS_harmony_categorical_int2.pdf", onefile = TRUE, width = 12)
plots_list
dev.off()

# Continuous
plots_list <- list()
continuous_vars <- c("nFeature_RNA","nCount_RNA","pmi","age","percent.mt")
for(reduction in reductions){
  plots_list[[reduction]] <- FeaturePlot(inh_sub, features = continuous_vars, reduction = reduction, ncol = 2)
}
pdf("/data/BBRF/Finalized_outputs/sketch_clusters/inh/4_UMAPS_PCS_harmony_continuous_int2.pdf", onefile = TRUE, width = 12, height = 16)
plots_list
dev.off()

# =====================
# 10. Marker Gene Analysis
# =====================
Idents(inh_sub) <- inh_sub$SCT_snn_res.0.1
mrks <- FindAllMarkers(inh_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dir.create("/data/BBRF/Finalized_outputs/sketch_clusters/inh/MarkerGenes/", recursive = TRUE)
saveRDS(mrks, "/data/BBRF/Finalized_outputs/sketch_clusters/inh/MarkerGenes/mrks_0.1_2.rds")

other_specific <- c("SST","PVALB","VIP","LAMP5","RELN","NPY","SNCG","CHODL","ADARB2","GRIK1","ZNF385D")
pdf("/data/BBRF/Finalized_outputs/sketch_clusters/inh/MarkerGenes/inh_sub_celltypes_1_2.pdf", width = 20)
VlnPlot(inh_sub, features = other_specific, pt.size = 0, split.by = "dataset")
dev.off()

# =====================
# 11. CTP Analysis
# =====================
dir.create("/data/BBRF/Finalized_outputs/sketch_clusters/inh/ctp/", recursive = TRUE)

ctp_batch <- propeller(x = inh_sub, clusters = inh_sub$subclasses, sample = inh_sub$sampleID, group = inh_sub$orBatch, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_batch, "/data/BBRF/Finalized_outputs/sketch_clusters/inh/ctp/ctp_subclasses_batch.csv")

ctp_tissue <- propeller(x = inh_sub, clusters = inh_sub$subclasses, sample = inh_sub$sampleID, group = inh_sub$tissue, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_tissue, "/data/BBRF/Finalized_outputs/sketch_clusters/inh/ctp/ctp_subclasses_tissue.csv")

ctp_version <- propeller(x = inh_sub, clusters = inh_sub$subclasses, sample = inh_sub$sampleID, group = inh_sub$version, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_version, "/data/BBRF/Finalized_outputs/sketch_clusters/inh/ctp/ctp_subclasses_version.csv")

ctp_sex <- propeller(x = inh_sub, clusters = inh_sub$subclasses, sample = inh_sub$sampleID, group = inh_sub$sex, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_sex, "/data/BBRF/Finalized_outputs/sketch_clusters/inh/ctp/ctp_subclasses_sex.csv")

ctp_dx <- propeller(x = inh_sub, clusters = inh_sub$subclasses, sample = inh_sub$sampleID, group = inh_sub$gDx, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_dx, "/data/BBRF/Finalized_outputs/sketch_clusters/inh/ctp/ctp_subclasses_dx.csv")

# =====================
# 12. Save Metadata
# =====================
meta <- data.frame(
  Cell_ID = rownames(inh_sub@meta.data),
  Classes = as.character(inh_sub$classes),
  Subclasses = as.character(inh_sub$subclasses),
  Cluster = inh_sub$SCT_snn_res.0.1
)
write.csv(meta, "/data/BBRF/Pheno_Demo/combined_meta_clusters_inh2.csv", row.names = FALSE)

# =====================
# 13. Save Final Seurat Object
# =====================
saveRDS(inh_sub, file = "/data/BBRF/Seurat_obj/Subclusters/inh_tissue_final2.rds")
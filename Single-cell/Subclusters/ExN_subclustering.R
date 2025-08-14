########################################
# Excitatory Neurons Subcluster Processing - oFC and dlPFC
# Author: Artemis Iatrou
# Date: 2025-08-14
# Description:
#   Subsets ex neurons by tissue, performs normalization, PCA, Harmony integration,
#   UMAP, clustering, Azimuth annotation, marker gene analysis, and saves outputs.
########################################

# =====================
# 0. Load Packages & Set Options
# =====================
library(Seurat)
library(future)
library(harmony)
library(Azimuth)
library(speckle)
options(future.globals.maxSize = 1000000 * 1024^2)
plan("multicore", workers = 96)

# =====================
# 1. Load Data & Subset by Tissue
# =====================
ex.sub <- readRDS("/data/Seurat_obj/Subclusters/ex.sub2.RDS")
ex.sub <- DietSeurat(ex.sub)

ex.sub_oFC <- subset(ex.sub, subset = tissue == "oFC")
ex.sub_dlpfc <- subset(ex.sub, subset = tissue == "dlPFC")

# =====================
# 2. Add Metadata & Factor Levels
# =====================
# Dataset identifiers
ex.sub_oFC$dataset <- "eb"
ex.sub_dlpfc$dataset <- "ai"

# Sample IDs
ex.sub_oFC$sampleID <- ex.sub_oFC$Donor
ex.sub_dlpfc$sampleID <- ex.sub_dlpfc$sample

# Tissue labels
ex.sub_oFC$tissue <- "oFC"
ex.sub_dlpfc$tissue <- "dlPFC"

# Other metadata
ex.sub_oFC$pr_clusters <- ex.sub_oFC$celltypes_final
ex.sub_oFC$pmi <- ex.sub_oFC$PMI
ex.sub_oFC$orBatch <- "Binder_oFC"
ex.sub_oFC$version <- "v3"

# Sex and diagnosis as factors
ex.sub_oFC$sex <- factor(ex.sub_oFC$Sex, levels = c("F", "M"), labels = c("female", "male"))
ex.sub_oFC$age <- ex.sub_oFC$Age
ex.sub_oFC$gDx <- factor(ex.sub_oFC$Classification, levels = c("Control", "MDD"))

ex.sub_dlpfc$gDx <- factor(ex.sub_dlpfc$primarydx, levels = c("Control", "MDD", "MDD"))

# =====================
# 3. Normalize, Find Variable Features, and Scale
# =====================
ex.sub_oFC <- NormalizeData(ex.sub_oFC)
ex.sub_oFC <- FindVariableFeatures(ex.sub_oFC)
ex.sub_oFC <- ScaleData(ex.sub_oFC, vars.to.regress = c("percent.mt", "nCount_RNA"))

ex.sub_dlpfc@assays$sketch <- NULL
ex.sub_dlpfc <- NormalizeData(ex.sub_dlpfc)
ex.sub_dlpfc <- FindVariableFeatures(ex.sub_dlpfc)
ex.sub_dlpfc <- ScaleData(ex.sub_dlpfc, vars.to.regress = c("percent.mt", "nCount_RNA"))

# =====================
# 4. PCA
# =====================
ex.sub_oFC <- RunPCA(ex.sub_oFC, npcs = 50)
ex.sub_dlpfc <- RunPCA(ex.sub_dlpfc, npcs = 50)

# =====================
# 5. Harmony Integration
# =====================
dir.create("/data/Finalized_outputs/sketch_clusters/ex/", recursive = TRUE)
set.seed(1234)

ex.sub_oFC <- RunHarmony(ex.sub_oFC, reduction.save = "Harmony_", 
                          group.by.vars = c("sampleID", "sex"), plot_convergence = FALSE)

ex.sub_dlpfc <- RunHarmony(ex.sub_dlpfc, reduction.save = "Harmony_", 
                            group.by.vars = c("orBatch", "sampleID", "version", "sex"), plot_convergence = FALSE)

# =====================
# 6. UMAP, Neighbors, and Clustering
# =====================
ex.sub_oFC <- RunUMAP(ex.sub_oFC, dims = 1:50, reduction = "Harmony_")
ex.sub_dlpfc <- RunUMAP(ex.sub_dlpfc, dims = 1:50, reduction = "Harmony_")

ex.sub_oFC <- FindNeighbors(ex.sub_oFC, reduction = "Harmony_", k.param = 50, nn.eps = 0)
ex.sub_dlpfc <- FindNeighbors(ex.sub_dlpfc, reduction = "Harmony_", k.param = 50, nn.eps = 0)

ex.sub_oFC <- FindClusters(ex.sub_oFC, algorithm = 1, resolution = 0.1)
ex.sub_dlpfc <- FindClusters(ex.sub_dlpfc, algorithm = 1, resolution = 0.1)

# =====================
# 7. Azimuth Annotation
# =====================
ex.sub_oFC <- RunAzimuth(ex.sub_oFC, reference = "/data/Azimuth_annotation")
ex.sub_dlpfc <- RunAzimuth(ex.sub_dlpfc, reference = "/data/Azimuth_annotation")

# =====================
# 8. Subset Specific Cell Types & Merge
# =====================
microglia_classes <- c("L2/3 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3", "L6b")

ex.sub_oFC_ss <- subset(ex.sub_oFC, subset = predicted.subclass %in% microglia_classes)
ex.sub_dlpfc_ss <- subset(ex.sub_dlpfc, subset = predicted.subclass %in% microglia_classes)

DefaultAssay(ex.sub_oFC_ss) <- "RNA"
DefaultAssay(ex.sub_dlpfc_ss) <- "RNA"

ex.sub_oFC_ss <- DietSeurat(ex.sub_oFC_ss)
ex.sub_dlpfc_ss <- DietSeurat(ex.sub_dlpfc_ss)

ex.sub. <- merge(x = ex.sub_oFC_ss, y = ex.sub_dlpfc_ss)
ex.sub. <- DietSeurat(ex.sub., layers = "counts")
ex.sub.@assays$sketch <- NULL

# =====================
# 9. Integration Across Datasets
# =====================
ex.sub.list <- SplitObject(ex.sub., split.by = "dataset")
features <- SelectIntegrationFeatures(object.list = ex.sub.list, nfeatures = 3000)

ex.sub. <- NormalizeData(ex.sub.)
VariableFeatures(ex.sub.) <- features
ex.sub. <- ScaleData(ex.sub., vars.to.regress = c("percent.mt", "nCount_RNA"))
ex.sub. <- RunPCA(ex.sub., npcs = 100)
ex.sub. <- RunHarmony(ex.sub., reduction.save = "Harmony_", 
                       group.by.vars = c("sampleID", "sex", "tissue"), plot_convergence = FALSE)

ex.sub. <- RunUMAP(ex.sub., dims = 1:50, reduction = "Harmony_")
ex.sub. <- FindNeighbors(ex.sub., reduction = "Harmony_", k.param = 50, nn.eps = 0)
ex.sub. <- FindClusters(ex.sub., algorithm = 1, resolution = c(0.1, 0.3))

ex.sub <- ex.sub.

# =====================
# 10. Marker Gene Analysis & Subclustering
# =====================
Idents(ex.sub) <- ex.sub$RNA_snn_res.0.1
mrks <- FindAllMarkers(ex.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dir.create("/data/Finalized_outputs/sketch_clusters/ex/MarkerGenes/", recursive = TRUE)
saveRDS(mrks, file = "/data/Finalized_outputs/sketch_clusters/ex/MarkerGenes/mrks_0.1.rds")

# Subset for specific subclasses
ex.sub.ss <- subset(ex.sub, subset = predicted.subclass %in% microglia_classes)
ex.sub.ss <- NormalizeData(ex.sub.ss)
ex.sub.ss <- FindVariableFeatures(ex.sub.ss)
ex.sub.ss <- ScaleData(ex.sub.ss, vars.to.regress = c("percent.mt", "nCount_RNA"))
ex.sub.ss <- RunPCA(ex.sub.ss, npcs = 50)
ex.sub.ss <- RunHarmony(ex.sub.ss, reduction.save = "Harmony_", 
                         group.by.vars = c("sampleID", "sex", "tissue"), plot_convergence = FALSE)
ex.sub.ss <- RunUMAP(ex.sub.ss, dims = 1:50, reduction = "Harmony_")
ex.sub.ss <- FindNeighbors(ex.sub.ss, reduction = "Harmony_", k.param = 50, nn.eps = 0)
ex.sub.ss <- FindClusters(ex.sub.ss, algorithm = 1, resolution = 0.2)

Idents(ex.sub.ss) <- ex.sub.ss$RNA_snn_res.0.1
mrks <- FindAllMarkers(ex.sub.ss, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(mrks, file = "/data/Finalized_outputs/sketch_clusters/ex/MarkerGenes/mrks_ss_0.1.rds")

ex.sub <- ex.sub.ss

# =====================
# 11. Save Outputs & Metadata
# =====================
dir.create("/data/Pheno_Demo/", recursive = TRUE)
meta <- data.frame(
  Cell_ID = rownames(ex.sub@meta.data),
  Classes = as.character(ex.sub$classes),
  Subclasses = as.character(ex.sub$subclasses),
  Cluster = ex.sub$L6_
)
write.csv(meta, "/data/Pheno_Demo/combined_meta_clusters_ex.csv", row.names = FALSE)

saveRDS(ex.sub, file = "/data/Seurat_obj/Subclusters/ex_tissue_final.rds")

# =====================
# 12. Propeller Analysis
# =====================
dir.create("/data/Finalized_outputs/sketch_clusters/ex/ctp/", recursive = TRUE)

ctp_1 <- propeller(x = ex.sub, clusters = ex.sub$subclasses, sample = ex.sub$sampleID,
                   group = ex.sub$orBatch, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_1, "/data/Finalized_outputs/sketch_clusters/ex/ctp/ctp_umap_2_batch_0.1.csv")

ctp_1.1 <- propeller(x = ex.sub, clusters = ex.sub$subclasses, sample = ex.sub$sampleID,
                     group = ex.sub$tissue, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_1.1, "/data/Finalized_outputs/sketch_clusters/ex/ctp/ctp_umap_2_tissue_0.1.csv")

ctp_1.2 <- propeller(x = ex.sub, clusters = ex.sub$subclasses, sample = ex.sub$sampleID,
                     group = ex.sub$version, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_1.2, "/data/Finalized_outputs/sketch_clusters/ex/ctp/ctp_umap_2_version_0.1.csv")

ctp_1.3 <- propeller(x = ex.sub, clusters = ex.sub$subclasses, sample = ex.sub$sampleID,
                     group = ex.sub$sex, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_1.3, "/data/Finalized_outputs/sketch_clusters/ex/ctp/ctp_umap_2_sex_0.1.csv")

ctp_1.4 <- propeller(x = ex.sub, clusters = ex.sub$subclasses, sample = ex.sub$sampleID,
                     group = ex.sub$gDx, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_1.4, "/data/Finalized_outputs/sketch_clusters/ex/ctp/ctp_umap_2_dx_0.1.csv")
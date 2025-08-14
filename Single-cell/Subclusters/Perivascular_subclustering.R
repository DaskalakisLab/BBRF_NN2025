########################################
# Perivascular cells Subcluster Processing - oFC and dlPFC
# Author: Artemis Iatrou
# Date: 2025-08-08
# Description:
#   Subsets perivascular cells by tissue, performs SCTransform, PCA, Harmony integration,
#   UMAP, clustering, Azimuth annotation, marker gene analysis, and saves outputs.
########################################

# =====================
# 0. Load Packages & Set Options
# =====================
library(Seurat)
library(future)
options(future.globals.maxSize = 1000000 * 1024^2)

# =====================
# 1. Load perivascular Seurat object
# =====================
peri_sub <- readRDS("/data/Seurat_obj/Subclusters/peri_sub.RDS")
peri_sub <- DietSeurat(peri_sub)

# =====================
# 2. Subset by tissue
# =====================
peri_sub_oFC <- subset(peri_sub, subset = tissue == "oFC")
peri_sub_dlpfc <- subset(peri_sub, subset = tissue == "dlPFC")

peri_sub_oFC$dataset <- "eb"
peri_sub_dlpfc$dataset <- "ai"

peri_sub_oFC$sampleID <- peri_sub_oFC$Donor
peri_sub_dlpfc$sampleID <- peri_sub_dlpfc$sample

peri_sub_oFC$tissue <- "oFC"
peri_sub_dlpfc$tissue <- "dlPFC"

peri_sub_oFC$pr_clusters <- peri_sub_oFC$celltypes_final
peri_sub_oFC$pmi <- peri_sub_oFC$PMI
peri_sub_oFC$orBatch <- "Binder_oFC"
peri_sub_oFC$version <- "v3"

peri_sub_oFC$sex <- as.factor(peri_sub_oFC$Sex)
levels(peri_sub_oFC$sex) <- c("female", "male")
peri_sub_oFC$age <- peri_sub_oFC$Age

peri_sub_oFC$gDx <- as.factor(peri_sub_oFC$Classification)
levels(peri_sub_oFC$gDx) <- c("Control", "MDD")

peri_sub_dlpfc$gDx <- as.factor(peri_sub_dlpfc$primarydx)
levels(peri_sub_dlpfc$gDx) <- c("Control", "MDD", "MDD")

# =====================
# 3. SCTransform and PCA
# =====================
peri_sub_oFC <- SCTransform(peri_sub_oFC, vars.to.regress = c("percent.mt", "nCount_RNA"))
peri_sub_dlpfc <- SCTransform(peri_sub_dlpfc, vars.to.regress = c("percent.mt", "nCount_RNA"))

peri_sub_oFC <- RunPCA(peri_sub_oFC, npcs = 50)
peri_sub_dlpfc <- RunPCA(peri_sub_dlpfc, npcs = 50)

dir.create("/data/Finalized_outputs/sketch_clusters/perivascular/", recursive = TRUE)

# =====================
# 4. Harmony Integration
# =====================
library(harmony)
set.seed(1234)

peri_sub_oFC <- RunHarmony(peri_sub_oFC, reduction.save= "Harmony_", group.by.vars = c("sampleID", "sex"), plot_convergence=FALSE)
peri_sub_dlpfc <- RunHarmony(peri_sub_dlpfc, reduction.save= "Harmony_", group.by.vars = c("orBatch", "sampleID", "version", "sex"), plot_convergence=FALSE)

# =====================
# 5. UMAP, Neighbors, Clustering
# =====================
peri_sub_oFC <- RunUMAP(peri_sub_oFC, dims = 1:30, reduction = "Harmony_")
peri_sub_dlpfc <- RunUMAP(peri_sub_dlpfc, dims = 1:30, reduction = "Harmony_")

peri_sub_oFC <- FindNeighbors(peri_sub_oFC, reduction = "Harmony_", k.param=30, nn.eps=0)
peri_sub_dlpfc <- FindNeighbors(peri_sub_dlpfc, reduction = "Harmony_", k.param=30, nn.eps=0)

peri_sub_oFC <- FindClusters(peri_sub_oFC, algorithm=1, resolution = c(0.1, 0.3))
peri_sub_dlpfc <- FindClusters(peri_sub_dlpfc, algorithm=1, resolution = c(0.1, 0.3))

# =====================
# 6. Azimuth Annotation
# =====================
library(Azimuth)
peri_sub_oFC <- RunAzimuth(peri_sub_oFC, reference = "/data/Azimuth_annotation")
peri_sub_dlpfc <- RunAzimuth(peri_sub_dlpfc, reference = "/data/Azimuth_annotation")

# =====================
# 7. Coray reference transfer
# =====================
coray <- readRDS("/data/PubData/Meninges_coray_coray_seurat_sct_pca_ready.rds")
coray <- SCTransform(coray)

anchors <- FindTransferAnchors(reference = coray, query = peri_sub_oFC, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = coray$cellID,
                                  weight.reduction = peri_sub_oFC[["Harmony_"]], dims = 1:30)
peri_sub_oFC[["prediction_ct_coray_all"]] <- predictions.assay$predicted.id

anchors <- FindTransferAnchors(reference = coray, query = peri_sub_dlpfc, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = coray$cellID,
                                  weight.reduction = peri_sub_dlpfc[["Harmony_"]], dims = 1:30)
peri_sub_dlpfc[["prediction_ct_coray_all"]] <- predictions.assay$predicted.id

# Subset to non-neuronal and non-glial cells
peri_sub_oFC_ss <- subset(peri_sub_oFC, subset = prediction_ct_coray_all %in% c("Neuron", "Oligo", "OPC", "Astrocyte"), invert=TRUE)
peri_sub_dlpfc_ss <- subset(peri_sub_dlpfc, subset = prediction_ct_coray_all %in% c("Neuron", "Oligo", "OPC", "Astrocyte"), invert=TRUE)

DefaultAssay(peri_sub_oFC_ss) <- "RNA"
DefaultAssay(peri_sub_dlpfc_ss) <- "RNA"

# =====================
# 8. Merge and integrated processing
# =====================
peri_sub. <- merge(x= peri_sub_oFC_ss, y = peri_sub_dlpfc_ss)
peri_sub. <- DietSeurat(peri_sub., layers = "counts")
peri_sub.@assays$SCT <- NULL

peri_sub.obj <- SplitObject(peri_sub., split.by = "dataset")
features <- SelectIntegrationFeatures(object.list = peri_sub.obj, nfeatures = 3000)

peri_sub. <- NormalizeData(peri_sub.)
VariableFeatures(peri_sub.) <- features
peri_sub. <- ScaleData(peri_sub., vars.to.regress = c("percent.mt", "nCount_RNA"))

peri_sub. <- RunPCA(peri_sub., npcs = 50)
peri_sub. <- RunHarmony(peri_sub., reduction.save= "Harmony_", group.by.vars = c("sampleID", "sex", "tissue"), plot_convergence=FALSE)

peri_sub. <- RunUMAP(peri_sub., dims = 1:10, reduction = "Harmony_")
peri_sub. <- FindNeighbors(peri_sub., reduction = "Harmony_", k.param=10, nn.eps=0)
peri_sub. <- FindClusters(peri_sub., algorithm=1, resolution = c(0.1, 0.3))

anchors <- FindTransferAnchors(reference = coray, query = peri_sub., normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = coray$cellID,
                                  weight.reduction = peri_sub.[["Harmony_"]], dims = 1:30)
peri_sub.[["prediction_ct_coray_all"]] <- predictions.assay$predicted.id

peri_sub. <- RunAzimuth(peri_sub., reference = "/data/Azimuth_annotation")

# =====================
# 9. Save final object
# =====================
saveRDS(peri_sub., file = "/data/Seurat_obj/Subclusters/peri_tissue_final.rds")

# =====================
# 10. Define markers and subclasses/classes
# =====================
detect.genes <- c("CLDN5", "PDGFRB", "ACTA2", "CEMIP", "MBP", "AQP4", "DSCAM", "SYT1")
peri.specific <- c("SLC12A7", "LAMC3", "CTNNA3", "ACTA2", "DLC1", "TAGLN")
fibro.specific <- c("KCNMA1", "SLC4A4", "PLXNA4", "NRG3", "MYRIP", "LAMA2", "FBLN1")

# Assign subclasses and classes (example)
peri_sub.$subclasses <- as.factor(peri_sub.$SCT_snn_res.0.1)
levels(peri_sub.$subclasses) <- c("Pericytes_1", "Pericytes_2", "Fibro_peri", "SMC", "Endo_cap", "Endo_veinous", "Fibro_men")

peri_sub.$classes <- as.factor(peri_sub.$SCT_snn_res.0.1)
levels(peri_sub.$classes) <- c("Mural", "Mural", "Fibro", "Mural", "Endo", "Endo", "Fibro")

# =====================
# 11. Propeller cluster proportion testing (CTP)
# =====================
library(speckle)
dir.create("/data/Finalized_outputs/sketch_clusters/perivascular/ctp/", recursive = TRUE)

# Subclass-level CTP
ctp_sub <- propeller(x = peri_sub., clusters = peri_sub.$subclasses,
                     sample = peri_sub.$sampleID,
                     group = peri_sub.$orBatch, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_sub, "/data/Finalized_outputs/sketch_clusters/perivascular/ctp/ctp_subclasses_batch.csv")

ctp_sub_tissue <- propeller(x = peri_sub., clusters = peri_sub.$subclasses,
                            sample = peri_sub.$sampleID,
                            group = peri_sub.$tissue, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_sub_tissue, "/data/Finalized_outputs/sketch_clusters/perivascular/ctp/ctp_subclasses_tissue.csv")

ctp_sub_sex <- propeller(x = peri_sub., clusters = peri_sub.$subclasses,
                         sample = peri_sub.$sampleID,
                         group = peri_sub.$sex, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_sub_sex, "/data/Finalized_outputs/sketch_clusters/perivascular/ctp/ctp_subclasses_sex.csv")

ctp_sub_dx <- propeller(x = peri_sub., clusters = peri_sub.$subclasses,
                        sample = peri_sub.$sampleID,
                        group = peri_sub.$gDx, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_sub_dx, "/data/Finalized_outputs/sketch_clusters/perivascular/ctp/ctp_subclasses_dx.csv")

# Class-level CTP
ctp_cls <- propeller(x = peri_sub., clusters = peri_sub.$classes,
                     sample = peri_sub.$sampleID,
                     group = peri_sub.$orBatch, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_cls, "/data/Finalized_outputs/sketch_clusters/perivascular/ctp/ctp_classes_batch.csv")

ctp_cls_tissue <- propeller(x = peri_sub., clusters = peri_sub.$classes,
                            sample = peri_sub.$sampleID,
                            group = peri_sub.$tissue, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_cls_tissue, "/data/Finalized_outputs/sketch_clusters/perivascular/ctp/ctp_classes_tissue.csv")

ctp_cls_sex <- propeller(x = peri_sub., clusters = peri_sub.$classes,
                         sample = peri_sub.$sampleID,
                         group = peri_sub.$sex, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_cls_sex, "/data/Finalized_outputs/sketch_clusters/perivascular/ctp/ctp_classes_sex.csv")

ctp_cls_dx <- propeller(x = peri_sub., clusters = peri_sub.$classes,
                        sample = peri_sub.$sampleID,
                        group = peri_sub.$gDx, trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_cls_dx, "/data/Finalized_outputs/sketch_clusters/perivascular/ctp/ctp_classes_dx.csv")

# =====================
# 12. Save metadata
# =====================
meta <- data.frame(
  Cell_ID = rownames(peri_sub.@meta.data),
  Classes = as.character(peri_sub.$classes),
  Subclasses = as.character(peri_sub.$subclasses),
  Cluster = peri_sub.$SCT_snn_res.0.1
)
write.csv(meta, "/data/Pheno_Demo/combined_meta_clusters_periv.csv", row.names = FALSE)
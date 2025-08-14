########################################
# Astrocytes Subcluster Processing - oFC and dlPFC + Das Reference
# Author: Artemis Iatrou
# Date: 2025-08-08
# Description:
#   Subsets astrocytes by tissue, performs SCTransform, PCA, Harmony integration,
#   UMAP, clustering, Azimuth annotation, marker gene analysis,
#   Das reference integration, module scoring, propeller analysis, and saves outputs.
########################################

# =====================
# 0. Load Packages & Set Options
# =====================

library(Seurat)
library(future)
library(harmony)
library(ggplot2)
library(Matrix)
library(dplyr)
library(readxl)
library(speckle)

options(future.globals.maxSize = 1000000 * 1024^2)

# =====================
# 1. Load Astrocyte Subcluster Object
# =====================

astro_sub <- readRDS("/data/Seurat_obj/Subclusters/astro_sub.RDS")
astro_sub <- DietSeurat(astro_sub)

astro_sub_oFC <- subset(astro_sub, subset = tissue == "oFC")
mito <- astro_sub_oFC$percent.mt
astro_sub_dlpfc <- subset(astro_sub, subset = tissue == "dlPFC")

astro_sub_oFC$dataset <- "eb"
astro_sub_dlpfc$dataset <- "ai"

astro_sub_oFC$sampleID <- astro_sub_oFC$Donor
astro_sub_dlpfc$sampleID <- astro_sub_dlpfc$sample

astro_sub_oFC$tissue <- "oFC"
astro_sub_dlpfc$tissue <- "dlPFC"

astro_sub_oFC$pr_clusters <- astro_sub_oFC$celltypes_final
astro_sub_oFC$pmi <- astro_sub_oFC$PMI
astro_sub_oFC$orBatch <- "Binder_oFC"
astro_sub_oFC$version <- "v3"

astro_sub_oFC$sex <- factor(astro_sub_oFC$Sex, levels = c("F", "M"), labels = c("female", "male"))
astro_sub_dlpfc$gDx <- factor(astro_sub_dlpfc$primarydx, levels = c("Control", "MDD", "MDD"))

astro_sub_oFC$gDx <- factor(astro_sub_oFC$Classification, levels = c("Control", "MDD"))
astro_sub_oFC$age <- astro_sub_oFC$Age

# =====================
# 2. SCTransform & PCA
# =====================

astro_sub_oFC <- SCTransform(astro_sub_oFC, vars.to.regress = c("percent.mt", "nCount_RNA"))
astro_sub_dlpfc <- SCTransform(astro_sub_dlpfc, vars.to.regress = c("percent.mt", "nCount_RNA"))

astro_sub_oFC <- RunPCA(astro_sub_oFC, npcs = 50)
astro_sub_dlpfc <- RunPCA(astro_sub_dlpfc, npcs = 50)

dir.create("/data/Finalized_outputs/sketch_clusters/astro/", recursive = TRUE)

set.seed(1234)

astro_sub_oFC <- RunHarmony(astro_sub_oFC, reduction.save = "Harmony_", 
                            group.by.vars = c("sampleID", "sex"), plot_convergence = FALSE)

astro_sub_dlpfc <- RunHarmony(astro_sub_dlpfc, reduction.save = "Harmony_", 
                              group.by.vars = c("orBatch", "sampleID", "version", "sex"), plot_convergence = FALSE)

# =====================
# 3. UMAP & Clustering
# =====================

astro_sub_oFC <- RunUMAP(astro_sub_oFC, dims = 1:30, reduction = "Harmony_")
astro_sub_dlpfc <- RunUMAP(astro_sub_dlpfc, dims = 1:30, reduction = "Harmony_")

astro_sub_oFC <- FindNeighbors(astro_sub_oFC, reduction = "Harmony_", k.param = 30, nn.eps = 0)
astro_sub_dlpfc <- FindNeighbors(astro_sub_dlpfc, reduction = "Harmony_", k.param = 30, nn.eps = 0)

astro_sub_oFC <- FindClusters(astro_sub_oFC, algorithm = 1, resolution = c(0.1, 0.3))
astro_sub_dlpfc <- FindClusters(astro_sub_dlpfc, algorithm = 1, resolution = c(0.1, 0.3))

# =====================
# 4. Azimuth Annotation
# =====================

astro_sub_oFC <- RunAzimuth(astro_sub_oFC, reference = "/data/Azimuth_annotation")
astro_sub_dlpfc <- RunAzimuth(astro_sub_dlpfc, reference = "/data/Azimuth_annotation")

astro_sub_oFC_ss <- subset(astro_sub_oFC, subset = predicted.subclass == "Astro")
astro_sub_dlpfc_ss <- subset(astro_sub_dlpfc, subset = predicted.subclass == "Astro")
astro_sub_oFC_ss$percent.mt <- mito

DefaultAssay(astro_sub_oFC_ss) <- "RNA"
DefaultAssay(astro_sub_dlpfc_ss) <- "RNA"

# =====================
# 5. Merge & Integrate Tissues
# =====================

astro_sub. <- merge(x = astro_sub_oFC_ss, y = astro_sub_dlpfc_ss)
astro_sub. <- DietSeurat(astro_sub., layers = "counts")
astro_sub.@assays$SCT <- NULL

astro_sub.obj <- SplitObject(astro_sub., split.by = "dataset")
features <- SelectIntegrationFeatures(object.list = astro_sub.obj, nfeatures = 3000)

astro_sub. <- NormalizeData(astro_sub.)
VariableFeatures(astro_sub.) <- features
astro_sub. <- ScaleData(astro_sub., vars.to.regress = c("percent.mt", "nCount_RNA"))
astro_sub. <- RunPCA(astro_sub., npcs = 50)
astro_sub. <- RunHarmony(astro_sub., reduction.save = "Harmony_", 
                         group.by.vars = c("sampleID", "sex", "tissue", "orBatch"), plot_convergence = FALSE)

astro_sub. <- RunAzimuth(astro_sub., reference = "/data/Azimuth_annotation")
astro_sub. <- RunUMAP(astro_sub., dims = 1:50, reduction = "Harmony_")
astro_sub. <- FindNeighbors(astro_sub., reduction = "Harmony_", k.param = 50, nn.eps = 0)
astro_sub. <- FindClusters(astro_sub., algorithm = 1, resolution = c(0.1, 0.3))

astro_sub <- astro_sub.

# =====================
# 6. Plot Categorical & Continuous Features
# =====================

# (Categorical and continuous plotting sections as in other cell types)

# =====================
# 7. Marker Gene Analysis
# =====================

Idents(astro_sub) <- astro_sub$RNA_snn_res.0.1
mrks <- FindAllMarkers(astro_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(mrks, "/data/Finalized_outputs/sketch_clusters/astro/MarkerGenes/mrks_0.1_2.rds")

other.specific <- c("AQP4", "GFAP", "CLU", "CHI3L1", "S100A10", "EMP1")
pdf("/data/Finalized_outputs/sketch_clusters/astro/MarkerGenes/astro_sub_celltypes_1_2.pdf")
VlnPlot(astro_sub, features = other.specific, pt.size = 0, split.by = "dataset")
dev.off()

# =====================
# 8. Subcluster Refinement & SCTransform
# =====================

astro_sub <- subset(astro_sub, idents = 4, invert = TRUE)
astro_sub_ss <- SCTransform(astro_sub, vars.to.regress = c("percent.mt", "nCount_RNA"))
astro_sub_ss <- RunPCA(astro_sub_ss, npcs = 50)
astro_sub_ss <- RunHarmony(astro_sub_ss, reduction.save = "Harmony_", 
                           group.by.vars = c("sampleID", "sex", "tissue"), plot_convergence = FALSE)
astro_sub_ss <- RunUMAP(astro_sub_ss, dims = 1:10, reduction = "Harmony_")
astro_sub_ss <- FindNeighbors(astro_sub_ss, reduction = "Harmony_", k.param = 10, nn.eps = 0)
astro_sub_ss <- FindClusters(astro_sub_ss, algorithm = 1, resolution = c(0.1, 0.3))

# =====================
# 9. Das Reference Integration
# =====================

########################################
# Load Counts, Metadata, and Create Seurat Object
########################################

library(Matrix)
library(Seurat)
library(readxl)
library(dplyr)
library(harmony)

# Load counts and annotations
counts <- readMM("/data/Das_astrocytes/ad_progression_download_ui_1-download_matrix_PFC_astro?w=")
key <- read.table("/data/Das_astrocytes/ad_progression_download_ui_1-download_row_annotation_PFC_astro?w=")
cluster <- read.table("/data/Das_astrocytes/ad_progression_download_ui_1-download_cell_annotation_PFC_astro?w=")

rownames(counts) <- key$V1
colnames(counts) <- cluster$V1

# Create initial Seurat object
seurat_obj <- CreateSeuratObject(counts = counts)

# Load metadata
meta <- readxl::read_excel("/data/Das_astrocytes/2025-04-24 _metadata.xlsx", sheet = "BA46")
metadata <- as.data.frame(meta[match(cluster$V1, meta$cell_annotation), ])
rownames(seurat_obj@meta.data) <- Cells(seurat_obj)
seurat_obj@meta.data <- metadata

# Normalize, find variable features, scale, and PCA
ast <- seurat_obj %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npсs = 20, verbose = FALSE)

# Harmony integration
ast <- RunHarmony(ast, group.by.vars = "SampleName", plot_convergence = TRUE, lambda = 1)

# Clustering and UMAP
ast <- ast %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.2) %>%
  RunUMAP(reduction = "harmony", dims = 1:20, n_neighbors = 50, seed.use = 42)

# Visualize clusters
DimPlot(ast, reduction = "umap", pt.size = 1, group.by = "seurat_clusters")

# Find cluster markers
cluster.markers <- FindAllMarkers(ast, min.pct = 0.1, only.pos = TRUE)

# DotPlot for selected genes
genes_unique <- c("SGCD","GRM3","SLC1A2","HPSE2","NRXN1","EGFR","GLUL","PTN","MERTK","GRIA2",
                  "KCNIP4","CNTNAP2","NRXN3","SYT1","NLGN1","FGF12","MEG3","FGF14","RYR2","GRM5",
                  "ETNPPL","DNAH7","SLC14A1","ADGRV1","TENM2","PCLO","GABRB1","NTRK2","RBMS3","GLUD1",
                  "VEGFA","HSPH1","ANGPTL4","DNAJB1","HSPB1","HSPA1A","HSP90AA1","HSPA1B","BAG3","MT2A",
                  "CHI3L1","ELOVL5","NRP1","DST","MAPRE2","F3","SCARA3","OSMR","SERPINA3",
                  "SPP1","TBXAS1","HDAC9","GRID2","MEF2A","FMN1","PTPRC","INPP5D","SYK","MEF2C",
                  "PLCE1","CREB5","COL21A1","LAMA2","FAP","MAP1B","SLC7A11","CACNA1C","C3",
                  "TNC","CD44","ADAMTSL3","GPC6","VCAN","SLC38A1","CRYAB","AQP1","CP","ANGPT1",
                  "DPP10","GRIA1","TJP2","SPARC","KCNN2","UBC","FOS","GFAP","HSPB8","AQP4")

DotPlot(ast, features = genes_unique, cols = c("grey90", "darkred"), cluster.idents = TRUE, dot.min = 1) +
  coord_flip()
ggsave("/data/Das_astrocytes/dotplot_genes_unique.pdf", width = 20, height = 20)

########################################
# Transfer Anchors and Predictions to astro_sub_ss_log
########################################

# Normalize and scale query
astro_sub_ss_log <- astro_sub_ss %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE)

DefaultAssay(astro_sub_ss_log) <- "RNA"

# Find transfer anchors
anchors <- FindTransferAnchors(
  reference = ast,
  query = astro_sub_ss_log,
  normalization.method = "LogNormalize",
  recompute.residuals = FALSE
)

# Transfer predicted cluster IDs
predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = ast$seurat_clusters,
  weight.reduction = astro_sub_ss_log[["Harmony_"]],
  dims = 1:30
)

astro_sub_ss_log[["prediction_astro"]] <- predictions.assay$predicted.id

# Compare predicted IDs to original subclasses
table(astro_sub_ss_log$subclasses, astro_sub_ss_log$prediction_astro)

# Continue downstream workflow with SCT assay, UMAP, neighbors, clusters
DefaultAssay(astro_sub_ss_log) <- "SCT"
astro_sub_ss_log <- astro_sub_ss_log %>%
  RunUMAP(dims = 1:10, reduction = "Harmony_") %>%
  FindNeighbors(reduction = "Harmony_", k.param = 10, nn.eps = 0) %>%
  FindClusters(algorithm = 1, resolution = c(0.08))

# Visualize transferred annotations
DimPlot(astro_sub_ss_log)

# =====================
# 10. Module Scoring
# =====================

# genes_split list remains unchanged
astro_sub_ss_log <- AddModuleScore(astro_sub_ss_log, genes_split, name = "ast_")

# =====================
# 11. Propeller Analysis
# =====================

ctp_2 <- propeller(x = astro_sub_ss_log, clusters = astro_sub_ss_log$subclasses,
                   sample = astro_sub_ss_log$sampleID, group = astro_sub_ss_log$orBatch,
                   trend = FALSE, robust = TRUE, transform = "logit")
write.csv(ctp_2, "/data/Finalized_outputs/sketch_clusters/astro/ctp/ctp_2_umap_batch_Das.csv")

# Repeat propeller for tissue, version, sex, gDx 

# =====================
# 12. Save Final Objects & Metadata
# =====================

saveRDS(astro_sub_ss, file = "/data/Seurat_obj/Subclusters/astro_tissue_final2.rds")

meta <- data.frame(
  Cell_ID = rownames(astro_sub_ss_log@meta.data),
  Classes = as.character(astro_sub_ss_log$classes),
  Subclasses = as.character(astro_sub_ss_log$subclasses),
  Cluster = astro_sub_ss_log$SCT_snn_res.0.08
)

write.csv(meta, "/data/Pheno_Demo/combined_meta_clusters_astro_DAS.csv", row.names = FALSE)
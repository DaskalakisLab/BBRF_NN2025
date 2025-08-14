########################################
# Across Brain Regions Bulk RNA-Seq Analysis
# Author: Artemis Iatrou
# Date: 2025-08-08
# Description:
#   This script combines mPFC, DLPFC, and dACC bulk RNA-Seq data,
#   performs normalization, surrogate variable analysis (SVA),
#   runs limma/voom with dream, visualizes results, and
#   conducts GSEA enrichment analyses (GO BP and positional).
########################################

# =====================
# Load Required Packages
# =====================
library(edgeR)
library(variancePartition)
library(BiocParallel)
library(sva)
library(bacon)
library(ggrepel)
library(ggplot2)
library(msigdbr)
library(dplyr)
library(tibble)
library(fgsea)

# =====================
# User Parameters
# =====================
# Define paths for input and output
data_dir <- "data/BBRF/Bulk"
output_dir <- "results/AcrossRegions"
pathway_dir <- file.path(output_dir, "Pathway")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)

regions <- c("mPFC", "DLPFC", "dACC")

# =====================
# Load & Annotate Metadata
# =====================
anno_list <- list()
for (region in regions) {
  anno <- readRDS(file.path(data_dir, region, paste0("anno_", region, ".RDS")))
  anno$R_ID <- paste0(anno$SampleID, "_", region)
  anno$Region <- region
  anno_list[[region]] <- anno
}

# Keep only common columns across all regions
common_cols <- Reduce(intersect, lapply(anno_list, colnames))
anno_list <- lapply(anno_list, function(df) df[, common_cols, drop = FALSE])
anno_combined <- do.call(rbind, anno_list)

# =====================
# Load Counts & Normalize
# =====================
counts_list <- list()
for (region in regions) {
  dge_obj <- readRDS(file.path(data_dir, region, paste0("DGE_", region, ".RDS")))
  counts <- dge_obj$counts
  colnames(counts) <- paste0(colnames(counts), "_", region)
  counts_list[[region]] <- counts
}

# Keep only common genes across all regions
common_genes <- Reduce(intersect, lapply(counts_list, rownames))
counts_list <- lapply(counts_list, function(mat) mat[common_genes, , drop = FALSE])
counts_combined <- do.call(cbind, counts_list)

# Filter low-expression genes
isexpr <- rowSums(cpm(counts_combined) > 1) >= 5
dge <- DGEList(counts_combined[isexpr, ])
dge <- calcNormFactors(dge)

# =====================
# Surrogate Variable Analysis (SVA)
# =====================
mm <- model.matrix(~ gDx, anno_combined)
mm0 <- model.matrix(~ 1, anno_combined)
fit <- svaseq(cpm(dge$counts), mod = mm, mod0 = mm0)

SVs <- fit$sv
colnames(SVs) <- paste0("sv", seq_len(ncol(SVs)))
anno_combined <- cbind(anno_combined, SVs)

# Save intermediate objects
save(anno_combined, counts_combined, dge,
     file = file.path(output_dir, "objects_DGE.rda"))

# =====================
# Covariate Correlations
# =====================
covariates <- c("gDx", "Age", "Sex", "bestpop", "PMI", "RIN", "mitoMapped", colnames(SVs))
form <- as.formula(paste("~", paste(covariates, collapse = "+")))
C <- canCorPairs(form, anno_combined)
plotCorrMatrix(C)

# =====================
# Dream Analysis
# =====================
param <- SnowParam(1, "SOCK", progressbar = TRUE)
form <- ~ gDx + Sex + Age + (1 | Region) + (1 | BrNum) + sv1 + sv2 + sv3

vobjDream <- voomWithDreamWeights(dge, form, anno_combined, BPPARAM = param)
fitmm <- dream(vobjDream, form, anno_combined)
fitmm <- eBayes(fitmm)
results <- topTable(fitmm, coef = "gDxMDD", number = Inf)

# Map gene symbols
RNAGeneMap <- readRDS("data/GeneMappings/RNAGeneMap.RDS")
results$symbol <- RNAGeneMap$symbol[match(rownames(results), RNAGeneMap$genes)]

# Bacon correction
z <- abs(qnorm(results$P.Value)) * sign(results$logFC)
results$z <- ifelse(is.na(z), 0, z)
bc <- bacon(results$z)
results$p_bacon <- as.numeric(pval(bc))
results$p_bacon_adj <- p.adjust(results$p_bacon, method = "fdr")

# Save results
saveRDS(results, file.path(output_dir, "results_gDxMDD_All.RDS"))

# =====================
# Volcano Plot
# =====================
results$lp <- -log10(results$P.Value)
results$Color <- 1
results[(results$P.Value < 0.05) & (results$logFC < 0), "Color"] <- 2
results[(results$P.Value < 0.05) & (results$logFC > 0), "Color"] <- 3
results[(results$adj.P.Val < 0.05) & (results$logFC < 0), "Color"] <- 4
results[(results$adj.P.Val < 0.05) & (results$logFC > 0), "Color"] <- 5

col_df <- data.frame(
  Vals = 1:5,
  Colors = c("grey", "lightblue", "pink", "darkblue", "darkred")
)
results <- results[order(results$P.Value), ]
results$genelabels <- FALSE
results$genelabels[1:10] <- TRUE
highlight_genes <- c("FKBP5", "NR3C1", "NR3C2", "CRH", "CRHR1",
                     "CRHR2", "CORT", "IL1RL1", "LINC00996")
results$genelabels[results$symbol %in% highlight_genes] <- TRUE

p <- ggplot(results) +
  geom_point(aes(x = logFC, y = lp, colour = factor(Color))) +
  geom_text_repel(aes(label = ifelse(genelabels, symbol, "")),
                  fontface = "italic", size = 3, max.overlaps = 1000) +
  scale_color_manual(values = col_df$Colors) +
  xlab("log2(FC)") +
  ylab("-log10(p-value)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 20, face = "bold"))
ggsave(file.path(output_dir, "volcano_plot.pdf"), p)

# =====================
# GSEA Enrichment
# =====================
run_fgsea <- function(pathways, ranks, min_size, max_size, out_file) {
  fgsea_sets <- pathways %>% split(x = .$gene_symbol, f = .$gs_name)
  fgseaRes <- fgseaMultilevel(fgsea_sets, stats = ranks,
                              minSize = min_size, maxSize = max_size, eps = 0)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  saveRDS(fgseaResTidy, out_file)
}

gsea_input <- results %>%
  arrange(desc(logFC)) %>%
  select(symbol, logFC)
ranks <- deframe(gsea_input)

# GO Biological Process
pathways_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
run_fgsea(pathways_bp, ranks, min_size = 10, max_size = 500,
          out_file = file.path(pathway_dir, "C5_BP_All.rds"))

# Positional Enrichment
pathways_pos <- msigdbr(species = "Homo sapiens", category = "C1")
run_fgsea(pathways_pos, ranks, min_size = 150, max_size = 500,
          out_file = file.path(pathway_dir, "C1_All.rds"))
########################################
# TWAS/SMR/DGE Overlap & Manhattan/Miami Plots
# Author: Artemis Iatrou
# Date: 2025-08-08
# Description:
#   This script reads TWAS, SMR, and DGE results,
#   identifies significant genes, computes overlaps,
#   plots Venn diagrams, and generates Manhattan/Miami plots
########################################

# =====================
# 0. Load Required Packages
# =====================
library(dplyr)
library(ggvenn)
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

# =====================
# 1. Load TWAS Data & Define Significant Genes
# =====================
TWAS_Cell2025 <- read_csv("/data/BBRF/SMR/TWAS_Cell2025.csv")

twas_dlpfc <- TWAS_Cell2025 %>% filter(Panel == "CMC DLPFC") %>% pull(ID)
twas_pec <- TWAS_Cell2025 %>% filter(Panel == "PsychENCODE") %>% pull(ID)
twas_fc <- TWAS_Cell2025 %>% filter(Panel == "GTEx Frontal Cortex") %>% pull(ID)

sig_twas <- c(twas_dlpfc, twas_fc, twas_pec)

# =====================
# 2. Load DGE and SMR Data
# =====================
dge <- readRDS("/data/BBRF/Bulk/AcrossRegions/results_gDxMDD_All.RDS")
smr <- read.table("/data/humgen/daskalakislab/GWAS/SMR/MDD_PGC_2025/eQTL/trait_eSMR.merged.tsv",
                  header = TRUE, fill = TRUE)

# Significant DGE
sig_dge <- dge %>% filter(adj.P.Val < 0.05)

# Significant SMR (BrainMeta)
sig_smr <- smr %>%
  filter(qtl_name == "eQTL_BrainMeta") %>%
  mutate(symbol = index,
         fdr = p.adjust(p_SMR, method = "BH")) %>%
  filter(fdr < 0.05 & p_HEIDI > 0.05)

# =====================
# 3. Plot Venn Diagram of Gene Overlap
# =====================
gene_lists <- list(
  "SMR" = sig_smr$symbol,
  "TWAS" = sig_twas,
  "DGE" = sig_dge$symbol
)
gene_lists <- lapply(gene_lists, function(x) x[!is.na(x)])

# Basic ggvenn
ggvenn(gene_lists)

# Identify common genes
common_genes <- Reduce(intersect, gene_lists)
print(common_genes)

# VennDiagram plot
myCol <- brewer.pal(3, "Accent")
venn.diagram(
  x = gene_lists,
  category.names = c("SMR", "TWAS", "DGE"),
  filename = '/data/BBRF/SMR/eSMR/Visualisations/PEC_TWAS.svg',
  output = TRUE,
  imagetype = "svg",
  height = 480,
  width = 480,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  col = c("#4BBF4B", "#8A7BC1", "#FDB65A"),
  fill = myCol,
  cex = 0.6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#4BBF4B", "#8A7BC1", "#FDB65A"),
  rotation = 1
)

# =====================
# 4. Prepare Data for Manhattan/Miami Plots
# =====================
load("/data/humgen/daskalakislab/jajoo/references/hg38/Homo_sapeins.GRCh38.109.rdata")
gtf_df <- as.data.frame(gtf)

# Map gene ranges
gtf_gene_ranges <- gtf_df %>%
  group_by(gene_name) %>%
  reframe(
    chr = unique(seqnames),
    start_min = min(start, na.rm = TRUE),
    end_max = max(end, na.rm = TRUE)
  )

# Merge with SMR data
smr_filtered <- smr %>%
  filter(qtl_name == "eQTL_BrainMeta") %>%
  mutate(fdr = p.adjust(p_SMR, method = "BH")) %>%
  select(symbol, p_SMR, p_HEIDI, fdr, b_SMR)

ds_comb <- merge(smr_filtered, gtf_gene_ranges, by.x = "symbol", by.y = "gene_name")

# Compute cumulative positions for Manhattan plot
data_cum <- ds_comb %>%
  group_by(chr) %>%
  summarise(max_bp = max(start_min)) %>%
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(chr, bp_add) %>%
  slice(1:22)

df <- ds_comb %>%
  inner_join(data_cum, by = "chr") %>%
  mutate(bp_cum = start_min + bp_add)

axis_set <- df %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum))

# Define top genes for labeling
df_ss_10 <- df %>%
  filter(fdr < 0.05 & p_HEIDI > 0.05) %>%
  arrange(fdr) %>%
  slice_head(n = 20)

df_fdr <- df %>% filter(fdr < 0.05)

# Intersect with DGE and TWAS
intersect_dge <- intersect(sig_dge$symbol, df_ss_10$symbol)
intersect_twas <- intersect(sig_twas, df_ss_10$symbol)

# Add plotting columns
df <- df %>%
  mutate(
    lp = -log10(p_SMR),
    Label = symbol %in% df_ss_10$symbol,
    Shape = case_when(
      symbol %in% intersect_dge & symbol %in% intersect_twas ~ "Both DGE and TWAS",
      symbol %in% intersect_dge ~ "DGE",
      symbol %in% intersect_twas ~ "TWAS",
      symbol %in% df_fdr$symbol ~ "FDR",
      TRUE ~ "Other"
    ),
    Color = case_when(
      Shape == "Both DGE and TWAS" ~ "DGE, TWAS, SMR",
      Shape == "DGE" ~ "DGE",
      Shape == "TWAS" ~ "TWAS",
      Shape == "FDR" ~ "FDR",
      as.numeric(chr) %% 2 == 0 ~ "Z2",
      TRUE ~ "Z1"
    ),
    Color = factor(Color, levels = c("Z1", "Z2", "FDR", "DGE, TWAS, SMR", "DGE", "TWAS"))
  )

# =====================
# 5. Manhattan Plot
# =====================
manhattan_plot <- ggplot(df) +
  geom_point(aes(x = bp_cum, y = lp, color = Color), size = 1) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(limits = c(0, 15), breaks = c(0, 5, 10, 15)) +
  scale_color_manual(values = c("grey70", "grey90", "#4BBF4B", "#2F81C7", "#FDB65A", "#8A7BC1")) +
  labs(x = "Chromosome", y = expression("-log"[10]*"p-value"), col = "FDR SMR") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(color = "white", fill = "white"),
    plot.background = element_rect(color = "white", fill = "white"),
    axis.text.x = element_text(size = 10)
  ) +
  geom_text_repel(
    data = subset(df, Label == TRUE & Color != "DGE, TWAS, SMR"),
    aes(label = symbol, x = bp_cum, y = lp),
    size = 3.5
  ) +
  geom_text_repel(
    data = subset(df, Color == "DGE, TWAS, SMR"),
    aes(label = symbol, x = bp_cum, y = lp),
    color = "#2F81C7",
    size = 4.5,
    fontface = "bold"
  )

ggsave("/data/BBRF/SMR/eSMR/Visualisations/ManhattanPlot_PEC.pdf",
       plot = manhattan_plot, width = 15, height = 5)

# =====================
# 6. Miami Plot (b_SMR)
# =====================
miami_plot <- ggplot(df) +
  geom_point(aes(x = bp_cum, y = b_SMR, color = Color), size = 1) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_color_manual(values = c("grey70", "grey90", "#4BBF4B", "#2F81C7", "#FDB65A", "#8A7BC1")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Chromosome", y = "b_SMR", col = "FDR SMR") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(color = "white", fill = "white"),
    plot.background = element_rect(color = "white", fill = "white"),
    axis.text.x = element_text(size = 10)
  ) +
  geom_text_repel(
    data = subset(df, Label == TRUE),
    aes(label = symbol, x = bp_cum, y = b_SMR),
    size = 6
  ) +
  geom_text_repel(
    data = subset(df, Color == "DGE, TWAS, SMR"),
    aes(label = symbol, x = bp_cum, y = b_SMR),
    color = "#2F81C7",
    size = 6
  )

ggsave("/data/BBRF/SMR/eSMR/Visualisations/MiamiPlot.pdf",
       plot = miami_plot, width = 15, height = 7)

print("TWAS/SMR/DGE analysis completed: Venn diagram, Manhattan, and Miami plots saved.")
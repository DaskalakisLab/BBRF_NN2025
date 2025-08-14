library(dplyr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)

# Subset eph/ephrin genes and label by cell type
eph_ephrin_genes_df <- imap_dfr(
  AD_data_list_PFC_cogdx,
  ~ .x %>%
    filter(grepl("^EPHA|^EPHB|^EFN", gene)) %>%
    mutate(cell_type = .y)
)

# Prepare AD volcano plot data
df <- eph_ephrin_genes_df %>%
  mutate(
    pval = 10^(-log10p_nm),
    fdr_p = p.adjust(pval, method = "BH"),
    pval_capped = pmin(pmax(log10p_nm, -100), 100),
    signif = ifelse(fdr_p < 0.05, "Significant", "Not Significant")
  )

n_celltypes <- length(unique(df$cell_type))
palette <- RColorBrewer::brewer.pal(min(n_celltypes, 8), "Dark2")

# Volcano plot AD
ggplot(df, aes(x = logFC_nb, y = pval_capped, color = cell_type, shape = signif)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(aes(label = gene), size = 2.5, max.overlaps = Inf) +
  scale_color_manual(values = palette) +
  scale_shape_manual(values = c("Not Significant" = 1, "Significant" = 19)) +
  theme_classic() +
  labs(x = "Log2 Fold Change", y = "-log10(p-value)", color = "Cell Type", shape = "Significance") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40")

# Subset eph/ephrin genes MDD
eph_ephrin_genes_mdd_df <- imap_dfr(
  MDD1_data_list,
  ~ .x %>%
    filter(grepl("^EPHA|^EPHB|^EFN", genes)) %>%
    mutate(cell_type = .y)
)

df2 <- eph_ephrin_genes_mdd_df %>%
  mutate(
    pval = P.Value,
    signif = ifelse(adj.P.Val < 0.05, "Significant", "Not Significant")
  )

# Volcano plot MDD
ggplot(df2, aes(x = logFC, y = -log10(pval), color = cell_type, shape = signif)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(aes(label = genes), size = 2.5, max.overlaps = Inf) +
  scale_shape_manual(values = c("Not Significant" = 1, "Significant" = 19)) +
  theme_classic() +
  labs(x = "Log2 Fold Change", y = "-log10(p-value)", color = "Cell Type", shape = "Significance") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40")

# Clean and merge
df_clean <- df %>% rename(logFC_df = logFC_nb)
df2_clean <- df2 %>% rename(gene = genes, logFC_df2 = logFC)
merged_df <- inner_join(df_clean, df2_clean, by = c("gene", "cell_type"))
merged_df$cell_type <- as.factor(merged_df$cell_type)
levels(merged_df$cell_type)[5] <- "Inh"
merged_df$broad <- as.factor(merged_df$cell_type)
levels(merged_df$broad) <- c("non-neuronal", "non-neuronal", "neuronal", "non-neuronal", "neuronal", "non-neuronal", "non-neuronal")

common_limits <- c(-0.6, 1)
common_breaks <- seq(-0.6, 1, by = 0.2)

colors <- RColorBrewer::brewer.pal(7, "Dark2")

# Scatterplot AD vs MDD
ggplot(merged_df, aes(x = logFC_df, y = logFC_df2, color = cell_type, label = gene)) +
  geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0), fill = "gray99", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0), fill = "gray95", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), fill = "gray95", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf), fill = "gray90", alpha = 0.3, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.4) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(max.overlaps = 15, size = 3, fontface = "italic") +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(limits = common_limits, breaks = common_breaks, labels = scales::label_number(accuracy = 0.1)) +
  scale_y_continuous(limits = common_limits, breaks = common_breaks, labels = scales::label_number(accuracy = 0.1)) +
  labs(x = "logFC (AD)", y = "logFC (MDD)", color = "Cell Type") +
  coord_fixed() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  facet_wrap(~broad)

# Blood protein heatmap
MDD_Combined <- read_csv("/data/SciencePaper/UKBBData/MDD_CombinedArrayI.csv")
bloodUKBB <- rbind.data.frame(cbind(def=MDD_Combined$DepressionDef,
                                    prot=MDD_Combined$UKBPPP_ProteinID_Panel,
                                    p=MDD_Combined$P,
                                    logFC=MDD_Combined$BETA,
                                    array=MDD_Combined$Panel,
                                    gene=MDD_Combined$gene))

blood_de_proteins <- bloodUKBB %>%
  group_by(def) %>%
  mutate(adj.P.Val = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  filter(grepl("^EPHA|^EPHB|^EFN", gene))

heatmap_mat <- blood_de_proteins %>%
  dplyr::select(gene, def, logFC) %>%
  mutate(logFC = as.numeric(logFC)) %>%
  pivot_wider(names_from = def, values_from = logFC) %>%
  column_to_rownames("gene") %>%
  as.matrix()
heatmap_mat[is.na(heatmap_mat)] <- 0

star_mat <- blood_de_proteins %>%
  mutate(star = ifelse(adj.P.Val < 0.05, "*", "")) %>%
  dplyr::select(gene, def, star) %>%
  pivot_wider(names_from = def, values_from = star, values_fill = "") %>%
  column_to_rownames("gene") %>%
  as.matrix()

col_fun <- colorRamp2(c(0,0.25), c("white",'#d12620'))
DepDx_colors <- c("a" = "#713dc9", "b" = "#c9713d", "c" = "#3dc971")

p <- Heatmap(
  heatmap_mat,
  name = "logFC",
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(star_mat[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))
  },
  row_dend_width = unit(0.25, "cm"),
  column_dend_height = unit(0.25, "cm"),
  top_annotation = columnAnnotation(
    DepDx = c("c", "c", "b", "b", "c", "a", "a"),
    col = list(DepDx = DepDx_colors),
    simple_anno_size = unit(1, "mm"),
    annotation_name_gp = grid::gpar(fontsize = 9),
    annotation_legend_param = list(title = "DepDx")
  ),
  cluster_column_slices = TRUE,
  column_split = c("c", "c", "b", "b", "c", "a", "a"),
  column_title = NULL,
  width = unit(7, "cm"),
  height = unit(7, "cm"),
  heatmap_legend_param = list(title = "logFC",
                              legend_height = unit(4, "cm"),
                              title_gp = grid::gpar(fontsize = 10),
                              labels_gp = grid::gpar(fontsize = 9),
                              direction = "horizontal")
)

pdf("/data/BBRF/Finalized_output/Blood/EPH_heatmap.pdf")
draw(p, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legen = TRUE)
dev.off()

# Heatmap of MDD single-cell logFC
celltypes <- names(MDD1_data_list_cc)
sig_gene_lists <- lapply(MDD1_data_list_cc, function(df) {
  df %>% filter(P.Value < 0.05) %>% pull(genes)
})
sig_genes_all <- unique(unlist(sig_gene_lists))

logfc_mat <- matrix(0, nrow = length(sig_genes_all), ncol = length(celltypes),
                    dimnames = list(sig_genes_all, celltypes))
for (ct in celltypes) {
  df <- MDD1_data_list_cc[[ct]]
  df_sig <- df %>% filter(genes %in% sig_genes_all)
  logfc_mat[df_sig$genes, ct] <- df_sig$logFC
}

row_order <- hclust(dist(logfc_mat))$order

Heatmap(logfc_mat[row_order, ],
        name = "logFC",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cluster_rows = FALSE, cluster_columns = TRUE,
        show_row_names = TRUE, show_column_names = TRUE)
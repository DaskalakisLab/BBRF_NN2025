###################################################################
### Cortical Subregions MOFA Analysis - Based on ARA | December ###
###################################################################

### BEFORE YOU START: RAW OMIC DATA HANDLING ###
# Refer to: https://github.com/bioFAM/MOFA2/blob/master/vignettes/getting_started_R.Rmd

# Set Working Directory
setwd('/data/BBRF/Bulk/MOFA')

# Load Libraries
library(bigreadr)
library(tidyr)
library(stringr)
library(reticulate)
library(MOFA2)
library(limma)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


# Ensure Output Directories Exist
output_dirs <- c(
  "./MOFA_INPUT", 
  "./MOFA_OUTPUT/run5", 
  "./MOFA_FIGURES/run5"
)

for (dir in output_dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

write.table("Multi-omics, 0.6 var, 30 Factors, convergence slow", "./run5_specs.txt")

########################### 2.1 - PREPARE DATA FOR MOFA ############################

# Load Annotation Data
anno_mpfc <- readRDS("/data/BBRF/Bulk/mPFC/CombinedResults_svs/all/annotation_gDxMDD_All.RDS")
anno_dlpfc <- readRDS("/data/BBRF/Bulk/DLPFC/CombinedResults_svs/all/annotation_gDxMDD_All.RDS")
anno_dacc <- readRDS("/data/BBRF/Bulk/dACC/CombinedResults_svs/all/annotation_gDxMDD_All.RDS")

# Identify Common Samples
common_samples <- Reduce(intersect, list(anno_mpfc$BrNum, anno_dlpfc$BrNum, anno_dacc$BrNum))

# Subset and Align Annotations
anno_mpfc <- anno_mpfc[match(common_samples, anno_mpfc$BrNum), ]
anno_dlpfc <- anno_dlpfc[match(common_samples, anno_dlpfc$BrNum), ]
anno_dacc <- anno_dacc[match(common_samples, anno_dacc$BrNum), ]

# Load Expression Data
exp_mpfc <- readRDS("/data/BBRF/Bulk/mPFC/CombinedResults_svs/all/expression_gDxMDD_All.RDS")
exp_dlpfc <- readRDS("/data/BBRF/Bulk/DLPFC/CombinedResults_svs/all/expression_gDxMDD_All.RDS")
exp_dacc <- readRDS("/data/BBRF/Bulk/dACC/CombinedResults_svs/all/expression_gDxMDD_All.RDS")

# Align Expression Data with Annotations
exp_mpfc <- exp_mpfc[, match(anno_mpfc$SampleID, colnames(exp_mpfc))]
exp_dlpfc <- exp_dlpfc[, match(anno_dlpfc$SampleID, colnames(exp_dlpfc))]
exp_dacc <- exp_dacc[, match(anno_dacc$SampleID, colnames(exp_dacc))]

colnames(exp_mpfc) <- anno_mpfc$BrNum
colnames(exp_dlpfc) <- anno_dlpfc$BrNum
colnames(exp_dacc) <- anno_dacc$BrNum

stopifnot(identical(colnames(exp_mpfc), colnames(exp_dacc)))

vars <- apply(exp_mpfc, 1, var)
exp_mpfc_F <- exp_mpfc[vars >= quantile(vars,0.6),]
dim(exp_mpfc_F)
#[1] 8974  138

vars <- apply(exp_dlpfc, 1, var)
exp_dlpfc_F <- exp_dlpfc[vars >= quantile(vars,0.6),]
dim(exp_dlpfc_F)
#[1] 9039  138

vars <- apply(exp_dacc, 1, var)
exp_dacc_F <- exp_dacc[vars >= quantile(vars,0.6),]
dim(exp_dacc_F)
#[1] 9164  138

####Add more omics from mPFC####
meth_mpfc <- readRDS("/data/MOFA/MOFA_INPUT/meth_mpfc.rds")
protein_mpfc <- readRDS("/data/MOFA/MOFA_INPUT/protein_mpfc.rds")

meth_mpfc <- meth_mpfc[, match(anno_mpfc$BrNum, colnames(meth_mpfc))]
protein_mpfc <- protein_mpfc[, match(anno_mpfc$BrNum, colnames(protein_mpfc))]

################################### 2.3 - PREPARE INPUT OBJECT ###################################

# Combine Expression Data into a List
input.list <- list(exp_mpfc_F, exp_dlpfc_F, exp_dacc_F, meth_mpfc, protein_mpfc)

# Save Input List
saveRDS(input.list, "./MOFA_INPUT/mpfc_multiomics_run5_input.rds")

# Convert Matrices to Numeric
input.list <- lapply(input.list, as.matrix)

################################ 3 - CREATE AND PREPARE MOFA OBJECT ################################

# Create MOFA Object
mofa.obj <- create_mofa_from_matrix(input.list)
views_names(mofa.obj) <- c("mPFC_rna", "dlPFC_rna", "dACC_rna", "mPFC_meth", "mPFC_protein")

# Save Data Overview Plot
pdf("./MOFA_FIGURES/run5/mpfc_omics_input.pdf", width = 5, height = 5)
plot_data_overview(mofa.obj, colors = c("orange2", "purple4", "darkgreen", "lightblue", "darkred"))
dev.off()

# Set Metadata Options
data_opts <- get_default_data_options(mofa.obj)

# Set Model Options
model_opts <- get_default_model_options(mofa.obj)
model_opts$num_factors <- 30

# Set Training Options
train_opts <- get_default_training_options(mofa.obj)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 1234

# Prepare MOFA Model
mofa.obj <- prepare_mofa(mofa.obj, data_opts, model_opts, train_opts)

# Save MOFA Object
saveRDS(mofa.obj, "./MOFA_OUTPUT/run5/MofaObj_1.rds")

################################ 4 - TRAIN MOFA MODEL ################################
library(reticulate)

# Train and Save Model
mofa.train <- run_mofa(mofa.obj, outfile = "./MOFA_OUTPUT/run5/subregion.hdf5", use_basilisk = T)

################################ 5 - EXTRACT AND SAVE RESULTS ################################

# Load Trained Model
mofa.train <- load_model("./MOFA_OUTPUT/run5/subregion.hdf5", remove_inactive_factors = FALSE)

# Save and Process Weights
weights <- get_weights(mofa.train, as.data.frame = TRUE)
saveRDS(weights, "./MOFA_OUTPUT/run5/subregion_WEIGHTS.rds")

weights <- reshape(weights, idvar = "feature", timevar = "factor", direction = "wide")

gene_map <- readRDS("/data/humgen/daskalakislab/dipietro/SciencePaper/Data/GeneMappings/RNAGeneMap.RDS")

processed_datasets <- list()

all_views <- unique(c(
  unique(weights$view.Factor1),
  unique(weights$view.Factor2),
  unique(weights$view.Factor3),
  unique(weights$view.Factor4),
  unique(weights$view.Factor5)
))


# Loop through each unique view
for (view in all_views) {
  # Filter rows for the current view across all Factor columns
  filtered_data <- weights[
    weights$view.Factor1 == view |
      weights$view.Factor2 == view |
      weights$view.Factor3 == view |
      weights$view.Factor4 == view |
      weights$view.Factor5 == view, 
  ]
  
  # Remove any column that starts with 'view.'
  filtered_data <- filtered_data[, !grepl("^view\\.", colnames(filtered_data))]
  
  # Replace 'value.' in column names with an empty string
  colnames(filtered_data) <- gsub("^value\\.", "", colnames(filtered_data))
  
  filtered_data$ensg <- sapply(filtered_data$feature, function(x) str_split(x, "_")[[1]][1])
  filtered_data <- left_join(filtered_data, gene_map, by = c("ensg" = "genes"))
  
  saveRDS(filtered_data, paste0("./MOFA_OUTPUT/run5/weights_", as.character(view), ".rds"))

  processed_datasets[[as.character(view)]] <- filtered_data
}

r2_list <- mofa.train@cache$variance_explained$r2_per_factor$group1
saveRDS(r2_list, "./MOFA_OUTPUT/run5/mofa_var_explained_r2_per_factor.rds")

# Save Variance Explained
varexp <- calculate_variance_explained(mofa.train)
saveRDS(varexp, "./MOFA_OUTPUT/run5/subregion_VAREXP.rds")

# Save Variance Explained Plot
pdf("./MOFA_FIGURES/run5/varexp1.pdf", width = 5, height = 5)
plot_variance_explained(mofa.train, plot_total = TRUE)[[2]]
dev.off()

source("/data/BBRF/Bulk/MOFA/Scripts/plot_variance_explained_sqrt.R")

pdf('./MOFA_FIGURES/run5/varexp2.pdf', width = 5, height = 5, onefile = T)
plot_variance_explained(mofa.train, max_r2=100)
plot_variance_explained(mofa.train, max_r2=20)
plot_variance_explained(mofa.train, max_r2=10)
plot_variance_explained_sqrt(mofa.train, max_r2=5, fill = c("sqrt_value"),legend.title = expression(sqrt(R^2)))
dev.off()



# Save Factor Correlation Plot
pdf("./MOFA_FIGURES/run5/factor_cor.pdf", width = 5, height = 5)
plot_factor_cor(mofa.train)
dev.off()

######################### 6 - EXPLORE RELATIONSHIP WITH DEMO DATA ################################
factors <- get_factors(mofa.train, as.data.frame = TRUE)
saveRDS(factors, './MOFA_OUTPUT/run5/subregion_FACTORS.rds')

factors = reshape(factors, idvar = 'sample', timevar = 'factor',direction = "wide")

# Remove "group.Factor" columns
factors <- factors[, !grepl("group\\.Factor", colnames(factors))]

# Remove the "value." prefix from column names
colnames(factors) <- gsub("^value\\.", "", colnames(factors))

# merge with anno
anno <- anno_mpfc
anno <- merge(x= anno, y = factors, by.x = "BrNum", by.y = "sample")

wtr <- readRDS("/data/BBRF/Bulk/MethAge/watermelon_age.rds")
anno <- merge(x= anno, y = wtr, by = "BrNum")


cor(anno$Age, anno$horvath)
cor(anno$Age, anno$hannum)
cor(anno$Age, anno$phenoage)

anno <- anno %>%
  mutate(female = case_when(
    Sex == "F" ~ 1,
    Sex == "M" ~ 0
  )) %>%
  mutate(psychosis = case_when(
    Psychosis == "Yes" ~ 1,
    Psychosis == "No" ~0
  )) %>%
  mutate(smoking = case_when(
    Smoking == "Yes" ~ 1,
    Smoking == "No" ~0
  ))

#check these variables
variables <- c("PMI", "Age","horvath", "BMI", "ancestryPC1", "ancestryPC2", "female", "psychosis", "smoking")

# Select only the numeric factor columns for correlation
factor_columns <- grep("^Factor", colnames(anno), value = TRUE)

# Initialize a matrix to store Spearman correlations
results <- matrix(NA, nrow = length(variables), ncol = length(factor_columns),
                  dimnames = list(variables, factor_columns))

# Loop over each variable and factor, calculate Spearman correlation
for (var in variables) {
  for (factor in factor_columns) {
    # Calculate Spearman correlation and store
      results[var, factor] <- cor(anno[[var]], anno[[factor]], method = "spearman")
  }
}

# Convert results to a data frame for easier viewing
results_df <- as.data.frame(results)
results_df <- tibble::rownames_to_column(results_df, "Variable")

#Visualize
# Required libraries
library(ggplot2)
library(reshape2)

# Convert the results to a data frame for plotting
results_df <- as.data.frame(results)
results_df$Variable <- rownames(results_df)
results_long <- melt(results_df, id.vars = "Variable", variable.name = "Factor", value.name = "Correlation")

##cluster if desired
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

library(tibble, lib.loc = "/opt/R/4.2.0/lib/R/library")
# Compute distance matrix and hierarchical clustering
results_wide <- results_long %>% 
  dplyr::select(Variable, Factor, Correlation) %>% 
  pivot_wider(names_from = Factor, values_from = Correlation) %>% 
  column_to_rownames(var = "Variable") 


row_order <- hclust(dist(results_wide))$order

# Reorder `Variable` based on clustering
results_long$Variable <- factor(results_long$Variable, 
                                levels = rownames(results_wide)[row_order])

results_long$Correlation <- abs(results_long$Correlation)
# Plot heatmap with clustered rows
heatmap_plot <- ggplot(results_long, aes(x = Factor, y = Variable, fill = Correlation)) +
  geom_tile(color = "grey50") +
  geom_text(aes(label = round(Correlation, 2)), size = 4, color = "black") +
  # scale_fill_gradient2(low = "#3d95c9", mid = "white", high = '#d12620', midpoint = 0, 
  #                      limits = c(-1, 1), name = "Spearman\nCorrelation") +
  scale_fill_gradient2(low = "white", high = '#d12620', midpoint = 0, 
                       limits = c(0, 1), name = "Spearman\nCorrelation") +
  labs(title = "", x = "", y = "") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())


ggsave(filename = "/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/abs_correlation_numeric_samplespec.pdf", w = 15, h = 8)

#Association with diagnosis

factors_long <- melt(factors, id.vars = "sample", variable.name = "Factor", value.name = "Value")

# Merge with gDx information
plot_data <- merge(factors_long, anno[, c("BrNum", "gDx")], by.x = "sample", by.y = "BrNum")

# Create the boxplot with points
boxplot_plot <- ggplot(plot_data, aes(x = gDx, y = Value, fill = gDx)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot with no outlier points
  geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.5, color = "black", size=0.2) + # Add individual points
  facet_wrap(~ Factor, scales = "free_y", ncol = 5) + # Create one boxplot per factor
  scale_fill_brewer(palette = "Dark2") + # Adjust colors for groups
  labs(title = "",
       x = "",
       y = "Factor Values") +
  theme_minimal() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(filename = "/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/factors_gDx.pdf", w = 10, h = 10)

#Association with sex

factors_long <- melt(factors, id.vars = "sample", variable.name = "Factor", value.name = "Value")

# Merge with sex information
plot_data <- merge(factors_long, anno[, c("BrNum", "Sex")], by.x = "sample", by.y = "BrNum")

# Create the boxplot with points
boxplot_plot <- ggplot(plot_data, aes(x = Sex, y = Value, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot with no outlier points
  geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.5, color = "black", size=0.2) + # Add individual points
  facet_wrap(~ Factor, scales = "free_y", ncol = 5) + # Create one boxplot per factor
  scale_fill_brewer(palette = "Set1") + # Adjust colors for groups
  labs(title = "",
       x = "",
       y = "Factor Values") +
  theme_minimal() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(filename = "/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/factors_Sex.pdf", w = 10, h = 10)

###
plot_data <- merge(factors_long, anno[, c("BrNum", "Batch")], by.x = "sample", by.y = "BrNum")

# Create the boxplot with points
boxplot_plot <- ggplot(plot_data, aes(x = Batch, y = Value, fill = Batch)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot with no outlier points
  geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.5, color = "black", size=0.2) + # Add individual points
  facet_wrap(~ Factor, scales = "free_y", ncol = 5) + # Create one boxplot per factor
  scale_fill_brewer(palette = "Set1") + # Adjust colors for groups
  labs(title = "",
       x = "",
       y = "Factor Values") +
  theme_minimal() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(filename = "/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/factors_Batch.pdf", w = 10, h = 10)


###Science like graphics####

mofa_heatmap<- as.matrix(r2_list)
col_info <- as.factor(sub("^.*?_", "", colnames(mofa_heatmap)))
levels(col_info) <- c("Methylation", "Protein", "RNA")
col_info <- factor(col_info, levels= c("RNA", "Protein", "Methylation"), )
capped_vector <- pmin(mofa_heatmap, 10)

modified_vector <- gsub("rna", "RNA", colnames(capped_vector))
modified_vector <- gsub("protein", "Protein", modified_vector)
modified_vector <- gsub("meth", "Methyl", modified_vector)
modified_vector -> colnames(capped_vector)

rowAnn <- rowAnnotation(omic = col_info,
                        col = list(omic = c("RNA" = "#94A46D", "Methylation"= "#5C667C", "Protein"= "#B07677")),
                        simple_anno_size = unit(1.5, "mm"),
                        show_annotation_name = TRUE,
                        annotation_label = NULL)

col_fun = colorRamp2(c(0, 5, 10), c("white", "#87c4ff", "#003B73"))
# "#FFF7F3" "#FDE0DD" "#FCC5C0" "#FA9FB5" "#F768A1" "#DD3497" "#AE017E" "#7A0177" "#49006A"
pdf("/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/var_explained_ch.pdf",w = 8.5, h = 3)
Heatmap(t(capped_vector), name = "Variance \nExplained (%)", cluster_rows = F, cluster_columns = F,
        col = col_fun,
        width = nrow(mofa_heatmap)*unit(5, "mm"), 
        height = ncol(mofa_heatmap)*unit(5, "mm"), 
        column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 10), 
        show_row_names = T,
        row_split = col_info,
        row_title = " ",
        rect_gp = gpar(col = "grey60", lwd = 1),
        right_annotation = rowAnn)
dev.off()

####Factor13xAge####
subregion_FACTORS <- readRDS("/data/BBRF/Bulk/MOFA/MOFA_OUTPUT/run5/subregion_FACTORS.rds")
library(tidyr)
subregion_FACTORS$value <- -subregion_FACTORS$value

# Widen the dataframe
widened_df <- subregion_FACTORS %>%
  pivot_wider(
    names_from = factor, # Column to use for new column names
    values_from = value  # Column to fill with values
  )

anno_ext <- merge(anno_mpfc, widened_df, by.x= "BrNum", by.y = "sample")
anno_ext <- merge(anno_ext, wtr, by = "BrNum")
anno_ext$Delta_Age <- anno_ext$Age - anno_ext$horvath

library(ggplot2)
library(grid)

# Scatterplot with regression lines

cor_values <- anno_ext %>%
  group_by(gDx) %>%
  summarise(correlation = cor(Age, Factor13, method = "spearman"))

cor_values_D <- anno_ext %>%
  group_by(gDx) %>%
  summarise(correlation = cor(Delta_Age, Factor13, method = "spearman"))

cor_values_ho <- anno_ext %>%
  group_by(gDx) %>%
  summarise(correlation = cor(horvath, Factor13, method = "spearman"))

cor_values_ha <- anno_ext %>%
  group_by(gDx) %>%
  summarise(correlation = cor(hannum, Factor13, method = "spearman"))

cor_values_pa <- anno_ext %>%
  group_by(gDx) %>%
  summarise(correlation = cor(phenoage, Factor13, method = "spearman"))


ggplot(anno_ext, aes(x = Age, y = Factor13, color = gDx)) +
  geom_point(size = 2, alpha = 0.7) +                     
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +    
  labs(
    x = "Age at Death",
    y = "Factor 13 score",
    color = "Dx",
    title = ""
  ) +
  scale_color_manual(values= c("#82d4a8", "#e86878"))+
  theme_minimal(base_size = 14) +                        
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  )+
  annotation_custom(
    grob = textGrob(expression(R[CTRL] == 0.87), gp = gpar(fontsize = 15, fontface = "bold")),
    xmin = 55, xmax = 0, ymin = 0.5, ymax = Inf
  ) +
  annotation_custom(
    grob = textGrob(expression(R[MDD] == 0.88), gp = gpar(fontsize = 15, fontface = "bold")),
    xmin = 55, xmax = 0, ymin = 0.5, ymax = 0.6
  )
ggsave(filename = "/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/abs_Age_Factor_scatter.pdf", device = "pdf", w = 6, h = 6)

ggplot(anno_ext, aes(x = horvath, y = Factor13, color = gDx)) +
  geom_point(size = 2, alpha = 0.7) +                     
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +    
  labs(
    x = "Biological Age",
    y = "Factor 13 score",
    color = "Dx",
    title = ""
  ) +
  scale_color_manual(values= c("#82d4a8", "#e86878"))+
  theme_minimal(base_size = 14) +                        
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  )+
  annotation_custom(
    grob = textGrob(expression(R[CTRL] == 0.91), gp = gpar(fontsize = 15, fontface = "bold")),
    xmin = 55, xmax = 0, ymin = 0.5, ymax = Inf
  ) +
  annotation_custom(
    grob = textGrob(expression(R[MDD] == 0.84), gp = gpar(fontsize = 15, fontface = "bold")),
    xmin = 55, xmax = 0, ymin = 0.5, ymax = 0.6
  )
ggsave(filename = "/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/abs_mAge_Factor_scatter.pdf", device = "pdf", w = 7, h = 7)


ggplot(anno_ext, aes(x = Age, y = Factor13, color = horvath, shape = gDx)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = gDx), color = "black") + 
  scale_color_viridis_c() + 
  labs(x = "Chronological Age", y = "Factor 13 Score", color = "Biological Age", linetype = "Diagnosis")+
  theme_minimal(base_size = 14) +                        
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  )+
  annotation_custom(
    grob = textGrob(expression(R[CTRL] == 0.87), gp = gpar(fontsize = 15, fontface = "bold")),
    xmin = 55, xmax = 0, ymin = 0.5, ymax = Inf
  ) +
  annotation_custom(
    grob = textGrob(expression(R[MDD] == 0.88), gp = gpar(fontsize = 15, fontface = "bold")),
    xmin = 55, xmax = 0, ymin = 0.5, ymax = 0.6
  )
ggsave(filename = "/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/abs_Age_mAge_Factor_scatter.pdf", device = "pdf", w = 6, h = 6)

##Age as factor
summary(anno_ext$Age)
anno_ext$Age_Group <- cut(anno_ext$Age, 
                          breaks = quantile(anno_ext$Age, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                          labels = c("Young", "Middle-Aged", "Older"), 
                          include.lowest = TRUE)
table(anno_ext$Age_Group)

ggplot(anno_ext, aes(x = Age_Group, y = Factor13, fill = gDx)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.4, position = position_dodge(width = 0.6)) +        
  geom_jitter(aes(color = horvath), size = 2, alpha = 0.7, position = position_dodge(width = 0.6)) +     
  labs(
    x = "Age Group",
    y = "Factor 4 score",
    color = "Biological Age",
    fill = "Diagnosis",
    title = ""
  ) +
  scale_color_viridis_c() + 
  scale_fill_manual(values= c("#82d4a8", "#e86878"))+
  theme_minimal(base_size = 14) +                        
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  )

anno_ext %>%
  group_by(Age_Group) %>%
  summarise(p_value = wilcox.test(Factor13 ~ gDx, data = pick())$p.value)

library(dplyr)

anno_ext %>%
  group_by(Age_Group) %>%
  summarise(p_value = t.test(Factor13 ~ gDx, data = cur_data())$p.value)

anno_ext %>%
  group_by(Age_Group) %>%
  summarise(p_value = wilcox.test(Factor13 ~ gDx, data = cur_data())$p.value)
####
ggplot(anno_ext, aes(x = Delta_Age, y = Factor13, shape = gDx)) +
  geom_point(aes(color = Age, size = 3, alpha = 0.7)) +
  geom_smooth(method = "lm", se = F, aes(linetype = gDx), color = "black") +
  scale_color_viridis_c() +
  labs(x = "Delta Age (Chronological - Biological)", y = "Factor 13 Score", color = "Chronological Age")+
  theme_minimal(base_size = 14) +                        
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  )+
  annotation_custom(
    grob = textGrob(expression(R[CTRL] == -0.35), gp = gpar(fontsize = 15, fontface = "bold")),
    xmin = 28, xmax = 0, ymin = 0.55, ymax = Inf
  ) +
  annotation_custom(
    grob = textGrob(expression(R[MDD] == -0.56), gp = gpar(fontsize = 15, fontface = "bold")),
    xmin = 28, xmax = 0, ymin = 0.55, ymax = 0.5
  )
ggsave(filename = "/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/abs_dAge_Factor_scatter.pdf", device = "pdf", w = 8, h = 8)

lm_model <- lm(Factor13 ~ Age : gDx, data = anno_ext)
summary(lm_model)

lm_model <- lm(Factor13 ~ horvath : gDx, data = anno_ext)
summary(lm_model)

lm_model <- lm(Factor13 ~ Delta_Age : gDx, data = anno_ext)
summary(lm_model)

lm_model <- lm(Factor13 ~ Delta_Age * gDx, data = anno_ext)
summary(lm_model)

library(smatr)
slope_test <- sma(Factor13 ~ horvath + gDx, data = anno_ext)
summary(slope_test)
####Factor4xDx####

# Widen the dataframe

# Boxplot with scatter

ggplot(anno_ext, aes(x = gDx, y = Factor4)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.2) +        
  geom_jitter(size = 2, alpha = 0.7, width = 0.2, aes(color = gDx)) +    
  labs(
    x = "Dx",
    y = "Factor 4 score",
    color = "",
    title = ""
  ) +
  scale_color_manual(values= c("#82d4a8", "#e86878"))+
  theme_minimal(base_size = 14) +                        
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  )
ggsave(filename = "/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/absDx_Factor_scatter.pdf", device = "pdf", w = 6, h = 6)

lm_model <- lm(Factor4 ~ gDx, data = anno_ext)
summary(lm_model)

t.test(Factor4 ~ gDx, data = anno_ext)

####Factor28xSex####

# Widen the dataframe

# Boxplot with scatter
anno_ext$Sex <- as.factor(anno_ext$Sex)
levels(anno_ext$Sex) <- c("Female", "Male")

ggplot(anno_ext, aes(x = Sex, y = Factor28)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.2) +        
  geom_jitter(size = 2, alpha = 0.7, width = 0.2, aes(color = gDx)) +    
  labs(
    x = "Sex",
    y = "Factor 28 score",
    color = "",
    title = ""
  ) +
  scale_color_manual(values= c("#82d4a8", "#e86878"))+
  theme_minimal(base_size = 14) +                        
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  )
ggsave(filename = "/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/absSex_Factor_scatter.pdf", device = "pdf", w = 6, h = 6)

lm_model <- lm(Factor28 ~ Sex, data = anno_ext)
summary(lm_model)

t.test(Factor28 ~ Sex, data = anno_ext)
#t = 3.0392, df = 135.38, p-value = 0.002847

lm_model <- lm(Factor28 ~ Sex*gDx, data = anno_ext)
summary(lm_model)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)    -0.01858    0.05269  -0.353    0.725
# SexMale        -0.01396    0.06418  -0.218    0.828
# gDxMDD          0.09832    0.06021   1.633    0.105
# SexMale:gDxMDD -0.12657    0.07705  -1.643    0.103
# 
# Residual standard error: 0.2041 on 134 degrees of freedom
# Multiple R-squared:  0.08352,	Adjusted R-squared:  0.063 
# F-statistic:  4.07 on 3 and 134 DF,  p-value: 0.008369

ggplot(anno_ext, aes(x = gDx, y = Factor28, fill = Sex)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values= c("#EBB8DD", "#43A5BE"))+
  theme_minimal(base_size = 14)+
  labs(
    x = "Dx",
    y = "Factor 28 score",
    color = "",
    title = ""
  ) 
ggsave(filename = "/data/BBRF/Bulk/MOFA/MOFA_FIGURES/run5/absSexDx_Factor_scatter.pdf", device = "pdf", w = 6, h = 6)

###GO####
library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(MOFAdata)
library(MOFA2)

dir.create("/data/BBRF/Bulk/MOFA/MOFA_OUTPUT/run5/GSEA", recursive = T)

factors <- paste0("Factor", 1:30)

regions <- c("dlPFC", "mPFC", "dACC")
pathways2 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
fgsea_sets <- pathways2 %>% split(x = .$gene_symbol, f = .$gs_name)

library(tibble)


for(region in regions){
  weights <- readRDS(paste0("/data/BBRF/Bulk/MOFA/MOFA_OUTPUT/run5/weights_", region, "_rna.rds"))
  
  for (factor in factors) {
    if (!factor %in% colnames(weights)) next
    
    results <- weights[,c(factor, "symbol")]
    
    results[,factor] <- as.numeric(results[,factor] )
    
    gsea_input <- results %>%
      arrange(desc(.data[[factor]])) %>%  # Use .data[[factor]] for dynamic column reference
      dplyr::select(symbol, !!sym(factor))
    
    ranks <- deframe(gsea_input)
    
    # Run GSEA
    fgseaRes <- fgsea::fgseaMultilevel(fgsea_sets, stats = ranks, 
                                       minSize = 10, maxSize = 500, eps = 0)
    
    # Prepare data to save and plot
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    fgseaResTidier <- fgseaResTidy %>% 
      dplyr::select(-leadingEdge, -ES, -log2err) %>% 
      arrange(padj) %>%
      mutate(lp = -log10(pval))
    
    
    saveRDS(fgseaResTidier, paste0("/data/BBRF/Bulk/MOFA/MOFA_OUTPUT/run5/GSEA/C5_BP_", factor, "_", region, ".rds"))
    save(fgseaResTidy, file = paste0("/data/BBRF/Bulk/MOFA/MOFA_OUTPUT/run5/GSEA/C5_BP_", factor, "_", region, ".rdata"))
  }}

  
  

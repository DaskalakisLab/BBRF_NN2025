#cellchat@net[["Control"]][["count"]]
ctrl_count <- cellchat@net[["Control"]][["count"]]
mdd_count <- cellchat@net[["MDD"]][["count"]]

chisq.test(matrix(c(ctrl_count, mdd_count), nrow = 2))


set.seed(123)

# Extract count matrices
mat_control <- cellchat@net[["Control"]][["count"]]
mat_mdd     <- cellchat@net[["MDD"]][["count"]]

# Observed difference in total communication
observed_diff <- sum(mat_mdd) - sum(mat_control)

# Combine into a single vector
combined <- c(as.vector(mat_control), as.vector(mat_mdd))

# Permutation test
n_perm <- 10000
perm_diffs <- numeric(n_perm)

for (i in 1:n_perm) {
  permuted <- sample(combined)
  perm_control <- matrix(permuted[1:length(mat_control)], nrow = nrow(mat_control))
  perm_mdd     <- matrix(permuted[(length(mat_control)+1):length(permuted)], nrow = nrow(mat_control))
  
  perm_diffs[i] <- sum(perm_mdd) - sum(perm_control)
}

# p-value
p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
hist(perm_diffs, breaks = 50, main = "Permutation Test: Communication Sum Difference")
abline(v = observed_diff, col = "red", lwd = 2)


#########
# Sum of outgoing signals (rows) and incoming (columns)
out_control <- rowSums(mat_control)
in_control  <- colSums(mat_control)

out_mdd <- rowSums(mat_mdd)
in_mdd  <- colSums(mat_mdd)

# Focus on immune cells
immune_cells <- c("Micro", "Macrophages", "Tcells")

data_summary <- data.frame(
  celltype = rep(immune_cells, each = 2),
  condition = rep(c("Control", "MDD"), times = length(immune_cells)),
  outgoing = c(out_control[immune_cells], out_mdd[immune_cells]),
  incoming = c(in_control[immune_cells], in_mdd[immune_cells])
)


library(ggplot2)

ggplot(data_summary, aes(x = celltype, y = outgoing, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Total Outgoing Signaling from Immune Cells", y = "Count", x = "") +
  theme_minimal()

ggplot(data_summary, aes(x = celltype, y = incoming, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Total Incoming Signaling to Immune Cells", y = "Count", x = "") +
  theme_minimal()


# Define immune-involved interactions
is_immune <- function(source, target) {
  source %in% immune_cells || target %in% immune_cells
}

# Flatten and tag interactions
get_tagged_df <- function(mat, condition) {
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("source", "target", "count")
  df$condition <- condition
  df$immune <- mapply(is_immune, df$source, df$target)
  return(df)
}

df_all <- rbind(get_tagged_df(mat_control, "Control"),
                get_tagged_df(mat_mdd, "MDD"))

# Collapse to immune vs. non-immune totals
summary_table <- df_all %>%
  group_by(condition, immune) %>%
  summarise(total = sum(count)) %>%
  tidyr::pivot_wider(names_from = condition, values_from = total) %>%
  column_to_rownames("immune")

# Perform test
fisher.test(as.matrix(summary_table))

library(CellChat)
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Micro")


set.seed(42)

# Step 1: Define immune-involved interactions
immune_cells <- c("Micro", "Macrophages", "Tcells")

is_immune <- function(source, target) {
  source %in% immune_cells || target %in% immune_cells
}

# Step 2: Flatten count matrices into long form
flatten_matrix <- function(mat, condition) {
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("source", "target", "count")
  df$condition <- condition
  df$immune <- mapply(is_immune, df$source, df$target)
  return(df)
}

df_control <- flatten_matrix(mat_control, "Control")
df_mdd     <- flatten_matrix(mat_mdd, "MDD")
df_all <- rbind(df_control, df_mdd)

# Step 3: Observed difference in immune communication
observed_diff <- sum(df_mdd$count[df_mdd$immune]) - sum(df_control$count[df_control$immune])

# Step 4: Permutation test
n_perm <- 10000
perm_diffs <- numeric(n_perm)

# Permute condition labels
for (i in 1:n_perm) {
  permuted_labels <- sample(df_all$condition)
  df_perm <- df_all
  df_perm$perm_condition <- permuted_labels
  
  sum_mdd <- sum(df_perm$count[df_perm$perm_condition == "MDD" & df_perm$immune])
  sum_ctrl <- sum(df_perm$count[df_perm$perm_condition == "Control" & df_perm$immune])
  
  perm_diffs[i] <- sum_mdd - sum_ctrl
}

# Step 5: Calculate p-value and plot
p_value <- mean(abs(perm_diffs) >= abs(observed_diff))

hist(perm_diffs, breaks = 50, main = "Permutation Test: Immune Communication Difference",
     xlab = "Difference (MDD - Control)", col = "skyblue")
abline(v = observed_diff, col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste0("Observed = ", round(observed_diff, 1), 
                                   "\nP = ", signif(p_value, 3)),
       bty = "n")

cat("Permutation test p-value:", p_value, "\n")


######
cell_groups <- list(
  Neurons = c("Ex", "Inh"),
  Glia = c("Oligo", "OPC", "Astro"),
  Immune = c("Micro", "Macrophages", "Tcells"),
  Perivascular = c("Endo", "Mural", "Fibro")
)

set.seed(42)

# Flatten matrix
flatten_matrix <- function(mat, condition) {
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("source", "target", "count")
  df$condition <- condition
  return(df)
}

df_control <- flatten_matrix(mat_control, "Control")
df_mdd     <- flatten_matrix(mat_mdd, "MDD")
df_all <- rbind(df_control, df_mdd)

# Assign group tags
assign_group <- function(cell, group_map) {
  for (grp in names(group_map)) {
    if (cell %in% group_map[[grp]]) return(grp)
  }
  return("Other")
}

df_all$source_group <- sapply(df_all$source, assign_group, group_map = cell_groups)
df_all$target_group <- sapply(df_all$target, assign_group, group_map = cell_groups)

# Mark if the interaction involves a group
group_names <- names(cell_groups)

# Create a long-format list of group comparisons
results <- list()

for (grp in group_names) {
  # Label interactions involving this group
  df_all$group_involved <- df_all$source_group == grp | df_all$target_group == grp
  
  # Observed difference
  obs_diff <- with(df_all, sum(count[condition == "MDD" & group_involved]) -
                     sum(count[condition == "Control" & group_involved]))
  
  # Permutation test
  n_perm <- 10000
  perm_diffs <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    perm_labels <- sample(df_all$condition)
    df_all$perm_condition <- perm_labels
    
    mdd_sum <- sum(df_all$count[df_all$perm_condition == "MDD" & df_all$group_involved])
    ctrl_sum <- sum(df_all$count[df_all$perm_condition == "Control" & df_all$group_involved])
    
    perm_diffs[i] <- mdd_sum - ctrl_sum
  }
  
  p_val <- mean(abs(perm_diffs) >= abs(obs_diff))
  
  # Save results
  results[[grp]] <- list(
    observed_diff = obs_diff,
    p_value = p_val,
    perm_diffs = perm_diffs
  )
  
  # Optional plot
  hist(perm_diffs, breaks = 50,
       main = paste("Permutation Test:", grp),
       xlab = paste("Difference (MDD - Control) for", grp),
       col = "lightgray", border = "white")
  abline(v = obs_diff, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste0("Obs = ", round(obs_diff, 1), "\nP = ", signif(p_val, 3)), bty = "n")
}


# Set up a 2x2 plotting area
par(mfrow = c(2, 2),  # 2 rows, 2 columns
    mar = c(4, 4, 3, 1))  # margins: bottom, left, top, right

# Loop through each group and plot
for (grp in names(results)) {
  grp_result <- results[[grp]]
  
  hist(grp_result$perm_diffs, breaks = 50,
       main = grp,
       xlab = "Difference (MDD - Control)",
       col = "lightgray", border = "white",
       xlim = range(c(grp_result$perm_diffs, grp_result$observed_diff)))
  
  abline(v = grp_result$observed_diff, col = "red", lwd = 2, lty = 2)
  
  legend("topright",
         legend = paste0("Obs = ", round(grp_result$observed_diff, 1), 
                         "\nP = ", signif(grp_result$p_value, 3)),
         bty = "n", text.col = "black")
}


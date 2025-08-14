object = cellchat
color.use = c("grey10", "#82d4a8", "#e86878")
comparison = c(1, 2)
signaling = NULL
signaling.label = NULL
top.label = 1
signaling.exclude = NULL
xlims = NULL
ylims = NULL
slot.name = "netP"
dot.size = 2.5
point.shape = c(21, 22, 24, 23)
label.size = 3
dot.alpha = 0.6
x.measure = "outdeg"
y.measure = "indeg"
xlabel = "Differential outgoing interaction strength"
ylabel = "Differential incoming interaction strength"
title = NULL
font.size = 10
font.size.title = 10
do.label = T
show.legend = T
show.axes = T

cell_type_order <- c("Ex", "Inh", "Oligo", "OPC", "Astro", "Micro", "Macrophages", "Tcells", "Endo", "Mural", "Fibro")


colors_dist <- c("#ff2d55", "#007aff", "#6a4ad3", "#b8a7d9","#ff9500", "#ffcc00", "#b38f00" ,"#675200", "#4cd964", "#2a8d3f", "#4cd9ab")
names(colors_dist) <- cell_type_order


# Initialize an empty list to store results
df_list <- list()

for (cell in cell_type_order) {
  idents.use <- cell
  dataset.name <- names(object@net)
  
  message(paste0("Visualizing differential outgoing and incoming signaling changes from ", 
                 dataset.name[comparison[1]], " to ", dataset.name[comparison[2]]))
  
  cell.levels <- levels(object@idents$joint)
  
  if (!(idents.use %in% cell.levels)) {
    warning(paste("Skipping", idents.use, "as it is not in cell levels"))
    next
  }
  
  signaling <- union(object@netP[[comparison[1]]]$pathways, 
                     object@netP[[comparison[2]]]$pathways)
  
  mat.all.merged <- list()
  
  for (ii in seq_along(comparison)) {
    if (length(slot(object, slot.name)[[comparison[ii]]]$centr) == 0) {
      stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores for each dataset separately!")
    }
    
    if (!all(c(x.measure, y.measure) %in% names(slot(object, slot.name)[[comparison[ii]]]$centr[[1]]))) {
      stop(paste0("`x.measure, y.measure` should be one of ", 
                  paste(names(slot(object, slot.name)[[comparison[ii]]]$centr[[1]]), collapse = ", ")))
    }
    
    centr <- slot(object, slot.name)[[comparison[ii]]]$centr
    outgoing <- matrix(0, nrow = length(cell.levels), ncol = length(centr))
    incoming <- matrix(0, nrow = length(cell.levels), ncol = length(centr))
    
    dimnames(outgoing) <- list(cell.levels, names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    
    for (i in seq_along(centr)) {
      outgoing[, i] <- centr[[i]][[x.measure]]
      incoming[, i] <- centr[[i]][[y.measure]]
    }
    
    mat.out <- t(outgoing)
    mat.in <- t(incoming)
    
    mat.all <- array(0, dim = c(length(signaling), ncol(mat.out), 2))
    mat.t <- list(mat.out, mat.in)
    
    for (i in seq_along(comparison)) {
      mat <- mat.t[[i]]
      mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
      mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
      idx <- match(rownames(mat1), signaling)
      mat[idx[!is.na(idx)], ] <- mat1
      dimnames(mat) <- list(signaling, colnames(mat1))
      mat.all[, , i] <- mat
    }
    
    dimnames(mat.all) <- list(dimnames(mat)[[1]], dimnames(mat)[[2]], c("outgoing", "incoming"))
    mat.all.merged[[ii]] <- mat.all
  }
  
  mat.all.merged.use <- list(mat.all.merged[[1]][, idents.use, ], mat.all.merged[[2]][, idents.use, ])
  idx.specific <- mat.all.merged.use[[1]] * mat.all.merged.use[[2]]
  mat.sum <- mat.all.merged.use[[2]] + mat.all.merged.use[[1]]
  
  out.specific.signaling <- rownames(idx.specific)[(mat.sum[, 1] != 0) & (idx.specific[, 1] == 0)]
  in.specific.signaling <- rownames(idx.specific)[(mat.sum[, 2] != 0) & (idx.specific[, 2] == 0)]
  
  mat.diff <- mat.all.merged.use[[2]] - mat.all.merged.use[[1]]
  idx <- rowSums(mat.diff) != 0
  mat.diff <- mat.diff[idx, , drop = FALSE]
  
  out.specific.signaling <- rownames(mat.diff) %in% out.specific.signaling
  in.specific.signaling <- rownames(mat.diff) %in% in.specific.signaling
  out.in.specific.signaling <- as.logical(out.specific.signaling * in.specific.signaling)
  
  specificity.out.in <- rep(0, nrow(mat.diff))
  specificity.out.in[out.in.specific.signaling] <- 2
  specificity.out.in[setdiff(which(out.specific.signaling), which(out.in.specific.signaling))] <- 1
  specificity.out.in[setdiff(which(in.specific.signaling), which(out.in.specific.signaling))] <- -1
  
  df <- as.data.frame(mat.diff)
  df$specificity.out.in <- specificity.out.in
  df$specificity <- 0
  df$specificity[(specificity.out.in != 0) & (rowSums(mat.diff >= 0) == 2)] <- 1
  df$specificity[(specificity.out.in != 0) & (rowSums(mat.diff <= 0) == 2)] <- -1
  
  out.in.category <- c("Shared", "Incoming specific", "Outgoing specific", "Incoming & Outgoing specific")
  specificity.category <- c("Shared", paste0(dataset.name[comparison[1]], " specific"), paste0(dataset.name[comparison[2]], " specific"))
  
  df$specificity.out.in <- plyr::mapvalues(df$specificity.out.in, from = c(0, -1, 1, 2), to = out.in.category)
  df$specificity.out.in <- factor(df$specificity.out.in, levels = out.in.category)
  df$specificity <- plyr::mapvalues(df$specificity, from = c(0, -1, 1), to = specificity.category)
  df$specificity <- factor(df$specificity, levels = specificity.category)
  
  point.shape.use <- point.shape[out.in.category %in% unique(df$specificity.out.in)]
  df$specificity.out.in <- droplevels(df$specificity.out.in, exclude = setdiff(out.in.category, unique(df$specificity.out.in)))
  color.use <- color.use[specificity.category %in% unique(df$specificity)]
  df$specificity <- droplevels(df$specificity, exclude = setdiff(specificity.category, unique(df$specificity)))
  
  df$labels <- rownames(df)
  
  # Add celltype column
  df$celltype <- idents.use
  df$color <- colors_dist[idents.use]
  
  # Store in list
  df_list[[idents.use]] <- df
}

# Combine all results into a single dataframe
df_aggregated <- do.call(rbind, df_list)

# Reset row names
rownames(df_aggregated) <- NULL


##plot
# df <- df_aggregated %>%
#   filter(celltype %in% c("Ex", "Inh"))
df_aggregated$broad_cell_group <- as.factor(df_aggregated$celltype)
levels(df_aggregated$broad_cell_group) <- c("Glia", "Perivascular", "Neuronal", "Perivascular", "Neuronal", "Immune", "Immune",
                                            "Perivascular", "Glia", "Glia", "Immune")
df_aggregated$broad_cell_group <- factor(df_aggregated$broad_cell_group, levels = c("Neuronal" , "Glia", "Immune", "Perivascular"))


thresh <- stats::quantile(abs(as.matrix(df_aggregated[, 1:2])), 
                          probs = 1 - top.label)
idx = abs(df_aggregated[, 1]) > thresh | abs(df_aggregated[, 2]) > thresh
data.label <- df_aggregated[idx, ]

gg <- ggplot(data = df_aggregated, aes(outgoing, incoming)) + 
  geom_point(aes(colour = specificity, fill = celltype, shape = specificity.out.in), 
             size = 4, stroke = 1.5) +
  theme_linedraw() + 
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", size = 0.25) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 0.25) + 
  theme(text = element_text(size = 18), 
        legend.key.height = grid::unit(0.15, "in"),
        axis.title = element_text(size = 20),  # Increase axis titles font size
        axis.text = element_text(size = 16),   # Increase axis labels font size
        strip.text = element_text(size = 18),  # Increase facet title font size
        legend.title = element_blank()) + 
  labs(title = "", x = xlabel, y = ylabel) + 
  theme(axis.line.x = element_line(size = 0.25), 
        axis.line.y = element_line(size = 0.25)) + 
  scale_fill_manual(values = colors_dist, drop = FALSE) + 
  scale_colour_manual(values = color.use, drop = FALSE) + 
  scale_shape_manual(values = point.shape.use) + 
  ggrepel::geom_text_repel(data = data.label, 
                           mapping = aes(label = labels, colour = specificity), 
                           size = 7, show.legend = F, segment.size = 0.5, 
                           segment.alpha = 0.5) +
  facet_wrap(~broad_cell_group, scales = "free", nrow = 2) + 
  guides(fill = "none") 

ggsave(filename = "/data/BBRF/Finalized_output/singlecell/CellChat/modified_in_out_cellgroup.pdf", h = 15, w = 15)


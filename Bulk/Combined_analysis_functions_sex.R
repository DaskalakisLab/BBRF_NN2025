inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

limma_wrapper <- function(anno, expr, covariates, out_dir, subset = "all", comparison_name = "gDxMDD") {
  
  # Subset the data based on the specified subset (all, male, female)
  if (subset == "male") {
    anno <- subset(anno, Sex == "M")
    expr <- expr[, match(anno$SampleID, colnames(expr))]
    comparison_name <- paste0(comparison_name, "_Males")
  } else if (subset == "female") {
    anno <- subset(anno, Sex == "F")
    expr <- expr[, match(anno$SampleID, colnames(expr))]
    comparison_name <- paste0(comparison_name, "_Females")
  } else {
    comparison_name <- paste0(comparison_name, "_All")
  }
  
  log_file <- file.path(out_dir, "analysis_log.txt")
  sink(log_file, append = TRUE)
  
  # Match expression matrix columns to annotation samples
  expr <- expr[, match(anno$SampleID, colnames(expr))]
  
  # Check if the order of columns matches between expr and anno
  if (!(identical(colnames(expr), anno$SampleID))) {
    stop("FATAL ERROR: orders do not match")
  } else {
    print("ORDER ALL GOOD")
  }
  
  # Create the design formula based on the covariates
  create_formula <- function(covariates) {
    form <- as.formula(paste0("~", paste(covariates, collapse = "+")))
    return(form)
  }
  
  form <- create_formula(covariates)
  print(form)
  
  # Create design matrix
  design <- model.matrix(form, data = anno)  # Fixed data = anno
  colnames(design) <- gsub(":", "_", colnames(design))
  print(colnames(design))
  print(dim(design))
  
  # Fit the linear model
  fit <- lmFit(expr, design)
  
  # Make contrast matrix
  myargs <- list("gDxMDD", levels=design)
  contrast.matrix<-do.call(makeContrasts, myargs)
  rownames(contrast.matrix)[1]<-"(Intercept)"
  # Apply contrasts to the fit
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  # Apply empirical Bayes
  fit2.ebayes <- eBayes(fit2)
  SE_table <- data.frame(s2.post = fit2.ebayes$s2.post, stdev.unscaled = fit2.ebayes$stdev.unscaled)
  SE_table$SE <- sqrt(SE_table$s2.post) * SE_table$gDxMDD
  
  # Run limma
  results <- limma::topTable(fit2.ebayes, number = nrow(fit2.ebayes), coef = 1, adjust = "BH")
  results[["ensID"]] <- rownames(results)
  
  results <- merge(results, SE_table, by = "row.names", all.x = TRUE)
  rownames(results) <- results$Row.names
  results <- results[, -1]
  
  # Bacon adjustment
  z <- abs(qnorm(results$P.Value)) * sign(results$logFC)
  results$z <- z
  results[is.na(results$z), "z"] <- 0
  
  bc <- bacon(results$z)
  results$p_bacon <- as.numeric(pval(bc))
  results$p_bacon_adj <- p.adjust(results$p_bacon, method = "fdr", n = nrow(results))
  
  lambda <- inflation(results$P.Value)
  lambda_bacon <- inflation(results$p_bacon)
  
  print(paste0("Lambda = ", round(lambda, 2)))
  print(paste0("Lambda corrected = ", round(lambda_bacon, 2)))
  
  results$lambda <- lambda
  results$lambda_bacon <- lambda_bacon
  
  # Load gene map and merge
  gene_map <- readRDS("/Volumes/humgen/daskalakislab/dipietro/SciencePaper/Data/RNA/Bulk/RNA_Gene_Map.RDS")
  results <- merge(results, gene_map, by.x = "ensID", by.y = "genes", all.x = TRUE)
  
  # Volcano plot preparation
  results$lp <- -log10(results$P.Value)
  results$Color <- 1
  
  results[(results$P.Value < 0.05) & (results$logFC < 0), "Color"] <- 2
  results[(results$P.Value < 0.05) & (results$logFC > 0), "Color"] <- 3
  results[(results$adj.P.Val < 0.05) & (results$logFC < 0), "Color"] <- 4
  results[(results$adj.P.Val < 0.05) & (results$logFC > 0), "Color"] <- 5
  
  col_df <- data.frame(Vals = c(1, 2, 3, 4, 5), 
                       Colors = c("grey", "lightblue", "pink", "darkblue", "darkred"))
  
  cols <- unique(results$Color)
  col_df <- col_df[col_df$Vals %in% cols, ]
  
  results <- results[order(results$P.Value), ]
  results$genelabels <- FALSE
  results$genelabels[1:30] <- TRUE
  
  results[results$symbol == "FKBP5", "genelabels"] <- TRUE
  results[results$symbol == "NR3C1", "genelabels"] <- TRUE
  results[results$symbol == "NR3C2", "genelabels"] <- TRUE
  results[results$symbol == "CRH", "genelabels"] <- TRUE
  results[results$symbol == "CRHR1", "genelabels"] <- TRUE
  results[results$symbol == "CRHR2", "genelabels"] <- TRUE
  
  # Plot the volcano plot
  p <- ggplot(results) +
    geom_point(aes(x = logFC, y = lp, colour = factor(Color))) +
    geom_text_repel(aes(x = logFC, y = lp, label = ifelse(genelabels == TRUE, symbol, "")), 
                    fontface = "italic", size = 3, max.overlaps = 1000) +
    scale_color_manual(values = col_df$Colors) +
    xlab("log2(FC)") +
    ylab("-log10(pvalue)") +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 20, face = "bold"))
  
  # Save the plot
  ggsave(filename = paste0(out_dir, "/Volcano_", comparison_name, ".pdf"),
         plot = p,
         device = "pdf",
         units = "in",
         width = 9,
         height = 12)
  
  saveRDS(anno, paste0(out_dir,"/annotation_", comparison_name, ".RDS"))
  saveRDS(expr, paste0(out_dir,"/expression_", comparison_name, ".RDS"))
  saveRDS(results, paste0(out_dir,"/results_", comparison_name, ".RDS"))
  
  limma_breakdown <- data.frame(
    n_ctrl = length(which(anno$gDx == "Control")),
    n_case = length(which(anno$gDx == "MDD")),
    n_genes = nrow(results),
    n_nominal = nrow(results[(results$P.Value < 0.05) & (results$adj.P.Val > 0.05), ]),
    n_fdr = nrow(results[(results$adj.P.Val < 0.05), ]),
    n_bacon = nrow(results[(results$p_bacon < 0.05), ]),
    n_bacon_adj = nrow(results[(results$p_bacon_adj < 0.05), ]),
    lambda = results$lambda[1], lambda_bacon = results$lambda_bacon[1]
  )
  write.csv(
    limma_breakdown,
    paste0(out_dir, "/breakdown_", comparison_name, ".csv"),
    quote = FALSE, row.names = FALSE
  )
  sink()
  print("Returning results")
  return(results)
}

limma_wrapper_sa <- function(anno, expr, covariates, out_dir, subset = "all", comparison_name = "gDxMDD") {
  
  # Subset the data based on the specified subset (all, male, female)
  if (subset == "male") {
    anno <- subset(anno, Sex == "M")
    expr <- expr[, match(anno$SampleID, colnames(expr))]
    comparison_name <- paste0(comparison_name, "_Males")
  } else if (subset == "female") {
    anno <- subset(anno, Sex == "F")
    expr <- expr[, match(anno$SampleID, colnames(expr))]
    comparison_name <- paste0(comparison_name, "_Females")
  } else {
    comparison_name <- paste0(comparison_name, "_All")
  }
  
  log_file <- file.path(out_dir, "analysis_log.txt")
  sink(log_file, append = TRUE)
  
  # Match expression matrix columns to annotation samples
  expr <- expr[, match(anno$SampleID, colnames(expr))]
  
  # Check if the order of columns matches between expr and anno
  if (!(identical(colnames(expr), anno$SampleID))) {
    stop("FATAL ERROR: orders do not match")
  } else {
    print("ORDER ALL GOOD")
  }
  
  # Create the design formula based on the covariates
  create_formula <- function(covariates) {
    form <- as.formula(paste0("~", paste(covariates, collapse = "+")))
    return(form)
  }
  
  form <- create_formula(covariates)
  print(form)
  
  # Create design matrix
  design <- model.matrix(form, data = anno)  # Fixed data = anno
  colnames(design) <- gsub(":", "_", colnames(design))
  print(colnames(design))
  print(dim(design))
  
  # Fit the linear model
  fit <- lmFit(expr, design)
  
  # Make contrast matrix
  myargs <- list("gDxMDD", levels=design)
  contrast.matrix<-do.call(makeContrasts, myargs)
  rownames(contrast.matrix)[1]<-"(Intercept)"
  # Apply contrasts to the fit
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  # Apply empirical Bayes
  fit2.ebayes <- eBayes(fit2)
  SE_table <- data.frame(s2.post = fit2.ebayes$s2.post, stdev.unscaled = fit2.ebayes$stdev.unscaled)
  SE_table$SE <- sqrt(SE_table$s2.post) * SE_table$gDxMDD
  
  # Run limma
  results <- limma::topTable(fit2.ebayes, number = nrow(fit2.ebayes), coef = 1, adjust = "BH")
  results[["ensID"]] <- rownames(results)
  
  results <- merge(results, SE_table, by = "row.names", all.x = TRUE)
  rownames(results) <- results$Row.names
  results <- results[, -1]
  
  # Bacon adjustment
  z <- abs(qnorm(results$P.Value)) * sign(results$logFC)
  results$z <- z
  results[is.na(results$z), "z"] <- 0
  
  bc <- bacon(results$z)
  results$p_bacon <- as.numeric(pval(bc))
  results$p_bacon_adj <- p.adjust(results$p_bacon, method = "fdr", n = nrow(results))
  
  lambda <- inflation(results$P.Value)
  lambda_bacon <- inflation(results$p_bacon)
  
  print(paste0("Lambda = ", round(lambda, 2)))
  print(paste0("Lambda corrected = ", round(lambda_bacon, 2)))
  
  results$lambda <- lambda
  results$lambda_bacon <- lambda_bacon
  
  # Load gene map and merge
  gene_map <- readRDS("/Volumes/humgen/daskalakislab/dipietro/SciencePaper/Data/RNA/Bulk/RNA_Gene_Map.RDS")
  results <- merge(results, gene_map, by.x = "ensID", by.y = "genes", all.x = TRUE)
  
  # Volcano plot preparation
  results$lp <- -log10(results$P.Value)
  results$Color <- 1
  
  results[(results$P.Value < 0.05) & (results$logFC < 0), "Color"] <- 2
  results[(results$P.Value < 0.05) & (results$logFC > 0), "Color"] <- 3
  results[(results$adj.P.Val < 0.05) & (results$logFC < 0), "Color"] <- 4
  results[(results$adj.P.Val < 0.05) & (results$logFC > 0), "Color"] <- 5
  
  col_df <- data.frame(Vals = c(1, 2, 3, 4, 5), 
                       Colors = c("grey", "lightblue", "pink", "darkblue", "darkred"))
  
  cols <- unique(results$Color)
  col_df <- col_df[col_df$Vals %in% cols, ]
  
  results <- results[order(results$P.Value), ]
  results$genelabels <- FALSE
  results$genelabels[1:30] <- TRUE
  
  results[results$symbol == "FKBP5", "genelabels"] <- TRUE
  results[results$symbol == "NR3C1", "genelabels"] <- TRUE
  results[results$symbol == "NR3C2", "genelabels"] <- TRUE
  results[results$symbol == "CRH", "genelabels"] <- TRUE
  results[results$symbol == "CRHR1", "genelabels"] <- TRUE
  results[results$symbol == "CRHR2", "genelabels"] <- TRUE
  
  # Plot the volcano plot
  p <- ggplot(results) +
    geom_point(aes(x = logFC, y = lp, colour = factor(Color))) +
    geom_text_repel(aes(x = logFC, y = lp, label = ifelse(genelabels == TRUE, symbol, "")), 
                    fontface = "italic", size = 3, max.overlaps = 1000) +
    scale_color_manual(values = col_df$Colors) +
    xlab("log2(FC)") +
    ylab("-log10(pvalue)") +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 20, face = "bold"))
  
  # Save the plot
  ggsave(filename = paste0(out_dir, "/Volcano_", comparison_name, ".pdf"),
         plot = p,
         device = "pdf",
         units = "in",
         width = 9,
         height = 12)
  
  saveRDS(anno, paste0(out_dir,"/annotation_", comparison_name, ".RDS"))
  saveRDS(expr, paste0(out_dir,"/expression_", comparison_name, ".RDS"))
  saveRDS(results, paste0(out_dir,"/results_", comparison_name, ".RDS"))
  
  limma_breakdown <- data.frame(
    n_ctrl = length(which(anno$gDx == "Control")),
    n_case = length(which(anno$gDx == "MDD")),
    n_genes = nrow(results),
    n_nominal = nrow(results[(results$P.Value < 0.05) & (results$adj.P.Val > 0.05), ]),
    n_fdr = nrow(results[(results$adj.P.Val < 0.05), ]),
    n_bacon = nrow(results[(results$p_bacon < 0.05), ]),
    n_bacon_adj = nrow(results[(results$p_bacon_adj < 0.05), ]),
    lambda = results$lambda[1], lambda_bacon = results$lambda_bacon[1]
  )
  write.csv(
    limma_breakdown,
    paste0(out_dir, "/breakdown_", comparison_name, ".csv"),
    quote = FALSE, row.names = FALSE
  )
  sink()
  print("Returning results")
  return(results)
}

limma_wrapper_ct <- function(anno, expr, covariates, out_dir, subset = "all", comparison_name = "ChildTrauma012") {
  
  # Subset the data based on the specified subset (all, male, female)
  if (subset == "male") {
    anno <- subset(anno, Sex == "M")
    expr <- expr[, match(anno$SampleID, colnames(expr))]
    comparison_name <- paste0(comparison_name, "_Males")
  } else if (subset == "female") {
    anno <- subset(anno, Sex == "F")
    expr <- expr[, match(anno$SampleID, colnames(expr))]
    comparison_name <- paste0(comparison_name, "_Females")
  } else {
    comparison_name <- paste0(comparison_name, "_All")
  }
  
  log_file <- file.path(out_dir, "analysis_log.txt")
  sink(log_file, append = TRUE)
  
  # Match expression matrix columns to annotation samples
  expr <- expr[, match(anno$SampleID, colnames(expr))]
  anno$ChildTrauma012 <- anno$`Child trauma_0/1/2` 
  
  # Check if the order of columns matches between expr and anno
  if (!(identical(colnames(expr), anno$SampleID))) {
    stop("FATAL ERROR: orders do not match")
  } else {
    print("ORDER ALL GOOD")
  }
  
  # Create the design formula based on the covariates
  covariates[1] <- "ChildTrauma012"
  
  create_formula <- function(covariates) {
    form <- as.formula(paste0("~", paste(covariates, collapse = "+")))
    return(form)
  }
  
  form <- create_formula(covariates)
  print(form)
  
  # Create design matrix
  design <- model.matrix(form, data = anno)  # Fixed data = anno
  colnames(design) <- gsub(":", "_", colnames(design))
  print(colnames(design))
  print(dim(design))
  
  # Fit the linear model
  fit <- lmFit(expr, design)
  
  # Make contrast matrix
  myargs <- list("ChildTrauma012", levels=design)
  contrast.matrix<-do.call(makeContrasts, myargs)
  rownames(contrast.matrix)[1]<-"(Intercept)"
  # Apply contrasts to the fit
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  # Apply empirical Bayes
  fit2.ebayes <- eBayes(fit2)
  SE_table <- data.frame(s2.post = fit2.ebayes$s2.post, stdev.unscaled = fit2.ebayes$stdev.unscaled)
  SE_table$SE <- sqrt(SE_table$s2.post) * SE_table$ChildTrauma012
  
  # Run limma
  results <- limma::topTable(fit2.ebayes, number = nrow(fit2.ebayes), coef = 1, adjust = "BH")
  results[["ensID"]] <- rownames(results)
  
  results <- merge(results, SE_table, by = "row.names", all.x = TRUE)
  rownames(results) <- results$Row.names
  results <- results[, -1]
  
  # Bacon adjustment
  z <- abs(qnorm(results$P.Value)) * sign(results$logFC)
  results$z <- z
  results[is.na(results$z), "z"] <- 0
  
  bc <- bacon(results$z)
  results$p_bacon <- as.numeric(pval(bc))
  results$p_bacon_adj <- p.adjust(results$p_bacon, method = "fdr", n = nrow(results))
  
  lambda <- inflation(results$P.Value)
  lambda_bacon <- inflation(results$p_bacon)
  
  print(paste0("Lambda = ", round(lambda, 2)))
  print(paste0("Lambda corrected = ", round(lambda_bacon, 2)))
  
  results$lambda <- lambda
  results$lambda_bacon <- lambda_bacon
  
  # Load gene map and merge
  gene_map <- readRDS("/Volumes/humgen/daskalakislab/dipietro/SciencePaper/Data/RNA/Bulk/RNA_Gene_Map.RDS")
  results <- merge(results, gene_map, by.x = "ensID", by.y = "genes", all.x = TRUE)
  
  # Volcano plot preparation
  results$lp <- -log10(results$P.Value)
  results$Color <- 1
  
  results[(results$P.Value < 0.05) & (results$logFC < 0), "Color"] <- 2
  results[(results$P.Value < 0.05) & (results$logFC > 0), "Color"] <- 3
  results[(results$adj.P.Val < 0.05) & (results$logFC < 0), "Color"] <- 4
  results[(results$adj.P.Val < 0.05) & (results$logFC > 0), "Color"] <- 5
  
  col_df <- data.frame(Vals = c(1, 2, 3, 4, 5), 
                       Colors = c("grey", "lightblue", "pink", "darkblue", "darkred"))
  
  cols <- unique(results$Color)
  col_df <- col_df[col_df$Vals %in% cols, ]
  
  results <- results[order(results$P.Value), ]
  results$genelabels <- FALSE
  results$genelabels[1:30] <- TRUE
  
  results[results$symbol == "FKBP5", "genelabels"] <- TRUE
  results[results$symbol == "NR3C1", "genelabels"] <- TRUE
  results[results$symbol == "NR3C2", "genelabels"] <- TRUE
  results[results$symbol == "CRH", "genelabels"] <- TRUE
  results[results$symbol == "CRHR1", "genelabels"] <- TRUE
  results[results$symbol == "CRHR2", "genelabels"] <- TRUE
  
  # Plot the volcano plot
  p <- ggplot(results) +
    geom_point(aes(x = logFC, y = lp, colour = factor(Color))) +
    geom_text_repel(aes(x = logFC, y = lp, label = ifelse(genelabels == TRUE, symbol, "")), 
                    fontface = "italic", size = 3, max.overlaps = 1000) +
    scale_color_manual(values = col_df$Colors) +
    xlab("log2(FC)") +
    ylab("-log10(pvalue)") +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 20, face = "bold"))
  
  # Save the plot
  ggsave(filename = paste0(out_dir, "/Volcano_", comparison_name, ".pdf"),
         plot = p,
         device = "pdf",
         units = "in",
         width = 9,
         height = 12)
  
  saveRDS(anno, paste0(out_dir,"/annotation_", comparison_name, ".RDS"))
  saveRDS(expr, paste0(out_dir,"/expression_", comparison_name, ".RDS"))
  saveRDS(results, paste0(out_dir,"/results_", comparison_name, ".RDS"))
  
  limma_breakdown <- data.frame(
    n_ctrl = length(which(anno$gDx == "Control")),
    n_case = length(which(anno$gDx == "MDD")),
    n_genes = nrow(results),
    n_nominal = nrow(results[(results$P.Value < 0.05) & (results$adj.P.Val > 0.05), ]),
    n_fdr = nrow(results[(results$adj.P.Val < 0.05), ]),
    n_bacon = nrow(results[(results$p_bacon < 0.05), ]),
    n_bacon_adj = nrow(results[(results$p_bacon_adj < 0.05), ]),
    lambda = results$lambda[1], lambda_bacon = results$lambda_bacon[1]
  )
  write.csv(
    limma_breakdown,
    paste0(out_dir, "/breakdown_", comparison_name, ".csv"),
    quote = FALSE, row.names = FALSE
  )
  sink()
  print("Returning results")
  return(results)
}



run_gsea_analysis <- function(results, category = "C2", subcategory = "REACTOME", species = "Homo sapiens", 
                              minSize = 10, maxSize = 500, out_dir = out_dir, comparison_name = "gDxMDD", subset = "all") {
  library(fgsea)
  library(ggplot2)
  library(msigdbr)
  library(dplyr)
  library(tibble)
  
  if (subset == "male") {
    comparison_name <- paste0(comparison_name, "_Males")
  } else if (subset == "female") {
    comparison_name <- paste0(comparison_name, "_Females")
  } else {
    comparison_name <- paste0(comparison_name, "_All")
  }
  
  # Load pathways based on specified category and subcategory
  pathways <- msigdbr(species = species, category = category, subcategory = subcategory)
  
  # Prepare fgsea gene sets
  fgsea_sets <- pathways %>% split(x = .$gene_symbol, f = .$gs_name)
  
  # Prepare the dataset for GSEA
  gsea_input <- results %>%
    arrange(desc(logFC)) %>%
    dplyr::select(symbol, logFC)
  
  ranks <- deframe(gsea_input)
  
  # Run GSEA
  fgseaRes <- fgseaMultilevel(fgsea_sets, stats = ranks, 
                              minSize = minSize, maxSize = maxSize, eps = 0)
  
  # Prepare data to save and plot
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidier <- fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -log2err) %>% 
    arrange(padj) 
  
  # Plot top pathways
  # Get the top 10 and bottom 10 pathways based on NES
  top_10_pathways <- fgseaResTidy %>% 
    filter(padj < 0.05) %>%
    arrange(desc(NES)) %>%
    head(n = 10)
  
  bottom_10_pathways <- fgseaResTidy %>% 
    filter(padj < 0.05) %>%
    arrange(NES) %>%
    head(n = 10)
  
  # Combine the top and bottom pathways
  top_bottom_pathways <- bind_rows(top_10_pathways, bottom_10_pathways)
  
  # Plot
  g <- ggplot(top_bottom_pathways, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = NES < 0)) +  # Fill based on whether NES is negative or positive
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = paste("Top and Bottom GSEA pathways -", category, subcategory)) + 
    theme_minimal()
  
  # Save the plot
  ggsave(filename = paste0(out_dir, "/Top_Bottom_GSEA_pathways_", category, "_", subcategory, ".pdf"),
         plot = g, width = 10, height = 12)
  
  # Save results as CSV and RData
  write.csv(fgseaResTidier, paste0(out_dir, "/gsea_", category, "_", subcategory, ".csv"), quote = F, row.names = F)
  save(fgseaResTidy, file = paste0(out_dir, "/gsea_", category, "_", subcategory, ".RData"))
  
  return(list(topPathways = fgseaResTidy, filePaths = list(csv = paste0(out_dir, "/gsea_", category, "_", subcategory, ".csv"), 
                                                           rdata = paste0(out_dir, "/gsea_", category, "_", subcategory, ".RData"),
                                                           pdf = paste0(out_dir, "/Top_GSEA_pathways_", category, "_", subcategory, ".pdf"))))
}

gene_set_enrichment_plot_complex_mod <- function(enrichment,
                                                 PThresh = 12,
                                                 ORcut = 3,
                                                 enrichOnly = FALSE,
                                                 gene_count_col = NULL,
                                                 gene_count_row = NULL,
                                                 anno_title_col = NULL,
                                                 anno_title_row = NULL,
                                                 column_order = NULL,
                                                 anno_add = NULL,
                                                 mypal = c(
                                                   "white",
                                                          grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(50)
                                                 )) {
  
  ## Check inputs
  stopifnot(is(enrichment, "data.frame"))
  stopifnot(all(c("ID", "test", "OR", "Pval", "padj") %in% colnames(enrichment)))
  stopifnot(ORcut <= PThresh)
  # stopifnot(length(xlabs) == length(unique(enrichment$ID)))
  
  ## Convert to -log10 scale and threshold the pvalues
  enrichment$log10_P_thresh <-
    round(-log10(enrichment$padj), 2)
  enrichment$log10_P_thresh[which(enrichment$log10_P_thresh > PThresh)] <-
    PThresh
  
  ## Change some values for the plot
  if (enrichOnly) {
    enrichment$log10_P_thresh[enrichment$OR < 1] <- 0
  }
  enrichment$OR_char <- as.character(round(enrichment$OR, 2))
  enrichment$OR_char[enrichment$log10_P_thresh < ORcut] <- ""
  
  ## Make into wide matrices
  make_wide <- function(var = "OR_char") {
    res <-
      reshape(
        enrichment,
        idvar = "ID",
        timevar = "test",
        direction = "wide",
        drop = colnames(enrichment)[!colnames(enrichment) %in% c("ID", "test", var)],
        sep = "_mypattern_"
      )[, -1, drop = FALSE]
    colnames(res) <-
      gsub(".*_mypattern_", "", colnames(res))
    rownames(res) <- unique(enrichment$ID)
    res <- res[, levels(as.factor(enrichment$test))]
    t(res)
  }
  
  ## Define matrix
  wide_or <- make_wide("OR_char")
  wide_p <- make_wide("log10_P_thresh")
  
  ## Reorder
  if (!is.null(column_order)) {
    stopifnot(setequal(column_order, colnames(wide_or)))
    wide_or <- wide_or[, column_order]
    wide_p <- wide_p[, column_order]
  }
  
  if (!is.null(anno_add)) {
    stopifnot(setequal(colnames(anno_add), colnames(wide_or)))
    stopifnot(setequal(rownames(anno_add), rownames(wide_or)))
    
    wide_or[] <- paste0(anno_add[rownames(wide_or), colnames(wide_or)], "\n", wide_or)
  }
  
  ## define annotations
  stopifnot(setequal(rownames(gene_count_col), colnames(wide_p)))
  stopifnot(setequal(rownames(gene_count_row), rownames(wide_p)))
  
  col_gene_anno <- ComplexHeatmap::columnAnnotation(
    `n genes` = ComplexHeatmap::anno_barplot(gene_count_col[colnames(wide_p), ]),
    annotation_label = anno_title_col
  )
  row_gene_anno <- ComplexHeatmap::rowAnnotation(
    `n genes` = ComplexHeatmap::anno_barplot(gene_count_row[rownames(wide_p), ]),
    annotation_label = anno_title_row
  )
  
  ComplexHeatmap::Heatmap(wide_p,
                          col = mypal,
                          name = "-log10(p-val)",
                          rect_gp = grid::gpar(col = "black", lwd = 1),
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          right_annotation = row_gene_anno,
                          top_annotation = col_gene_anno,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid::grid.text(wide_or[i, j], x, y, gp = grid::gpar(fontsize = 10))
                          }
  )
}






run_pipeline <- function(demo_data, expr_data, covariates, covariates2, out_dir, bayesSpace_registration, bayes_anno, spatial_reg_enrich_threshold, gsea_categories) {
  
  subsets <- c("all", "male", "female")
  
  for (subset in subsets) {
    # Dynamically set covariates based on subset
    if (subset == "all") {
      current_covariates <- covariates
    } else {
      current_covariates <- covariates2
    }
    
    # Step 0: Load libraries
    required_packages <- c("dplyr", "limma", "fgsea", "msigdbr", "spatialLIBD", "ggplot2", "jaffelab", "here", "purrr", "bacon", "ggrepel")
    
    lapply(required_packages, function(pkg) {
      if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
      library(pkg, character.only = TRUE)
    })
    
    print("Packages loaded!")
    
    # Create output directory based on subset
    subset_out_dir <- file.path(out_dir, subset)
    
    print(paste("Attempting to create directory:", subset_out_dir))
    dir.create(subset_out_dir, showWarnings = TRUE, recursive = TRUE)
    
    if (!dir.exists(subset_out_dir)) {
      stop("Failed to create directory: ", subset_out_dir)
    } else {
      print("Directory created successfully!")
    } 
    
    # Log file for each subset
    log_file <- paste0(subset_out_dir, "/log_", subset, "_covariates_", paste(current_covariates, collapse = "_"), ".txt")
    write(paste("Subset:", subset), log_file)
    write(paste("Covariates used:", paste(current_covariates, collapse = ", ")), log_file, append = TRUE)
    
    # Step 1: Differential Expression Analysis using Limma
    print(paste("Running limma for subset:", subset))
    results <- limma_wrapper(
      anno = demo_data,
      expr = expr_data,
      covariates = current_covariates,
      out_dir = subset_out_dir,
      subset = subset,
      comparison_name = "gDxMDD"
    )
    
    print("Limma done!")
    write("Limma done!", log_file, append = TRUE)
    
    # Step 2: GSEA Analysis
    gsea_results <- list()
    for (category in gsea_categories) {
      cat <- category[1]
      subcat <- category[2]
      gsea_results[[paste(cat, subcat, sep = "_")]] <- run_gsea_analysis(
        results = results,
        category = cat,
        subcategory = subcat,
        out_dir = subset_out_dir,
        subset = subset
      )
    }
    
    print("GSEA done!")
    write("GSEA done!", log_file, append = TRUE)
    
    # Step 3: Spatial Registration Enrichment Analysis
    # Prepare significant gene list for MDD
    sigGeneSign <- cbind(
      results$ensID[which(results$adj.P.Val <= spatial_reg_enrich_threshold)],
      sign(results$logFC[which(results$adj.P.Val <= spatial_reg_enrich_threshold)])
    )
    
    if (nrow(sigGeneSign) <= 2) {
      # No significant genes
      print("No significant genes found for spatial enrichment.")
      write("No significant genes found for spatial enrichment.", log_file, append = TRUE)
      next
    }
    
    sigGeneSign[,1] <- as.character(sigGeneSign[,1])
    sigGeneSign[,2] <- as.numeric(sigGeneSign[,2])
    sigGeneSign <- data.frame(sigGeneSign)
    colnames(sigGeneSign) <- c("gene", "change")
    cleaned_gene_id <- sub("\\..*", "", sigGeneSign$gene)
    sigGeneSign$gene <- cleaned_gene_id
    mdd_gene_list <- list(MDD = sigGeneSign)
    
    # Run enrichment analysis
    enriched <- map(mdd_gene_list, function(gl) {
      map(bayesSpace_registration, ~
            gene_set_enrichment(gene_list = gl, modeling_results = .x, model_type = "enrichment") %>% 
            left_join(bayes_anno, by = "test") %>% 
            mutate(
              test = factor(layer_combo, levels = bayes_anno$layer_combo[bayes_anno$layer_combo %in% layer_combo])
            ) %>% 
            dplyr::select(-c(layer_combo, Annotation, fdr_cut, model_type))
      )
    })
    
    enriched <- map_depth(enriched, 2, ~ {
      .x$padj <- p.adjust(.x$Pval, method = "BH")
      return(.x)
    })
    
    saveRDS(enriched, file = paste0(subset_out_dir, "/enriched_MDD.rds"))
    
    # Summarize results
    gene_enrichment_count <- map(bayesSpace_registration, function(r) {
      en_count <- get_gene_enrichment_count(r)
      rownames(en_count) <- bayes_anno$layer_combo[match(rownames(en_count), bayes_anno$test)]
      layer_order <- bayes_anno$layer_combo[bayes_anno$layer_combo %in% rownames(en_count)]
      return(en_count[layer_order, , drop = FALSE])
    })
    
    gene_list_count <- map(mdd_gene_list, get_gene_list_count)
    
    pdf(paste0(subset_out_dir, "/Enrich_MDD_k09_adj.pdf"), height = 8, width = 10)
    gene_set_enrichment_plot_complex_mod(enriched$MDD$k09,
                                         gene_count_col = gene_list_count[["MDD"]],
                                         gene_count_row = gene_enrichment_count[["k09"]],
                                         anno_title_col = "n DE Genes",
                                         anno_title_row = "n Domain\nGenes"
    )
    dev.off()
    
    walk2(enriched, names(enriched), function(enriched, ds_name) {
      pdf(paste0(subset_out_dir,"/Enrich_MDD_", ds_name, ".pdf"), height = 8, width = 9)
      map2(enriched, names(enriched), function(x, k) {
        message(ds_name, " - ", k)
        
        print(gene_set_enrichment_plot_complex(x,
                                               gene_count_col = gene_list_count[[ds_name]],
                                               gene_count_row = gene_enrichment_count[[k]],
                                               anno_title_col = "n DE Genes",
                                               anno_title_row = "n Domain\nGenes"
        ))
      })
      dev.off()
    })
    
    print("Spatial enrichment done!")
    write("Spatial enrichment done!", log_file, append = TRUE)
  }
}


limma_analysis_split <- function(anno, expr, covariates, out_dir, subset = "all") {
  
  # Subset for MDD vs Control
  anno_MDD <- anno[anno$PrimaryDx %in% c("MDD", "Control"), ]
  expr_MDD <- expr[, match(anno_MDD$SampleID, colnames(expr))]
  print("Running limma for MDD vs Control...")
  results_MDD <- limma_wrapper_sa(anno_MDD, expr_MDD, covariates, out_dir, subset, comparison_name = "PrimaryDxMDD")
  
  # Subset for PTSD vs Control
  anno_PTSD <- anno[anno$PrimaryDx %in% c("PTSD", "Control"), ]
  expr_PTSD <- expr[, match(anno_PTSD$SampleID, colnames(expr))]
  print("Running limma for PTSD vs Control...")
  results_PTSD <- limma_wrapper_sa(anno_PTSD, expr_PTSD, covariates, out_dir, subset, comparison_name = "gDxPTSD")
  
  # Subset for PTSD vs Control
  anno_CT <- anno[complete.cases(demo$`Child trauma_0/1/2`), ]
  expr_CT <- expr[, match(anno_CT$SampleID, colnames(expr))]
  print("Running limma for CT...")
  results_CT <- limma_wrapper_ct(anno_CT, expr_CT, covariates, out_dir, subset, comparison_name = "ChildTrauma012")
  
  return(list(MDD_results = results_MDD, PTSD_results = results_PTSD, CT_results = results_CT))
}
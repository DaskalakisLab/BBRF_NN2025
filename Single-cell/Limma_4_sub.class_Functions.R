plot_volcano <- function(results, dx, contrast, region, ptsd_vs_mdd=F, ptsd_mdd=F, file ){
  
  if (ptsd_vs_mdd){
    dx <- "PTSD vs MDD"
  }else if (ptsd_mdd){
    dx <- "PTSD+MDD"
  }
  title <- paste0("Dx: ", dx, ", Contrast: ", contrast, ", Region: ", region)
  results$lp <- -log10(results$P.Value)
  results$Color <- 1
  
  results[(results$P.Value < 0.05) & (results$logFC < 0), "Color"] <- 2
  results[(results$P.Value < 0.05) & (results$logFC > 0), "Color"] <- 3
  results[(results$adj.P.Val < 0.05) & (results$logFC < 0), "Color"] <- 4
  results[(results$adj.P.Val < 0.05) & (results$logFC > 0), "Color"] <- 5
  
  col_df <- data.frame(Vals=c(1,2,3,4,5), 
                       Colors=c("grey","lightblue","pink","#2c7891","#ff0f39"))
  
  cols <- unique(results$Color)
  col_df <- col_df[col_df$Vals %in% cols,]
  
  results <- results[order(results$P.Value),]
  results$genelabels <- F
  results$genelabels[1:10] <- T
  results[results$genes=="ARL17B", "genelabels"] <- T
  results[results$genes=="FKBP5", "genelabels"] <- T
  results[results$genes=="CORT", "genelabels"] <- T
  results[results$genes=="CRH", "genelabels"] <- T
  
  p<-ggplot(results) +
    geom_point(aes(x = logFC, y = lp, colour = factor(Color))) +
    geom_text_repel(aes(x = logFC, y = lp, label = ifelse(genelabels == T, genes,"")),fontface="italic",size=3,max.overlaps=1000) +
    scale_color_manual(values=col_df$Colors)+
    xlab("log2(FC)") +
    ylab("-log10(pvalue)") +
    ggtitle(title)+
    theme_bw()+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 10,face="bold")
    )
  
  ggsave(filename=file,
         plot = p,
         device = "pdf",
         units = "in",
         width = 6,
         height = 8)
  
  
}

inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

DGE_limma <- function(group_name, covariates, output_dir, min.count = min.count, min.prop = min.prop, sva=F){
  sink(file = paste0(output_dir, "log_", group_name, ".txt"), append = TRUE)
  
  for(i in 1:length(class)){
    for(j in 1:length(subdir_list[[class[i]]])){
      cat("Starting with ", class[i],"/", subdir_list[[class[i]]][j], "\n")
      pg=paste0(pwd11,class[i],"/", subdir_list[[class[i]]][j] )
      
      
      pheno <- readRDS(paste0("/data/humgen/daskalakislab/aiatrou/BBRF/Pheno_Demo/Pseudobulk_Dec/phenofile_class_", class[i], ".rds"))
      
      allfiles <- list.files(path = pg)
      dge <- readDGE(allfiles, path=pg,columns = c(1,2))
      dge <- DGEList(dge)
      
      # Calculate the library size (column-wise sum of counts)
      library_sizes <- colSums(dge$counts)
      
      # Exclude samples with library size of zero
      dge <- dge[, library_sizes > 0]
      
      if (ncol(dge) == 0) {
        cat("Skipping ", class[i], " ", subdir_list[[class[i]]][j], " - no valid samples left after filtering.\n")
        next
      }
      
      dge$samples$group <- as.factor(pheno$gDx[match(rownames(dge$samples), pheno$sampleID)])
      keep <- filterByExpr(dge, min.count = min.count, min.prop = min.prop, group = dge$samples$group)
      
      dge <- dge[keep, keep.lib.sizes=FALSE]
      cat("Data dimensions after filtering: ", dim(dge), "\n")
      
      sample_names=colnames(dge)
      pheno <- pheno[match(sample_names,pheno$sampleID),]
      nsamples <- ncol(dge)
      dge <- calcNormFactors(dge, method = "TMM") 
      
      
      
      if(sva == T){
        mm <- model.matrix(~ as.factor(gDx) + as.factor(sex) + as.numeric(age) +
                             as.factor(version)+ as.factor(tissue), pheno)
        mm0 <- model.matrix(~ as.factor(sex) +age +
                              as.factor(version)+ as.factor(tissue), pheno)
        fit <- svaseq(cpm(dge), mod=mm, mod0=mm0)
        SVs <- fit$sv
        num_svs <- ncol(SVs)
        
        colnames(SVs) <- paste0("sv", seq_len(num_svs)) 
        
        pheno <- cbind(pheno, SVs)
        
        covariates_all <- c(covariates, paste0("sv", seq_len(num_svs)))
        
      } else {
        covariates_all <- c(covariates)
        mod = model.matrix(formula(paste0("~", paste(covariates_all, collapse = "+"))), data = pheno)
      }
      
      #
      form <- as.formula(paste("~", paste(covariates_all, collapse = " + ")))
      C <- canCorPairs(form, pheno)
      
      pdf(paste0(output_dir, "corrmat_", class[i],"_", subdir_list[[class[i]]][j] , ".pdf"))
      plotCorrMatrix(C)
      dev.off()
      
      #
      
      #Find SVs that correlate highly (> 0.7) with the covariates of interest
      high_corr_svs <- c()
      
      for (cov in covariates) {
        sv_corr <- C[cov, grep("^sv", colnames(C))]
        
        # Find SVs with correlation > 0.7
        correlated_svs <- names(sv_corr[abs(sv_corr) > 0.7])
        
        # Store the names of these highly correlated SVs
        high_corr_svs <- c(high_corr_svs, correlated_svs)
      }
      
      covariates_filtered <- covariates_all[!covariates_all %in% high_corr_svs]
      #
      
      mod = model.matrix(formula(paste0("~", paste(covariates_filtered, collapse = "+"))), data = pheno)
      
      #
      
      colnames(mod) <- make.names(colnames(mod))
      colnames(mod)[1] = "Int"
      
      vGene <- voom(dge,mod, plot = F)
      
      fitGene = lmFit(vGene, mod)
      eBGene = eBayes(fitGene)
      SE_table <- data.frame(s2.post=eBGene$s2.post, stdev.unscaled=eBGene$stdev.unscaled[,2])
      names(SE_table) <- c("s2.post","stdev.unscaled")
      SE_table$SE <- sqrt(SE_table$s2.post)*SE_table$stdev.unscaled
      
      results = limma::topTable(eBGene,coef=2, number = Inf, sort.by = "P")
      results$genes=row.names(results)
      
      results <- merge(results, SE_table, by="row.names", all.x=T)
      rownames(results) <- results$Row.names
      results <- results[,2:ncol(results)]
      
      p1=apply(cbind(results$P.Value,10^-320),1,max)
      z=qnorm((p1/2),low=F)*sign(results$logFC)
      results$z=z
      results[is.na(results$z),"z"] <- 0
      
      bc<-bacon(results$z)
      p_new=as.numeric(pval(bc))
      results$p_bacon=p_new
      results$p_bacon_adj=p.adjust(results$p_bacon, method = "fdr", n = nrow(results))
      
      lambda <- inflation(results$P.Value)
      lambda_bacon <- inflation(results$p_bacon)
      cat("lambda = ", lambda , "\n")
      cat("lambda_b = ", lambda_bacon, "\n")
      
      results$lambda <- lambda
      results$lambda_bacon <- lambda_bacon
      
      output_file <- paste0(output_dir, "/DEGs_", group_name, "_", class[i], "_", subdir_list[[class[i]]][j] )
      
      saveRDS(results, file = paste0(output_file, ".RDS"))
      filename <- paste0(output_file, ".txt")
      write.table(results, filename, sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE, append = FALSE)
      
      limma_breakdown <- data.frame(
        n_ctrl = length(which(pheno$gDx == "Control")),
        n_case = length(which(pheno$gDx == "MDD")),
        n_genes = nrow(results),
        n_svs = num_svs,
        n_nominal = nrow(results[(results$P.Value < 0.05) & (results$adj.P.Val > 0.05), ]),
        n_fdr = nrow(results[(results$adj.P.Val < 0.05), ]),
        n_bacon = nrow(results[(results$p_bacon < 0.05), ]),
        n_bacon_adj = nrow(results[(results$p_bacon_adj < 0.05), ]),
        lambda = results$lambda[1],
        lambda_bacon = results$lambda_bacon[1]
      )
      write.csv(
        limma_breakdown,
        paste0(output_file, "_breakdown.csv"),
        quote = FALSE, row.names = FALSE
      )
      
      plot_volcano(
        results,
        dx = group_name,
        contrast = group_name,
        region = class[i],
        file = paste0(output_file, "_volcano.pdf")
      )
      
      cat(group_name, " - ", class[i], "/", subdir_list[[class[i]]][j] , " done!", "\n")
    }}
  sink()
}

DGE_limma_sub <- function(group_name, covariates, output_dir, min.count = min.count, min.prop = min.prop, sva=F){
  sink(file = paste0(output_dir, "log_", group_name, ".txt"), append = TRUE)
  
  for(i in 1:length(class)){
    for(j in 1:length(subdir_list[[class[i]]])){
      cat("Starting with ", class[i],"/", subdir_list[[class[i]]][j], "\n")
      pg=paste0(pwd11,class[i],"/", subdir_list[[class[i]]][j] )
      
      
      pheno <- readRDS(paste0("/data/humgen/daskalakislab/aiatrou/BBRF/Pheno_Demo/Pseudobulk_Dec/phenofile_subclass_", class[i], ".rds"))
      
      allfiles <- list.files(path = pg)
      dge <- readDGE(allfiles, path=pg,columns = c(1,2))
      dge <- DGEList(dge)
      
      # Calculate the library size (column-wise sum of counts)
      library_sizes <- colSums(dge$counts)
      
      # Exclude samples with library size of zero
      dge <- dge[, library_sizes > 0]
      
      if (ncol(dge) == 0) {
        cat("Skipping ", class[i], " ", subdir_list[[class[i]]][j], " - no valid samples left after filtering.\n")
        next
      }
      
      dge$samples$group <- as.factor(pheno$gDx[match(rownames(dge$samples), pheno$sampleID)])
      keep <- filterByExpr(dge, min.count = min.count, min.prop = min.prop, group = dge$samples$group)
      
      dge <- dge[keep, keep.lib.sizes=FALSE]
      cat("Data dimensions after filtering: ", dim(dge), "\n")
      
      sample_names=colnames(dge)
      pheno <- pheno[match(sample_names,pheno$sampleID),]
      nsamples <- ncol(dge)
      dge <- calcNormFactors(dge, method = "TMM") 
      
      
      
      if(sva == T){
        mm <- model.matrix(~ as.factor(gDx) + as.factor(sex) + as.numeric(age) +
                             as.factor(version)+ as.factor(tissue), pheno)
        mm0 <- model.matrix(~ as.factor(sex) +age +
                              as.factor(version)+ as.factor(tissue), pheno)
        fit <- svaseq(cpm(dge), mod=mm, mod0=mm0)
        SVs <- fit$sv
        num_svs <- ncol(SVs)
        
        colnames(SVs) <- paste0("sv", seq_len(num_svs)) 
        
        pheno <- cbind(pheno, SVs)
        
        covariates_all <- c(covariates, paste0("sv", seq_len(num_svs)))
        
      } else {
        covariates_all <- c(covariates)
        mod = model.matrix(formula(paste0("~", paste(covariates_all, collapse = "+"))), data = pheno)
      }
      
      #
      form <- as.formula(paste("~", paste(covariates_all, collapse = " + ")))
      C <- canCorPairs(form, pheno)
      
      pdf(paste0(output_dir, "corrmat_", class[i],"_", subdir_list[[class[i]]][j] , ".pdf"))
      plotCorrMatrix(C)
      dev.off()
      
      #
      
      #Find SVs that correlate highly (> 0.7) with the covariates of interest
      high_corr_svs <- c()
      
      for (cov in covariates) {
        sv_corr <- C[cov, grep("^sv", colnames(C))]
        
        # Find SVs with correlation > 0.7
        correlated_svs <- names(sv_corr[abs(sv_corr) > 0.7])
        
        # Store the names of these highly correlated SVs
        high_corr_svs <- c(high_corr_svs, correlated_svs)
      }
      
      covariates_filtered <- covariates_all[!covariates_all %in% high_corr_svs]
      #
      
      mod = model.matrix(formula(paste0("~", paste(covariates_filtered, collapse = "+"))), data = pheno)
      
      #
      
      colnames(mod) <- make.names(colnames(mod))
      colnames(mod)[1] = "Int"
      
      vGene <- voom(dge,mod, plot = F)
      
      fitGene = lmFit(vGene, mod)
      eBGene = eBayes(fitGene)
      SE_table <- data.frame(s2.post=eBGene$s2.post, stdev.unscaled=eBGene$stdev.unscaled[,2])
      names(SE_table) <- c("s2.post","stdev.unscaled")
      SE_table$SE <- sqrt(SE_table$s2.post)*SE_table$stdev.unscaled
      
      results = limma::topTable(eBGene,coef=2, number = Inf, sort.by = "P")
      results$genes=row.names(results)
      
      results <- merge(results, SE_table, by="row.names", all.x=T)
      rownames(results) <- results$Row.names
      results <- results[,2:ncol(results)]
      
      p1=apply(cbind(results$P.Value,10^-320),1,max)
      z=qnorm((p1/2),low=F)*sign(results$logFC)
      results$z=z
      results[is.na(results$z),"z"] <- 0
      
      bc<-bacon(results$z)
      p_new=as.numeric(pval(bc))
      results$p_bacon=p_new
      results$p_bacon_adj=p.adjust(results$p_bacon, method = "fdr", n = nrow(results))
      
      lambda <- inflation(results$P.Value)
      lambda_bacon <- inflation(results$p_bacon)
      cat("lambda = ", lambda , "\n")
      cat("lambda_b = ", lambda_bacon, "\n")
      
      results$lambda <- lambda
      results$lambda_bacon <- lambda_bacon
      
      output_file <- paste0(output_dir, "/DEGs_", group_name, "_", class[i], "_", subdir_list[[class[i]]][j] )
      
      saveRDS(results, file = paste0(output_file, ".RDS"))
      filename <- paste0(output_file, ".txt")
      write.table(results, filename, sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE, append = FALSE)
      
      limma_breakdown <- data.frame(
        n_ctrl = length(which(pheno$gDx == "Control")),
        n_case = length(which(pheno$gDx == "MDD")),
        n_genes = nrow(results),
        n_svs = num_svs,
        n_nominal = nrow(results[(results$P.Value < 0.05) & (results$adj.P.Val > 0.05), ]),
        n_fdr = nrow(results[(results$adj.P.Val < 0.05), ]),
        n_bacon = nrow(results[(results$p_bacon < 0.05), ]),
        n_bacon_adj = nrow(results[(results$p_bacon_adj < 0.05), ]),
        lambda = results$lambda[1],
        lambda_bacon = results$lambda_bacon[1]
      )
      write.csv(
        limma_breakdown,
        paste0(output_file, "_breakdown.csv"),
        quote = FALSE, row.names = FALSE
      )
      
      plot_volcano(
        results,
        dx = group_name,
        contrast = group_name,
        region = class[i],
        file = paste0(output_file, "_volcano.pdf")
      )
      
      cat(group_name, " - ", class[i], "/", subdir_list[[class[i]]][j] , " done!", "\n")
    }}
  sink()
}

DGE_limma_sex_specific_subset <- function(group_name, covariates, output_dir, min.count = 10, min.prop = 0.7, sex = c("male", "female"), sva=F) {
  sink(file = paste0(output_dir, "log_sex_subset_", group_name, "_", sex, ".txt"), append = TRUE)
  
  for (i in 1:length(class)) {
    for (j in 1:length(subdir_list[[class[i]]])) {
      cat("Starting with ", class[i], "/", subdir_list[[class[i]]][j], " for sex: ", sex, "\n")
      pg <- paste0(pwd11, class[i], "/", subdir_list[[class[i]]][j])
      
      pheno <- readRDS(paste0("/data/humgen/daskalakislab/aiatrou/BBRF/Pheno_Demo/Pseudobulk_Dec/phenofile_class_", class[i], ".rds"))
      
      # Filter phenotype data for the specified sex
      pheno <- pheno[pheno$sex == sex, ]
      
      allfiles <- list.files(path = pg)
      dge <- readDGE(allfiles, path = pg, columns = c(1, 2))
      dge <- DGEList(dge)
      
      # Calculate the library size (column-wise sum of counts)
      library_sizes <- colSums(dge$counts)
      
      # Exclude samples with library size of zero
      dge <- dge[, library_sizes > 0]
      
      # Subset the DGE object to include only samples in the filtered phenotype data
      dge <- dge[, colnames(dge) %in% pheno$sampleID]
      
      if (ncol(dge) == 0) {
        cat("Skipping ", class[i], " ", subdir_list[[class[i]]][j], " - no valid samples left after filtering.\n")
        next
      }
      
      dge$samples$group <- as.factor(pheno$gDx[match(rownames(dge$samples), pheno$sampleID)])
      keep <- filterByExpr(dge, min.count = min.count, min.prop = min.prop, group = dge$samples$group)
      
      dge <- dge[keep, keep.lib.sizes = FALSE]
      cat("Data dimensions after filtering: ", dim(dge), "\n")
      
      sample_names <- colnames(dge)
      pheno <- pheno[match(sample_names, pheno$sampleID), ]
      dge <- calcNormFactors(dge, method = "TMM")
      
      if(sva == T){
        mm <- model.matrix(~ as.factor(gDx) + as.numeric(age) +
                             as.factor(version)+ as.factor(tissue), pheno)
        mm0 <- model.matrix(~ as.numeric(age) + as.factor(version)+ as.factor(tissue), pheno)
        fit <- svaseq(cpm(dge), mod=mm, mod0=mm0)
        SVs <- fit$sv
        num_svs <- ncol(SVs)
        
        colnames(SVs) <- paste0("sv", seq_len(num_svs)) 
        
        pheno <- cbind(pheno, SVs)
        
        covariates_all <- c(covariates, paste0("sv", seq_len(num_svs)))
        
      } else {
        covariates_all <- c(covariates)
        mod = model.matrix(formula(paste0("~", paste(covariates_all, collapse = "+"))), data = pheno)
      }
      
      
      #
      form <- as.formula(paste("~", paste(covariates_all, collapse = " + ")))
      C <- canCorPairs(form, pheno)
      
      pdf(paste0(output_dir, "corrmat_", class[i],"_", subdir_list[[class[i]]][j] , ".pdf"))
      plotCorrMatrix(C)
      dev.off()
      
      #
      
      #Find SVs that correlate highly (> 0.7) with the covariates of interest
      high_corr_svs <- c()
      
      for (cov in covariates) {
        sv_corr <- C[cov, grep("^sv", colnames(C))]
        
        # Find SVs with correlation > 0.7
        correlated_svs <- names(sv_corr[abs(sv_corr) > 0.7])
        
        # Store the names of these highly correlated SVs
        high_corr_svs <- c(high_corr_svs, correlated_svs)
      }
      
      covariates_filtered <- covariates_all[!covariates_all %in% high_corr_svs]
      #
      
      mod = model.matrix(formula(paste0("~", paste(covariates_filtered, collapse = "+"))), data = pheno)
      
      #
      
      colnames(mod) <- make.names(colnames(mod))
      colnames(mod)[1] = "Int"
      
      vGene <- voom(dge,mod, plot = F)
      
      fitGene = lmFit(vGene, mod)
      eBGene = eBayes(fitGene)
      SE_table <- data.frame(s2.post=eBGene$s2.post, stdev.unscaled=eBGene$stdev.unscaled[,2])
      names(SE_table) <- c("s2.post","stdev.unscaled")
      SE_table$SE <- sqrt(SE_table$s2.post)*SE_table$stdev.unscaled
      
      results = limma::topTable(eBGene,coef=2, number = Inf, sort.by = "P")
      results$genes=row.names(results)
      
      results <- merge(results, SE_table, by="row.names", all.x=T)
      rownames(results) <- results$Row.names
      results <- results[,2:ncol(results)]
      
      p1=apply(cbind(results$P.Value,10^-320),1,max)
      z=qnorm((p1/2),low=F)*sign(results$logFC)
      results$z=z
      results[is.na(results$z),"z"] <- 0
      
      bc<-bacon(results$z)
      p_new=as.numeric(pval(bc))
      results$p_bacon=p_new
      results$p_bacon_adj=p.adjust(results$p_bacon, method = "fdr", n = nrow(results))
      
      lambda <- inflation(results$P.Value)
      lambda_bacon <- inflation(results$p_bacon)
      cat("lambda = ", lambda , "\n")
      cat("lambda_b = ", lambda_bacon, "\n")
      
      results$lambda <- lambda
      results$lambda_bacon <- lambda_bacon
      
      # Save results
      output_file <- paste0(output_dir, "/SexSubset_DEGs_", group_name, "_", class[i], "_", subdir_list[[class[i]]][j], "_", sex)
      saveRDS(results, file = paste0(output_file, ".RDS"))
      write.table(results, paste0(output_file, ".txt"), sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE, append = FALSE)
      
      limma_breakdown <- data.frame(
        n_ctrl = length(which(pheno$gDx == "Control")),
        n_case = length(which(pheno$gDx == "MDD")),
        n_genes = nrow(results),
        n_svs = num_svs,
        n_nominal = nrow(results[(results$P.Value < 0.05) & (results$adj.P.Val > 0.05), ]),
        n_fdr = nrow(results[(results$adj.P.Val < 0.05), ]),
        n_bacon = nrow(results[(results$p_bacon < 0.05), ]),
        n_bacon_adj = nrow(results[(results$p_bacon_adj < 0.05), ]),
        lambda = results$lambda[1],
        lambda_bacon = results$lambda_bacon[1]
      )
      write.csv(
        limma_breakdown,
        paste0(output_file, "_breakdown.csv"),
        quote = FALSE, row.names = FALSE
      )
      
      plot_volcano(
        results,
        dx = group_name,
        contrast = group_name,
        region = class[i],
        file = paste0(output_file, "_volcano.pdf")
      )
      
      cat(group_name, " - ", class[i], "/", subdir_list[[class[i]]][j], " for sex: ", sex, " done!", "\n")
    }}
  sink()
  
}

run_RE_limma_analysis <- function(group, output_dir, covariates, celltype ) {
  # Subset pheno_mod
  pheno_modX <- subset(pheno_mod, pheno_mod$gDx %in% c("Control", group))
  #pheno_modX$gDx <- droplevels(pheno_modX$gDx)
  dge <- dge[, colnames(dge) %in% pheno_modX$SampleID_cl]
  pheno_modX <- pheno_modX[pheno_modX$SampleID_cl%in% colnames(dge), ]
  
  mod = model.matrix(formula(paste0("~", paste0(covariates[1:2], collapse = "*"), "+", paste(covariates[3:length(covariates)], collapse = "+"))), data = pheno_modX)
  
  colnames(mod) <- make.names(colnames(mod))
  colnames(mod)[1] = "Int"
  
  vGene_temp = voom(dge,mod, plot=F)
  dupCorr = duplicateCorrelation(vGene_temp$E, mod, block=pheno_modX$sampleID)
  print(paste("dupCorr2 done"))
  
  vGene = voom(dge, mod, plot=FALSE, block=pheno_modX$sampleID, correlation=dupCorr$consensus)
  dupcor <- duplicateCorrelation(vGene, mod, block=pheno_modX$sampleID)
  print(paste("dupCorr2 done"))
  
  fitGene <- lmFit(vGene, mod, block=pheno_modX$sampleID, correlation=dupcor$consensus)
  
  eBGene = eBayes(fitGene)
  
  SE_table <- data.frame(s2.post=eBGene$s2.post, stdev.unscaled=eBGene$stdev.unscaled[,2])
  names(SE_table) <- c("s2.post","stdev.unscaled")
  SE_table$SE <- sqrt(SE_table$s2.post)*SE_table$stdev.unscaled
  
  results = limma::topTable(eBGene,coef=2, number = Inf, sort.by = "P")
  results$genes=row.names(results)
  
  results <- merge(results, SE_table, by="row.names", all.x=T)
  rownames(results) <- results$Row.names
  results <- results[,2:ncol(results)]
  
  p1=apply(cbind(results$P.Value,10^-320),1,max)
  z=qnorm((p1/2),low=F)*sign(results$logFC)
  results$z=z
  results[is.na(results$z),"z"] <- 0
  
  bc<-bacon(results$z)
  #bc<-bacon(z,na.exclude = T)
  p_new=as.numeric(pval(bc))
  results$p_bacon=p_new
  results$p_bacon_adj=p.adjust(results$p_bacon, method = "fdr", n = nrow(results))
  
  lambda <- inflation(results$P.Value)
  lambda_bacon <- inflation(results$p_bacon)
  
  results$lambda <- lambda
  results$lambda_bacon <- lambda_bacon
  
  # Save results
  saveRDS(results, file = paste0(output_dir, "/DEGs_", group, "_", celltype,"_RE.RDS"))
  filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_RE.txt")
  write.table(results, filename, sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE, append = FALSE)
  print("DGE done")
  
  # Create limma breakdown data frame and write to CSV
  limma_breakdown <- data.frame(
    n_genes = nrow(results),
    n_nominal = nrow(results[(results$P.Value < 0.05) & (results$adj.P.Val > 0.05),]),
    n_fdr = nrow(results[(results$adj.P.Val < 0.05),]),
    n_bacon = nrow(results[(results$p_bacon < 0.05),]),
    n_bacon_adj = nrow(results[(results$p_bacon_adj < 0.05),]),
    lambda = results$lambda[1], lambda_bacon = results$lambda_bacon[1]
  )
  breakdown_filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_RE_breakdown.csv")
  write.csv(limma_breakdown, breakdown_filename, quote = FALSE, row.names = FALSE)
  
  # Plot volcano plot
  volcano_filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_RE_volcano.pdf")
  plot_volcano(results, dx = group, contrast = group, region = celltype, file = volcano_filename)
  print("Volcano done")
}

run_RE_limma_analysis_sex_specific <- function(group, output_dir, covariates, celltype, sex = c("male", "female")) {
  # Subset pheno_mod for the specified group and sex
  pheno_modX <- subset(pheno_mod, gDx %in% c("Control", group))
  
  subset <- which(pheno_modX$sex == sex)
  pheno_modX <- pheno_modX[subset, ]
  
  # Subset DGE to include only samples in pheno_modX
  dge <- dge[, colnames(dge) %in% pheno_modX$SampleID_cl]
  pheno_modX <- pheno_modX[pheno_modX$SampleID_cl %in% colnames(dge), ]
  
  if (nrow(pheno_modX) == 0) {
    cat("No samples available for sex:", sex, "and group:", group, ". Skipping.\n")
    return(NULL)
  }
  
  # Create model matrix with specified covariates
  mod <- model.matrix(
    formula(paste0("~", paste0(covariates[1:2], collapse = "*"), "+", paste(covariates[3:length(covariates)], collapse = "+"))),
    data = pheno_modX
  )
  colnames(mod) <- make.names(colnames(mod))
  colnames(mod)[1] <- "Int"
  
  # Perform voom transformation and duplicate correlation estimation
  vGene_temp <- voom(dge, mod, plot = FALSE)
  dupCorr <- duplicateCorrelation(vGene_temp$E, mod, block = pheno_modX$sampleID)
  print(paste("dupCorr1 done"))
  
  vGene <- voom(dge, mod, plot = FALSE, block = pheno_modX$sampleID, correlation = dupCorr$consensus)
  dupcor <- duplicateCorrelation(vGene, mod, block = pheno_modX$sampleID)
  print(paste("dupCorr2 done"))
  
  # Fit linear model with duplicate correlation
  fitGene <- lmFit(vGene, mod, block = pheno_modX$sampleID, correlation = dupcor$consensus)
  eBGene <- eBayes(fitGene)
  
  # Extract results
  SE_table <- data.frame(
    s2.post = eBGene$s2.post,
    stdev.unscaled = eBGene$stdev.unscaled[, 2]
  )
  SE_table$SE <- sqrt(SE_table$s2.post) * SE_table$stdev.unscaled
  
  results <- limma::topTable(eBGene, coef = 2, number = Inf, sort.by = "P")
  results$genes <- rownames(results)
  results <- merge(results, SE_table, by = "row.names", all.x = TRUE)
  rownames(results) <- results$Row.names
  results <- results[, -1]
  
  # Compute z-scores and apply BACON for inflation correction
  results$z <- qnorm(results$P.Value / 2, lower.tail = FALSE) * sign(results$logFC)
  results$z[is.na(results$z)] <- 0
  
  bc <- bacon(results$z)
  results$p_bacon <- as.numeric(pval(bc))
  results$p_bacon_adj <- p.adjust(results$p_bacon, method = "fdr")
  
  results$lambda <- inflation(results$P.Value)
  results$lambda_bacon <- inflation(results$p_bacon)
  
  # Save results
  prefix <- paste0(output_dir, "/DEGs_", group, "_", celltype, "_", sex, "_RE")
  saveRDS(results, file = paste0(prefix, ".RDS"))
  write.table(results, paste0(prefix, ".txt"), sep = " ", quote = FALSE, row.names = FALSE)
  
  # Create and save breakdown summary
  limma_breakdown <- data.frame(
    n_genes = nrow(results),
    n_nominal = sum(results$P.Value < 0.05 & results$adj.P.Val > 0.05, na.rm = TRUE),
    n_fdr = sum(results$adj.P.Val < 0.05, na.rm = TRUE),
    n_bacon = sum(results$p_bacon < 0.05, na.rm = TRUE),
    n_bacon_adj = sum(results$p_bacon_adj < 0.05, na.rm = TRUE),
    lambda = results$lambda[1],
    lambda_bacon = results$lambda_bacon[1]
  )
  write.csv(limma_breakdown, paste0(prefix, "_breakdown.csv"), quote = FALSE, row.names = FALSE)
  
  # Plot volcano
  volcano_filename <- paste0(prefix, "_volcano.pdf")
  plot_volcano(results, dx = group, contrast = group, region = celltype, file = volcano_filename)
  
  print(paste("Analysis completed for sex:", sex, "and celltype:", celltype))
}

run_RE_CT_limma_analysis <- function(group, output_dir, covariates, celltype ) {
  # Subset pheno_mod
  pheno_modX <- pheno_mod
  pheno_modX$CT <- as.numeric(as.character(pheno_modX$CT))
  
  mod = model.matrix(formula(paste0("~", paste0(covariates[1:2], collapse = "*"), "+", paste(covariates[3:length(covariates)], collapse = "+"))), data = pheno_modX)
  
  colnames(mod) <- make.names(colnames(mod))
  colnames(mod)[1] = "Int"
  
  vGene_temp = voom(dge,mod, plot=F)
  dupCorr = duplicateCorrelation(vGene_temp$E, mod, block=pheno_modX$sampleID)
  print(paste("dupCorr1 done"))
  
  vGene = voom(dge, mod, plot=FALSE, block=pheno_modX$sampleID, correlation=dupCorr$consensus)
  dupcor <- duplicateCorrelation(vGene, mod, block=pheno_modX$sampleID)
  print(paste("dupCorr2 done"))
  
  fitGene <- lmFit(vGene, mod, block=pheno_modX$sampleID, correlation=dupcor$consensus)
  
  eBGene = eBayes(fitGene)
  
  SE_table <- data.frame(s2.post=eBGene$s2.post, stdev.unscaled=eBGene$stdev.unscaled[,2])
  names(SE_table) <- c("s2.post","stdev.unscaled")
  SE_table$SE <- sqrt(SE_table$s2.post)*SE_table$stdev.unscaled
  
  results = limma::topTable(eBGene,coef=2, number = Inf, sort.by = "P")
  results$genes=row.names(results)
  
  results <- merge(results, SE_table, by="row.names", all.x=T)
  rownames(results) <- results$Row.names
  results <- results[,2:ncol(results)]
  
  p1=apply(cbind(results$P.Value,10^-320),1,max)
  z=qnorm((p1/2),low=F)*sign(results$logFC)
  results$z=z
  results[is.na(results$z),"z"] <- 0
  
  bc<-bacon(results$z)
  #bc<-bacon(z,na.exclude = T)
  p_new=as.numeric(pval(bc))
  results$p_bacon=p_new
  results$p_bacon_adj=p.adjust(results$p_bacon, method = "fdr", n = nrow(results))
  
  lambda <- inflation(results$P.Value)
  lambda_bacon <- inflation(results$p_bacon)
  
  results$lambda <- lambda
  results$lambda_bacon <- lambda_bacon
  
  # Save results
  saveRDS(results, file = paste0(output_dir, "/DEGs_", group, "_", celltype,"_RE.RDS"))
  filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_RE.txt")
  write.table(results, filename, sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE, append = FALSE)
  print("DGE done")
  
  # Create limma breakdown data frame and write to CSV
  limma_breakdown <- data.frame(
    n_genes = nrow(results),
    n_nominal = nrow(results[(results$P.Value < 0.05) & (results$adj.P.Val > 0.05),]),
    n_fdr = nrow(results[(results$adj.P.Val < 0.05),]),
    n_bacon = nrow(results[(results$p_bacon < 0.05),]),
    n_bacon_adj = nrow(results[(results$p_bacon_adj < 0.05),]),
    lambda = results$lambda[1], lambda_bacon = results$lambda_bacon[1]
  )
  breakdown_filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_RE_breakdown.csv")
  write.csv(limma_breakdown, breakdown_filename, quote = FALSE, row.names = FALSE)
  
  # Plot volcano plot
  volcano_filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_RE_volcano.pdf")
  plot_volcano(results, dx = group, contrast = group, region = celltype, file = volcano_filename)
  print("Volcano done")
}


run_RE_CT_contrast_limma_analysis <- function(group, output_dir, covariates, celltype ) {
  # Subset pheno_mod
  pheno_modX <- pheno_mod
  pheno_modX$CT <- factor(pheno_modX$CT, levels = c(2, 0, 1))  
  
  mod = model.matrix(formula(paste0("~", paste0(covariates[1:2], collapse = "*"), "+", paste(covariates[3:length(covariates)], collapse = "+"))), data = pheno_modX)
  
  colnames(mod) <- make.names(colnames(mod))
  colnames(mod)[1] = "Int"
  
  vGene_temp = voom(dge,mod, plot=F)
  dupCorr = duplicateCorrelation(vGene_temp$E, mod, block=pheno_modX$sampleID)
  print(paste("dupCorr1 done"))
  
  vGene = voom(dge, mod, plot=FALSE, block=pheno_modX$sampleID, correlation=dupCorr$consensus)
  dupcor <- duplicateCorrelation(vGene, mod, block=pheno_modX$sampleID)
  print(paste("dupCorr2 done"))
  
  fitGene <- lmFit(vGene, mod, block=pheno_modX$sampleID, correlation=dupcor$consensus)
  
  eBGene = eBayes(fitGene)
  
  #CT0
  SE_table <- data.frame(s2.post=eBGene$s2.post, stdev.unscaled=eBGene$stdev.unscaled[,"CT0"])
  names(SE_table) <- c("s2.post","stdev.unscaled")
  SE_table$SE <- sqrt(SE_table$s2.post)*SE_table$stdev.unscaled
  
  results = limma::topTable(eBGene,coef="CT0", number = Inf, sort.by = "P")
  results$genes=row.names(results)
  
  results <- merge(results, SE_table, by="row.names", all.x=T)
  rownames(results) <- results$Row.names
  results <- results[,2:ncol(results)]
  
  p1=apply(cbind(results$P.Value,10^-320),1,max)
  z=qnorm((p1/2),low=F)*sign(results$logFC)
  results$z=z
  results[is.na(results$z),"z"] <- 0
  
  bc<-bacon(results$z)
  p_new=as.numeric(pval(bc))
  results$p_bacon=p_new
  results$p_bacon_adj=p.adjust(results$p_bacon, method = "fdr", n = nrow(results))
  
  lambda <- inflation(results$P.Value)
  lambda_bacon <- inflation(results$p_bacon)
  
  results$lambda <- lambda
  results$lambda_bacon <- lambda_bacon
  
  results$logFC <- -results$logFC
  
  # Save results
  saveRDS(results, file = paste0(output_dir, "/DEGs_", group, "_", celltype,"_CT0_RE.RDS"))
  filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_CT0_RE.txt")
  write.table(results, filename, sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE, append = FALSE)
  print("DGE done")
  
  # Create limma breakdown data frame and write to CSV
  limma_breakdown <- data.frame(
    n_genes = nrow(results),
    n_nominal = nrow(results[(results$P.Value < 0.05) & (results$adj.P.Val > 0.05),]),
    n_fdr = nrow(results[(results$adj.P.Val < 0.05),]),
    n_bacon = nrow(results[(results$p_bacon < 0.05),]),
    n_bacon_adj = nrow(results[(results$p_bacon_adj < 0.05),]),
    lambda = results$lambda[1], lambda_bacon = results$lambda_bacon[1]
  )
  breakdown_filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_CT0_RE_breakdown.csv")
  write.csv(limma_breakdown, breakdown_filename, quote = FALSE, row.names = FALSE)
  
  # Plot volcano plot
  volcano_filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_CT0_RE_volcano.pdf")
  plot_volcano(results, dx = group, contrast = group, region = celltype, file = volcano_filename)
  print("Volcano done")
  
  #CT1
  SE_table <- data.frame(s2.post=eBGene$s2.post, stdev.unscaled=eBGene$stdev.unscaled[,"CT1"])
  names(SE_table) <- c("s2.post","stdev.unscaled")
  SE_table$SE <- sqrt(SE_table$s2.post)*SE_table$stdev.unscaled
  
  results = limma::topTable(eBGene,coef="CT1", number = Inf, sort.by = "P")
  results$genes=row.names(results)
  
  results <- merge(results, SE_table, by="row.names", all.x=T)
  rownames(results) <- results$Row.names
  results <- results[,2:ncol(results)]
  
  p1=apply(cbind(results$P.Value,10^-320),1,max)
  z=qnorm((p1/2),low=F)*sign(results$logFC)
  results$z=z
  results[is.na(results$z),"z"] <- 0
  
  bc<-bacon(results$z)
  p_new=as.numeric(pval(bc))
  results$p_bacon=p_new
  results$p_bacon_adj=p.adjust(results$p_bacon, method = "fdr", n = nrow(results))
  
  lambda <- inflation(results$P.Value)
  lambda_bacon <- inflation(results$p_bacon)
  
  results$lambda <- lambda
  results$lambda_bacon <- lambda_bacon
  
  results$logFC <- -results$logFC
  # Save results
  saveRDS(results, file = paste0(output_dir, "/DEGs_", group, "_", celltype,"_CT1_RE.RDS"))
  filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_CT1_RE.txt")
  write.table(results, filename, sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE, append = FALSE)
  print("DGE done")
  
  # Create limma breakdown data frame and write to CSV
  limma_breakdown <- data.frame(
    n_genes = nrow(results),
    n_nominal = nrow(results[(results$P.Value < 0.05) & (results$adj.P.Val > 0.05),]),
    n_fdr = nrow(results[(results$adj.P.Val < 0.05),]),
    n_bacon = nrow(results[(results$p_bacon < 0.05),]),
    n_bacon_adj = nrow(results[(results$p_bacon_adj < 0.05),]),
    lambda = results$lambda[1], lambda_bacon = results$lambda_bacon[1]
  )
  breakdown_filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_CT1_RE_breakdown.csv")
  write.csv(limma_breakdown, breakdown_filename, quote = FALSE, row.names = FALSE)
  
  # Plot volcano plot
  volcano_filename <- paste0(output_dir, "/DEGs_", group, "_", celltype,"_CT1_RE_volcano.pdf")
  plot_volcano(results, dx = group, contrast = group, region = celltype, file = volcano_filename)
  print("Volcano done")
}


########################################
# Enhanced Pathway Analysis for MDD
# Author: Artemis Iatrou
# Date: 2025-08-08
# Description:
#   This script reads differential expression pathway results
#   for multiple brain regions, identifies significant pathways,
#   computes shared and region-specific pathways relative to RE,
#   and performs semantic similarity clustering of GO terms.
########################################

# =====================
# 0. Load Required Packages
# =====================
library(tidyverse)
library(GO.db)
library(AnnotationDbi)
library(igraph)
library(CePa)

# =====================
# 1. Define Parameters
# =====================
brain_regions <- c("RE", "mPFC", "dACC", "DLPFC")
subcat <- "C5_BP"  # Pathway category
descriptive_info <- list()
data_dir <- "data/Bulk"   
output_dir <- "output/Pathways"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =====================
# 2. Load Pathway Data & Identify Significant Pathways
# =====================
for (region in brain_regions) {
  
  # Load region-specific pathway data
  df <- if (region == "RE") {
    readRDS(file.path(data_dir, "AcrossRegions", "Pathway", paste0(subcat, "_All.rds")))
  } else {
    read_csv(file.path(data_dir, region, "CombinedResults_svs", "all", paste0("gsea_", subcat, ".csv")),
             show_col_types = FALSE)
  }
  
  # Proceed only if data exists
  if (nrow(df) > 0) {
    
    df <- df %>%
      mutate(
        ss = str_extract(pathway, "(?<=GO).*?(?=_)"),      
        sig = ifelse(padj < 0.05, 1, 0),                  
        region = region,                                 
        NES_sign = if_else(NES > 0, "positive", "negative") 
      )
    
    if (region == "RE") {
      # Store RE significant pathways
      re_significant <- df %>% filter(sig == 1)
      top_re_pathways <- re_significant %>%
        mutate(lp = -log10(padj)) %>%
        group_by(NES_sign) %>%
        arrange(NES_sign, desc(lp)) %>%
        slice_head(n = 10) %>%
        pull(pathway)
      
      descriptive_info[[region]] <- list(
        all_data = df,
        re_significant = re_significant,
        re_pos = re_significant %>% filter(NES > 0),
        re_neg = re_significant %>% filter(NES < 0),
        top_pathways = top_re_pathways
      )
      
    } else if (!is.null(re_significant)) {
      # Compute shared, RE-specific, and region-specific significant pathways
      shared <- df %>% filter(pathway %in% re_significant$pathway, sig == 1)
      re_specific <- re_significant %>% filter(!pathway %in% df$pathway)
      region_specific <- df %>% filter(!pathway %in% re_significant$pathway, sig == 1)
      
      descriptive_info[[region]] <- list(
        all_data = df,
        shared_significant = shared,
        shared_positive = shared %>% filter(NES > 0),
        shared_negative = shared %>% filter(NES < 0),
        re_specific = re_specific,
        re_specific_positive = re_specific %>% filter(NES > 0),
        re_specific_negative = re_specific %>% filter(NES < 0),
        region_specific = region_specific,
        region_specific_positive = region_specific %>% filter(NES > 0),
        region_specific_negative = region_specific %>% filter(NES < 0)
      )
    }
  }
}

# =====================
# 3. Process Region-Specific Pathways
# =====================
clean_pathways <- function(paths) {
  paths %>%
    sub("GOBP_", "", .) %>%
    gsub("_", " ", .) %>%
    tolower()
}

region_specific_terms <- lapply(descriptive_info, function(info) {
  if (!is.null(info$region_specific)) {
    cleaned <- clean_pathways(info$region_specific$pathway)
    matched_terms <- Term(GOTERM)[match(cleaned, Term(GOTERM))]
    matched_terms[complete.cases(matched_terms)]
  } else {
    NULL
  }
})

# =====================
# 4. Semantic Similarity Clustering
# =====================
process_go_similarity <- function(go_terms) {
  go_ids <- names(go_terms)
  sm <- GO_similarity(go_ids, measure = "Rel", ont = "BP")
  diag(sm) <- 1
  sm[rowSums(sm > 0) > 1, colSums(sm > 0) > 1]
}

semantic_similarity <- lapply(region_specific_terms, process_go_similarity)

# Extract hierarchical GO relationships
lt <- as.list(GOBPCHILDREN)
df_edges <- do.call(rbind, lapply(names(lt), function(nm) data.frame(from = nm, to = lt[[nm]])))
df_edges <- df_edges[!is.na(df_edges$to), ]
g <- graph.edgelist(as.matrix(df_edges), directed = TRUE)

# Cluster GO terms per region
go_clusters <- lapply(semantic_similarity, function(sm) {
  cn <- intersect(rownames(sm), unlist(df_edges))
  sm2 <- sm[cn, cn]
  g_sub <- induced_subgraph(g, cn)
  cl <- cluster_terms(sm2, method = "louvain")
  gl <- tapply(rownames(sm2), cl, function(go_id) induced_subgraph(g_sub, go_id))
  lapply(gl, function(g) {
    membership <- components(g)$membership
    induced_subgraph(g, membership == which.max(table(membership)))
  })
})

# Extract meaningful GO terms
go_list <- lapply(go_clusters, function(gl2) {
  lapply(gl2, function(g) {
    s1 <- degree(g, mode = "in")
    s2 <- spread(g, mode = "out")
    go_id <- V(g)$name[which(s1 == 0)]
    s <- s2[s1 == 0]
    go_id <- go_id[s >= mean(s)]
    term <- suppressMessages(AnnotationDbi::select(GO.db, keys = go_id, columns = "TERM", keytype = "GOID")$TERM)
    paste0(go_id, ": ", term)
  })
})

# =====================
# 5. Save Plots for Each Region
# =====================
pdf_files <- paste0(output_dir, "/GO_clusters_smaller_", names(go_clusters), ".pdf")

for (i in seq_along(go_clusters)) {
  pdf(pdf_files[i], width = 2.6, height = 1.3)
  ht_clusters(semantic_similarity[[i]],
              cl = cluster_terms(semantic_similarity[[i]], method = "louvain"),
              exclude_words = c("process", "regulation", "response"),
              draw_word_cloud = FALSE)
  dev.off()
}

print("GO term cluster plots saved as PDFs.")
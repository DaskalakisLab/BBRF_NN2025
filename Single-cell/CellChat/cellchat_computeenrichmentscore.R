computeEnrichmentScore_mod <- function (df, measure = c("ligand", "signaling", "LR-pair"), 
          species = c("mouse", "human"), color.use = NULL, color.name = "Dark2", 
          n.color = 8, scale = c(4, 0.8), min.freq = 0, max.words = 200, 
          random.order = FALSE, rot.per = 0, return.data = FALSE, seed = 1, 
          ...) 
{measure <- match.arg(measure)
  species <- match.arg(species)
  LRpairs <- as.character(unique(df$interaction_name))
  ES <- vector(length = length(LRpairs))
  for (i in 1:length(LRpairs)) {
    df.i <- subset(df, interaction_name == LRpairs[i])
    if (length(which(rowSums(is.na(df.i)) > 0)) > 0) {
      df.i <- df.i[-which(rowSums(is.na(df.i)) > 0), , 
                   drop = FALSE]
    }
    ES[i] = mean(abs(df.i$ligand.logFC) * abs(df.i$receptor.logFC) * 
                   abs(df.i$ligand.pct.2 - df.i$ligand.pct.1) * abs(df.i$receptor.pct.2 - 
                                                                      df.i$receptor.pct.1))
  }
  if (species == "mouse") {
    CellChatDB <- CellChatDB.mouse
  }
  else if (species == "human") {
    CellChatDB <- CellChatDB.human
  }
  df.es <- CellChatDB$interaction[LRpairs, c("ligand", "receptor", 
                                             "pathway_name")]
  df.es$score <- ES
  df.es <- df.es[complete.cases(df.es$score),]
  
  df.es.ensemble <- df.es %>% group_by(ligand) %>% summarize(total = sum(score))
  set.seed(seed)
  if (is.null(color.use)) {
    color.use <- RColorBrewer::brewer.pal(n.color, color.name)
  }
  wordcloud::wordcloud(words = df.es.ensemble$ligand, freq = df.es.ensemble$total, 
                       min.freq = min.freq, max.words = max.words, scale = scale, 
                       random.order = random.order, rot.per = rot.per, colors = color.use)
  if (return.data) {
    return(df.es.ensemble)
  }
}

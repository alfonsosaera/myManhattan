myManhattan <- function(df, graph.title = "", highlight = NULL, highlight.col = "green",
                        col = c("lightblue", "blue"), even.facet = FALSE, chrom.lab = NULL,
                        suggestiveline = 1e-05, genomewideline = 5e-08, font.size = 12,
                        axis.size = 0.5, significance = NULL, report = FALSE){
  require(ggplot2)
  require(stats)
  y.title <- expression(-log[10](italic(p)))
  if (length(col) > length(unique(df$CHR))){
    chrom.col <- col[1:length(unique(df$CHR))]
  } else if (!(length(col) > length(unique(df$CHR)))){
    chrom.col <- rep(col, length(unique(df$CHR))/length(col))
    if (length(chrom.col) < length(unique(df$CHR))){
      dif <- length(unique(df$CHR)) - length(chrom.col)
      chrom.col <- c(chrom.col, col[1:dif])
    }
  }
  y.max <- floor(max(-log10(df$P))) + 1
  if (y.max %% 2 != 0){
    y.max <- y.max + 1
  }
  if (!is.null(chrom.lab)){
    if (length(unique(df$CHR)) != length(chrom.lab)){
      warning("Number of chrom.lab different of number of chromosomes in dataset, argument ignored.")
    } else {
      df$CHR <- factor(df$CHR, levels = unique(df$CHR), labels=chrom.lab)
    }
  }
  g <- ggplot(df) +
    geom_point(aes(BP, -log10(P), colour = as.factor(CHR)), size = 1)
  if (!is.null(significance)){
    if (is.numeric(significance)){
      genomewideline <- significance
      suggestiveline <- genomewideline / 0.005
    } else if (significance == "Bonferroni"){
      BFlevel <- 0.05 / length(df$SNP)
      cat("Bonferroni correction significance level:", BFlevel, "\n")
      genomewideline <- BFlevel
      suggestiveline <- BFlevel / 0.005
    } else if (significance == "FDR"){
      df$fdr <- p.adjust(df$P, "fdr")
      genomewideline <- 0.05
      suggestiveline <- FALSE
      y.title <- expression(-log[10](italic(q)))
      g <- ggplot(df) +
        geom_point(aes(BP, -log10(fdr), colour = as.factor(CHR)), size = 1)
      if (!is.null(highlight)) {
        if (is.numeric(highlight)){
          highlight <- as.character(df$SNP[df$P < highlight])
        }
        if (any(!(highlight %in% df$SNP))){
          warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
        } else {
          g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                              aes(BP, -log10(fdr), group=SNP, colour=SNP),
                              color = highlight.col, size = 1)
          highlight <- NULL
          y.max <- floor(max(-log10(df$fdr))) + 1
          if (y.max %% 2 != 0){
            y.max <- y.max + 1
          }
        }
      }
    }
  }
  if (even.facet){
    g <- g + facet_grid(.~CHR, scale = "free_x", switch = "x")
  } else {
    g <- g + facet_grid(.~CHR, scale = "free_x", space = "free_x", switch = "x")
  }
  g <- g + scale_colour_manual(values = chrom.col) +
    scale_y_continuous(expand = c(0, 0), limit = c(0, y.max),
                       breaks = seq(from = 0, to = y.max, by = 2)) +
    scale_x_continuous() +
    theme(strip.background = element_blank(), legend.position = "none",
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.x=unit(0.1, "lines"),
          axis.line.y = element_line(size = axis.size, color = "black"),
          axis.ticks.y = element_line(size = axis.size, color = "black"),
          axis.ticks.length = unit(axis.size * 10, "points"),
          plot.title = element_text(hjust = (0.5), size = font.size + 8),
          axis.title.y = element_text(size = font.size + 5),
          axis.title.x = element_text(size = font.size + 5),
          axis.text = element_text(size = font.size),
          strip.text.x = element_text(size = font.size))+
    labs(title = graph.title, x = "Chromosome", y = y.title)
  if (!is.null(highlight)) {
    if (is.numeric(highlight)){
      highlight <- as.character(df$SNP[df$P < highlight])
    }
    if (any(!(highlight %in% df$SNP))){
      warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
    } else {
      g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                          aes(BP, -log10(P), group=SNP, colour=SNP),
                          color = highlight.col, size = 1)
    }
  }
  if (suggestiveline){
    g <- g + geom_hline(yintercept = -log10(suggestiveline), color = "blue")
  }
  if (genomewideline){
    g <- g + geom_hline(yintercept = -log10(genomewideline), color = "red")
  }
  if (report){
    if (significance == "FDR"){
      rep <- df[df$fdr < 0.05, ]
    } else if (significance == "Bonferroni"){
      rep <- df[df$P < BFlevel, ]
    } else if (is.numeric(significance)){
      rep <- df[df$P < significance, ]
    } else {
      cat("using default significance level, 5e-8")
      rep <- df[df$P < 5e-8, ]
    }
    print(rep)
  }
  return(g)
}

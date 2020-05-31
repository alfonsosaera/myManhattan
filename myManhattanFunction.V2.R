# https://www.r-graph-gallery.com/101_Manhattan_plot.html

myManhattan <- function(df, graph.title = "", chrom.lab = NULL, 
                        col = c("lightblue", "navy"),
                        suggestiveline = 1e-05, suggestivecolor = "blue",
                        genomewideline = 5e-08, genomewidecolor = "red",
                        highlight = NULL, highlight.names = FALSE, 
                        highlight.col = "green", name.col = highlight.col,
                        significance = NULL, report = FALSE,
                        sep = 0.05, back.panels = TRUE, back.panels.col = "#ebebeb",
                        font.size = 12, axis.size = 0.5, inf.corr = 0.95, y.step = 2, point.size = 1){
  # libraries
  require(ggplot2)
  require(stats)
  require(dplyr)
  require(tidyr)
  require(ggrepel)
  myMin <- min(df$P[df$P != 0]) * inf.corr
  df$P[df$P == 0] <- myMin
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
  for (j in (1:y.step)){
    if (y.max %% y.step != 0){
      y.max <- y.max + 1
    } else {
      break
    }
  }
  separator <- as.integer(df %>%
                            group_by(CHR) %>%
                            summarise(CHR_size = max(BP)) %>% 
                            summarise(floor(max(CHR_size) * sep)))
  tmp.df <- df %>%
    group_by(CHR) %>%
    mutate(new_BP = BP - min(BP) + 1) %>%
    summarise(CHR_size = max(new_BP) + separator) %>%
    mutate(cum_size = cumsum(as.numeric(CHR_size )) - CHR_size)
  df <- df %>%
    group_by(CHR) %>%
    mutate(new_BP = BP - min(BP) + 1) %>%
    ungroup() %>%
    inner_join(tmp.df, by = "CHR") %>%
    mutate(plot_BP = new_BP + cum_size)
  if (highlight.names){ 
    if (is.numeric(highlight)){
      highlight <- as.character(df$SNP[df$P < highlight])
    }
    df$label <- ""
    df$label[which(df$SNP %in% highlight)] <- as.character(df$SNP[which(df$SNP %in% highlight)])
  }
  axis.info.df <- df %>%
    group_by(CHR) %>%
    summarise(x_min=min(plot_BP), 
              x_max=max(plot_BP), 
              center=(min(plot_BP) + max(plot_BP))/2,
              y_min="A", y_max="B")
  g <- ggplot(df) 
  if (back.panels){
    polygon.df <- axis.info.df %>%
      select(-center) %>%
      pivot_longer(cols=c(x_min,x_max), names_to="x_name",values_to="x") %>%
      pivot_longer(cols=c(y_min,y_max), names_to="y_name",values_to="y") %>%
      group_by(CHR) %>%
      mutate(y=c(0,Inf,Inf,0))
    g <- g + 
      geom_polygon(data = polygon.df,aes(x=x,y=y), fill = back.panels.col)
  }
  g <- g +
    geom_point(aes(x=plot_BP, y=-log10(P), color = as.factor(CHR)), size = point.size)
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
        geom_point(aes(BP, -log10(fdr), colour = as.factor(CHR)), size = point.size)
      if (!is.null(highlight)) {
        if (is.numeric(highlight)){
          highlight <- as.character(df$SNP[df$P < highlight])
        }
        if (any(!(highlight %in% df$SNP))){
          warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
        } else {
          g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                              aes(BP, -log10(fdr), group=SNP, colour=SNP),
                              color = highlight.col, size = point.size)
          highlight <- NULL
          y.max <- floor(max(-log10(df$fdr))) + 1
          for (j in (1:y.step)){
            if (y.max %% y.step != 0){
              y.max <- y.max + 1
            } else {
              break
            }
          }
        }
      }
    }
  }
  g <- g + scale_colour_manual(values = chrom.col) +
    scale_y_continuous(expand = c(0, 0), limit = c(0, y.max),
                       breaks = seq(from = 0, to = y.max, by = y.step)) +
    theme(strip.background = element_blank(), legend.position = "none",
          panel.background = element_rect(fill = "white", 
                                          colour = NA),
          axis.line.y = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title = graph.title, x = "Chromosome", y = y.title)
  if (!is.null(chrom.lab)){
    if (length(unique(df$CHR)) != length(chrom.lab)){
      warning("Number of chrom.lab different of number of chromosomes in dataset, argument ignored.")
      g <- g +
        scale_x_continuous(expand = c(0.01, 0),
                           breaks = axis.info.df$center, labels = axis.info.df$CHR)
    } else {
      g <- g +
        scale_x_continuous(expand = c(0.01, 0),
                           breaks = axis.info.df$center, labels = axis.info.df$CHR)
    }
  } else {
    g <- g +
      scale_x_continuous(expand = c(0.01, 0),
                         breaks = axis.info.df$center, labels = axis.info.df$CHR)
  }
  if (suggestiveline){
    g <- g + 
      annotate("segment", 
               x=min(df$plot_BP), y=-log10(suggestiveline), 
               xend = max(df$plot_BP), yend=-log10(suggestiveline), 
               color = suggestivecolor)
  }
  if (genomewideline){
    g <- g + annotate("segment", 
                      x=min(df$plot_BP), y=-log10(genomewideline), 
                      xend = max(df$plot_BP), yend=-log10(genomewideline), 
                      color = genomewidecolor)
  }
  if (!is.null(highlight)) {
    if (is.numeric(highlight)){
      highlight <- as.character(df$SNP[df$P < highlight])
    }
    if (any(!(highlight %in% df$SNP))){
      warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
    } else {
      g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                          aes(plot_BP, -log10(P), group=SNP, colour=SNP),
                          color = highlight.col, size = point.size)
      if (highlight.names){
        if (highlight.col != name.col){
          highlight.col <- name.col
        }
        g <- g + geom_label_repel(aes(plot_BP, -log10(P), label = label),
                                  color = highlight.col)
      }
    }
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

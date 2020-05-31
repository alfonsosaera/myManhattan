# https://www.r-graph-gallery.com/101_Manhattan_plot.html

setwd("/media/alfonso/data/GIT/myManhattan")

# install.packages("qqman")
# library(qqman)
# manhattan(ex, chr="CHR", bp="BP", snp="SNP", p="P" )

library(dplyr)
library(tidylog)
library(ggplot2)

ex <- read.table("example.txt")

ex$CHR <- factor(ex$CHR, levels = unique(ex$CHR), labels = c(1:22,"X","Y","MT")) # labels cannot be repetated

hola.df <- ex %>%
  group_by(CHR) %>%
  summarise(first=min(BP))

separator <- as.integer(ex %>%
                          group_by(CHR) %>%
                          summarise(CHR_size = max(BP)) %>% 
                          summarise(floor(max(CHR_size) / 20)))

tmp.df <- ex %>%
  group_by(CHR) %>%
  mutate(new_BP = BP - min(BP) + 1) %>%
  summarise(CHR_size = max(new_BP) + separator) %>%
  mutate(cum_size = cumsum(as.numeric(CHR_size )) - CHR_size)

ex <- ex %>%
  filter(CHR != 25) %>%
  group_by(CHR) %>%
  mutate(new_BP = BP - min(BP) + 1) %>%
  ungroup() %>%
  inner_join(tmp.df, by = "CHR") %>%
  mutate(plot_BP = new_BP + cum_size)

decoration.df <- ex %>%
  group_by(CHR) %>%
  summarise(x_min=min(plot_BP), 
            x_max=max(plot_BP), 
            center=(min(plot_BP) + max(plot_BP))/2,
            y_min="A", y_max="B")

polygon.df <- decoration.df %>%
  select(-center) %>%
  pivot_longer(cols=c(x_min,x_max), names_to="x_name",values_to="x") %>%
  pivot_longer(cols=c(y_min,y_max), names_to="y_name",values_to="y") %>%
  group_by(CHR) %>%
  mutate(y=c(0,Inf,Inf,0))

max_y <- max(-log10(ex$P)) * 1.01

y.max <- floor(max(-log10(ex$P))) + 1
if (y.max %% 2 != 0){
  y.max <- y.max + 1
}

ggplot(ex) +
  geom_polygon(data = polygon.df,aes(x=x,y=y), fill = "#ebebeb") + # panel.background = element_rect(fill = "grey92", colour = NA)
  geom_point(aes(x=plot_BP, y=-log10(P), color = as.factor(CHR))) +
  scale_y_continuous(limits = c(0, y.max), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0.01, 0),
                     breaks = decoration.df$center, labels = decoration.df$CHR) +
  labs(x="", y= expression(-log[10](italic(p)))) +
  scale_color_manual(values = rep(c("lightblue", "navy"), 13)) +
  theme(strip.background = element_blank(), legend.position = "none",
        panel.background = element_rect(fill = "white", 
                                        colour = NA),
        axis.line.y = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())



suggestiveline = 1e-05
suggestivecolor = "blue"
genomewideline = 5e-08
genomewidecolor = "red"

p <- ggplot(ex) +
  geom_polygon(data = polygon.df,aes(x=x,y=y), fill = "#ebebeb") + # panel.background = element_rect(fill = "grey92", colour = NA)
  geom_point(aes(x=plot_BP, y=-log10(P), color = as.factor(CHR))) +
  scale_y_continuous(limits = c(0, y.max), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0.01, 0),
                     breaks = decoration.df$center, labels = c(1:15,"",17,"",19,"",21,"","X","Y","")) +
  labs(x="", y= expression(-log[10](italic(p)))) +
  scale_color_manual(values = rep(c("lightblue", "navy"), 13)) +
  theme(strip.background = element_blank(), legend.position = "none",
        panel.background = element_rect(fill = "white", 
                                        colour = NA),
        axis.line.y = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())  

print(p)

p + guides(x = guide_axis(n.dodge = 2)) + 
  annotate("segment", x=min(ex$plot_BP), y=-log10(suggestiveline), xend = max(ex$plot_BP), yend=-log10(suggestiveline), color = suggestivecolor) + 
  annotate("segment", x=min(ex$plot_BP), y=-log10(genomewideline), xend = max(ex$plot_BP), yend=-log10(genomewideline), color = genomewidecolor)

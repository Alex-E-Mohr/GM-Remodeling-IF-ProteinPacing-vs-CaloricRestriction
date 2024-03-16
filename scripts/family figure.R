# load packages
library(patchwork); packageVersion("patchwork")     
library(ggplot2); packageVersion("ggplot2")


# Import data
DA <- read.csv(file = './combined_family_group_results.csv')

level_order1 <- factor(DA$feature, level = c('Butyricicoccaceae',
                                             'Streptococcaceae',
                                             'Ruminococcaceae',
                                             'Rikenellaceae',
                                             'Tannerellaceae',
                                             'Christensenellaceae'))

level_order2 <- factor(DA$value, level=c('IF_WK4', 'IF_WK8', 'HH_WK4', 'HH_WK8'),
                       label=c('IF-P (WK0-4)', 'IF-P (WK0-8)', 'CR (WK0-4)', 'CR (WK0-8)'))

p1 <- ggplot(DA, aes(x= forcats::fct_rev(level_order2), y = forcats::fct_rev(level_order1), fill = coef)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red") +
  geom_text(aes(label = label, color = coef), col = "black", nudge_y = 0, 
            size = 2,
            fontface = "bold",
            show.legend = TRUE) +
  scale_y_discrete(position = "right", guide = guide_axis(angle = -45)) +
  #ylab("Family") +
  theme_bw()  +
  theme(axis.title = element_blank()) +
  theme(text = element_text(size = 6, color = "black")) +
  theme(axis.text=element_text(colour="black"),
        panel.border = element_rect(colour = NA, fill=NA)) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  guides(fill = guide_colourbar(title = "Beta\nCoefficient", title.position = "top",
                                title.hjust = 0.5)) +
  coord_flip() +
  theme(legend.key.size = unit(0.2, 'cm'))

p1

# import interactions
DA1 <- read.csv(file = './combined_family_interaction_results.csv')

level_order3 <- factor(DA1$feature, level = c('Butyricicoccaceae',
                                              'Streptococcaceae',
                                              'Ruminococcaceae',
                                              'Rikenellaceae',
                                              'Tannerellaceae',
                                              'Christensenellaceae'))

level_order4 <- factor(DA1$value, level=c('WK4', 'WK8'),
                       label=c('IF-P vs. CR\n(WK0-4)', 'IF-P vs. CR\n(WK0-8)'))

colors <- c("#252525", "#636363", "#969696", "#cccccc", "#f7f7f7")

p2<- ggplot(DA1, aes(x= forcats::fct_rev(level_order4), y = forcats::fct_rev(level_order3), fill = label)) + 
  geom_tile(color = "black") +
  scale_fill_manual(values=colors, name="p.adj") +
  scale_y_discrete(guide = guide_axis(angle = -45)) +
  theme_bw()  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(text = element_text(size = 6, color = "black")) +
  theme(axis.text=element_text(colour="black"),
        panel.border = element_rect(colour = NA, fill=NA)) +
  theme(legend.title.align=0.5) +
  #theme(legend.position="bottom", legend.box = "horizontal") +
  coord_flip() +
  theme(legend.key.size = unit(0.2, 'cm'))

p2


# combine plots

p1 + plot_spacer() + p2 +
  plot_layout(heights = c(3.5, -1.2, 2), guides = "collect")

#png
ggsave("./Fig1_h.png", height = 5.5, width = 6.5, units = "cm")

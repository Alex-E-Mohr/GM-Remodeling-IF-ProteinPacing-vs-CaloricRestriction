# load packages
library(ggplot2); packageVersion("ggplot2")

# Import IFP data
IFP_pathway <- read.csv(file = './IFP_pathway.csv')

level_orderIFP <- factor(IFP_pathway$Metabolite, level = c('Citrate cycle (TCA cycle)',
                                                           'Alanine, aspartate and glutamate metabolism',
                                                           'Fatty acid degradation',
                                                           'Fatty acid elongation',
                                                           'Tyrosine metabolism',
                                                           'Butanoate metabolism',
                                                           'Glycerophospholipid metabolism',
                                                           'Inositol phosphate metabolism',
                                                           'Ascorbate and aldarate metabolism',
                                                           'Amino sugar and nucleotide sugar metabolism',
                                                           'Nicotinate and nicotinamide metabolism',
                                                           'Glycine, serine and threonine metabolism',
                                                           'Arginine and proline metabolism',
                                                           'Pentose and glucuronate interconversions'))

IFP_plot <- ggplot(IFP_pathway, aes(x = level_orderIFP, y = p, fill = Impact)) +
  geom_segment( aes(x=level_orderIFP, xend=level_orderIFP, y=0, yend=p), size = 0.5, color="grey") +
  geom_point(color = "black", shape=21, stroke=0.5) +
  coord_flip() +
  labs(title="IF-P Pathway Enrichment", y="-Log10(P-value)") +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 6)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_blank()) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.key.size = unit(0.2, 'cm')) +
  theme(plot.title = element_text(size=6))

IFP_plot

# Import CR data
CR_pathway <- read.csv(file = './CR_pathway.csv')

level_orderCR <- factor(CR_pathway$Metabolite, level = c('Butanoate metabolism',
                                                       'Glutathione metabolism',
                                                       'Tyrosine metabolism',
                                                       'Phenylalanine metabolism',
                                                       'Tryptophan metabolism',
                                                       'Purine metabolism',
                                                       'Arginine biosynthesis',
                                                       'Glycine, serine and threonine metabolism',
                                                       'Pyrimidine metabolism',
                                                       'Glycerolipid metabolism',
                                                       'Taurine and hypotaurine metabolism',
                                                       'Amino sugar and nucleotide sugar metabolism',
                                                       'D-Glutamine and D-glutamate metabolism',
                                                       'Citrate cycle (TCA cycle)',
                                                       'Primary bile acid biosynthesis',
                                                       'Pantothenate and CoA biosynthesis',
                                                       'Pyruvate metabolism',
                                                       'Glyoxylate and dicarboxylate metabolism',
                                                       'Glycolysis / Gluconeogenesis',
                                                       'Valine, leucine and isoleucine degradation',
                                                       'Phenylalanine, tyrosine and tryptophan biosynthesis',
                                                       'Alanine, aspartate and glutamate metabolism',
                                                       'Cysteine and methionine metabolism',
                                                       'Nicotinate and nicotinamide metabolism'))

CR_plot <- ggplot(CR_pathway, aes(x = level_orderCR, y = p, fill = Impact)) +
  geom_segment( aes(x=level_orderCR, xend=level_orderCR, y=0, yend=p), size = 0.5, color="grey") +
  geom_point(color = "black", shape=21, stroke=0.5) +
  coord_flip() +
  labs(title="CR Pathway Enrichment", y="-Log10(P-value)") +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 6)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_blank()) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.key.size = unit(0.2, 'cm')) +
  theme(plot.title = element_text(size=6))

CR_plot

# combine plots

IFP_plot / CR_plot +
  plot_layout(heights = c(1.1, 2))

#png
ggsave("./figures/Fig2c.png", height = 11, width = 9, units = "cm")

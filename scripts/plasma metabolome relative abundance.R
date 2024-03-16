# load packages
library(patchwork); packageVersion("patchwork")     
library(ggplot2); packageVersion("ggplot2")


# Import data
DA <- read.csv(file = './plasma_metabolome_DA.csv')

level_order1 <- factor(DA$Metabolite, level = c('2,3-Dihydroxybenzoic acid',
                                                'Protocatechuic acid',
                                                'Myoinositol',
                                                'Dulcitol',
                                                'D-Mannitol',
                                                'Agmatine',
                                                'Sorbitol',
                                                'Choline',
                                                'N-Acetylglutamine',
                                                'Oxaloacetic acid',
                                                'Cytidine',
                                                'Decanoylcarnitine',
                                                'L-(+)-Arabinose',
                                                'Xylitol',
                                                'Kynurenine',
                                                'Asparagine',
                                                'D-Galacturonic acid',
                                                'Pyroglutamic acid',
                                                '5-Hydroxyindoleacetic acid',
                                                '2-Hydroxyglutarate',
                                                'Dimethylglycine',
                                                'Tryptophan',
                                                'Urate',
                                                'Indole-3-lactic acid',
                                                'Asymmetric dimethylarginine',
                                                'Alanine',
                                                'Sarcosine',
                                                'Palmitic acid',
                                                'Valeric acid',
                                                '9-Octadecynoic acid',
                                                'Acetylcarnitine',
                                                'Malonic acid'))

# New facet label names for supp variable
DA$Group[DA$Group == 'IF_P'] <- 'IF-P'

p1 <- ggplot(DA, aes(x= Time, y = level_order1, fill = Abundance)) + 
  geom_tile() +
  facet_grid(~factor(Group, levels=c('IF-P', 'CR'))) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red") +
  scale_x_discrete(labels=c("WK0" = "0", "WK4" = "4", "WK8" = "8")) +
  theme_classic()  +
  labs(x="Time (weeks)") +
  theme(axis.title.y = element_blank()) +
  theme(text = element_text(size = 6, color = "black")) +
  theme(axis.text=element_text(colour="black"),
        panel.border = element_rect(colour = NA, fill=NA)) +
  theme(strip.background = element_rect(colour = "black", fill = "black", size=1)) +
  theme(strip.text = element_text(colour = 'white', size=5, face="bold")) +
  theme(legend.position="top", legend.box = "horizontal") +
  guides(fill = guide_colourbar(title = "Log10 Abundance", title.position = "top",
                                title.hjust = 0.5, barheight = 0.5)) +
  theme(strip.text.x = element_text(size = 6))

p1
  

# import interactions
GLM <- read.csv(file = './plasma_metabolome_GLM.csv')

level_order <- factor(GLM$Metabolite, level = c('2,3-Dihydroxybenzoic acid',
                                                'Protocatechuic acid',
                                                'Myoinositol',
                                                'Dulcitol',
                                                'D-Mannitol',
                                                'Agmatine',
                                                'Sorbitol',
                                                'Choline',
                                                'N-Acetylglutamine',
                                                'Oxaloacetic acid',
                                                'Cytidine',
                                                'Decanoylcarnitine',
                                                'L-(+)-Arabinose',
                                                'Xylitol',
                                                'Kynurenine',
                                                'Asparagine',
                                                'D-Galacturonic acid',
                                                'Pyroglutamic acid',
                                                '5-Hydroxyindoleacetic acid',
                                                '2-Hydroxyglutarate',
                                                'Dimethylglycine',
                                                'Tryptophan',
                                                'Urate',
                                                'Indole-3-lactic acid',
                                                'Asymmetric dimethylarginine',
                                                'Alanine',
                                                'Sarcosine',
                                                'Palmitic acid',
                                                'Valeric acid',
                                                '9-Octadecynoic acid',
                                                'Acetylcarnitine',
                                                'Malonic acid'))


p2 <- ggplot(GLM, aes(x = level_order, y = logFC, fill = sig)) +
  geom_point(color = "black", size = 1, shape=21, stroke=0.5) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype="dotted", color="red", linewidth = .25) +
  labs(y="LogFC (IF-P/CR)") +
  theme_classic() +
  theme(text = element_text(size = 6, color = "black")) +
  theme(axis.text.x = element_text(size = 5, color = "black")) +
  theme(axis.text.y=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_blank()) +
  scale_fill_manual(values=c("white", "Black")) +
  guides(fill = guide_legend(title = "FDR < 0.10", title.position = "top",
                                title.hjust = 0.5)) +
  theme(legend.position="top")

p2


# combine plots
p1 + plot_spacer() + p2 +
  plot_layout(widths = c(1.5, -.15, 1))


#png
ggsave("./Fig2b.png", height = 11, width = 9, units = "cm")

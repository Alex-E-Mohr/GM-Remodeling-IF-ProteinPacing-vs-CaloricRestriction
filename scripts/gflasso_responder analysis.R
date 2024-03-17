# Load packages
library(tidyverse)
library(corrplot)
library(pheatmap)
library(janitor)
library(viridis)
library(gflasso)
library(phyloseq)

### Prepare species dataframe

# Import species data
species <- read.csv(file = './data/species_GFLASSO.csv')

# Create new dataframe
species_ld2 <- data.frame(
  "id" = species$Sample,
  "OTU" = species$OTU,
  "Abundance" = species$Abundance)
head(species_ld2)

# convert to wide dataframe and prepare for GFLASSO
species_wide <- reshape(species_ld2, idvar = "id", timevar = "OTU", direction = "wide")
species_wide <- species_wide[order(species_wide$id),]
species_wide <- species_wide %>% `row.names<-`(., NULL) %>% column_to_rownames(var = "id") %>% clean_names()

# remove "log10_" prior to the taxa name
for ( col in 1:ncol(species_wide)){
  colnames(species_wide)[col] <-  sub("abundance_", "", colnames(species_wide)[col])}

# correlation plot
DS <- cor(species_wide, method = "pearson")
corrplot(DS)

################################################################################

### prepare fecal metabolomic dataframe
fecal_wide <- read.csv(file = './fecal_GFLASSO.csv')

fecal_wide <- fecal_wide[order(fecal_wide$Sample.ID),]
fecal_wide <- fecal_wide %>% `row.names<-`(., NULL) %>% column_to_rownames(var = "Sample.ID") %>% clean_names()

################################################################################

# check to see if any NAs
scale(species_wide)
scale(fecal_wide)

### Perform gflasso ###
set.seed(999)

#CV <- readRDS(file = "CV_gflasso.rds")
system.time(CV <- cv_gflasso(X = scale(fecal_wide), Y = scale(species_wide), R = DS, nCores = 8, 
                             additionalOpts = list(delta_conv = 1e-5, iter_max = 1e5)))

cv_plot_gflasso(CV)

gfMod <- gflasso(X = scale(fecal_wide), Y = scale(species_wide), R = DS, opts = list(lambda = CV$optimal$lambda,
                                                                                   gamma = CV$optimal$gamma, 
                                                                                   delta_conv = 1e-5,
                                                                                   iter_max = 1e5))

gfMod
#saveRDS(gfMod, "ps_objects/GFLASSO.rds")


colnames(gfMod$B) <- colnames(species_wide)
Lasso_data <- gfMod$B[abs(rowSums(gfMod$B)) > 0.2, ]
#rownames(Lasso_data) <- sapply(strsplit(rownames(Lasso_data), '_g_'), getElement, 2)
pheatmap(Lasso_data, show_rownames = T, fontsize = 6, col=viridis(40))
#pheatmap(Lasso_data, scale = "row", col=viridis(40), margins = c(4,15),
#         legend = TRUE, fontsize = 8)


### make better heatmap ###
corr_plot <- as.data.frame(Lasso_data)
corr_plot <- cbind(Metabolite = rownames(corr_plot), corr_plot)
rownames(corr_plot) <- 1:nrow(corr_plot)
corr_plot <- gather(corr_plot, Species, coefficient, collinsella_sgb14861:lachnospiraceae_bacterium_nsj_29)

level_order1 <- factor(corr_plot$Species, level = c('blautia_hydrogenotrophica',
                                                    'clostridium_leptum',
                                                    'collinsella_sgb14861',
                                                    'ggb3511_sgb4688',
                                                    'ggb38744_sgb14842',
                                                    'lachnospiraceae_bacterium_nsj_29',
                                                    'massiliimalia_timonensis',
                                                    'phascolarctobacterium_sgb4573',
                                                    'anaerostipes_hadrus',
                                                    'blautia_massiliensis',
                                                    'eubacterium_rectale',
                                                    'eubacterium_ventriosum',
                                                    'faecalicatena_contorta',
                                                    'mediterraneibacter_glycyrrhizinilyticus',
                                                    'roseburia_inulinivorans',
                                                    'streptococcus_salivarius'),
                      label = c('Blautia hydrogenotrophica',
                                'Clostridium leptum',
                                'Collinsella SGB14861',
                                'GGB3511 SGB4688',
                                'GGB38744 SGB14842',
                                'Lachnospiraceae bacterium NSJ 29',
                                'Massiliimalia timonensis',
                                'Phascolarctobacterium SGB4573',
                                'Anaerostipes hadrus',
                                'Blautia massiliensis',
                                'Eubacterium rectale',
                                'Eubacterium ventriosum',
                                'Faecalicatena contorta',
                                'Mediterraneibacter glycyrrhizinilyticus',
                                'Roseburia inulinivorans',
                                'Streptococcus salivarius'))

level_order2 <- factor(corr_plot$Metabolite, level = c('x_2e_4e_2_7_dimethyl_2_4_octadienedioic_acid',
                                                       'x3_7_dimethyluric_acid',
                                                       'methyl_2_3_dihydro_3_hydroxy_2_oxo_1h_indole_3_acetate',
                                                       'x4_guanidinobutanal',
                                                       'x7_ketocholesterol',
                                                       'acesulfame',
                                                       'l_aspartic_acid',
                                                       'azelaic_acid',
                                                       'cadaverine',
                                                       'famotidine',
                                                       'histamine',
                                                       'homocarnosine',
                                                       'itaconic_acid',
                                                       'leukotriene_b4_dimethylamide',
                                                       'linoleic_acid',
                                                       'malonic_acid',
                                                       'methyl_sorbate',
                                                       'n_nonanoylglycine',
                                                       'octadecanamide',
                                                       'palmitic_acid',
                                                       'phytosphingosine',
                                                       'ser_val',
                                                       'sorbitan_laurate',
                                                       'succinic_anhydride'),
                       label = c('(2E,4E)-2,7-Dimethyl-2,4-octadienedioic acid',
                                 '3,7-Dimethyluric acid',
                                 '3-Hydroxy-2-oxo-1H-indole-3-acetic acid',
                                 '4-Guanidinobutanal',
                                 '7-Ketocholesterol',
                                 'Acesulfame',
                                 'Aspartic acid',
                                 'Azelaic',
                                 'Cadaverine',
                                 'Famotidine',
                                 'Histamine',
                                 'Homocarnosine',
                                 'Itaconic',
                                 'Leukotriene B4 dimethylamide',
                                 'Linoleic acid',
                                 'Malonic acid',
                                 'Methyl sorbate',
                                 'N-Nonanoylglycine',
                                 'Octadecanamide',
                                 'Palmitic acid',
                                 'Phytosphingosine',
                                 'Ser-Val',
                                 'Sorbitan laurate',
                                 'Succinic anhydride'
                       ))

GFLASSO <- ggplot(corr_plot, aes(level_order2, forcats::fct_rev(level_order1), fill = coefficient)) +
  geom_tile(aes(fill=coefficient), color = "grey") +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  theme(axis.title = element_blank()) +
  theme(text = element_text(size = 5, color = "black")) +
  theme(axis.text = element_text(size = 5, color = "black")) +
  guides(fill = guide_colourbar(title = "GFLASSO\nCoefficients", title.position = "top",
                                title.hjust = 0.5, barwidth = 0.75, barheight = 4.5)) +
  theme(legend.title.align=0.5) +
  scale_fill_viridis() +
  theme(axis.text.y = element_text(face = "italic")) 

GFLASSO

#png
ggsave("./Fig4_j.png", height = 6, width = 16, units = "cm")

# Spearman correlations
cyto_asv <- merge(inflam_ld2, genus_wide, by.x = "row.names", by.y = "row.names")

# Corr with IL-4
library(corrr)
IL4_plot <-cyto_asv %>% select(., il_4, 16:85) %>%
  corrr::correlate(method = "spearman") %>%
  corrr::focus(il_4) %>%
  filter(abs(il_4) > 0.1) %>%
  #mutate(rowname = factor(rowname, levels = rowname[order(Fiber)])) %>%
  mutate(name = fct_reorder(term, desc(il_4))) %>%
  ggplot(aes(x = name, y = il_4)) +
  geom_bar(stat = "identity") +
  ylab("Correlation (rho) with IL-4\n (Spearman)") +
  xlab("Genus") + 
  theme_bw() +
  theme(text = element_text(size = 6, color = "black")) +
  theme(axis.text = element_text(size = 5, color = "black")) +
  coord_flip()

IL4_plot

library(rstatix)
il4_cor <- cyto_asv %>% select(., il_4, 16:85) %>%
  cor_test(method = "spearman", vars = "il_4")
il4_cor <- il4_cor %>% filter(., var1 == "il_4") 
il4_cor1 <- as.data.frame(p.adjust(il4_cor$p, method = "fdr"))
il4_cor <- cbind(il4_cor, il4_cor1)

# Corr with IL-6
IL6_plot <- cyto_asv %>% select(., il_6, 16:85) %>%
  corrr::correlate(method = "spearman") %>%
  corrr::focus(il_6) %>%
  filter(abs(il_6) > 0.1) %>%
  #mutate(rowname = factor(rowname, levels = rowname[order(Fiber)])) %>%
  mutate(name = fct_reorder(term, desc(il_6))) %>%
  ggplot(aes(x = name, y = il_6)) +
  geom_bar(stat = "identity") +
  ylab("Correlation (rho) with IL-6\n (Spearman)") +
  xlab("Genus") + 
  theme_bw() + 
  theme(text = element_text(size = 6, color = "black")) +
  theme(axis.text = element_text(size = 6, color = "black")) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

IL6_plot

il6_cor <- cyto_asv %>% select(., il_6, 16:85) %>%
  cor_test(method = "spearman", vars = "il_6")
il6_cor <- il6_cor %>% filter(., var1 == "il_6") 
il6_cor1 <- as.data.frame(p.adjust(il6_cor$p, method = "fdr"))
il6_cor <- cbind(il6_cor, il6_cor1)


# Corr with IL-8
IL8_plot <- cyto_asv %>% select(., il_8, 16:85) %>%
  corrr::correlate(method = "spearman") %>%
  corrr::focus(il_8) %>%
  filter(abs(il_8) > 0.1) %>%
  #mutate(rowname = factor(rowname, levels = rowname[order(Fiber)])) %>%
  mutate(name = fct_reorder(term, desc(il_8))) %>%
  ggplot(aes(x = name, y = il_8)) +
  geom_bar(stat = "identity") +
  ylab("Correlation (rho) with IL-8\n (Spearman)") +
  xlab("Genus") + 
  theme_bw() + 
  theme(text = element_text(size = 6, color = "black")) +
  theme(axis.text = element_text(size = 6, color = "black")) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

IL8_plot

il8_cor <- cyto_asv %>% select(., il_8, 16:85) %>%
  cor_test(method = "spearman", vars = "il_8")
il8_cor <- il8_cor %>% filter(., var1 == "il_8") 
il8_cor1 <- as.data.frame(p.adjust(il8_cor$p, method = "fdr"))
il8_cor <- cbind(il8_cor, il8_cor1)


# Corr with IL-13
IL13_plot <- cyto_asv %>% select(., il_13, 16:85) %>%
  corrr::correlate(method = "spearman") %>%
  corrr::focus(il_13) %>%
  filter(abs(il_13) > 0.1) %>%
  #mutate(rowname = factor(rowname, levels = rowname[order(Fiber)])) %>%
  mutate(name = fct_reorder(term, desc(il_13))) %>%
  ggplot(aes(x = name, y = il_13)) +
  geom_bar(stat = "identity") +
  ylab("Correlation (rho) with IL-13\n (Spearman)") +
  xlab("Genus") + 
  theme_bw() + 
  theme(text = element_text(size = 6, color = "black")) +
  theme(axis.text = element_text(size = 5, color = "black")) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

IL13_plot

il13_cor <- cyto_asv %>% select(., il_13, 16:85) %>%
  cor_test(method = "spearman", vars = "il_13")
il13_cor <- il13_cor %>% filter(., var1 == "il_13") 
il13_cor1 <- as.data.frame(p.adjust(il13_cor$p, method = "fdr"))
il13_cor <- cbind(il13_cor, il13_cor1)

# plot together
library(patchwork)

S2b <- IL4_plot | IL13_plot
S2b

#png
ggsave("./FigS2_b.png", height = 8, width = 18, units = "cm")

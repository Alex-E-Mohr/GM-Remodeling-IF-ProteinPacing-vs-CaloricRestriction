# Load packages
library(tidyverse)
library(corrplot)
library(pheatmap)
library(janitor)
library(viridis)
library(gflasso)
library(phyloseq)

### Prepare cytokine dataframe

# Import data
inflammation <- read.csv(file = './inflammation.csv')

# Import metadata
metadata <- read.table(file = "./metadata.txt", sep = "\t", header = TRUE)

# Join dataframes 
inflam <- merge(metadata, inflammation)

# Log2 transformation
inflam <- inflam %>%
  mutate(across(starts_with("GM_"), ~log2(.)))

# Remove missing samples
inflam <- inflam %>%
  filter(Sample.ID != "SM24_01" & Sample.ID != "SM27_01" & Sample.ID != "SM37_02" & PP_HH != "HH" & Time != "WK0" & subject != "SM36")

# Create new dataframe
inflam_ld2 <- inflam %>%
  select(Sample.ID, starts_with("log2")) %>%
  column_to_rownames(var = "Sample.ID") %>%
  clean_names()

# Calculate correlations
DS <- cor(inflam_ld2, method = "spearman")
corrplot(DS)

### Prepare genus dataframe

# Load phyloseq object
ps <- readRDS("./phyloseq_object.rds")

# Clean up data
ps <- ps %>%
  filter_taxa(function(x) sum(x > 0) > 2, prune = TRUE) %>%
  filter_taxa(function(x) sum(x > 3) > (0.30 * length(x)), TRUE) %>%
  tax_glom("Genus") %>%
  transform("compositional")

# Create a dataframe of agglomerated data at Genus level
genus_df <- psmelt(ps) %>%
  filter(Time != "WK0" & PP_HH != "HH" & subject != "SM36") %>%
  select(Sample, Abundance, Genus)

# Convert to wide dataframe
genus_wide <- genus_df %>%
  pivot_wider(names_from = Genus, values_from = Abundance) %>%
  column_to_rownames(var = "Sample")

# Remove "abundance_" prior to the taxa name
colnames(genus_wide) <- sub("abundance_", "", colnames(genus_wide))

# Perform gflasso
set.seed(999)
CV <- cv_gflasso(X = scale(genus_wide), Y = scale(inflam_ld2), R = DS, nCores = 2, 
                 additionalOpts = list(delta_conv = 1e-5, iter_max = 1e5))

cv_plot_gflasso(CV)

gfMod <- gflasso(X = scale(genus_wide), Y = scale(inflam_ld2), R = DS, 
                 opts = list(lambda = CV$optimal$lambda, gamma = CV$optimal$gamma, 
                             delta_conv = 1e-5, iter_max = 1e5))

colnames(gfMod$B) <- colnames(inflam_ld2)
Lasso_data <- gfMod$B[abs(rowSums(gfMod$B)) > 0.3, ]
rownames(Lasso_data) <- sapply(strsplit(rownames(Lasso_data), '_g_'), getElement, 2)

# heatmap
corr_plot <- as.data.frame(Lasso_data) %>%
  rownames_to_column(var = "genus") %>%
  pivot_longer(-genus, names_to = "cytokine", values_to = "coefficient")

level_order <- c('GM_CSF', 'IFNy', 'IL_1B', 'IL_2', 'IL_4', 'IL_5', 'IL_6', 
                 'IL_8', 'IL_10', 'IL_12p70', 'IL_13', 'IL_17A', 'IL_23', 'TNFa')

GFLASSO <- ggplot(corr_plot, aes(genus, fct_rev(fct_relevel(cytokine, level_order)), fill = coefficient)) +
  geom_tile() +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  theme(axis.title = element_blank(),
        text = element_text(size = 6, color = "black"),
        axis.text = element_text(size = 6, color = "black")) +
  guides(fill = guide_colourbar(title = "GFLASSO\nCoefficients", title.position = "top",
                                title.hjust = 0.5, barwidth = 1, barheight = 5)) +
  theme(legend.title.align = 0.5) +
  scale_fill_viridis()

# Spearman correlations
cyto_asv <- merge(inflam_ld2, genus_wide, by.x = "row.names", by.y = "row.names")

IL4_plot <- cyto_asv %>%
  select(il_4, everything()) %>%
  corrr::correlate(method = "spearman") %>%
  corrr::focus(il_4) %>%
  filter(abs(il_4) > 0.1) %>%
  mutate(name = fct_reorder(term, desc(il_4))) %>%
  ggplot(aes(x = name, y = il_4)) +
  geom_bar(stat = "identity") +
  ylab("Correlation (rho) with IL-4\n (Spearman)") +
  xlab("Genus") + 
  theme_bw() +
  theme(text = element_text(size = 6, color = "black"),
        axis.text = element_text(size = 5, color = "black")) +
  coord_flip()

il4_cor <- cyto_asv %>%
  select(il_4, everything()) %>%
  cor_test(method = "spearman", vars = "il_4") %>%
  filter(var1 == "il_4") %>%
  mutate(il4_cor1 = p.adjust(p, method = "fdr"))

# Combine plots
S2b <- IL4_plot | IL13_plot

# Save plots as PNG
ggsave("./figures/FigS2_a.png", GFLASSO, height = 6, width = 18, units = "cm")
ggsave("./figures/FigS2_b.png", S2b, height = 8, width = 18, units = "cm")

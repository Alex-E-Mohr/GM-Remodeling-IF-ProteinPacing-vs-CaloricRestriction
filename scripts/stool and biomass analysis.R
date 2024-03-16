# Load necessary packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(emmeans)
library(lmtest)
library(sandwich)
library(rcompanion)
library(nlme)
library(ggplot2)
library(plyr)

# Load data
metadata <- read.table(file = "./metadata.txt", sep = "\t", header = TRUE, row.names = 1)

# Normality Test
# Histograms
plotNormalHistogram(metadata$biomass, main = "biomass")
plotNormalHistogram(metadata$BSS, main = "BSS")
plotNormalHistogram(metadata$pH, main = "pH")
plotNormalHistogram(metadata$pH, main = "fecal_WT_g")

# Test for normality via Shapiro-Wilk normality test
normality_tests <- list(
  shapiro_biomass = shapiro.test(metadata$biomass),
  shapiro_biomass_log = shapiro.test(log10(metadata$biomass)),
  shapiro_BSS = shapiro.test(metadata$BSS),
  shapiro_BSS_log = shapiro.test(log(metadata$BSS)),
  shapiro_pH = shapiro.test(metadata$pH),
  shapiro_pH_log = shapiro.test(log(metadata$pH)),
  shapiro_fecal_WT = shapiro.test(metadata$fecal_WT_g),
  shapiro_fecal_WT_log = shapiro.test(log(metadata$fecal_WT_g))
)

normality_tests

# Transform pH and fecal weight.
metadata$log10biomass <- log10(metadata$biomass)
metadata$logPH <- log(metadata$pH)
metadata$logFecalWT <- log(metadata$fecal_WT_g)

### LME - Biomass
LME_biomass <- lme(log10biomass ~ Time_num * PP_HH + Age + Sex, data = metadata,
                   random = ~ 1 | as.factor(subject), correlation = corAR1())
LME_biomass
anova(LME_biomass)

### LME - BSS
LME_BSS <- lme(BSS ~ Time_num * PP_HH + Age + Sex, data = metadata,
               random = ~ 1 | as.factor(subject), correlation = corAR1())
LME_BSS
anova(LME_BSS)

### LME - pH
LME_pH <- lme(logPH ~ Time_num * PP_HH + Age + Sex, data = metadata,
              random = ~ 1 | as.factor(subject), correlation = corAR1())
LME_pH
anova(LME_pH)

### LME - FecalWT
LME_FecalWT <- lme(fecal_WT_g ~ Time_num * PP_HH + Age + Sex, data = metadata,
                   random = ~ 1 | as.factor(subject), correlation = corAR1())
LME_FecalWT
anova(LME_FecalWT)

# Plot Biomass
metadata$PP_HH <- revalue(metadata$PP_HH, c("HH" = "CR", "PP" = "IF-P"))

level_order <- factor(metadata$PP_HH, level = c("CR", "IF-P"))

biomass_plot <- ggplot(data = metadata, aes(x = Time, y = biomass, fill = level_order, color = level_order)) +
  geom_boxplot(lwd = 0.5, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 1, stroke = 0.5,
              position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.2)) +
  geom_line(aes(group = subject), size = 0.5, alpha = 0.2) +
  facet_grid(~factor(PP_HH, levels = c('IF-P', 'CR'))) +
  scale_fill_manual(name = NULL, values = c("magenta3", "darkturquoise")) +
  scale_color_manual(name = NULL, values = c("black", "black")) +
  labs(y = "16S rRNA (gene copies/g)", x = "Time (weeks)") +
  scale_x_discrete(labels = c("WK0" = "0", "WK4" = "4", "WK8" = "8")) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "black", size = 1),
        strip.text = element_text(colour = 'white', size = 6, face = "bold"),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 6),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

biomass_plot

saveRDS(biomass_plot, "./figures/biomass.rds")
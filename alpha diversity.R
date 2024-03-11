# Load necessary packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(emmeans)
library(lmtest)
library(sandwich)
library(picante)
library(rcompanion)
library(multcomp)
library(plyr)

# Load phyloseq object
phylo <- readRDS("./phyloseq_object.rds")

# Remove singletons
phylo <- filter_taxa(phylo, function(x) {sum(x > 0) > 1}, prune = TRUE)

# Print phyloseq object
print(phylo)

# Look at sequencing depth
summary(sample_sums(phylo))

# Plot rarefaction curve
asv_tab <- t(abundances(phylo))
plot_rare <- vegan::rarecurve(asv_tab, step = 50, label = FALSE, sample = min(rowSums(asv_tab), col = "blue", cex = 0.6))

# Check association of ASVs with number of reads
ggplot(data = data.frame("total_reads" = sample_sums(phylo),
                         "observed" = estimate_richness(phylo, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")

# Normalize samples to equal sampling depth
set.seed(111) # Keep result reproducible
(ps_rare <- rarefy_even_depth(phylo, rngseed = 243, replace = FALSE))

head(sample_sums(ps_rare))

### Phylogenetic diversity
asvtab <- as.data.frame(ps_rare@otu_table)
tree <- phylo@phy_tree

# Check if the tree is rooted
phylo@phy_tree

# Create data frame with Observed ASVs and PD
df.pd <- pd(t(asvtab), tree, include.root = TRUE) 
head(df.pd)

adiv <- data.frame(
  "Observed" = estimate_richness(ps_rare, measures = "Observed"),
  "PD" = sample_data(df.pd)$PD,
  "id" = sample_data(ps_rare)$subject,
  "time" = sample_data(ps_rare)$Time,
  "group" = sample_data(ps_rare)$PP_HH,
  "sex" = sample_data(ps_rare)$Sex,
  "bmi" = sample_data(ps_rare)$BMI,
  "BSS" = sample_data(ps_rare)$BSS,
  "pH" = sample_data(ps_rare)$pH,
  "age" = sample_data(ps_rare)$Age)

head(adiv)

# Normality Test
# Histograms
plotNormalHistogram(adiv$Observed, main = "Observed")
plotNormalHistogram(adiv$PD, main = "Faith's PD")

# Test for normality via Shapiro-Wilks normality test
shapiro.test(adiv$Observed)
shapiro.test(log(adiv$Observed))
shapiro.test(adiv$PD)
shapiro.test(log(adiv$PD)) 

# Transformations needed. Use logarithmic transformation
adiv$logObs <- log(adiv$Observed)
adiv$logPD <- log(adiv$PD)
head(adiv)

# Summary statistics Observed
adiv %>% 
  group_by(group, time) %>% 
  mutate(
    avg = mean(Observed), 
    min = min(Observed), 
    max = max(Observed), 
    SD = sd(Observed))

# Summary statistics PD
adiv %>% 
  group_by(group, time) %>% 
  mutate(
    avg = mean(PD), 
    min = min(PD), 
    max = max(PD), 
    SD = sd(PD))

### LME - Observed ASVs
LME.Observed <- lme(logObs ~ time * group + age + sex, data = adiv, 
                    random = ~ 1 | as.factor(id), correlation = corAR1())
LME.Observed
anova(LME.Observed) 

# Test pairwise comparison using the 'emmeans' package
emmeans(LME.Observed, pairwise ~ time)

### LME - Faith's PD
LME.PD <- lme(logPD ~ time * group + age + bmi + sex, data = adiv, 
              random = ~ 1 | as.factor(id), correlation = corAR1())
LME.PD
anova(LME.PD)

# Test pairwise comparison using the 'emmeans' package
emmeans(LME.PD, pairwise ~ time)

# Plots
adiv$group <- revalue(adiv$group, c("HH" = "CR", "PP" = "IF-P"))

# Observed ASVs
ObsASV.plot <- ggplot(data = adiv, aes(x = time, y = Observed, fill = group, color = group)) +
  geom_boxplot(lwd = 0.5, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 1, stroke = 0.5,
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
  geom_line(aes(group = id), size = 0.5, alpha = 0.2) +
  facet_grid(~ factor(group, levels = c('IF-P', 'CR'))) +
  scale_fill_manual(name = NULL, values = c("magenta3", "darkturquoise")) +
  scale_color_manual(name = NULL, values = c("black", "black")) +
  labs(y = "Observed ASVs", x = "Time (weeks)") +
  scale_x_discrete(labels = c("WK0" = "0", "WK4" = "4", "WK8" = "8")) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "black", size = 1),
        strip.text = element_text(colour = 'white', size = 6, face = "bold"),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 6),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

ObsASV.plot

# Plot Faith's PD
PD.plot <- ggplot(data = adiv, aes(x = time, y = PD, fill = group, color = group)) +
  geom_boxplot(lwd = 0.5, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 1, stroke = 0.5,
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
  geom_line(aes(group = id), size = 0.5, alpha = 0.2) +
  facet_grid(~ factor(group, levels = c('IF-P', 'CR'))) +
  scale_fill_manual(name = NULL, values = c("magenta3", "darkturquoise")) +
  scale_color_manual(name = NULL, values = c("black", "black")) +
  labs(y = "Phylogenetic Diversity", x = "Time (weeks)") +
  scale_x_discrete(labels = c("WK0" = "0", "WK4" = "4", "WK8" = "8")) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "black", size = 1),
        strip.text = element_text(colour = 'white', size = 6, face = "bold"),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 6),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

PD.plot

### Correlation and multiple regression analyses
# Load data
clinical <- read.table(file = "SPSS data/clinical data.txt", sep = "\t", header = TRUE)

# Prepare dataframes for joining alpha diversity data with clinical data
adiv_corr <- subset(adiv, select = -c(id, time, group, sex, bmi, age))
adiv_corr <- adiv_corr %>% rownames_to_column(var = "Sample.ID")

# Merge dataframes
merged_adiv <- merge(x = clinical, y = adiv_corr, by = "Sample.ID")

# Remove WK4 and WK8 data only looking at baseline
merged_adiv <- subset(merged_adiv, Time != "WK4")
merged_adiv <- subset(merged_adiv, Time != "WK8")

# Perform correlation analyses
cor.test(merged_adiv$logObs, merged_adiv$PostWtkg_PerChange, method = "pearson")
cor.test(merged_adiv$logObs, merged_adiv$LBM_PerChange, method = "pearson")
cor.test(merged_adiv$logObs, merged_adiv$BF_PerChange, method = "pearson")
cor.test(merged_adiv$logObs, merged_adiv$VAF_PerChange, method = "pearson")
cor.test(merged_adiv$logPD, merged_adiv$PostWtkg_PerChange, method = "pearson")
cor.test(merged_adiv$logPD, merged_adiv$LBM_PerChange, method = "pearson")
cor.test(merged_adiv$logPD, merged_adiv$BF_PerChange, method = "pearson")
cor.test(merged_adiv$logPD, merged_adiv$VAF_PerChange, method = "pearson")



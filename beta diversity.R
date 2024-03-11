# Load libraries
library(vegan)
library(labdsv)
library(ggplot2)
library(grid)
library(dplyr)
library(reshape2)
library(microbiome) # data analysis and visualization
library(tidyverse)
library(emmeans)
library(lmtest)
library(sandwich)
library(phyloseq) # also the basis of data object. Data analysis and visualization
library(microbiomeutilities) # some utility tools 
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(plyr)
library(patchwork)
library(ggsignif)
library(ggtext)
library(ggplot2)
library(nlme)
library(car)

# Load phyloseq object
phylo <- readRDS("./phyloseq_object.rds")

# Normalize the samples to equal sampling depth (7300) 
set.seed(111) # keep result reproductive
(ps_rare <- rarefy_even_depth(phylo, rngseed = 243, replace = FALSE))

head(sample_sums(ps_rare))
ps_rare

# Calculate Bray Curtis distance matrix
bray.dm <- phyloseq::distance(ps_rare, method = "bray")
bray.dm

# Make a data frame from the sample_data
bray.df <- data.frame(sample_data(ps_rare))
bray.df <- rownames_to_column(bray.df, var = "sample")
bray.df

### PERMANOVAs
set.seed(247)

# Nested PERMANOVA
nested_interaction <- adonis2(bray.dm ~ PP_HH:Time_num*as.character(subject),
                              data = bray.df, parallel = 4)
nested_interaction
str(nested_interaction)


### Assess intra-participant distances

# Start with the same dataframe
bray.df
bray.dm

# Create separate data frames
meta.samples <- bray.df[c("sample", "subject","Time_num", "PP_HH", "BMI", "Sex", "Age")] 
head(meta.samples)

meta.name <- bray.df[c("sample", "subject","Time_num", "PP_HH", "BMI", "Sex", "Age")] 
names(meta.name)[names(meta.name) == "sample"] <- "name"
head(meta.name)

# Convert bray curtis distances into a matrix
matrix_bray <- as.matrix(bray.dm) %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name)

matrix_bray
matrix_bray <- left_join(meta.samples, matrix_bray, by = "sample")
matrix_bray
matrix_bray <- left_join(meta.name, matrix_bray, by = "name") 
matrix_bray

# Check if number of time counts are accurate, again since we joined distances from all time points
table(matrix_bray$Time_num.y) # Yes, there are only prepost values

# Structure of data? Make sure 'Time_num' is numeric!
str(matrix_bray)

matrix_bray$Time_num.x <- as.numeric(as.character(matrix_bray$Time_num.x))
matrix_bray$Time_num.y <- as.numeric(as.character(matrix_bray$Time_num.y))

matrix_bray

# Next we want to identify 3 different categories:
# 1. Different participants (by group) at the same time point (i.e., pre)
# 2. Different participants (by group) at the same time point (i.e., post)
# 3. Same participant (by group) at pre and post times (i.e., same)
# GET RID OF ANYTHING ELSE!

dist.intra <- matrix_bray %>%
  mutate(comparison = case_when(
    subject.x == subject.y & Time_num.x == 2 & Time_num.y == 1 ~ "1",         # same week 4
    subject.x == subject.y & Time_num.x == 3 & Time_num.y == 1 ~ "2",         # same week 8
    TRUE ~ NA_character_)) 

dist.intra
dist.intra %>%
  count("comparison")

# Get rid of the NAs
dist.intra <- dist.intra %>%
  drop_na()

dist.intra %>%
  count("comparison")

dist.intra

# plot

order <- factor(dist.intra$PP_HH.y, level=c('IF-P', 'CR'))

beta_plot <- dist.intra %>%
  ggplot(aes(x=comparison, y=value, fill=order)) +
  geom_boxplot(lwd=0.5, outlier.shape = NA) +
  geom_jitter(shape=21, size=1, stroke = 0.5, alpha = 0.8,
              position=position_jitterdodge(dodge.width=0.8,
                                            jitter.width=0.2)) +
  scale_fill_manual(name=NULL,
                    values=c("darkturquoise", "magenta3")) +
  scale_color_manual(name=NULL,
                     values=c("black", "black")) +
  labs(y="Bray Curtis distance", x="Time (weeks)") +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "black", size=1)) +
  theme(strip.text = element_text(colour = 'white', size=12, face="bold")) +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 6)) +
  theme(panel.border = element_rect(color = "black", fill=NA, size=1)) +
  theme(legend.position = c(0.16, 0.885),
        legend.background = element_rect(fill = "white", color = "black")) +
  scale_x_discrete(breaks=c("1", "2"),
                   labels=c("4", "8")) +
  scale_y_continuous(limits=c(0.2,1.0)) +
  stat_compare_means(label = "p.format", method = "wilcox", size = 2,
                     label.y = c(0.78, 0.84)) +
  theme(legend.key.size = unit(0.2, 'cm'))

saveRDS(beta_plot, "./figures/Beta_diversity.rds")

# PNG
ggsave("./figures/Fig1_g.png", height = 6, width = 5, units = "cm")

### Linear mixed-effects models

# Assess group*time intra-individual distances

# Check distribution and normality
hist(dist.intra$value)
qqPlot(dist.intra$value)
shapiro.test(dist.intra$value) # Shapiro-Wilks normality test
shapiro.test(log(dist.intra$value)) # Best kept untransformed

intra.lme <- lme(value ~ comparison*PP_HH.y+Age.y+Sex.y, data = dist.intra, 
                 random = ~ 1|as.factor(subject.y), correlation = corAR1())
summary(intra.lme)
anova(intra.lme)

# Test pairwise comparison using the 'emmeans' package
emmeans(intra.lme, pairwise ~ comparison*PP_HH.y)

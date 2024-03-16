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
library(dplyr) # data handling 
library(patchwork)
library(ggsignif)
library(ggtext)
library(ggplot2)
library(nlme); packageVersion(nlme)
library(car)

################################################################################

# Load phyloseq object
phyloMETAB <- readRDS("./phyloseq_metabolome.rds")

# Calculate canberra distance matrix
can.dm <- phyloseq::distance(phyloMETAB, method = "canberra")
can.dm

# make a data frame from the sample_data
can.df <- data.frame(sample_data(phyloMETAB))
can.df <- rownames_to_column(phyloMETAB, var = "sample")
can.df


### Assess intra-participant distances

# start with the same dataframe
can.df
can.df <- rownames_to_column(can.df, var = "sample")
can.dm

# Create separate data frames
meta.samples <- can.df[c("sample", "subject","Time_num", "PP_HH", "BMI", "Sex", "Age")] 
head(meta.samples)

meta.name <- can.df[c("sample", "subject","Time_num", "PP_HH", "BMI", "Sex", "Age")] 
names(meta.name)[names(meta.name) == "sample"] <- "name"
head(meta.name)

# convert canberra distances into a matrix
matrix_can <- as.matrix(can.dm) %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name)

matrix_can
matrix_can <- left_join(meta.samples, matrix_can, by = "sample")
matrix_can
matrix_can <- left_join(meta.name, matrix_can, by = "name") 
matrix_can

# check if number of time counts are accurate, again since we joined distances from all time points
table(matrix_can$Time_num.y) # yes there are only prepost values

# structure of data? Make sure 'Time_num' is numeric!
str(matrix_can)

matrix_can$Time_num.x <- as.numeric(as.character(matrix_can$Time_num.x))
matrix_can$Time_num.y <- as.numeric(as.character(matrix_can$Time_num.y))

matrix_can

# Next we want to identify 3 different categories:
# 1. Different participants (by group) at the same time point (i.e., pre)
# 2. Different participants (by group) at the same time point (i.e., post)
# 3. Same participant (by group) at pre and post times (i.e., same)
# GET RID OF ANYTHING ELSE!

dist.intra <- matrix_can %>%
  mutate(comparison = case_when(
    subject.x == subject.y & Time_num.x == 2 & Time_num.y == 1 ~ "1",         # same week 4
    subject.x == subject.y & Time_num.x == 3 & Time_num.y == 1 ~ "2",         # same week 8
    TRUE ~ NA_character_)) 

dist.intra
dist.intra %>%
  count("comparison")

# get rid of the NAs
dist.intra <- dist.intra %>%
  drop_na()

dist.intra %>%
  count("comparison")

dist.intra


# create data frames with just comparison 1 and then just 2 and 3.

# New facet label names for supp variable
labs <- c("CR", "IF-P")
names(labs) <- c("HH", "PP")
dist.intra$PP_HH.y[dist.intra$PP_HH.y == 'PP'] <- 'IF-P'
dist.intra$PP_HH.y[dist.intra$PP_HH.y == 'HH'] <- 'CR'


order <- factor(dist.intra$PP_HH.y, level=c('IF-P', 'CR'))

dist.intra %>%
  ggplot(aes(x=comparison, y=value, fill=order)) +
  geom_boxplot(lwd=0.5, outlier.shape = NA) +
  geom_jitter(shape=21, size=1, stroke = 0.5, alpha = 0.8,
              position=position_jitterdodge(dodge.width=0.8,
                                            jitter.width=0.2)) +
  scale_fill_manual(name=NULL,
                    values=c("darkturquoise", "magenta3")) +
  scale_color_manual(name=NULL,
                     values=c("black", "black")) +
  labs(y="Canberra distance", x="Time (weeks)") +
  theme_bw() +
  #theme(legend.position="none") +
  theme(strip.background = element_rect(colour = "black", fill = "black", size=1)) +
  theme(strip.text = element_text(colour = 'white', size=12, face="bold")) +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 6)) +
  theme(panel.border = element_rect(color = "black", fill=NA, size=1)) +
  theme(legend.position = c(0.09, 0.905),
        legend.background = element_rect(fill = "white", color = "black")) +
  scale_x_discrete(breaks=c("1", "2"),
                   labels=c("4", "8")) +
  scale_y_continuous(limits=c(0.2,1.0)) +
  stat_compare_means(label = "p.format", method = "wilcox", size = 2,
                     label.y = c(0.84, 0.82)) +
  theme(legend.key.size = unit(0.2, 'cm'))

#png
ggsave("./Fig.png", height = 6, width = 7, units = "cm")

### Linear mixed-effects models

# Assess group*time intra-individual distances

# check distribution and normality
hist(dist.intra$value)
qqPlot(dist.intra$value)
shapiro.test(dist.intra$value) # Shapiro-Wilks normality test
shapiro.test(log(dist.intra$value)) # best kept untransformed

intra.lme <- lme(value ~ comparison*PP_HH.y+Age.y+Sex.y, data = dist.intra, 
                 random = ~ 1|as.factor(subject.y), correlation = corAR1())
summary(intra.lme)
anova(intra.lme)

# test pairwise comparison using the 'emmeans' package
emmeans(intra.lme, pairwise ~ comparison:PP_HH.y)



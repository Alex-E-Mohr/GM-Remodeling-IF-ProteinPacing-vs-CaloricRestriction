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
library(plyr)



################################################################################

# Load phyloseq object
phylo <- readRDS("./phyloseq_responder.rds")


# Calculate bray curtis distance matrix
bray.dm <- phyloseq::distance(phylo, method = "bray")
bray.dm

# make a data frame from the sample_data
bray.df <- data.frame(sample_data(phylo))
bray.df <- rownames_to_column(bray.df, var = "sample")
bray.df

### PERMANOVAs
set.seed(247)

# Nested PERMANOVA
nested_interaction <- adonis2(bray.dm ~ group:time+as.character(subject),
                              data = bray.df, parallel = 4)
nested_interaction
str(nested_interaction)


### Assess intra-participant distances

# start with the same dataframe
bray.df
bray.dm

# Create separate data frames
meta.samples <- bray.df[c("sample", "subject","time", "group", "preBMI", "sex", "age")] 
head(meta.samples)

meta.name <- bray.df[c("sample", "subject","time", "group", "preBMI", "sex", "age")] 
names(meta.name)[names(meta.name) == "sample"] <- "name"
head(meta.name)

# convert bray curtis distances into a matrix
matrix_bray <- as.matrix(bray.dm) %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name)

matrix_bray
matrix_bray <- left_join(meta.samples, matrix_bray, by = "sample")
matrix_bray
matrix_bray <- left_join(meta.name, matrix_bray, by = "name") 
matrix_bray

# check if number of time counts are accurate, again since we joined distances from all time points
table(matrix_bray$time.y) # yes there are only prepost values

# structure of data? Make sure 'Time_num' is numeric!
str(matrix_bray)

matrix_bray$time.x <- as.numeric(as.character(matrix_bray$time.x))
matrix_bray$time.y <- as.numeric(as.character(matrix_bray$time.y))

matrix_bray

# Next we want to identify 3 different categories:
# 1. Different participants (by group) at the same time point (i.e., pre)
# 2. Different participants (by group) at the same time point (i.e., post)
# 3. Same participant (by group) at pre and post times (i.e., same)
# GET RID OF ANYTHING ELSE!

dist.intra <- matrix_bray %>%
  mutate(comparison = case_when(
    subject.x == subject.y ~ "1",
    TRUE ~ NA_character_)) 

dist.intra
dist.intra %>%
  count("comparison")

# get rid of the NAs
dist.intra <- dist.intra %>%
  filter(comparison == 1)

dist.intra %>%
  count("comparison")

dist.intra


# create data frames with just comparison 1 and then just 2 and 3.

# plot
dist.intra %>%
  ggplot(aes(x=group.y, y=value, fill=group.y, color=group.y)) +
  geom_boxplot(lwd=1.5) +
  geom_jitter(shape=21, size=3, stroke = 1.5,
              position=position_jitterdodge(dodge.width=0.8,
                                            jitter.width=0.2)) +
  scale_fill_manual(name=NULL,
                    values=c("#E69F00", "#56B4E9")) +
  scale_color_manual(name=NULL,
                     values=c("black", "black")) +
  labs(y="Bray Curtis dissimilarity", x="Time (weeks)") +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 6)) +
  theme(panel.border = element_rect(color = "black", fill=NA, size=1)) +
  scale_y_continuous(limits=c(0.2,1.0)) +
  stat_compare_means(label = "p.format", method = "wilcox", size = 2,
                     label.y = c(0.78, 0.84))




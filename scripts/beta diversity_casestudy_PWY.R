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


# Load phyloseq object
phylo <- readRDS("./humann_PWY_casestudy.rds")

# Calculate bray Curtis distance matrix
bray.dm <- phyloseq::distance(phylo, method = "bray")
bray.dm

# make a data frame from the sample_data
bray.df <- data.frame(sample_data(phylo))
bray.df <- rownames_to_column(bray.df, var = "sample")
bray.df


### Assess intra-participant distances

# start with the same dataframe
bray.df
bray.dm

# Create separate data frames
meta.samples <- bray.df[c("sample", "subject","time")] 
head(meta.samples)

meta.name <- bray.df[c("sample", "subject","time")] 
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

dist.intra <- data.frame (Subject = c("SM8", "SM8", "SM8", "SM8", "SM8", "SM8", "SM8"),
                          Time  = c("WK0", "WK4", "WK8", "WK12", "WK16", "WK32", "WK52"),
                  Bray = c("0", "0.053098173", "0.023657260", "0.030001398", "0.076141247", "0.051620950", "0.007214996"))


dist.intra$Bray <- as.numeric(as.character(dist.intra$Bray))

str(dist.intra)

# Plot
# Ensure time is properly ordered
dist.intra$Time<-factor(dist.intra$Time, c("WK0", "WK4", "WK8", "WK12", "WK16", "WK32", "WK52"))

dist.intra <- dist.intra %>%
  ggplot(aes(x=Time, y=Bray)) +
  geom_line(linetype = "solid", aes(group = Subject), color="darkgray") +
  geom_point(shape=21, size=1, stroke = 1) +
  labs(y="Bray Curtis dissimilarity (Pathways)", x="Time (weeks)") +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.y = element_text(colour = "black", size = 6), 
        axis.text.x = element_text(colour = "black", size = 6), 
        axis.title = element_text(colour = "black", size = 6)) +
  theme(axis.text=element_text(colour="black")) +
  theme(panel.border = element_rect(color = "black", fill=NA, size=1)) +
  scale_y_continuous(limits=c(0,1.0)) +
  scale_x_discrete(labels=c("WK0" = "0", "WK12" = "12", "WK16" = "16", "WK32" = "32", 
                            "WK4" = "4", "WK52" = "52", "WK8" = "8"))

dist.intra

saveRDS(dist.intra, "./BC_PWY_casestudy.rds")


################################################################################

### Pathway PERMANOVA

# remove taxa not seen >5 in 30% of samples
(phylo <- filter_taxa(phylo, function(x) sum(x > 5) > (0.3*length(x)), TRUE)) 

# transform counts to percent relative abundance
ps.rel = transform_sample_counts(phylo, function(x) x/sum(x)*100)

# dataframes for permanova
otu <- abundances(ps.rel)
meta <- meta(ps.rel)

permanova <- adonis(t(otu) ~ time_num, data = meta, permutations=999)

# P-value
print(as.data.frame(permanova$aov.tab)["time", "Pr(>F)"])

# check which taxa contribute most to the community differences
coef <- coefficients(permanova)["time_num",]
sort(coef)
top.coef <- coef[rev(order(abs(coef)))[1:20]]

# quick visual check
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

# Create dataframe with top coefficents
top.coef.df <- as.data.frame(top.coef)
top.coef.df <- rownames_to_column(top.coef.df, var = "OTU")

# clean up dataframe for plotting
Species.coef = subset(top.coef.df, select = c(OTU, top.coef)) %>%
  mutate(OTU = str_replace_all(OTU, "_", " ")) %>% #remove _
  mutate(OTU = str_replace_all(OTU, "\\[|\\]", "")) %>% #remove brackets around taxa
  mutate(OTU = str_replace_all(OTU, " group", "")) %>% #remove group
  mutate(pos = top.coef >= 0) # create a new column called pos, which indicates whether the value is positive or negative

write.csv(Species.coef, "./PWY_COEF.csv", row.names=FALSE)

Species.coef <- read.csv("./PWY_COEF.csv", header = T)

# plot
casestudy_coef <- ggplot(Species.coef, aes(reorder(OTU, +top.coef),top.coef, fill=pos, alpha = 0.7)) +
  geom_bar(stat = "identity", color ="black") +
  scale_fill_manual(values=c("gray", "black")) +
  theme_bw() +
  theme(legend.position="none") +
  ylab("Model Coefficents") +
  theme(axis.title.y = element_blank()) +
  theme(axis.text=element_text(colour="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(text = element_text(size = 6)) +
  coord_flip() 

casestudy_coef

saveRDS(casestudy_coef, "./casestudy_PWY_coef.rds")

### combine plots
dist.intra + casestudy_coef 

#png
ggsave("./Fig5_fg.png", height = 6, width = 12, units = "cm")


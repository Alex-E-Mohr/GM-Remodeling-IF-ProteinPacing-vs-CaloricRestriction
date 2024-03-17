#load packages
library(microbiome)
library(phyloseq)
library(RColorBrewer)
library(ggpubr)
library(dplyr)  
library(phyloseq)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(patchwork)
library(metagMisc)
library(ggplot2)
library(plyr); packageVersion("plyr")
library(forcats); packageVersion("forcats")
library(ggalluvial)
library(reshape2)

# Load phyloseq object
phylo <- readRDS("./phyloseq_casestudy.rds")

# taxa at Phylum level
table(tax_table(phylo)[, 2])

# remove unassigned features
phylo <- subset_taxa(phylo, Phylum !="Bacteria_unclassified")

# filter very low abundance ASVs. Now have 622 ASVs
phylo<-filter_taxa(phylo, function(x) mean(x) >5, TRUE)
ntaxa(phylo)

# Number of taxa at each phylogenic level
(phylo_phylum <- tax_glom(phylo, "Phylum"))
(phylo_Class <- tax_glom(phylo, "Class"))
(phylo_Order <- tax_glom(phylo, "Order"))
(phylo_Family <- tax_glom(phylo, "Family"))
(phylo_Genus <- tax_glom(phylo, "Genus"))
(phylo_Genus <- tax_glom(phylo, "Species"))

# transform counts to percent relative abundance
ps.rel = transform_sample_counts(phylo, function(x) x/sum(x)*100)

# agglomerate taxa to Species level
glom <- tax_glom(ps.rel, taxrank = "Species", NArm = FALSE)
glom

# create a dataframe of agglomerated data at species level.
dat <- psmelt(glom) 
dat2 <- reshape(dat, idvar = "OTU", timevar = "time_num", direction = "wide")
write.csv(dat2, "./casestudy_abundances.csv", row.names=FALSE)
str(dat)
means <- ddply(dat, ~Species, function(x) c(mean=mean(x$Abundance))) #find the mean count per species
means # Bacteroides has the greatest coverage

# find genera whose mean relative abundance is less than 1%
other <- means[means$mean <= 1.34,]$Species
str(other)
other
# change their name to "Other" to make plot less complex
dat[dat$Species %in% other,]$Species <- 'Other'

# clean up dataframe for plotting
dat = subset(dat, select = c(Sample, time, Species, Abundance)) %>%
  mutate(Species = str_replace_all(Species, "_", " ")) %>% #remove _
  mutate(Species = str_replace_all(Species, " group", "")) %>% #remove group
  mutate(Species = str_replace_all(Species, "\\[|\\]", "")) #remove brackets around taxa
dat
  
unique(dat$Species) # check out the unique taxa names

# Remove "Other" from the plot
df <- dat[!dat$Species == "Other", ]

# Ensure time is properly ordered
df$time<-factor(df$time, c("WK0", "WK4", "WK8", "WK12", "WK16", "WK32", "WK52"))

# Generate color palette
colorCount <- length(unique(df$Species))
colorCount # need 20 colors
pal<- c("black", "cyan", "aquamarine1", "darkorange", "chartreuse","darkviolet", "azure3", "darkseagreen3","darkred",
        "grey", "plum", "gold", "pink", "lightslategray", "lightyellow","wheat", "navy", "steelblue",
        "coral", "darkgreen")
getPalette <- colorRampPalette(pal)


rel_abun = ggplot(df, aes( x = time, y = Abundance, alluvium = Species)) + 
    geom_alluvium(aes(fill = Species), colour = "black", alpha = 0.6, decreasing = FALSE) + 
    labs(x = "Time (weeks)", y = "Species Relative Abundance (%)", fill = "", colour = "") +
    scale_fill_manual(values=getPalette(colorCount)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,80)) + 
    scale_x_discrete(expand = c(0.02,0.02)) + 
    theme(axis.text.y = element_text(colour = "black", size = 5), 
          axis.text.x = element_text(colour = "black", size = 5), 
          axis.title = element_text(colour = "black", size = 5), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
          legend.text = element_text(size = 5, face = "italic")) +
  #guides(fill=guide_legend(ncol=1,byrow=TRUE)) +
  theme(legend.key.size = unit(.25, "cm")) +
  scale_x_discrete(labels=c("WK0" = "0", "WK12" = "12", "WK16" = "16", "WK32" = "32", 
                            "WK4" = "4", "WK52" = "52", "WK8" = "8"))
  
rel_abun


#png
ggsave("./Fig5_e.png", height = 6, width = 18, units = "cm")


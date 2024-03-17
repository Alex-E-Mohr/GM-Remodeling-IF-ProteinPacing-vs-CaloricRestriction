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
phylo <- readRDS("./humann_PWY_casestudy.rds")

# filter 
phylo<-filter_taxa(phylo, function(x) mean(x) >5, TRUE)
ntaxa(phylo)

# transform counts to percent relative abundance
ps.rel = transform_sample_counts(phylo, function(x) x/sum(x)*100)

# create a dataframe of agglomerated data
dat <- psmelt(ps.rel) 
dat2 <- reshape(dat, idvar = "OTU", timevar = "time_num", direction = "wide")
write.csv(dat2, "./casestudy_PWY_abundances.csv", row.names=FALSE)
str(dat)
means <- ddply(dat, ~OTU, function(x) c(mean=mean(x$Abundance))) 
means # Bacteroides has the greatest coverage

# find features whose mean relative abundance is less than 1%
other <- means[means$mean <= 1,]$OTU
str(other)
other
# change their name to "Other" to make plot less complex
dat[dat$OTU %in% other,]$OTU <- 'Other'

# clean up dataframe for plotting
dat = subset(dat, select = c(Sample, time, OTU, Abundance)) %>%
  mutate(OTU = str_replace_all(OTU, "_", " ")) %>% #remove _
  mutate(OTU = str_replace_all(OTU, " group", "")) %>% #remove group
  mutate(OTU = str_replace_all(OTU, "\\[|\\]", "")) #remove brackets around taxa
dat
  
unique(dat$OTU) # check out the unique taxa names

# Remove "Other" from the plot
df <- dat[!dat$OTU == "Other", ]

# Ensure time is properly ordered
df$time<-factor(df$time, c("WK0", "WK4", "WK8", "WK12", "WK16", "WK32", "WK52"))

# Generate color palette
colorCount <- length(unique(df$OTU))
colorCount # need 20 colors
pal<- c("black", "cyan", "aquamarine1", "darkorange", "chartreuse","darkviolet", "azure3", "darkseagreen3","darkred",
        "grey", "plum", "gold", "pink", "lightslategray", "lightyellow","wheat", "navy", "steelblue",
        "coral", "darkgreen")
getPalette <- colorRampPalette(pal)


rel_abun = ggplot(df, aes( x = time, y = Abundance, alluvium = OTU)) + 
    geom_alluvium(aes(fill = OTU), colour = "black", alpha = 0.6, decreasing = FALSE) + 
    labs(x = "Time (weeks)", y = "Pathway Relative Abundance (%)", fill = "", colour = "") +
    scale_fill_manual(values=getPalette(colorCount)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,80)) + 
    scale_x_discrete(expand = c(0.02,0.02)) + 
    theme(axis.text.y = element_text(colour = "black", size = 6), 
          axis.text.x = element_text(colour = "black", size = 6), 
          axis.title = element_text(colour = "black", size = 6), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
          legend.text = element_text(size = 6, face = "italic")) +
  theme(legend.key.size = unit(.3, "cm")) +
  scale_x_discrete(labels=c("WK0" = "0", "WK12" = "12", "WK16" = "16", "WK32" = "32", 
                            "WK4" = "4", "WK52" = "52", "WK8" = "8")) 
  #guides(fill=guide_legend(ncol=2,byrow=TRUE))
  
rel_abun


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
phylo <- readRDS("./phyloseq_subgroup_metab.rds")

# taxa at Main Class level
table(tax_table(phylo)[, 3])

# filter very low abundance ASVs. Now have 622 ASVs
phylo<-filter_taxa(phylo, function(x) mean(x) >5, TRUE)
ntaxa(phylo)

# transform counts to percent relative abundance
ps.rel = transform_sample_counts(phylo, function(x) x/sum(x)*100)

# agglomerate
glom <- tax_glom(ps.rel, taxrank = "Sub_Class", NArm = FALSE)
glom

# create a dataframe of agglomerated
dat <- psmelt(glom) 
#dat2 <- reshape(dat, idvar = "OTU", timevar = "time_num", direction = "wide")
#write.csv(dat, "./casestudy_abundances_metab.csv", row.names=FALSE)
#str(dat)
means <- ddply(dat, ~Sub_Class*group*time, function(x) c(mean=mean(x$Abundance))) #find the mean count per species
means # Bacteroides has the greatest coverage


# find genera whose mean relative abundance is less than 1%
other <- means[means$mean <= 1,]$Sub_Class
str(other)
other
# change their name to "Other" to make plot less complex
dat[dat$Sub_Class %in% other,]$Sub_Class <- 'Other'

# Remove "Other" from the plot
dat <- dat[!dat$Sub_Class == "Other", ]

# clean up dataframe for plotting
dat = subset(dat, select = c(group, time, Sub_Class, Abundance)) %>%
  mutate(Sub_Class = str_replace_all(Sub_Class, "_", " ")) %>% #remove _
  mutate(Sub_Class = str_replace_all(Sub_Class, " group", "")) %>% #remove group
  mutate(Sub_Class = str_replace_all(Sub_Class, "\\[|\\]", "")) #remove brackets around taxa
dat
  
unique(dat$Sub_Class) # check out the unique taxa names

# Remove "Other" from the plot
#df <- dat[!dat$Species == "Other", ]

# Ensure time is properly ordered
dat$time<-factor(dat$time, c("WK0", "WK8"))

# Generate color palette
colorCount <- length(unique(dat$Sub_Class))
colorCount # need 14 colors
pal<- c("lightslategray", "cyan", "aquamarine1", "darkorange", "chartreuse","red", "azure3", "darkseagreen3","darkred",
        "grey", "plum", "gold", "pink", "lightslategray")
getPalette <- colorRampPalette(pal)

str(dat)

rel_abun = ggplot(means, aes( x = time, y = mean, alluvium = Sub_Class)) +
    geom_alluvium(aes(fill = Sub_Class), colour = "black", alpha = 0.6, decreasing = FALSE) + 
    facet_wrap(~factor(group, levels=c('High', 'Low'))) +
    labs(x = "Time (weeks)", y = "Metabolite Relative Abundance (%)", fill = "", colour = "") +
    scale_fill_manual(values=getPalette(colorCount)) +
    #scale_y_continuous(expand = c(0,0), limits = c(0,80)) + 
    scale_x_discrete(expand = c(0.02,0.02)) + 
    theme(axis.text.y = element_text(colour = "black", size = 5), 
          axis.text.x = element_text(colour = "black", size = 5), 
          axis.title = element_text(colour = "black", size = 5), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
          legend.text = element_text(size = 5)) +
  theme(strip.background = element_rect(colour = "black", fill = "black", size=1)) +
  theme(strip.text = element_text(colour = 'white', size=6, face="bold")) +
  #guides(fill=guide_legend(ncol=1,byrow=TRUE)) +
  theme(legend.key.size = unit(.25, "cm")) #+
  #scale_x_discrete(labels=c("WK0" = "0", "WK8" = "8"))
  
rel_abun

#png
ggsave("./Fig4g.png", height = 6, width = 8, units = "cm")


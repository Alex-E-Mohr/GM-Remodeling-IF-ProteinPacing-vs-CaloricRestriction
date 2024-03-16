# Load packages
library(tidyverse); packageVersion("tidyverse")    
library(phyloseq); packageVersion("phyloseq")          
library(microbiome); packageVersion("microbiome")     
library(emmeans); packageVersion("emmeans")
library(lmtest); packageVersion("lmtest")    
library(sandwich); packageVersion("sandwich")
library(picante); packageVersion("picante")
library(rcompanion); packageVersion("rcompanion")
library(multcomp); packageVersion("multcomp")
library(plyr); packageVersion("plyr")
library(patchwork)

# Load phyloseq object
phylo <- readRDS("./phyloseq_responder.rds")

# use print option to see the data saved as phyloseq object.
print(phylo)

# First, look at sequencing depth (i.e., how much of the taxa have been sampled)
# for each sample.
summary(sample_sums(phylo))


# Create data frame with Observed ASVs and Shannon.
adiv <- data.frame(
  "Observed" = estimate_richness(phylo, measures = "Observed"),
  "Shannon" = estimate_richness(phylo, measures = "Shannon"),
  "id" = sample_data(phylo)$subject,
  "time" = sample_data(phylo)$time,
  "group" = sample_data(phylo)$group,
  "sex" = sample_data(phylo)$sex,
  "bmi" = sample_data(phylo)$preBMI,
  "age" = sample_data(phylo)$age)
head(adiv)
adiv

# Normality Test

# Histograms
plotNormalHistogram(adiv$Observed, main= "Observed")
plotNormalHistogram(adiv$Shannon, main= "Shannon")

# Test for normality via Shapiro-Wilks normality test
shapiro.test(adiv$Observed)
shapiro.test(log(adiv$Observed))

shapiro.test(adiv$Shannon)
shapiro.test(log(adiv$Shannon)) 

# Transformations needed. Use logarithmic transformation
adiv$logObs=log(adiv$Observed)
adiv$logShannon=log(adiv$Shannon)
head(adiv)


### LME - Observed ASVs
LME.Observed <- lme(logObs ~ time*group, data = adiv, 
                    random = ~ 1|as.factor(id), correlation = corAR1())
LME.Observed
anova(LME.Observed)


### LME - Shannon
LME.Shannon <- lme(logShannon ~ time*group, data = adiv, 
              random = ~ 1|as.factor(id), correlation = corAR1())
LME.Shannon
anova(LME.Shannon)


################################################################################

### Plots
head(adiv)
adiv$group <- revalue(x = adiv$group, 
                      c("HH" = "CR", "PP" = "IF-P"))

# Observed ASVs
ObsASV.plot <- ggplot(data=adiv, aes(x=time, y=Observed, fill=group, color=group)) +
  geom_boxplot(lwd=0.5, outlier.shape = NA) +
  geom_jitter(shape=21, size=0.5, stroke = 0.5,
              position=position_jitterdodge(dodge.width=0.8,
                                            jitter.width=0.2)) +
  scale_fill_manual(name=NULL,
                    values=c("#E69F00", "#56B4E9")) +
  scale_color_manual(name=NULL,
                     values=c("black", "black")) +
  labs(y="Observed Species", x="Time\n(weeks)") +
  scale_x_discrete(labels=c("WK0" = "0", "WK8" = "8")) +
  theme_bw() +
  theme(legend.position="none") +
  theme(strip.background = element_rect(colour = "black", fill = "black", size=1)) +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 5)) +
  theme(axis.text = element_text(size = 5)) +
  theme(panel.border = element_rect(color = "black", fill=NA, size=1))

ObsASV.plot

#saveRDS(ObsASV.plot, "./ObsASV.plot.rds")


# Plot Shannon
Shannon.plot <- ggplot(data=adiv, aes(x=time, y=Shannon, fill=group, color=group)) +
  geom_boxplot(lwd=0.5, outlier.shape = NA) +
  geom_jitter(shape=21, size=0.5, stroke = 0.5,
              position=position_jitterdodge(dodge.width=0.8,
                                            jitter.width=0.2)) +
  scale_fill_manual(name=NULL,
                    values=c("#E69F00", "#56B4E9")) +
  scale_color_manual(name=NULL,
                     values=c("black", "black")) +
  labs(y="Shannon Index", x="Time\n(weeks)") +
  scale_x_discrete(labels=c("WK0" = "0", "WK8" = "8")) +
  theme_bw() +
  theme(legend.position="none") +
  theme(strip.background = element_rect(colour = "black", fill = "black", size=1)) +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 5)) +
  theme(axis.text = element_text(size = 5)) +
  theme(panel.border = element_rect(color = "black", fill=NA, size=1))

Shannon.plot

#saveRDS(PD.plot, "./PD.plot.rds")

### combine plots
# Load biomass figure

ObsASV.plot + Shannon.plot

#png
ggsave("./Fig4_cd.png", height = 4, width = 4, units = "cm")

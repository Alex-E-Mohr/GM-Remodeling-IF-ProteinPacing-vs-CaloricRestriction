# Load packages
library(tidyverse); packageVersion("tidyverse")    
library(phyloseq); packageVersion("phyloseq")          
library(microbiome); packageVersion("microbiome")     
library(sandwich); packageVersion("sandwich")
library(rcompanion); packageVersion("rcompanion")
library(multcomp); packageVersion("multcomp")
library(plyr); packageVersion("plyr")

# Load phyloseq object
phylo <- readRDS("./phyloseq_casestudy.rds")

# use print option to see the data saved as phyloseq object.
print(phylo)

# Create data frame with Observed ASVs and Shannon.
adiv <- data.frame(
  "Observed" = estimate_richness(phylo, measures = "Observed"),
  "Shannon" = estimate_richness(phylo, measures = "Shannon"),
  "id" = sample_data(phylo)$subject,
  "time" = sample_data(phylo)$time,
  "time_num" = sample_data(phylo)$time_num,
  "group" = sample_data(phylo)$group,
  "sex" = sample_data(phylo)$sex,
  "Wt" = sample_data(phylo)$Wt,
  "Wt_ch" = sample_data(phylo)$Wt_ch)
head(adiv)
adiv

################################################################################

### Plots

library(hrbrthemes)

level_order <- c("WK0", "WK4", "WK8", "WK12", "WK16", "WK32", "WK52") 

# Observed ASVs
scaled_adiv <- adiv %>%
  mutate(Observed_tr = (Observed - min(Observed))/ (max(Observed) - min(Observed)),
         Observed_min = min(Observed),
         Observed_max = max(Observed),
         Wt_ch_tr = (Wt_ch - min(Wt_ch))/ (max(Wt_ch) - min(Wt_ch)),
         Wt_ch_min = min(Wt_ch),
         Wt_ch_max = max(Wt_ch))
         
plot <- scaled_adiv %>%
  ggplot(aes(x = factor(time, level = level_order), y = Observed_tr)) +
  geom_point(color = "red", size = 0.5) +
  geom_point(aes(y=Wt_ch_tr), color = "blue", size = 0.5) +
  scale_y_continuous(labels = seq(80, 260, 20),
                     breaks = (seq(80, 260, 20) - 185)/(240-185),
                     limits = (c(80, 260) - 185)/(240-185),
                     name = "Observed Species",
                     sec.axis = sec_axis(trans = ~.,
                                         labels = seq(0, 110, 10),
                                         breaks = (seq(0, 110, 10) - 66.56423)/(100-66.56423),
                                         name = "% Baseline BW")
                     ) +
  theme_bw() +
  theme(axis.title.y.left = element_text(color = "red"),
        axis.title.y.right = element_text(color = "blue")) +
  theme(text = element_text(size = 6, color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.text=element_text(colour="black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

ObsASV.plot <- plot +
  scale_x_discrete(labels=c("WK0" = "0", "WK12" = "12", "WK16" = "16", "WK32" = "32", 
                            "WK4" = "4", "WK52" = "52", "WK8" = "8")) 

ObsASV.plot

#saveRDS(ObsASV.plot, "./ObsASV_casestudy.rds")


# Plot Shannon
scaled_adiv <- adiv %>%
  mutate(Shannon_tr = (Shannon - min(Shannon))/ (max(Shannon) - min(Shannon)),
         Shannon_min = min(Shannon),
         Shannon_max = max(Shannon),
         Wt_ch_tr = (Wt_ch - min(Wt_ch))/ (max(Wt_ch) - min(Wt_ch)),
         Wt_ch_min = min(Wt_ch),
         Wt_ch_max = max(Wt_ch))

plot <- scaled_adiv %>%
  ggplot(aes(x = factor(time, level = level_order), y = Shannon_tr)) +
  geom_point(color = "red", size = 0.5) +
  geom_point(aes(y=Wt_ch_tr), color = "blue", size = 0.5) +
  scale_y_continuous(labels = seq(2.4, 4.2, 0.2),
                     breaks = (seq(2.4, 4.2, 0.2) - 3.523862)/(4.093061-3.523862),
                     limits = (c(2.4, 4.2) - 3.523862)/(4.093061-3.523862),
                     name = "Shannon Index",
                     sec.axis = sec_axis(trans = ~.,
                                         labels = seq(0, 110, 10),
                                         breaks = (seq(0, 110, 10) - 66.56423)/(100-66.56423),
                                         name = "% Baseline BW")
  ) +
  xlab("Time (weeks)") +
  theme_bw() +
  theme(axis.title.y.left = element_text(color = "red"),
        axis.title.y.right = element_text(color = "blue")) +
  theme(text = element_text(size = 6, color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.text=element_text(colour="black"))

Shannon.plot <- plot +
  scale_x_discrete(labels=c("WK0" = "0", "WK12" = "12", "WK16" = "16", "WK32" = "32", 
                            "WK4" = "4", "WK52" = "52", "WK8" = "8")) 

Shannon.plot

#saveRDS(Shannon.plot, "./Shannon_casestudy.rds")

### combine plots
ObsASV.plot / Shannon.plot 

#png
ggsave("./Fig5_ab.png", height = 6, width = 6, units = "cm")

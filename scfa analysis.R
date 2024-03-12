# Load packages
library(emmeans); packageVersion("emmeans")
library(nlme); packageVersion("nlme")  
library(tidyverse); packageVersion("tidyverse") 
library(rcompanion); packageVersion("rcompanion")
library(ggplot2); packageVersion("ggplot2")
library(IMIFA); packageVersion("IMIFA")

# Import scfa data
scfa <- read.csv(file = './scfa_data.csv')

# Import metadata
metadata <- read.table(file = "./metadata.txt", sep = "\t", header = T)

# Join dataframes 
scfa <- merge(metadata, scfa)

### Normality Tests ###

# Histograms
plotNormalHistogram(scfa$Acetic_acid, main= "Acetic_acid")
plotNormalHistogram(scfa$Propionic_acid, main= "Propionic_acid")
plotNormalHistogram(scfa$Butyric_acid, main= "Butyric_acid")
plotNormalHistogram(scfa$Valeric_acid, main= "Valeric_acid")

# Test for normality via Shapiro-Wilks normality test
shapiro.test(scfa$Acetic_acid)
shapiro.test(log2(scfa$Acetic_acid))

shapiro.test(scfa$Propionic_acid)
shapiro.test(log2(scfa$Propionic_acid))

shapiro.test(scfa$Butyric_acid)
shapiro.test(log2(scfa$Butyric_acid))

shapiro.test(scfa$Valeric_acid)
shapiro.test(log2(scfa$Valeric_acid))


# Transformations
scfa$log2_Acetic_acid=log2(scfa$Acetic_acid)
scfa$log2_Propionic_acid=log2(scfa$Propionic_acid)
scfa$log2_Butyric_acid=log2(scfa$Butyric_acid)
scfa$log2_Valeric_acid=log2(scfa$Valeric_acid)
head(scfa)


### LMEs ###

# Acetic_acid
LME.Acetic_acid <- lme(log2_Acetic_acid ~ Time*PP_HH+Age+Sex, data = scfa, 
                       random = ~ 1|as.factor(subject), correlation = corAR1())
LME.Acetic_acid
anova(LME.Acetic_acid)

# test pairwise comparison using the 'emmeans' package
emmeans(LME.Acetic_acid, pairwise ~ Time:PP_HH) 


# Propionic_acid
LME.Propionic_acid <- lme(log2_Propionic_acid ~ Time*PP_HH+Age+Sex, data = scfa, 
                          random = ~ 1|as.factor(subject), correlation = corAR1())
LME.Propionic_acid
anova(LME.Propionic_acid)

# test pairwise comparison using the 'emmeans' package
emmeans(LME.Propionic_acid, pairwise ~ Time:PP_HH) 


# Butyric_acid
LME.Butyric_acid <- lme(log2_Butyric_acid ~ Time*PP_HH+Age+Sex, data = scfa, 
                        random = ~ 1|as.factor(subject), correlation = corAR1())
LME.Butyric_acid
anova(LME.Butyric_acid)

# test pairwise comparison using the 'emmeans' package
emmeans(LME.Butyric_acid, pairwise ~ Time:PP_HH) 

# Valeric_acid
LME.Valeric_acid <- lme(log2_Valeric_acid ~ Time*PP_HH+Age+Sex, data = scfa, 
                        random = ~ 1|as.factor(subject), correlation = corAR1())
LME.Valeric_acid
anova(LME.Valeric_acid)

# test pairwise comparison using the 'emmeans' package
emmeans(LME.Valeric_acid, pairwise ~ Time:PP_HH)


# Create data frame for plotting
scfa_plot <- data.frame(
  "group" = scfa$PP_HH,
  "id" = scfa$subject,
  "time" = scfa$Time,
  "Acetic_acid" = scfa$log2_Acetic_acid,
  "Propionic_acid" = scfa$log2_Propionic_acid,
  "Butyric_acid" = scfa$log2_Butyric_acid,
  "Valeric_acid" = scfa$log2_Valeric_acid
  )
head(scfa_plot)


# transform to long format
data_long <- gather(scfa_plot, scfa, measurement, Acetic_acid:Valeric_acid, factor_key=TRUE)
data_long

# New facet label names
labs <- c("Acetate", "Propionate", "Butyrate", "Valerate")
names(labs) <- c("Acetic_acid", "Propionic_acid", "Butyric_acid", "Valeric_acid")

level_order <- factor(data_long$group, level = c("IF-P","CR"))

ggplot(data=data_long, aes(x=time, y=measurement, fill=level_order, color=level_order)) +
  geom_boxplot(lwd=.2, outlier.shape = NA) +
  geom_jitter(shape=21, size=0.8, stroke = .2,
              position=position_jitterdodge(dodge.width=0.8,
                                            jitter.width=0.2)) +
  facet_grid(~scfa, labeller = labeller(scfa = labs)) + 
  scale_fill_manual(name=NULL,
                    values=c("darkturquoise", "magenta3")) +
  scale_color_manual(name=NULL,
                     values=c("black", "black")) +
  labs(y="Log2(SCFA Concentration)", x="Time (weeks)") +
  scale_x_discrete(labels=c("WK0" = "0", "WK4" = "4",
                            "WK8" = "8")) +
  theme_bw() + 
  #theme(legend.position="none") +
  theme(strip.background = element_rect(colour = "black", fill = "black", size=1)) +
  theme(strip.text = element_text(colour = 'white', size=6, face="bold")) +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 6, color = "black")) +
  theme(panel.border = element_rect(color = "black", fill=NA, size=1)) +
  theme(legend.position = c(0.968, 0.92),
        legend.background = element_rect(fill = "white", color = "black"))  +
  theme(legend.key.size = unit(0.2, 'cm'))


#png
ggsave("./figures/FigS1_c.png", height = 7, width = 18, units = "cm")


### Correlations with cytokines ####

### Prepare cytokine dataframe

# Import data
inflammation <- read.csv(file = './inflammation.csv')

# Import metadata
metadata <- read.table(file = "./metadata.txt", sep = "\t", header = T)

# Join dataframes 
inflam <- merge(metadata, inflammation)

# Log2 transformation
inflam$log2_GM_CSF=log2(inflam$GM_CSF)
inflam$log2_IFNy=log2(inflam$IFNy)
inflam$log2_IL_1B=log2(inflam$IL_1B)
inflam$log2_IL_2=log2(inflam$IL_2)
inflam$log2_IL_4=log2(inflam$IL_4)
inflam$log2_IL_5=log2(inflam$IL_5)
inflam$log2_IL_6=log2(inflam$IL_6)
inflam$log2_IL_8=log2(inflam$IL_8)
inflam$log2_IL_10=log2(inflam$IL_10)
inflam$log2_IL_12p70=log2(inflam$IL_12p70)
inflam$log2_IL_13=log2(inflam$IL_13)
inflam$log2_IL_17A=log2(inflam$IL_17A)
inflam$log2_IL_23=log2(inflam$IL_23)
inflam$log2_TNFa=log2(inflam$TNFa)
head(inflam)

# Remove missing samples
inflam_1 <- inflam[!(inflam$Sample.ID=="SM24_01"),]
inflam_2 <- inflam_1[!(inflam_1$Sample.ID=="SM27_01"),]
inflam_ld <- inflam_2[!(inflam_2$Sample.ID=="SM37_02"),]


# Create new dataframe
inflam_ld2 <- data.frame(
  "Sample.ID" = inflam$Sample.ID,
  "GM_CSF" = inflam$log2_GM_CSF,
  "IFNy" = inflam$log2_IFNy,
  "IL_1B" = inflam$log2_IL_1B,
  "IL_2" = inflam$log2_IL_2,
  "IL_4" = inflam$log2_IL_4,
  "IL_5" = inflam$log2_IL_5,
  "IL_6" = inflam$log2_IL_6,
  "IL_8" = inflam$log2_IL_8,
  "IL_10" = inflam$log2_IL_10,
  "IL_12p70" = inflam$log2_IL_12p70,
  "IL_13" = inflam$log2_IL_13,
  "IL_17A" = inflam$log2_IL_17A,
  "IL_23" = inflam$log2_IL_23,
  "TNFa" = inflam$log2_TNFa)
head(inflam_ld2)


# merge dataframes
merged <- merge(inflam_ld2, scfa)

# Corr with Acetic_acid
library(corrr)
merged %>% select(., log2_Acetic_acid, 2:15) %>%
  corrr::correlate(method = "spearman") %>%
  corrr::focus(log2_Acetic_acid) %>%
  #filter(abs(log2_Acetic_acid) > 0.1) %>%
  #mutate(rowname = factor(rowname, levels = rowname[order(Fiber)])) %>%
  mutate(name = fct_reorder(term, desc(log2_Acetic_acid))) %>%
  ggplot(aes(x = name, y = log2_Acetic_acid)) +
  geom_bar(stat = "identity") +
  ylab("Correlation (rho) with Acetic acid\n (Spearman)") +
  xlab("Cytokines") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

library(rstatix)
Acetic_acid_cor <- merged %>% select(., log2_Acetic_acid, 2:15) %>%
  cor_test(method = "spearman", vars = "log2_Acetic_acid")
Acetic_acid_cor <- Acetic_acid_cor %>% filter(., var1 == "log2_Acetic_acid") 
Acetic_acid_cor1 <- as.data.frame(p.adjust(Acetic_acid_cor$p, method = "fdr"))
Acetic_acid_cor <- cbind(Acetic_acid_cor, Acetic_acid_cor1)


# Corr with Propionic_acid
merged %>% select(., log2_Propionic_acid, 2:15) %>%
  corrr::correlate(method = "spearman") %>%
  corrr::focus(log2_Propionic_acid) %>%
  #filter(abs(log2_Propionic_acid) > 0.1) %>%
  #mutate(rowname = factor(rowname, levels = rowname[order(Fiber)])) %>%
  mutate(name = fct_reorder(term, desc(log2_Propionic_acid))) %>%
  ggplot(aes(x = name, y = log2_Propionic_acid)) +
  geom_bar(stat = "identity") +
  ylab("Correlation (rho) with Propionic acid\n (Spearman)") +
  xlab("Cytokines") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

Propionic_acid_cor <- merged %>% select(., log2_Propionic_acid, 2:15) %>%
  cor_test(method = "spearman", vars = "log2_Propionic_acid")
Propionic_acid_cor <- Propionic_acid_cor %>% filter(., var1 == "log2_Propionic_acid") 
Propionic_acid_cor1 <- as.data.frame(p.adjust(Propionic_acid_cor$p, method = "fdr"))
Propionic_acid_cor <- cbind(Propionic_acid_cor, Propionic_acid_cor1)


# Corr with Butyric_acid
merged %>% select(., log2_Butyric_acid, 2:15) %>%
  corrr::correlate(method = "spearman") %>%
  corrr::focus(log2_Butyric_acid) %>%
  #filter(abs(log2_Butyric_acid) > 0.1) %>%
  #mutate(rowname = factor(rowname, levels = rowname[order(Fiber)])) %>%
  mutate(name = fct_reorder(term, desc(log2_Butyric_acid))) %>%
  ggplot(aes(x = name, y = log2_Butyric_acid)) +
  geom_bar(stat = "identity") +
  ylab("Correlation (rho) with Butyric acid\n (Spearman)") +
  xlab("Cytokines") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

Butyric_acid_cor <- merged %>% select(., log2_Butyric_acid, 2:15) %>%
  cor_test(method = "spearman", vars = "log2_Butyric_acid")
Butyric_acid_cor <- Butyric_acid_cor %>% filter(., var1 == "log2_Butyric_acid") 
Butyric_acid_cor1 <- as.data.frame(p.adjust(Butyric_acid_cor$p, method = "fdr"))
Butyric_acid_cor <- cbind(Butyric_acid_cor, Butyric_acid_cor1)


# Corr with Valeric_acid
merged %>% select(., log2_Valeric_acid, 2:15) %>%
  corrr::correlate(method = "spearman") %>%
  corrr::focus(log2_Valeric_acid) %>%
  #filter(abs(log2_Valeric_acid) > 0.1) %>%
  #mutate(rowname = factor(rowname, levels = rowname[order(Fiber)])) %>%
  mutate(name = fct_reorder(term, desc(log2_Valeric_acid))) %>%
  ggplot(aes(x = name, y = log2_Valeric_acid)) +
  geom_bar(stat = "identity") +
  ylab("Correlation (rho) with Valeric acid\n (Spearman)") +
  xlab("Cytokines") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

Valeric_acid_cor <- merged %>% select(., log2_Valeric_acid, 2:15) %>%
  cor_test(method = "spearman", vars = "log2_Valeric_acid")
Valeric_acid_cor <- Valeric_acid_cor %>% filter(., var1 == "log2_Valeric_acid") 
Valeric_acid_cor1 <- as.data.frame(p.adjust(Valeric_acid_cor$p, method = "fdr"))
Valeric_acid_cor <- cbind(Valeric_acid_cor, Valeric_acid_cor1)





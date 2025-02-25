# load packages
library(Maaslin2); packageVersion("Maaslin2")     
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(ggrepel); packageVersion("ggrepel")
library(dplyr)

# Load phyloseq object
phylo <- readRDS("./phyloseq_responder.rds")

# remove taxa not seen >5 in 30% of samples
(phylo <- filter_taxa(phylo, function(x) sum(x > 5) > (0.3*length(x)), TRUE))   

# Agglomerate to Family level
(ps_family <- tax_glom(phylo, "Family")) #detected 49 taxa

# Here we look at what we have at the Family level
table(tax_table(ps_family)[, 5])

# create OTU dataframe with annotations
OTU <- data.frame(otu_table(ps_family))
TAX <- data.frame(tax_table(ps_family))

library(tibble)
OTU <- tibble::rownames_to_column(OTU, "OTU")

TAX <- TAX %>%
  dplyr::select(Family)

TAX <- select(TAX, c('Family'))
TAX <- tibble::rownames_to_column(TAX, "OTU")

OTUTAX <- merge(TAX, OTU, by="OTU")
OTUTAX = subset(OTUTAX, select = -c(OTU))
OTUTAX <- data.frame(OTUTAX[,-1], row.names = OTUTAX[,1])

str(OTUTAX)

# check metadata
head(sample_data(ps_family))

# create dataframe
metadata <- data.frame(sample_data(ps_family))

# create interactions
metadata$High = (metadata$group == "High") *
  metadata$time_num
metadata$Low = (metadata$group == "Low") *
  metadata$time_num

str(metadata)

# interactions at week 4
res <- Maaslin2(
  input_data = OTUTAX,
  input_metadata = metadata,
  output = "./Maaslin2_default_output",
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  fixed_effects = c("High", "Low", "age", "sex"), 
  reference = c("High, 1", "Low, 1","Sex, Male"), 
  random_effects = ("subject"),
  min_prevalence = 0.0,
  min_abundance = 0.0,
  plot_heatmap = TRUE,
  plot_scatter = TRUE)

str(res)

res <- res$results %>%
  dplyr::filter(metadata == "High") %>%
  dplyr::select(-qval) %>%
  dplyr::arrange(pval) %>%
  dplyr::mutate(qval = p.adjust(pval, method = "fdr"))

head(res)


# add a column of NAs
res$diffexpressed <- "NO"
# if coef > 2.0 and q-value < 0.1, set as "UP" 
res$diffexpressed[res$coef > 2 & res$qval < 0.1] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res$diffexpressed[res$coef < -2 & res$qval < 0.1] <- "DOWN"
res$delabel <- NA
res$delabel[res$diffexpressed != "NO"] <- res$feature[res$diffexpressed != "NO"]


# Volcano plot
DA_WG <- ggplot(data=res, aes(x=coef, y=-log10(qval))) + 
  geom_point(color="darkgray", alpha = 0.5) +
  geom_vline(xintercept=c(-2, 2), col="black", linetype="dashed", size=.5) +
  geom_hline(yintercept=-log10(0.1), col="black", linetype="dashed", size=.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 6))

DA_WG


DA_WL <- ggplot(data=res, aes(x=coef, y=-log10(qval), 
                                       col=diffexpressed,
                                       label=delabel)) + 
  geom_point(alpha = 0.5) +
  geom_text_repel(size = 2, fontface = "italic") +
  geom_vline(xintercept=c(-2, 2), col="black", linetype="dashed", size=.5) +
  geom_hline(yintercept=-log10(0.1), col="black", linetype="dashed", size=.5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 6)) + 
  theme(legend.position="none") +
  xlim(-4, 4) +
  scale_color_manual(values=c("#999999"))

DA_WL



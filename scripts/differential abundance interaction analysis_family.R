# load packages
library(Maaslin2); packageVersion("Maaslin2")     
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse") 

# Load phyloseq object
phylo <- readRDS("./phyloseq_object.rds")

# remove taxa not seen > 30% of samples
(phylo <- filter_taxa(phylo, function(x) sum(x > 5) > (0.3*length(x)), TRUE))   

# Agglomerate to Family level
(ps_family <- tax_glom(phylo, "Family")) #detected 28 taxa

# Here we look at what we have at the family level
table(tax_table(ps_genus)[, 5])

# create OTU dataframe with annotations
OTU <- data.frame(otu_table(ps_family))
TAX <- data.frame(tax_table(ps_family))

library(tibble)
OTU <- tibble::rownames_to_column(OTU, "OTU")

TAX <- select(TAX, c('Family'))
TAX <- tibble::rownames_to_column(TAX, "OTU")

OTUTAX <- merge(TAX, OTU, by="OTU")
OTUTAX = subset(OTUTAX, select = -c(OTU))
OTUTAX <- data.frame(OTUTAX[,-1], row.names = OTUTAX[,1])

# check metadata
head(sample_data(ps_family))

# create dataframe
metadata <- data.frame(sample_data(ps_family))

### Week 4 ###

# remove week 8 samples
metadata_wk4 <-metadata[!(metadata$Time=="WK8"),]

# create interactions for week 4
metadata_wk4$PP = (metadata_wk4$PP_HH == "PP") *
  metadata_wk4$Time_num
metadata_wk4$HH = (metadata_wk4$PP_HH == "HH") *
  metadata_wk4$Time_num

# interactions at week 4
wk4_res <- Maaslin2(
  input_data = OTUTAX,
  input_metadata = metadata_wk4,
  output = "./Maaslin2_output",
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  fixed_effects = c("PP", "HH", "Age", "Sex"), 
  reference = c("PP, 0", "HH, 0","Sex, Male"), 
  random_effects = ("subject"),
  min_prevalence = 0.0,
  min_abundance = 0.0,
  plot_heatmap = TRUE,
  plot_scatter = TRUE)

str(wk4_res)

wk4_res <- wk4_res$results %>%
  dplyr::filter(metadata == "PP") %>%
  dplyr::select(-qval) %>%
  dplyr::arrange(pval) %>%
  dplyr::mutate(qval = p.adjust(pval, method = "fdr"))

head(wk4_res)

# save results
write.csv(wk4_res,"Maaslin2_output\\WK4_family_results.csv", row.names = FALSE)

### Week 8 ###

# remove week 4 samples
metadata_wk8 <-metadata[!(metadata$Time=="WK4"),]

# create interactions for week 8
metadata_wk8$PP = (metadata_wk8$PP_HH == "PP") *
  metadata_wk8$Time_num
metadata_wk8$HH = (metadata_wk8$PP_HH == "HH") *
  metadata_wk8$Time_num

# interactions at week 8
wk8_res <- Maaslin2(
  input_data = OTUTAX,
  input_metadata = metadata_wk8,
  output = "./Maaslin2_output",
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  fixed_effects = c("PP", "HH", "Age", "Sex"), 
  reference = c("PP, 0", "HH, 0","Sex, Male"), 
  random_effects = ("subject"),
  min_prevalence = 0.0,
  min_abundance = 0.0,
  plot_heatmap = TRUE,
  plot_scatter = TRUE)

str(wk8_res)

wk8_res <- wk8_res$results %>%
  dplyr::filter(metadata == "PP") %>%
  dplyr::select(-qval) %>%
  dplyr::arrange(pval) %>%
  dplyr::mutate(qval = p.adjust(pval, method = "fdr"))

head(wk8_res)

# save results
write.csv(wk8_res,"Maaslin2_output\\WK8_family_results.csv", row.names = FALSE)

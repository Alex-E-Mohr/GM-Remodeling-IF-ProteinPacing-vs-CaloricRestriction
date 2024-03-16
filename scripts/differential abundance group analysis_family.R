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
table(tax_table(ps_family)[, 5])


### IF-P (PP) Group ###

# subset for PP group
PP_family <- subset_samples(ps_family, PP_HH == "PP")

# create OTU dataframe with annotations
PP_OTU <- data.frame(otu_table(PP_family))
PP_tax <- data.frame(tax_table(PP_family))

library(tibble)
PP_OTU <- tibble::rownames_to_column(PP_OTU, "OTU")

PP_tax <- select(PP_tax, c('Family'))
PP_tax <- tibble::rownames_to_column(PP_tax, "OTU")

PP_OTUTAX <- merge(PP_tax, PP_OTU, by="OTU")
PP_OTUTAX = subset(PP_OTUTAX, select = -c(OTU))
PP_OTUTAX <- data.frame(PP_OTUTAX[,-1], row.names = PP_OTUTAX[,1])

# check metadata
head(sample_data(PP_family))

# create dataframe
PP_metadata <- data.frame(sample_data(PP_family))

PP_res <- Maaslin2(
  input_data = PP_OTUTAX,
  input_metadata = PP_metadata,
  output = "./Maaslin2_output",
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  fixed_effects = c("Time", "Age", "Sex"), 
  reference = c("Time, WK0", "Sex,Male"), 
  random_effects = ("subject"),
  min_prevalence = 0.0,
  min_abundance = 0.0,
  plot_heatmap = TRUE,
  plot_scatter = TRUE)

PP_res_df <- PP_res$results %>%
  dplyr::filter(metadata == "Time") %>%
  dplyr::select(-qval) %>%
  dplyr::arrange(pval) %>%
  dplyr::mutate(qval = p.adjust(pval, method = "fdr"))

head(PP_res_df)

# save results
write.csv(PP_res_df,"Maaslin2_output\\IFP_family_results.csv", row.names = FALSE)


### CR (HH) Group ###

# subset for HH group
HH_family <- subset_samples(ps_family, PP_HH == "HH")

# create OTU dataframe with annotations
HH_OTU <- data.frame(otu_table(HH_family))
HH_tax <- data.frame(tax_table(HH_family))

library(tibble)
HH_OTU <- tibble::rownames_to_column(HH_OTU, "OTU")

HH_tax <- select(HH_tax, c('Family'))
HH_tax <- tibble::rownames_to_column(HH_tax, "OTU")

HH_OTUTAX <- merge(HH_tax, HH_OTU, by="OTU")
HH_OTUTAX = subset(HH_OTUTAX, select = -c(OTU))
HH_OTUTAX <- data.frame(HH_OTUTAX[,-1], row.names = HH_OTUTAX[,1])

# check metadata
head(sample_data(HH_family))

# create dataframe
HH_metadata <- data.frame(sample_data(HH_family))

HH_res <- Maaslin2(
  input_data = HH_OTUTAX,
  input_metadata = HH_metadata,
  output = "./Maaslin2_output",
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  fixed_effects = c("Time", "Age", "Sex"),
  reference = c("Time, WK0", "Sex,Male"), 
  random_effects = ("subject"),
  min_prevalence = 0.0,
  min_abundance = 0.0,
  plot_heatmap = TRUE,
  plot_scatter = TRUE)

HH_res_df <- HH_res$results %>%
  dplyr::filter(metadata == "Time") %>%
  dplyr::select(-qval) %>%
  dplyr::arrange(pval) %>%
  dplyr::mutate(qval = p.adjust(pval, method = "fdr"))

head(HH_res_df)

# save results
write.csv(HH_res_df,"Maaslin2_output\\CR_family_results.csv", row.names = FALSE)
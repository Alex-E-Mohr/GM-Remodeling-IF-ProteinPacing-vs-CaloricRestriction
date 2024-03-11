# load packages
library(Maaslin2); packageVersion("Maaslin2")     
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse") 

# Load phyloseq object
phylo <- readRDS("./phyloseq objects/phyloseq_object.rds")

# remove taxa not seen >5 in 30% of samples
(phylo <- filter_taxa(phylo, function(x) sum(x > 5) > (0.3*length(x)), TRUE))   

# Agglomerate to Genus level
(ps_genus <- tax_glom(phylo, "Genus")) #detected 127 taxa

# Here we look at what we have at the genus level
table(tax_table(ps_genus)[, 6])

# remove 56 "uncultured"
ps_genus <- subset_taxa(ps_genus, Genus != "uncultured")


### IF-P (PP) Group ###

# subset for PP group
PP_genus <- subset_samples(ps_genus, PP_HH == "PP")

# create OTU dataframe with annotations
PP_OTU <- data.frame(otu_table(PP_genus))
PP_tax <- data.frame(tax_table(PP_genus))

library(tibble)
PP_OTU <- tibble::rownames_to_column(PP_OTU, "OTU")

PP_tax <- select(PP_tax, c('Genus'))
PP_tax <- tibble::rownames_to_column(PP_tax, "OTU")

PP_OTUTAX <- merge(PP_tax, PP_OTU, by="OTU")
PP_OTUTAX = subset(PP_OTUTAX, select = -c(OTU))
PP_OTUTAX <- data.frame(PP_OTUTAX[,-1], row.names = PP_OTUTAX[,1])

# check metadata
head(sample_data(PP_genus))

# create dataframe
PP_metadata <- data.frame(sample_data(PP_genus))

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
# write.csv(PP_res_df,"Maaslin2_output\\IFP_genus_results.csv", row.names = FALSE)


### CR (HH) Group ###

# subset for HH group
HH_genus <- subset_samples(ps_genus, PP_HH == "HH")

# create OTU dataframe with annotations
HH_OTU <- data.frame(otu_table(HH_genus))
HH_tax <- data.frame(tax_table(HH_genus))

library(tibble)
HH_OTU <- tibble::rownames_to_column(HH_OTU, "OTU")

HH_tax <- select(HH_tax, c('Genus'))
HH_tax <- tibble::rownames_to_column(HH_tax, "OTU")

HH_OTUTAX <- merge(HH_tax, HH_OTU, by="OTU")
HH_OTUTAX = subset(HH_OTUTAX, select = -c(OTU))
HH_OTUTAX <- data.frame(HH_OTUTAX[,-1], row.names = HH_OTUTAX[,1])

# check metadata
head(sample_data(HH_genus))

# create dataframe
HH_metadata <- data.frame(sample_data(HH_genus))

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
# write.csv(HH_res_df,"Maaslin2_output\\CR_genus_results.csv", row.names = FALSE)

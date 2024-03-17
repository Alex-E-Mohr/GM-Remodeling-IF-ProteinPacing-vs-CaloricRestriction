# load packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")


#Building MetaPhlAn species abundance ps object
s_abund <- read.table("./species_relab.txt", sep = "\t", header = T, 
                    skip = 0, comment.char = "")

s_tax_tab <- s_abund %>%
  dplyr::select(taxonomy) %>%
  tidyr::separate(taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|") %>%
  dplyr::mutate(spec_row = Species) %>%
  tibble::column_to_rownames(var = "spec_row")

s_tax_tab <- data.frame(row.names = row.names(s_tax_tab),
                        Kingdom = str_replace(s_tax_tab[,1], "k__",""),
                        Phylum = str_replace(s_tax_tab[,2], "p__",""),
                        Class = str_replace(s_tax_tab[,3], "c__",""),
                        Order = str_replace(s_tax_tab[,4], "o__",""),
                        Family = str_replace(s_tax_tab[,5], "f__",""),
                        Genus = str_replace(s_tax_tab[,6], "g__",""),
                        Species = str_replace(s_tax_tab[,7], "s__",""),
                        stringsAsFactors = FALSE)

s_otu_tab <- s_abund %>%
  tidyr::separate(taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|") %>%
  dplyr::select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>%
  tibble::column_to_rownames(var = "Species")

names(s_otu_tab) <- gsub(names(s_otu_tab), pattern = "_taxonomic_profile", replacement = "") 

head(colSums(s_otu_tab))
s_otu_tab <- s_otu_tab / 100  #convert to proportion with unit sum of 1
head(colSums(s_otu_tab))

s_meta <- read.table(file = "./metadata_casestudy.txt", sep = "\t", header = T, row.names = 1)

(ps_species <- phyloseq(sample_data(s_meta),
                             otu_table(s_otu_tab, taxa_are_rows = TRUE),
                             tax_table(as.matrix(s_tax_tab))))

taxa_names(ps_species) <- gsub("s__", "", taxa_names(ps_species))    #cleaning up names



read_count_df <- read.table("./read_count_df.txt", sep = "\t", header = T, 
                            skip = 0, comment.char = "")
read_count_df <- read_count_df %>%
  dplyr::select(Sample, `number_reads`)


otu_tab <- data.frame(otu_table(ps_species))
otu_tab <- otu_tab %>%
  t(.) %>%
  data.frame(.) %>%
  rownames_to_column(var = "Sample")

otu_tab %>% select(-Sample) %>% rowSums(.) #confirm each row sums as expected

otu_tab <- left_join(otu_tab, read_count_df)

otu_tab <- otu_tab %>% 
  mutate_at(vars(-number_reads, -Sample), funs(. * number_reads)) %>%   #convert to counts and round up
  mutate_at(vars(-number_reads, -Sample), funs(ceiling(.))) %>%
  select(-number_reads) %>%
  column_to_rownames(var = "Sample")

otu_tab <- t(otu_tab)
otu_table(ps_species) <- otu_table(otu_tab, taxa_are_rows = TRUE)

read_count_df
sample_sums(ps_species)

saveRDS(ps_species, "./phyloseq_casestudy.rds")
rm(list=ls())

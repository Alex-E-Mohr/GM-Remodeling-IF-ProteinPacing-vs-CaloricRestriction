# load packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")


#Building MetaPhlAn species abundance ps object
otu <- read.table(file = "./metabolite_OTU_responder.txt", sep = "\t", header = T, 
                  skip = 0, comment.char = "")
otu <- data.frame(otu[,-1], row.names = otu[,1])


head(colSums(otu))
#otu_tab <- otu / 100  #convert to proportion with unit sum of 1
#head(colSums(otu_tab))

meta <- read.table(file = "./metadata_casestudy.txt", sep = "\t", header = T, row.names = 1)

tax <- read.table(file = "./metabolite_TAX_responder.txt", header = T, row.names = 1)


(ps_metabolite <- phyloseq(sample_data(meta),
                             otu_table(otu, taxa_are_rows = TRUE),
                             tax_table(as.matrix(tax))))

taxa_names(ps_metabolite) <- gsub("_", " ", taxa_names(ps_metabolite))    #cleaning up names


table(tax_table(ps_metabolite)[, 2])


saveRDS(ps_metabolite, "./phyloseq_casestudy_metab.rds")
rm(list=ls())

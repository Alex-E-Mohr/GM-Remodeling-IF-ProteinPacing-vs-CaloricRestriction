#Load packages
library(tidyverse); packageVersion("tidyverse")
library(BiocManager); packageVersion("BiocManager")
library(phyloseq); packageVersion("phyloseq")
library(microbiome); packageVersion("microbiome")
library(ggplot2); packageVersion("ggplot2")

#Build phyloseq project

otu <- read.table(file = "./Metabolite_table.tsv", sep = "\t", header = T, row.names = 1, 
                  skip = 0, comment.char = "")

metadata <- read.table(file = "./metadata.txt", sep = "\t", header = T, row.names = 1)
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata)

# merge the data
phyloseq_object <- phyloseq(OTU, SAMPLE)
phyloseq_object

# check sample names
sample_names(phyloseq_object)

# check variable names
sample_variables(phyloseq_object)

# check the the group variable in more detail
sample_data(phyloseq_object)$Group
table(sample_data(phyloseq_object)$Group)

# check the the PP_HH variable in more detail
sample_data(phyloseq_object)$PP_HH
table(sample_data(phyloseq_object)$PP_HH)

# check the appearance of the metadata
metadata <- data.frame(sample_data(phyloseq_object))
head(metadata)

##### Save the phyloseq object #####

saveRDS(phyloseq_object, "./phyloseq_metabolome.rds")


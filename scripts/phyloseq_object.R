#Load packages
library(tidyverse); packageVersion("tidyverse")
library(BiocManager); packageVersion("BiocManager")
library(phyloseq); packageVersion("phyloseq")
library(microbiome); packageVersion("microbiome")
library(ggplot2); packageVersion("ggplot2")

#Build phyloseq project

otu <- read.table(file = "./table.tsv", sep = "\t", header = T, row.names = 1, 
                  skip = 1, comment.char = "")

taxonomy <- read.table(file = "./taxonomy.tsv", sep = "\t", header = T ,row.names = 1)

# clean the taxonomy, Silva format
tax <- taxonomy %>%
  select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ")

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "d__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)


tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}

metadata <- read.table(file = "./metadata.txt", sep = "\t", header = T, row.names = 1)
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax.clean))
SAMPLE <- sample_data(metadata)
TREE = read_tree("./tree.nwk")

# merge the data
phyloseq_object <- phyloseq(OTU, TAX, SAMPLE, TREE)
phyloseq_object

# check number of samples
nsamples(phyloseq_object) # 117

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

# unfiltered ASVs = 3,406
ntaxa(phyloseq_object)

# taxa at kingdom level
table(tax_table(phyloseq_object)[, 1])

# remove 10 "Archaea", 1 "Eukaryota", and 3 "Unassigned features
phyloseq_object <- subset_taxa(phyloseq_object, Kingdom !="Archaea")
phyloseq_object <- subset_taxa(phyloseq_object, Kingdom !="Eukaryota")
phyloseq_object <- subset_taxa(phyloseq_object, Kingdom !="Unassigned")

# taxa at phyla level
table(tax_table(phyloseq_object)[, 2])

# remove 36 "unclassified Bacteria"
phyloseq_object <- subset_taxa(phyloseq_object, Phylum !="Unclassified Bacteria")

# taxa at class level
table(tax_table(phyloseq_object)[, 3])

# taxa at order level
table(tax_table(phyloseq_object)[, 4])

# remove 20 "Chloroplast"
phyloseq_object <- subset_taxa(phyloseq_object, Order != "Chloroplast")

# Here we look at what we have at the family level
table(tax_table(phyloseq_object)[, 5])

# remove 10 "Mitochondria"
phyloseq_object <- subset_taxa(phyloseq_object, Family != "Mitochondria")

# ASVs = 3,326
ntaxa(phyloseq_object)

# examine number of reads for each sample to identify potentially problematic 
# samples and plot their distribution.

sample_sums(phyloseq_object)

sort(sample_sums(phyloseq_object))

min(sample_sums(phyloseq_object))

mean(sample_sums(phyloseq_object))

max(sample_sums(phyloseq_object))

# The number of reads per sample ranges from 6,586 to 219,065, with a mean of 38,886.83.

hist(sample_sums(phyloseq_object), main="Histogram: Read Counts", xlab="Total Reads", 
     border="black", col="grey", las=1, breaks=12)

# add read count
reads_sample <- microbiome::readcount(phyloseq_object)
sample_data(phyloseq_object)$reads_sample <- reads_sample
sample_data(phyloseq_object)

# examine total reads per sample and the overall distribution
readsumsdf = data.frame(nreads = sort(taxa_sums(phyloseq_object), TRUE), sorted = 1:ntaxa(phyloseq_object), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(phyloseq_object), 
                                                        TRUE), sorted = 1:nsamples(phyloseq_object), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

# The distribution looks typical for the distribution from an amplicon-based microbiome census

##### Save the phyloseq object #####

saveRDS(phyloseq_object, "./phyloseq_object.rds")
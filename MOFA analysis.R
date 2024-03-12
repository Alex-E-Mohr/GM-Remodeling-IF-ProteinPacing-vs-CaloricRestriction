# load packages
library(mia)
library(MOFA2)
library(phyloseq)
library(basilisk)
library(ggplot2)
library(patchwork)
library(rstatix)
library(ggExtra)
library(reshape)
library(plyr)
library(ggpubr)
library(scales)
library(ggrepel)


### load microbiome data ###

# Load phyloseq object
phylo <- readRDS("./phyloseq_object.rds")

# remove taxa not seen >2 times in 10% of samples
(phylo <- filter_taxa(phylo, function(x) sum(x > 2) > (0.1*length(x)), TRUE))

# convert phyloseq to TSE object
microbiome_TSE <- makeTreeSummarizedExperimentFromPhyloseq(phylo) 
microbiome_TSE

# check metadata
colData(microbiome_TSE)


### load plasma metabolome data ###
# note: log10 transformed and Pareto scaled already
counts <- read.csv(file = './Metabolite_Table.csv', row.names = 1)
# samples <- read.table(file = './metadata.txt', sep = "\t", header = T, row.names = 1)
samples <- sample_data(phylo)

# make sure the number of participants match for counts and samples

# Convert into right format
counts <- as.matrix(counts)
assays <-  SimpleList(concs = counts)
colData <- DataFrame(colData)

# Create a TreeSE
metabolite_TSE <- TreeSummarizedExperiment(assays = assays, colData = samples)
metabolite_TSE

# check metadata
colData(metabolite_TSE)


### load plasma cytokine data ###
# counts2 <- read.csv(file = '/inflammation.csv', row.names = 1)

# Convert into right format
# counts2 <- as.matrix(counts2)
# assays2 <-  SimpleList(concs = counts2)

# Create a TreeSE
# cytokines_TSE <- TreeSummarizedExperiment(assays = assays2, colData = samples)
# cytokines_TSE


### combine these experiments into a MultiAssayExperiment (MAE)

# Create an ExperimentList that includes experiments
experiments <- ExperimentList(microbiome = microbiome_TSE, 
                              metabolome = metabolite_TSE)

# Create a MAE
mae <- MultiAssayExperiment(experiments = experiments,
                            colData = samples)
mae

experiments(mae)
colData(mae)

# Microbiome data
mae[[1]]

# Metabolite data
mae[[2]]

# Cytokine data
# mae[[3]]


### clean/transform data where necessary

# Removing duplicates at the microbiome data
# which are also in form e.g. "Ambiguous" and "uncultured" taxa
mae[[1]] <- mae[[1]][!duplicated(rownames(assay(mae[[1]]))), ]

# Agglomerate microbiome data at genus level
mae[[1]] <- agglomerateByPrevalence(mae[[1]], rank = "Genus")

# Transforming microbiome data with rclr
mae[[1]] <- transformCounts(mae[[1]], method = "relabundance")
mae[[1]] <- transformCounts(mae[[1]], assay_name = "relabundance", method = "rclr")

# Transforming cytokine data with z-transform
# mae[[3]] <- transformFeatures(mae[[3]], assay_name = "concs", method = "z")

# Removing assays no longer needed
assay(mae[[1]], "counts") <- NULL
assay(mae[[1]], "relabundance") <- NULL
# assay(mae[[3]], "concs") <- NULL

colData(mae)$PP_HH


### Build MOFA model

model <- create_mofa_from_MultiAssayExperiment(mae,
                                               extract_metadata = TRUE)
model

# Multi-group mode requested.
# Two important remarks:
# - The aim of the multi-group framework is to identify the sources of variability 
# *within* the groups. If your aim is to find a factor that 'separates' the groups, 
# you DO NOT want to use the multi-group framework. Please see the FAQ on the MOFA2 webpage.
#- It is important to account for the group effect before selecting highly variable 
# features (HVFs). We suggest that either you calculate HVFs per group and then 
# take the union, or regress out the group effect before HVF selection


# Visualize data structure
plot_data_overview(model)

### Prepare MOFA object
?prepare_mofa

# specify model options
model_opts <- get_default_model_options(model)
head(model_opts) # all likelihoods are gaussian

# specify training options
train_opts <- get_default_training_options(model)
train_opts$drop_factor_threshold <- 0.02 # drop factors from the model that describe <2% of the variance
head(train_opts)

# prepare the model
model_prepared <- prepare_mofa(model, model_options = model_opts, training_options = train_opts)


### Build MOFA model


# Train the model
model_trained <- run_mofa(model_prepared, use_basilisk = TRUE)
model_trained

# Add sample metadata to the model
clinical <- read.csv(file = './clinical_correlations.csv')
metadata <- cbind(sample = rownames(samples), samples)
metadata <- merge(clinical, metadata)
samples_metadata(model_trained) <- metadata

# Save mofa file for future use
# saveRDS(model_trained, "./MOFA_overall.rds")

# Load MOFA object
model_trained <- readRDS("./MOFA_overall.rds")


### explore data


### Correlation between factors ###
# make sure factors are largely uncorrelated. Greater correlation between factors
# suggests a poor model fit.
plot_factor_cor(model_trained)

### Total variance explained per view ###
# in this data set using k=8 factors the metabolome explains ~2x the variance of
# the microbiome.
plot_variance_explained(model_trained, plot_total = TRUE)

# Total variance explained per view
head(get_variance_explained(model_trained)$r2_total[[1]])

### Variance decompostion by Factor ###
# here we have formed 8 factors base on our >2% cutoff for explained variance.
# Factors 1 and 6 capture variability that is present across both the microbiome
# and metabolome. It's etiology is likely important for the weight loss intervention.
# Factors 2 and 5 capture a very strong source of variation that is exclusive to the metabolome
# Factors 3 and 4 capture some source of variation that is exclusive to the microbiome
plot_variance_explained(model_trained, max_r2=15) + coord_flip()

# Variance explained for every factor in per view
(get_variance_explained(model_trained)$r2_per_factor[[1]])

# plot each factor by group. Test by Wilcox
plot_factor(model_trained, 
            factor = 1:8, 
            color_by = "PP_HH", 
            dot_size = 3,
            dodge = TRUE,
            stroke = 0.4,
            add_violin = T,
            violin_alpha = 0.25,
            add_boxplot = T,
            boxplot_alpha = 0.25) +
  scale_fill_manual(values=category.colors) +
  stat_compare_means(size=3, aes(label = paste0("p = ", after_stat(p.format))), label.y = 3.2) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

### Association analysis
correlate_factors_with_covariates(model_trained, 
                                  covariates = c("Sex","BMI","Age", "PP_HH", "Time_num"), 
                                  plot="log_pval")

cormatrix <- correlate_factors_with_covariates(model_trained,
                                               covariates = c("PA","kcal","cho","fiber",
                                                              "sugar","fat","pro","sodium",
                                                              "wtkg","waist","bmi","lbm",
                                                              "lbm_bw","ffm","ffm_bw",
                                                              "fatmass","perc_bf",
                                                              "vaf","sat","android"),
                                               plot="r")

cormatrix

### Plot feature weights ###
# weights provide a score for each feature on each factor. The sign of the weights
# indicates the direction of the effect: a positive weight indicated that the feature
# has higher levels in the cells with positive factor values, and vice-versa. Features 
# with large positive values will be more expressed in the HH samples, whereas 
# features with large negative values will be more expressed in the PP samples.

category.colors <- c(
  "PP" = "darkturquoise", 
  "HH" = "magenta3")

# Factor 1: Microbiome
F1_micro <- plot_weights(model_trained,
                         view = 1,
                         factor = 1,
                         text_size = 1,
                         dot_size = 1,
                         nfeatures = 10,     # Top number of features to highlight
                         scale = T) +           # Scale weights from -1 to 1
  theme(text=element_text(color="black", size = 5))

F1_micro 

plot_top_weights(model_trained,
                 view = 1,
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T)           # Scale weights from -1 to 1

plot_data_scatter(model_trained, 
                  view = 1,
                  factor = 1,  
                  features = 5,
                  sign = "positive",
                  color_by = "PP_HH") +
  labs(y="CLR abundance") +
  scale_fill_manual(values=category.colors)# most significant Faecalibacterium, Fusicatenibacter


# Factor 1: Metabolome
F1_metab <- plot_weights(model_trained,
                         view = 2,
                         factor = 1,
                         text_size = 1,
                         dot_size = 1,
                         nfeatures = 10,     # Top number of features to highlight
                         scale = T) +           # Scale weights from -1 to 1
  theme(text=element_text(color="black", size = 5))

F1_metab

plot_top_weights(model_trained,
                 view = 2,
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T)           # Scale weights from -1 to 1

plot_data_scatter(model_trained, 
                  view = 2,
                  factor = 1,  
                  features = 8,
                  sign = "positive",
                  color_by = "PP_HH") +
  labs(y="Log10 abundance") + # most significant Faecalibacterium, Fusicatenibacter
  scale_fill_manual(values=category.colors)


# Factor 6: Microbiome
plot_weights(model_trained,
             view = 1,
             factor = 6,
             text_size = 0,
             dot_size = 2,
             nfeatures = 10,     # Top number of features to highlight
             scale = T)           # Scale weights from -1 to 1

plot_top_weights(model_trained,
                 view = 1,
                 factor = 6,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T)           # Scale weights from -1 to 1

plot_data_scatter(model_trained, 
                  view = 1,
                  factor = 6,  
                  features = 6,
                  sign = "positive",
                  color_by = "PP_HH") + 
  labs(y="CLR abundance") +
  scale_fill_manual(values=category.colors)


# Factor 6: Metabolome
plot_weights(model_trained,
             view = 2,
             factor = 6,
             text_size = 0,
             dot_size = 2,
             nfeatures = 10,     # Top number of features to highlight
             scale = T)           # Scale weights from -1 to 1

plot_top_weights(model_trained,
                 view = 2,
                 factor = 6,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T)           # Scale weights from -1 to 1

plot_data_scatter(model_trained, 
                  view = 2,
                  factor = 6,  
                  features = 10,
                  sign = "positive",
                  color_by = "PP_HH") + 
  labs(y="Log10 abundance") +
  scale_fill_manual(values=category.colors)


F1_micro + F1_metab

#png
ggsave("./Fig3de.png", height = 6, width = 9, units = "cm")


### Inspection of combinaions of Factors ###
# here we inspect the two factors that both explained variance in the microbiome
# and metabolome. This plot appears to show separation between the two groups. I.e.,
# difference in (multi-omic) molecular profile between the two dietary exposures.

plot_factors(model_trained, 
             factors = c(1,6), 
             color_by = "PP_HH",
             shape_by = "Time",
             dot_size = 2.5,
             show_missing = T) + 
  scale_fill_manual(values=category.colors) +
  geom_hline(yintercept=-1.1, linetype="dashed") +
  geom_vline(xintercept=(-0.25), linetype="dashed")


plot_factors(model_trained, 
             factors = c(1,6), 
             color_by = "PP_HH", 
             dot_size = 4) + 
  scale_fill_manual(values=category.colors) +
  stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) +
  scale_color_manual(values=category.colors)

### Prediction analysis
suppressPackageStartupMessages(library(randomForest))

# Prepare data
df <- as.data.frame(get_factors(model_trained, factors=c(1,6))[[1]])
df


# Train the model for diet intervention
df$group <- as.factor(model_trained@samples_metadata$PP_HH)
model.DIET <- randomForest(group ~ ., data=df[!is.na(df$group),], ntree=500)
model.DIET
varImpPlot(model.DIET)
df$group <- NULL

# Do predictions
model_trained@samples_metadata$DIET.pred <- stats::predict(model.DIET, df)

# Plot the predictions
model_trained@samples_metadata$DIET.pred_logical <- c("True",
                                                      "Predicted")[as.numeric(is.na(model_trained@samples_metadata$PP_HH))+1]

plot_factors(model_trained, 
             factors = c(1,6), 
             color_by = "DIET.pred",
             shape_by = "DIET.pred_logical",
             dot_size = 2.5,
             show_missing = T) + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")


### Extracting data for downstream analysis and plotting


str(model_trained)
### variance ###
variance <- get_variance_explained(model_trained,
                                   as.data.frame = TRUE)


# total variance
total_var <- variance$r2_total

t_var <- ggplot(total_var, aes(y=value, x=view)) +
  geom_bar(stat="identity") +
  scale_y_continuous(limits = c(0,40), expand = c(0, 0)) +
  xlab("") +
  ylab("Variance explained (%)") +
  theme_classic() +
  theme(text = element_text(size = 6, color = "black"),
        axis.text = element_text(color = "black")) +
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

t_var

# factor variance 
factor_var <- variance$r2_per_factor

level_order1 <- factor(factor_var$factor, level = c("Factor8","Factor7","Factor6",
                                                    "Factor5","Factor4","Factor3",
                                                    "Factor2","Factor1"),
                       label = c("Factor 8","Factor 7","Factor 6",
                                 "Factor 5","Factor 4","Factor 3",
                                 "Factor 2","Factor 1"))

library(viridis)

f_var <- ggplot(factor_var, aes(x=view, y=level_order1, fill=value)) +
  geom_tile() +
  geom_text(aes(label=scales::percent(value*0.01, accuracy=0.01)), color = "white",
            size = 1.5, fontface = "bold") +
  #scale_fill_gradient(low="white", high = "blue") +
  scale_fill_viridis() +
  theme_bw() +
  theme(text = element_text(size = 6),
        axis.text = element_text(color = "black")) +
  #theme(legend.title.align=0.5, legend.position="bottom") +
  #guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5,
  #                              title="Variance (%)"),
  #       size = guide_legend(title.position="top", title.hjust = 0.5)) +
  theme(axis.title = element_blank()) +
  scale_x_discrete(labels=c("microbiome" = "Microbiome", "metabolome" = "Metabolome"),
                   guide = guide_axis(angle = 45)) +
  #theme(legend.key.size = unit(0.3, 'cm'))
  theme(legend.position = "none")

f_var


# combine plots
library(patchwork)

t_var / plot_spacer() / f_var +
  plot_layout(heights = c(1, -0.3, 2))

#png
ggsave("./Fig3a.png", height = 9, width = 3, units = "cm")

### weights ###
weights <- get_weights(model_trained, 
                       views = "all", 
                       factors = "all", 
                       as.data.frame = TRUE)
head(weights)

### factors ###
factors <- get_factors(model_trained,
                       factors = "all",
                       as.data.frame = TRUE)
head(factors)

factors_16 <- factors %>%
  subset(factor!= "Factor2" & 
           factor!="Factor3" &
           factor!="Factor4" &
           factor!="Factor5" &
           factor!="Factor7" &
           factor!="Factor8")

factor16_df <- cast(factors_16, sample ~ factor)
factors16_df <- merge(factor16_df, metadata)

factors16_df$PP_HH <- revalue(x = factors16_df$PP_HH, 
                              c("HH" = "CR", "PP" = "IF-P"))

order <- factor(factors16_df$PP_HH, level=c('IF-P', 'CR'))

plot_16 <- ggplot(factors16_df, aes(x = Factor1, y = Factor6, 
                                    fill = order)) + 
  stat_ellipse(geom = "polygon",
               aes(fill = order), 
               alpha = 0.25) +
  geom_point(size = 2, shape = 21, stroke = 0.5,
             alpha = 1, color = "black") +
  scale_fill_manual(name=NULL,
                    values=c("darkturquoise", "magenta3")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", 
                                         color = "black",
                                         linewidth = 0.25),
        legend.position = c(0.913, 0.077)) +
  theme(text = element_text(size = 6, color = "black"),
        axis.text = element_text(size = 6, color = "black")) +
  xlab("Factor 1") + ylab("Factor 6") +
  theme(legend.key.size = unit(0.2, 'cm'))

plot_16

#png
ggsave("./Fig3c.png", height = 7, width = 7, units = "cm")


barplot_1 <- ggplot(factors16_df, aes(x = order, y = Factor1, fill = PP_HH)) +
  geom_boxplot(lwd=0.25, outlier.shape = NA, color = "black") +
  scale_fill_manual(name=NULL,
                    values=c("magenta3", "darkturquoise")) +
  scale_color_manual(name=NULL,
                     values=c("black", "black")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(text = element_blank(),
        title = element_blank(), 
        axis.ticks = element_blank(),
        line = element_blank())

barplot_1

ggsave("./Fig3c_boxplot1.png", height = 6, width = 1.5, units = "cm")


barplot_6 <- ggplot(factors16_df, aes(x = order, y = Factor6, fill = PP_HH)) +
  geom_boxplot(lwd=0.25, outlier.shape = NA, color = "black") +
  scale_fill_manual(name=NULL,
                    values=c("magenta3", "darkturquoise")) +
  scale_color_manual(name=NULL,
                     values=c("black", "black")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(text = element_blank(),
        title = element_blank(), 
        axis.ticks = element_blank(),
        line = element_blank())

barplot_6

ggsave("./Fig3c_boxplot6.png", height = 6, width = 1.5, units = "cm")


### data ###
data <- get_data(model_trained, 
                 views = "all", 
                 as.data.frame = TRUE)
head(data)



### correlation of factors with covariates ###
factor_df <- cast(factors, sample ~ factor)
factor_cor <- merge(factor_df, metadata)
factor_cor = subset(factor_cor, select = -c(sample, subject, PP_HH, time,
                                            Group, Time, Time_num, Race,
                                            BMI, BMI_Class, Sex, Age, ffm, lbm,
                                            BSS, pH, biomass, Acetate_mMol.g,
                                            Propionate_mMol.g, Isobutyrate_mMol.g,
                                            Butyrate_mMol.g, Valerate_mMol.g,
                                            reads_sample))

cor.mat <- factor_cor %>% cor_mat(method = "spearman")
p_values <- cor.mat %>% cor_get_pval()

cor.mat %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = F)

# Gather/collapse correlation matrix into long format
select_cor<- cor.mat %>% 
  cor_gather() %>%
  subset(var2!= "Factor1" & 
           var2!="Factor2" &
           var2!="Factor3" &
           var2!="Factor4" &
           var2!="Factor5" &
           var2!="Factor6" &
           var2!="Factor7" &
           var2!="Factor8" &
           var1!="wtkg" &
           var1!="waist" &
           var1!="bmi" &
           var1!="lbm" &
           var1!="lbm_bw" &
           var1!="ffm" &
           var1!="ffm_bw" &
           var1!="fatmass" &
           var1!="perc_bf" &
           var1!="vaf" &
           var1!="sat" &
           var1!="android" &
           var1!="PA" &
           var1!="kcal" &
           var1!="cho" &
           var1!="fiber" &
           var1!="sugar" &
           var1!="fat" &
           var1!="pro" &
           var1!="sodium")

cor_df <- select_cor %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj",
                   cutpoints = c(0, 1e-04, 0.001, 0.01, 0.1, 1),
                   symbols = c("****", "***", "**", "*", ""))


# Plot

level_order1 <- factor(cor_df$var1, level = c("Factor8","Factor7","Factor6",
                                              "Factor5","Factor4","Factor3",
                                              "Factor2","Factor1"),
                       label = c("Factor 8","Factor 7","Factor 6",
                                 "Factor 5","Factor 4","Factor 3",
                                 "Factor 2","Factor 1"))

level_order2 <- factor(cor_df$var2, level = c("wtkg","bmi","waist", "perc_bf","fatmass",
                                              "vaf","sat","android","ffm_bw","lbm_bw",
                                              "PA","kcal","cho","sugar","fiber","fat",
                                              "pro","sodium"),
                       label = c("BW","BMI","WC", "BF%", "FM", "VAF", "SAT", "Android",
                                 "FFM/BW","LBM/BW","PA", "KCAL", "CHO", "Sugar", "Fiber",
                                 "Fat","PRO", "Sodium"))

plot_corr <- cor_df %>%  
  ggplot(aes(x=level_order2, y = level_order1, color = cor, size = -log10(p))) + 
  geom_point() +
  labs(size = "-log10"~(italic(P)~"-value")) +
  #scale_size("p-values", trans="log10", range=c(10, 1)) +
  scale_color_gradient2(low="blue", mid="white",
                        high="red", space = "Lab", name = "Spearman rho") +
  theme_classic()  +
  theme(text = element_text(size = 6, color = "black")) +
  theme(axis.text=element_text(colour="black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
  theme(axis.title = element_blank()) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label = p.adj.signif, color = cor), col = "black", nudge_y = 0, 
            size = 3,
            fontface = "bold",
            show.legend = TRUE) +
  theme(legend.title.align=0.5, legend.position="top") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))

plot_corr

#png
ggsave("./Fig3b.png", height = 8, width = 9, units = "cm")


################################################################################

### Factor 1 plot ###
F1_weights <- get_weights(model_trained, 
                          views = "all", 
                          factors = c(1), 
                          as.data.frame = TRUE)
head(F1_weights)

#long to wide format
# cast(micro_weights, feature ~ factor)

F1_weights <-reshape(F1_weights, idvar = "feature", timevar = "factor", direction = "wide")
row.names(F1_weights) <- 1:191

# scale the weights
F1_WTscaled <- rescale(F1_weights$value.Factor1, to = c(-1, 1))
F1_WTscaled <- data.frame(F1_WTscaled)

F1_WTplot <- cbind(F1_weights, F1_WTscaled)

F1_WTplot

write.csv(F1_WTplot, "./MOFA_Factor1_weights.csv", row.names=FALSE)

F1_WTplot <- read.csv(file = './MOFA_Factor1_weights.csv')

ggplot(F1_WTplot, aes(x = view.Factor1, y = F1_WTscaled, fill = F1_WTscaled)) +
  geom_hline(yintercept = 0, linetype="solid", color="gray", linewidth = .5) +
  geom_jitter(shape = 21, size=1.5, stroke = 0.5, position=position_jitter(0.25)) +
  theme_bw() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 6)) +
  scale_fill_viridis() +
  #geom_text_repel(aes(label=ifelse(Annotation>0,as.character(feature),'')), 
  #                color = "black", 
  #                size=1, max.overlaps = Inf, box.padding = 1, 
  #                nudge_y = 0.1, min.segment.length = 0,
  #                segment.size = 0) +
  coord_flip() +
  ggtitle("Factor 1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_blank()) +
  ylab("Scaled Feature Weights") +
  ylim(-1.2, 1.2) +
  guides(fill = guide_colourbar(barheight = 0.3, barwidth = 13)) +
  theme(legend.title = element_blank(), legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

#png
ggsave("./Fig3d.png", height = 7, width = 9, units = "cm")


### Factor 6 plot ###
F6_weights <- get_weights(model_trained, 
                          views = "all", 
                          factors = c(6), 
                          as.data.frame = TRUE)
head(F6_weights)

#long to wide format
# cast(micro_weights, feature ~ factor)

F6_weights <-reshape(F6_weights, idvar = "feature", timevar = "factor", direction = "wide")
row.names(F6_weights) <- 1:191

# scale the weights
F6_WTscaled <- rescale(F6_weights$value.Factor6, to = c(-1, 1))
F6_WTscaled <- data.frame(F6_WTscaled)

F6_WTplot <- cbind(F6_weights, F6_WTscaled)

F6_WTplot

write.csv(F6_WTplot, "./MOFA_Factor6_weights.csv", row.names=FALSE)

F6_WTplot <- read.csv(file = './MOFA_Factor6_weights.csv')

ggplot(F6_WTplot, aes(x = view.Factor6, y = F6_WTscaled, fill = F6_WTscaled)) +
  geom_hline(yintercept = 0, linetype="solid", color="gray", linewidth = .5) +
  geom_jitter(shape = 21, size=1.5, stroke = 0.5, position=position_jitter(0.25)) +
  theme_bw() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 6)) +
  scale_fill_viridis() +
  #geom_text_repel(aes(label=ifelse(Annotation>0,as.character(feature),'')), 
  #                color = "black", 
  #                size=1, max.overlaps = Inf, box.padding = 1, 
  #                nudge_y = 0.1, min.segment.length = 0,
  #                segment.size = 0) +
  coord_flip() +
  ggtitle("Factor 6") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_blank()) +
  ylab("Scaled Feature Weights") +
  ylim(-1.2, 1.2) +
  guides(fill = guide_colourbar(barheight = 0.3, barwidth = 13)) +
  theme(legend.title = element_blank(), legend.position="bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

#png
ggsave("./Fig3d.png", height = 7, width = 9, units = "cm")

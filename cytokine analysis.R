# Load packages
library(emmeans)
library(nlme)  
library(tidyverse) 
library(rcompanion)
library(ggplot2)

# Import data
inflammation <- read.csv(file = './inflammation.csv')
metadata <- read.table(file = "./metadata.txt", sep = "\t", header = TRUE)
inflam <- merge(metadata, inflammation)

# Normality Tests
plotNormalHistograms <- function(variable) {
  plotNormalHistogram(variable, main = deparse(substitute(variable)))
}

variables <- c("GM_CSF", "IFNy", "IL_1B", "IL_2", "IL_4", "IL_5", 
               "IL_6", "IL_8", "IL_10", "IL_12p70", "IL_13", 
               "IL_17A", "IL_23", "TNFa")

lapply(variables, function(variable) {
  plotNormalHistograms(inflam[[variable]])
  shapiro.test(inflam[[variable]])
  shapiro.test(log2(inflam[[variable]]))
})

# Logarithmic Transformations
log_variables <- paste0("log2_", variables)
inflam[, log_variables] <- lapply(variables, function(variable) log2(inflam[[variable]]))

# Data Cleaning
inflam_clean <- inflam[!(inflam$Sample.ID %in% c("SM24_01", "SM27_01", "SM37_02")),]

# Linear Mixed-Effects Models (LMEs)
lmes <- lapply(log_variables, function(log_variable) {
  lme_formula <- as.formula(paste(log_variable, "~ Time * PP_HH + Age + Sex"))
  lme(lme_formula, data = inflam_clean, random = ~ 1|as.factor(subject), correlation = corAR1())
})

# Pairwise Comparisons
emmeans_list <- lapply(lmes, function(lme_model) {
  emmeans(lme_model, pairwise ~ Time:PP_HH)
})

# Data Visualization
plot_variable <- function(variable, data) {
  ggplot(data, aes(x = time, y = .data[[variable]], fill = group, color = group)) +
    geom_boxplot(lwd = 0.5, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 1, stroke = 0.5,
                position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
    facet_wrap(~factor(group, levels = c('IF-P', 'CR'))) +
    geom_line(aes(group = id), size = 0.5, alpha = 0.2) +
    scale_fill_manual(name = NULL, values = c("magenta3", "darkturquoise")) +
    scale_color_manual(name = NULL, values = c("black", "black")) +
    labs(y = paste("Log2(", variable, " Concentration)", sep = ""), x = "Time (weeks)") +
    scale_x_discrete(labels = c("WK0" = "0", "WK4" = "4", "WK8" = "8")) +
    theme_bw() +
    theme(legend.position = "none", strip.background = element_rect(colour = "black", fill = "black", size = 1),
          strip.text = element_text(colour = 'white', size = 6, face = "bold"),
          axis.text = element_text(size = 6, color = "black"),
          axis.title = element_text(size = 6, color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    ylim(NA, max(data[[variable]], na.rm = TRUE) * 1.1)
}

plots <- lapply(log_variables[c(5, 7, 8, 12)], function(variable) {
  plot_variable(variable, inflam_clean)
})

# Save plots as PNG
ggsave(filename = "./figures/Fig2a.png", plot = do.call("grid.arrange", c(plots, ncol = 4)), 
       height = 6, width = 18, units = "cm")



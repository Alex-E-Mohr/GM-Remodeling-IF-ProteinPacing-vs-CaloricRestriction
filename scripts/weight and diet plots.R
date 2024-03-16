# Load necessary packages
library(tidyverse)
library(emmeans)
library(lmtest)
library(sandwich)
library(rcompanion)
library(ggpubr)
library(nlme)
library(ggrepel)
library(Rmisc)
library(rstatix)
library(reshape2)

# Load data
clinical_DF <- read.table(file = "./clinical_data.txt", sep = "\t", header = TRUE)

# Define order for factor levels
level_order <- factor(clinical_DF$PP_HH, level = c("IF-P","CR"))

# Daily kcal plot
kcal_daily <- ggplot(data = clinical_DF, aes(x = Time, y = KCAL_d, fill = level_order, color = level_order)) +
  geom_boxplot(lwd = 0.5, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 1, stroke = 0.5,
              position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.2)) +
  scale_fill_manual(name = NULL, values = c("darkturquoise", "magenta3")) +
  scale_color_manual(name = NULL, values = c("black", "black")) +
  labs(y = "Adjusted daily consumption (kcals)", x = "Time (weeks)") +
  scale_x_discrete(labels = c("WK0" = "0", "WK4" = "4", "WK8" = "8")) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"),
        text = element_text(size = 6),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(0.885, 0.925),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title.align = 0) +
  scale_y_continuous(limits = c(500, 4500)) +
  stat_compare_means(aes(group = PP_HH), label = "p.format",
                     method = "t.test", label.y = c(4250, 2300, 2700), size = 2) +
  theme(legend.key.size = unit(0.2, 'cm'))

# % body weight loss
bw_change <- ggplot(data = clinical_DF, aes(x = Time, y = wt_loss_per, group = subject, fill = PP_HH)) +
  geom_point(shape = 21, size = 1, stroke = 0.5) +
  geom_line(size = 0.5, alpha = 0.2) +
  facet_grid(~factor(PP_HH, levels = c('IF-P', 'CR'))) +
  labs(x = "Time (weeks)", y = expression(Delta *" % Body Weight")) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "black", size = 1),
        strip.text = element_text(colour = 'white', size = 6, face = "bold"),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        text = element_text(size = 6)) +
  scale_fill_manual(values = c("magenta3", "darkturquoise")) +
  scale_x_discrete(labels = c("WK0" = "0", "WK4" = "4", "WK8" = "8"))

# Combine plots
combined_plots <- kcal_daily + bw_change +
  plot_layout(widths = c(1, 1))

# Save as PNG
ggsave("./Figures/Fig1_bc.png", height = 8, width = 12, units = "cm")

# Reshape dataframe for macronutrient plot
macro.df <- melt(clinical_DF, id.vars = c('Sample.ID', 'PP_HH', 'Time'), measure.vars = c('CHO', 'SUG', 'FIB', 'FAT', 'PRO'))

# Rename Time labels
macro.df$Time <- recode(macro.df$Time, "WK0" = "Baseline", "WK4" = "Week 4", "WK8" = "Week 8")

# Define factor levels for macro nutrients
level_order <- factor(macro.df$PP_HH, level = c("IF-P", "CR"))

# Define function for daily macronutrient intake plot
macro_daily_plot <- function(data) {
  ggplot(data = data, aes(x = variable, y = value, fill = level_order, color = level_order)) +
    geom_boxplot(lwd = 0.5, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 1, stroke = 0.5,
                position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
    facet_grid(~Time) +
    scale_fill_manual(name = NULL, values = c("darkturquoise", "magenta3")) +
    scale_color_manual(name = NULL, values = c("black", "black")) +
    labs(y = "Daily Intake (g)", x = "Macronutrient") +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_rect(colour = "black", fill = "black", size = 1),
          strip.text = element_text(colour = 'white', size = 6, face = "bold"),
          axis.text = element_text(colour = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 6),
          axis.title.x = element_blank()) +
    scale_y_continuous(limits = c(0, 460)) +
    stat_compare_means(label = "p.signif", method = "t.test", size = 2, label.y = c(450, 210, 70, 210, 210))
}

# Generate plots for each time point
macro_daily_WK0 <- macro_daily_plot(macro.df.WK0)
macro_daily_WK4 <- macro_daily_plot(macro.df.WK4)
macro_daily_WK8 <- macro_daily_plot(macro.df.WK8)

# Combine macronutrient plots
macro_daily_combined <- macro_daily_WK0 + macro_daily_WK4 + macro_daily_WK8

# Save combined macronutrient plot as PNG
ggsave("./Figures/FigS1_b.png", height = 8, width = 12, units = "cm")
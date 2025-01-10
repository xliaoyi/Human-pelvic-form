library(ggplot2)

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")

df <- read.csv("key_results/error_23kps_anno_vs_pred_longer.csv")

ggplot(df, aes(x = hip_region, y = error, fill = group)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75)) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(x = "", y = "Euclidean distance (pixel)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top")


# over all
library(tidyverse)

# Pivot the data so that each group becomes a separate column
df <- read.csv("key_results/error_23kps_anno_vs_pred_wider.csv")

# Plot the data
ggplot(df, aes(x = error_anno_vs_pred, y = error_pred_vs_pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2) +
  xlim(c(0, 10)) + 
  ylim(c(0, 10)) + 
  coord_equal() +
  labs(x = "Euclidean distance Annotation vs Predicted (pixel)", 
       y = "Euclidean distance Predicted vs Predicted (pixel)") +
  theme_bw()
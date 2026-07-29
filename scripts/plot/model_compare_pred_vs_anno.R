rm(list = ls())

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(ggpubr)

# selected_landmarks = c('iliac_spine_left', 'iliac_spine_right', 'iliopubic_eminence_left', 
#                        'iliopubic_eminence_right', 'inferior_pubic_ramus_left', 
#                        'inferior_pubic_ramus_right', 'pubic_arch', 'sciatic_notch_left', 
#                        'sciatic_notch_right', 'sacrum', 'pubic_tubercle',
#                        'inferior_iliac_spine_left', 'inferior_iliac_spine_right',
#                        'acetabular_inferior_left', 'acetabular_inferior_right')

df <- read.csv("key_results/error_23kps_anno_vs_pred_longer.csv")
# df <- df %>% filter(hip_region %in% selected_landmarks)

### For all landmarks
df$group <- factor(df$group, levels = c("anno_vs_pred", "pred_vs_pred"), 
                   labels = c("Error between annotation and first prediction", "Error between first prediction and second prediction"))
p <- 
  ggbarplot(df, x = "group", y = "error", 
           fill = "group", color = "black",
           width = 0.9,
           position = position_dodge(0.9), 
           add = "mean_se") +
  # geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA) +
  # geom_boxplot(outlier.shape=4,
  #              outlier.size=1,
  #              outlier.alpha = 0.3,
  #              size = 0.3) +
  labs(x = "",
       y = "Euclidean Distance (Pixel)") +
  scale_fill_brewer(palette="Paired", name=NULL) + 
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",            # Position of legend outside the plot, at the top
    legend.justification = "center",    # Center justification for the legend
    legend.title = element_blank(),     # Remove legend title
    axis.text.x = element_blank()       # Remove x-axis text
  ) +
  guides(fill = guide_legend(order=1)) 

pdf("out_fig/model_compare_pred_vs_anno_all_together.pdf",
    width = 4, 
    height = 5)
p
dev.off()
### For single landmarks

df <- read.csv("key_results/error_23kps_anno_vs_pred_longer.csv")
# Assuming you want to change the legend content for "group" from "anno_vs_pred" to, say, "Annotation vs Prediction"
df$group <- factor(df$group, levels = c("anno_vs_pred", "pred_vs_pred"), 
                   labels = c("Error between annotation and first prediction", "Error between first prediction and second prediction"))

p <- ggplot(df, aes(x = hip_region, y = error, fill = group)) +
  geom_boxplot(outlier.shape=4,
               outlier.size=1,
               outlier.alpha = 0.3,
               size = 0.3) +
  labs(x = "",
       y = "Euclidean distance (Pixel)") +
  scale_fill_brewer(palette="Paired", name=NULL) + 
  theme_bw() +
  theme(
    legend.position = c(0.95, 0.95),              # Position of legend (centered horizontally, at the top)
    legend.justification = c(0.95, 0.95),         # Justification at the top-center of the legend
    legend.title = element_blank(),           # Remove legend title
    axis.text.x = element_text(angle=45, hjust=1), # Rotate x-axis labels by 90 degrees
    plot.margin = margin(t=10, r=10, b=10, l=40, unit="pt")
  ) +
  guides(fill = guide_legend(order=1)) 

pdf("out_fig/model_compare_pred_vs_anno.pdf",
    width = 8, 
    height = 6)
p
dev.off()
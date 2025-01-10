rm(list = ls())

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape")

library(ggplot2)
library(tidyverse)

mvsf_cor <- read.csv("key_results/genetic_correlation/rg_hip_female_vs_male_240112.csv")

# rename phenotypes
rename_map_pheno <- list(
  "acetabular_diameter" = "Acetabular diameter",
  "ear_left2ear_right" = "Head width",
  "iliac_isthmus_breadth" = "Iliac isthmus breadth",
  "oblique_pelvic_inlet_length" = "Oblique pelvic inlet length",
  "pelvic_height" = "Pelvic height",
  "pelvic_inlet_width" = "Pelvic inlet width",
  "pelvic_width" = "Pelvic width",
  "shoulder_width" = "Shoulder width",
  "subpubic_angle" = "Subpubic angle",
  "pelvic_inlet_area" = "Pelvic inlet area",
  "bi_acetabular_width" = "Biacetabular width"
)
# Apply renaming to the df dataframe
mvsf_cor <- mvsf_cor %>%
  mutate(
    pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno)
  )

mvsf_cor <- mvsf_cor %>% 
  filter(!pheno %in% c("Pelvic inlet area", "Biacetabular width"))

# Create a vector of colors, one for each phenotype
colors <- c("Pelvic height" = "#76B7B2", "Pelvic width" = "#B07AA1", "Pelvic inlet width" = "#4E79A7",
            "Oblique pelvic inlet length" = "#FF9DA7", "Iliac isthmus breadth" = "#59A14F", 
            "Acetabular diameter" = "#E15759", "Subpubic angle" = "#EDC948", #"bi_acetabular_width" = "#F28E2B", 
            "Humerus" = "#d3d4d3", "Forearm" = "#d3d4d3", "Torso length" = "#d3d4d3", "Femur" = "#d3d4d3", "Shoulder width" = "#d3d4d3", "Tibia" = "#d3d4d3")
# Sort the data frame by 'cor' values in ascending order
df <- mvsf_cor[order(mvsf_cor$cor), ]

df$pheno <- factor(df$pheno, levels = unique(df$pheno))

# Create the ggplot
p <- ggplot(df, aes(x = pheno, y = cor, color = pheno)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = cor - se, ymax = cor + se), width = 0, size = 1) +
  scale_color_manual(values = colors) +
  labs(x = "", y = "Genetic correlation") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.x = element_line(size = 0.2, linetype = "dashed", colour = "lightgrey"), 
        panel.grid.major.y = element_line(size = 0.2, linetype = "dashed", colour = "lightgrey"),
        legend.position = "none",
        panel.grid.minor = element_blank()) 

ggsave(filename = "out_fig/genetic_correlation.pdf", plot = p,
       width = 4, height = 4, units = "in")

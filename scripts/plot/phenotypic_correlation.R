rm(list = ls())

library(ggplot2)
library(dplyr)
library(pheatmap)
library(Hmisc)

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")

select_pheno <- c("pelvic_height", "pelvic_width", "pelvic_inlet_width", "subpubic_angle",
                  "oblique_pelvic_inlet_length", "ear_left2ear_right", "arm_devide_torso",
                  "shoulder_width", "iliac_isthmus_breadth", "acetabular_diameter")
# both
both_df <- read.csv("key_results/hip_select_pheno_cm_residual.csv")
both_df <- both_df %>% select(c("eid", select_pheno))
both_cor_mat <- cor(both_df[, 1:ncol(both_df)], method = "pearson")

# Make the pheatmap
pdf("out_fig/phenotypic_cor_both.pdf", width = 8, height = 7)
pheatmap(both_cor_mat, 
         display_numbers = round(both_cor_mat, 3),  # Display rounded-off correlation coefficients in the cells
         number_format = "%.3f",  # Control the number of decimal places
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         main = "Phenotypic Correlation",
         color = colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(201),
         breaks = seq(-1, 1, by = 0.01)) 
dev.off()

# male
male_df <- read.csv('key_results/hip_select_pheno_cm_male_residual.csv')
male_cor_mat <- cor(male_df[, 9:ncol(male_df)], method = "pearson")

# Make the pheatmap
pheatmap(male_cor_mat, 
         display_numbers = round(both_cor_mat, 3),  # Display rounded-off correlation coefficients in the cells
         number_format = "%.3f",  # Control the number of decimal places
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         color = colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(21),
         breaks = seq(-1, 1, by = 0.1)) 

# female
female_df <- read.csv('key_results/hip_select_pheno_cm_female_residual.csv')
female_cor_mat <- cor(female_df[, 9:ncol(female_df)], method = "pearson")

# Make the pheatmap
pheatmap(female_cor_mat, 
         display_numbers = round(both_cor_mat, 3),  # Display rounded-off correlation coefficients in the cells
         number_format = "%.3f",  # Control the number of decimal places
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         color = colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(21),
         breaks = seq(-1, 1, by = 0.1)) 

####################################################
#########        Both male and female      #########
####################################################

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")

select_pheno <- c("pelvic_height", "pelvic_width", "pelvic_inlet_width", "subpubic_angle",
                  "oblique_pelvic_inlet_length", "ear_left2ear_right", "arm_devide_torso",
                  "shoulder_width", "iliac_isthmus_breadth", "acetabular_diameter")
# both
both_df <- read.csv("key_results/hip_select_pheno_cm_residual.csv")
both_data <- both_df %>% select(select_pheno)
both_cor_mat <- rcorr(as.matrix(both_data), type = "pearson")

both_cor_mat$P[is.na(both_cor_mat$P)] <- 0

hclust_out <- hclust(dist(both_cor_mat$r), method = "complete")
both_cluster_order <- order.dendrogram(as.dendrogram(hclust_out))

both_cor_mat_ordered_r <- both_cor_mat$r[both_cluster_order, both_cluster_order]
both_cor_mat_ordered_p <- both_cor_mat$P[both_cluster_order, both_cluster_order]

cor_df <- as.data.frame(both_cor_mat_ordered_r)
pvalue_df <- as.data.frame(both_cor_mat_ordered_p)

# Add variable names
cor_df$pheno <- rownames(cor_df)
pvalue_df$pheno <- rownames(pvalue_df)

# Reshape dataframes from wide to long format
cor_df <- reshape2::melt(cor_df, id.vars = "pheno")
pvalue_df <- reshape2::melt(pvalue_df, id.vars = "pheno")

# Merge correlation and p-value dataframes
df <- merge(cor_df, pvalue_df, by = c("variable", "pheno"))

# Rename columns
colnames(df) <- c("pheno1", "pheno2", "correlation", "p_value")

# Order phenotypes
df$pheno1 <- factor(df$pheno1, levels = colnames(both_cor_mat_ordered_r))
df$pheno2 <- factor(df$pheno2, levels = colnames(both_cor_mat_ordered_r))

# Bonferroni-corrected p-value threshold
n <- dim(df)[1]
alpha <- 0.05
bonferroni_threshold <- alpha / n

# Add column indicating whether Bonferroni-corrected p-value is significant
df$significant_bonferroni <- df$p_value < bonferroni_threshold

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = pheno1, y = pheno2)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 20), fill = df$correlation), shape = 21) +
  geom_text(aes(label=round(correlation,2), size = 2.5), show.legend = FALSE) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(-1, -0.5, 0, 0.5, 1), 
                       labels=c(-1, -0.5, 0, 0.5, 1), 
                       limits = c(-1, 1), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  # geom_point(data = subset(df, significant_bonferroni), 
  #            aes(x = pheno1, y = pheno2), 
  #            color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Correlation", size = "-log10(p-value)") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  ggtitle("Both sex phenotypic correlation")

ggsave(filename = "out_fig/phenotypic_cor_both.pdf", plot = p,
       width = 8, height = 7, units = "in")

####################################################
######### Combine male and female together #########
####################################################

# Function to replace upper triangle of a matrix
replace_upper_triangle <- function(mat1, mat2) {
  mat1[upper.tri(mat1)] <- mat2[upper.tri(mat2)]
  return(mat1)
}

# # regress separately
# # male
# male_df <- read.csv('key_results/hip_select_pheno_cm_male_wek_residual.csv')
# male_data <- male_df[, 9:ncol(male_df)]
# male_cor_mat <- rcorr(as.matrix(male_data), type = "pearson")
# 
# # female
# female_df <- read.csv('key_results/hip_select_pheno_cm_female_wek_residual.csv')
# female_data <- female_df[, 9:ncol(female_df)]
# female_cor_mat <- rcorr(as.matrix(female_data), type = "pearson")

# regress together
male_female_df <- read.csv("key_results/hip_select_pheno_cm_residual_regress_on_female_and_male.csv")

# filter
select_pheno <- c("pelvic_height", "pelvic_width", "pelvic_inlet_width", "subpubic_angle",
                  "oblique_pelvic_inlet_length", "ear_left2ear_right", "arm_devide_torso",
                  "shoulder_width", "iliac_isthmus_breadth", "acetabular_diameter")
male_female_df <- male_female_df %>% select(c("eid", "sex", select_pheno))

rename_vector <- c(
  "acetabular_diameter" = "Acetabular diameter",
  "arm_devide_torso" = "Arm:torso",
  "ear_left2ear_right" = "Head width",
  "iliac_isthmus_breadth" = "Iliac isthmus breadth",
  "oblique_pelvic_inlet_length" = "Oblique pelvic inlet length",
  "pelvic_height" = "Pelvic height",
  "pelvic_inlet_width" = "Pelvic inlet width",
  "pelvic_width" = "Pelvic width",
  "shoulder_width" = "Shoulder width",
  "subpubic_angle" = "Subpubic angle"
)
current_colnames <- colnames(male_female_df)
new_colnames <- ifelse(current_colnames %in% names(rename_vector), rename_vector[current_colnames], current_colnames)
colnames(male_female_df) <- new_colnames

# male
male_df <- male_female_df %>% filter(sex == 1)
male_data <- male_df[, 3:ncol(male_df)]
male_cor_mat <- rcorr(as.matrix(male_data), type = "pearson")

# female
female_df <- male_female_df %>% filter(sex == 0)
female_data <- female_df[, 3:ncol(female_df)]
female_cor_mat <- rcorr(as.matrix(female_data), type = "pearson")

# Get the column order from the female correlation matrix clustering
hclust_out <- hclust(dist(female_cor_mat$r), method = "complete")
female_cluster_order <- order.dendrogram(as.dendrogram(hclust_out))

# Reorder the male and female correlation matrices to match the female order
male_cor_mat_ordered_r <- male_cor_mat$r[female_cluster_order, female_cluster_order]
male_cor_mat_ordered_p <- male_cor_mat$P[female_cluster_order, female_cluster_order]

female_cor_mat_ordered_r <- female_cor_mat$r[female_cluster_order, female_cluster_order]
female_cor_mat_ordered_p <- female_cor_mat$P[female_cluster_order, female_cluster_order]

# Create a combined correlation matrix
combined_cor_mat_r <- replace_upper_triangle(male_cor_mat_ordered_r, female_cor_mat_ordered_r)
combined_cor_mat_p <- replace_upper_triangle(male_cor_mat_ordered_p, female_cor_mat_ordered_p)

# Convert matrices to dataframes
cor_df <- as.data.frame(combined_cor_mat_r)
pvalue_df <- as.data.frame(combined_cor_mat_p)

# Add variable names
cor_df$pheno <- rownames(cor_df)
pvalue_df$pheno <- rownames(pvalue_df)

# Reshape dataframes from wide to long format
cor_df <- reshape2::melt(cor_df, id.vars = "pheno")
pvalue_df <- reshape2::melt(pvalue_df, id.vars = "pheno")

# Merge correlation and p-value dataframes
df <- merge(cor_df, pvalue_df, by = c("variable", "pheno"))

# remove NA
df[is.na(df)] <- 0

# Rename columns
colnames(df) <- c("pheno1", "pheno2", "correlation", "p_value")

# Order phenotypes
df$pheno1 <- factor(df$pheno1, levels = colnames(combined_cor_mat_r))
df$pheno2 <- factor(df$pheno2, levels = colnames(combined_cor_mat_r))

# Bonferroni-corrected p-value threshold
n <- dim(df)[1] - dim(df)[1]**0.5
alpha <- 0.05
bonferroni_threshold <- alpha / n

# Add column indicating whether Bonferroni-corrected p-value is significant
df$significant_bonferroni <- df$p_value < bonferroni_threshold

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = pheno1, y = pheno2)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 25), fill = df$correlation), shape = 22) +
  # geom_text(aes(label=round(correlation,2), size = 2.5), show.legend = FALSE) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(-1, -0.5, 0, 0.5, 1), 
                       labels=c(-1, -0.5, 0, 0.5, 1), 
                       limits = c(-1, 1), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 10)) +
  # geom_point(data = subset(df, significant_bonferroni), 
  #            aes(x = pheno1, y = pheno2), 
  #            color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Female Phenotypic Correlation", y = "Male Phenotypic Correlation", fill = "Correlation", size = "-log10(p-value)") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  ggtitle("Phenotypic correlation")

ggsave(filename = "out_fig/phenotypic_cor_male_female_reg_separate.pdf", plot = p,
       width = 7, height = 6, units = "in")


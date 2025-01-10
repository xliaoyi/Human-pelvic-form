rm(list = ls())

library(tidyverse)
library(patchwork)

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")


select_pheno <- c("pelvic_height", "pelvic_width", "pelvic_inlet_width", "subpubic_angle",
                  "oblique_pelvic_inlet_length", "ear_left2ear_right", "arm_devide_torso",
                  "shoulder_width", "iliac_isthmus_breadth", "acetabular_diameter")

# Function to replace upper triangle of a matrix
replace_upper_triangle <- function(mat1, mat2) {
  mat1[upper.tri(mat1)] <- mat2[upper.tri(mat2)]
  return(mat1)
}

#######################################################
#               Phenotypic correlation
#######################################################

# regress together
male_female_df <- read.csv("key_results/hip_select_pheno_cm_residual_regress_on_female_and_male.csv")

# filter
male_female_df <- male_female_df %>% select(c("eid", "sex", select_pheno))

# rename phenotypes
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

# order
order <- colnames(combined_cor_mat_r)

# Convert matrices to dataframes
cor_df <- as.data.frame(combined_cor_mat_r)
pvalue_df <- as.data.frame(combined_cor_mat_p)

# Add variable names
cor_df$pheno1 <- rownames(cor_df)
pvalue_df$pheno1 <- rownames(pvalue_df)

# Reshape dataframes from wide to long format
cor_df <- reshape2::melt(cor_df, id.vars = "pheno1", variable.name = "pheno2", value.name = "correlation")
pvalue_df <- reshape2::melt(pvalue_df, id.vars = "pheno1", variable.name = "pheno2", value.name = "p_value")

# Merge correlation and p-value dataframes
df <- merge(cor_df, pvalue_df, by = c("pheno1", "pheno2"))

# remove NA
df[is.na(df)] <- 0

# Order phenotypes
df$pheno1 <- factor(df$pheno1, levels = colnames(combined_cor_mat_r))
df$pheno2 <- factor(df$pheno2, levels = colnames(combined_cor_mat_r))

# Bonferroni-corrected p-value threshold
n <- dim(df)[1] - dim(df)[1]**0.5
alpha <- 0.05
bonferroni_threshold <- alpha / n

# Add column indicating whether Bonferroni-corrected p-value is significant
df$significant_bonferroni <- df$p_value < bonferroni_threshold

pheno_df <- df

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

pheno <- ggplot(pheno_df, aes(x = pheno1, y = pheno2)) +
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
  labs(x = "Male phenotypic correlation", y = "Female phenotypic correlation", fill = "Correlation", size = "-log10(p-value)") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  ggtitle("Phenotypic correlation")
ggsave(filename = "out_fig/phenotypic_cor_male_female.pdf", plot = pheno,
       width = 7, height = 6, units = "in")

#######################################################
#               Genetic correlation
#######################################################

# read data
male_df <- read.csv("key_results/genetic_correlation/rg_20230909/rg_hip_male.csv")
male_df <- male_df %>% filter(pheno1 %in% select_pheno) %>% filter(pheno2 %in% select_pheno)

female_df <- read.csv("key_results/genetic_correlation/rg_20230909/rg_hip_female.csv")
female_df <- female_df %>% filter(pheno1 %in% select_pheno) %>% filter(pheno2 %in% select_pheno)

# rename phenotypes
rename_map <- list(
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

# Apply renaming to the df dataframe
male_df <- male_df %>%
  mutate(
    pheno1 = ifelse(pheno1 %in% names(rename_map), rename_map[pheno1], pheno1),
    pheno2 = ifelse(pheno2 %in% names(rename_map), rename_map[pheno2], pheno2)
  )
female_df <- female_df %>%
  mutate(
    pheno1 = ifelse(pheno1 %in% names(rename_map), rename_map[pheno1], pheno1),
    pheno2 = ifelse(pheno2 %in% names(rename_map), rename_map[pheno2], pheno2)
  )

#######################################################
# get male cor matrix
male_cor_mat_r <- male_df %>% select(c("pheno1", "pheno2", "cor"))
phenos <- unique(c(male_df$pheno1, male_df$pheno2))
for (i in phenos) {
  new_row <- data.frame(
    pheno1 = i,
    pheno2 = i,
    cor = 1
  )
  male_cor_mat_r <- rbind(new_row, male_cor_mat_r)
}

male_cor_mat_r <- male_cor_mat_r %>%
  spread(key = pheno2, value = cor)

rownames(male_cor_mat_r) <- male_cor_mat_r$pheno1
male_cor_mat_r$pheno1 <- NULL

male_cor_mat_r[upper.tri(male_cor_mat_r)] <- t(male_cor_mat_r)[upper.tri(t(male_cor_mat_r))]
male_cor_mat_ordered_r <- male_cor_mat_r[order, order]

#######################################################
# get male p matrix
male_cor_mat_p <- male_df %>% select(c("pheno1", "pheno2", "p_value"))
phenos <- unique(c(male_df$pheno1, male_df$pheno2))
for (i in phenos) {
  new_row <- data.frame(
    pheno1 = i,
    pheno2 = i,
    p_value = 0
  )
  male_cor_mat_p <- rbind(new_row, male_cor_mat_p)
}

male_cor_mat_p <- male_cor_mat_p %>%
  spread(key = pheno2, value = p_value)

rownames(male_cor_mat_p) <- male_cor_mat_p$pheno1
male_cor_mat_p$pheno1 <- NULL

male_cor_mat_p[upper.tri(male_cor_mat_p)] <- t(male_cor_mat_p)[upper.tri(t(male_cor_mat_p))]
male_cor_mat_ordered_p <- male_cor_mat_p[order, order]

#######################################################
# get female cor matrix
female_cor_mat_r <- female_df %>% select(c("pheno1", "pheno2", "cor"))
phenos <- unique(c(female_df$pheno1, female_df$pheno2))
for (i in phenos) {
  new_row <- data.frame(
    pheno1 = i,
    pheno2 = i,
    cor = 1
  )
  female_cor_mat_r <- rbind(new_row, female_cor_mat_r)
}

female_cor_mat_r <- female_cor_mat_r %>%
  spread(key = pheno2, value = cor)

rownames(female_cor_mat_r) <- female_cor_mat_r$pheno1
female_cor_mat_r$pheno1 <- NULL

female_cor_mat_r[upper.tri(female_cor_mat_r)] <- t(female_cor_mat_r)[upper.tri(t(female_cor_mat_r))]
female_cor_mat_ordered_r <- female_cor_mat_r[order, order]

#######################################################
# get female p matrix
female_cor_mat_p <- female_df %>% select(c("pheno1", "pheno2", "p_value"))
phenos <- unique(c(female_df$pheno1, female_df$pheno2))
for (i in phenos) {
  new_row <- data.frame(
    pheno1 = i,
    pheno2 = i,
    p_value = 0
  )
  female_cor_mat_p <- rbind(new_row, female_cor_mat_p)
}

female_cor_mat_p <- female_cor_mat_p %>%
  spread(key = pheno2, value = p_value)

rownames(female_cor_mat_p) <- female_cor_mat_p$pheno1
female_cor_mat_p$pheno1 <- NULL

female_cor_mat_p[upper.tri(female_cor_mat_p)] <- t(female_cor_mat_p)[upper.tri(t(female_cor_mat_p))]
female_cor_mat_ordered_p <- female_cor_mat_p[order, order]



# Create a combined correlation matrix
combined_cor_mat_r <- replace_upper_triangle(male_cor_mat_ordered_r, female_cor_mat_ordered_r)
combined_cor_mat_p <- replace_upper_triangle(male_cor_mat_ordered_p, female_cor_mat_ordered_p)

# Convert matrices to dataframes
cor_df <- as.data.frame(combined_cor_mat_r)
pvalue_df <- as.data.frame(combined_cor_mat_p)

# Add variable names
cor_df$pheno1 <- rownames(cor_df)
pvalue_df$pheno1 <- rownames(pvalue_df)

# Reshape dataframes from wide to long format
cor_df <- reshape2::melt(cor_df, id.vars = "pheno1", variable.name = "pheno2", value.name = "correlation")
pvalue_df <- reshape2::melt(pvalue_df, id.vars = "pheno1", variable.name = "pheno2", value.name = "p_value")

# Merge correlation and p-value dataframes
df <- merge(cor_df, pvalue_df, by = c("pheno1", "pheno2"))

# remove NA
df[is.na(df)] <- 0

# Order phenotypes
df$pheno1 <- factor(df$pheno1, levels = colnames(combined_cor_mat_r))
df$pheno2 <- factor(df$pheno2, levels = colnames(combined_cor_mat_r))


# Bonferroni-corrected p-value threshold
n <- dim(df)[1] - dim(df)[1]**0.5
alpha <- 0.05
bonferroni_threshold <- alpha / n

# Add column indicating whether Bonferroni-corrected p-value is significant
df$significant_bonferroni <- df$p_value < bonferroni_threshold

genome_df <- df

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)
genome <- ggplot(genome_df, aes(x = pheno1, y = pheno2)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 6), fill = df$correlation), shape = 22) +
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
  labs(x = "Male genetic correlation", y = "Female genetic correlation", fill = "Correlation", size = "-log10(p-value)") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  ggtitle("Genetic correlation")
ggsave(filename = "out_fig/genetic_cor_male_female.pdf", plot = genome,
       width = 7, height = 6, units = "in")



#######################################################
#      Both Phenotypic and Genetic correlation
#######################################################

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/key_results/genetic_correlation/")

df <- read.csv("genetic_and_phenotype_cor.csv")
df <- df %>% filter(pheno1 %in% select_pheno) %>% filter(pheno2 %in% select_pheno)

n <- nrow(df)

# Bonferroni-corrected p-value threshold
alpha <- 0.05
bonferroni_threshold <- alpha / n

# Add column indicating whether Bonferroni-corrected p-value is significant
df$significant_bonferroni <- df$p_value < bonferroni_threshold

color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = pheno1, y = pheno2)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 10), fill = cor), shape = 22) +
  scale_fill_gradientn(colors = color_palette, limits = c(-1, 1)) +
  scale_size_continuous(range = c(0, 10)) +
  # geom_point(data = subset(df, significant_bonferroni), color = "black", shape = 8, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Genetic correlation", y = "Phenotypic correlation", fill = "Correlation", size = "-log10(p-value)")

ggsave(filename = "../../out_fig/genetic_phenotypic_cor.pdf", plot = p, width = 7, height = 6, units = "in")

print(p)
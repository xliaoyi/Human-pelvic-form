rm(list = ls())

library(ggplot2)
library(dplyr)
library(pheatmap)
library(Hmisc)

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")


####################################################
#########       Both male and female       #########
####################################################

# both
both_df <- read.csv("key_results/genetic_correlation/rg_20230909/rg_hip_both.csv")

###################################################
# get cor dataframe

both_df_cor <- both_df %>% select(c("pheno1", "pheno2", "cor"))
both_df_cor_swap <- both_df_cor %>% rename(pheno1 = pheno2, pheno2 = pheno1)
both_df_cor <- rbind(both_df_cor, both_df_cor_swap)

phenos <- unique(c(both_df_cor$pheno1, both_df_cor$pheno2))
for (i in phenos) {
  new_row <- data.frame(
    pheno1 = i,
    pheno2 = i,
    cor = 1
  )
  both_df_cor <- rbind(both_df_cor, new_row)
}

both_df_cor <- both_df_cor %>%
  spread(key = pheno2, value = cor)

# Set the 'pheno1' column as row names
row.names(both_df_cor) <- both_df_cor$pheno1
# Remove the 'pheno1' column
both_df_cor$pheno1 <- NULL

###################################################
# get p value dataframe

both_df_p <- both_df %>% select(c("pheno1", "pheno2", "p_value"))
both_df_p_swap <- both_df_p %>% rename(pheno1 = pheno2, pheno2 = pheno1)
both_df_p <- rbind(both_df_p, both_df_p_swap)

phenos <- unique(c(both_df_p$pheno1, both_df_p$pheno2))
for (i in phenos) {
  new_row <- data.frame(
    pheno1 = i,
    pheno2 = i,
    p_value = 0
  )
  both_df_p <- rbind(both_df_p, new_row)
}

both_df_p <- both_df_p %>%
  spread(key = pheno2, value = p_value)

# Set the 'pheno1' column as row names
row.names(both_df_p) <- both_df_p$pheno1
# Remove the 'pheno1' column
both_df_p$pheno1 <- NULL

# cluster
hclust_out <- hclust(dist(both_df_cor), method = "complete")
cluster_order <- order.dendrogram(as.dendrogram(hclust_out))

# Reorder the male and female correlation matrices to match the female order
cor_mat_ordered_r <- both_df_cor[cluster_order, cluster_order]
cor_mat_ordered_p <- both_df_p[cluster_order, cluster_order]

# Convert matrices to dataframes
cor_df <- as.data.frame(cor_mat_ordered_r)
pvalue_df <- as.data.frame(cor_mat_ordered_p)

# Add variable names
cor_df$pheno <- rownames(cor_df)
pvalue_df$pheno <- rownames(pvalue_df)

# Reshape dataframes from wide to long format
cor_df <- reshape2::melt(cor_df, id.vars = "pheno") %>% rename(cor = value)
pvalue_df <- reshape2::melt(pvalue_df, id.vars = "pheno") %>% rename(p = value)

# Merge correlation and p-value dataframes
df <- merge(cor_df, pvalue_df, by = c("variable", "pheno"))

# remove NA
df[is.na(df)] <- 0

# Rename columns
colnames(df) <- c("pheno1", "pheno2", "correlation", "p_value")

# Order phenotypes
df$pheno1 <- factor(df$pheno1, levels = colnames(cor_mat_ordered_r))
df$pheno2 <- factor(df$pheno2, levels = colnames(cor_mat_ordered_r))

# Bonferroni-corrected p-value threshold
n <- (dim(df)[1] / 2) - 6
alpha <- 0.05
bonferroni_threshold <- alpha / n

# Add column indicating whether Bonferroni-corrected p-value is significant
df$significant_bonferroni <- df$p_value < bonferroni_threshold

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = pheno1, y = pheno2)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 20), fill = df$correlation), shape = 21) +
  # geom_text(aes(label=round(correlation,2), size = 2.5), show.legend = FALSE) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(-1, -0.5, 0, 0.5, 1), 
                       labels=c(-1, -0.5, 0, 0.5, 1), 
                       limits = c(-1, 1), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, significant_bonferroni),
             aes(x = pheno1, y = pheno2),
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Correlation", size = "-log10(p-value)") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  ggtitle("Genetic correlation of both male and female")

ggsave(filename = "out_fig/genetic_correlation_both.pdf", plot = p,
       width = 8, height = 7, units = "in")

####################################################
######### Combine male and female together #########
####################################################

select_pheno <- c("pelvic_height", "pelvic_width", "pelvic_inlet_width", "subpubic_angle",
                  "oblique_pelvic_inlet_length", "ear_left2ear_right", "arm_devide_torso",
                  "shoulder_width", "iliac_isthmus_breadth", "acetabular_diameter")

# read data
male_df <- read.csv("key_results/genetic_correlation/rg_20230909/rg_hip_male.csv")
male_df <- male_df %>% filter(pheno1 %in% select_pheno) %>% filter(pheno2 %in% select_pheno)

female_df <- read.csv("key_results/genetic_correlation/rg_20230909/rg_hip_female.csv")
female_df <- female_df %>% filter(pheno1 %in% select_pheno) %>% filter(pheno2 %in% select_pheno)

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

###################################################
# get cor dataframe

male_df_cor <- male_df %>% select(c("pheno1", "pheno2", "cor"))
female_df_cor <- female_df %>% select(c("pheno1", "pheno2", "cor"))
female_df_cor <- female_df_cor %>% rename(pheno1 = pheno2, pheno2 = pheno1)
mf_df_cor <- rbind(male_df_cor, female_df_cor)

phenos <- unique(mf_df_cor$pheno1)
for (i in phenos) {
  new_row <- data.frame(
    pheno1 = i,
    pheno2 = i,
    cor = 1
  )
  mf_df_cor <- rbind(mf_df_cor, new_row)
}

mf_df_cor <- mf_df_cor %>%
  spread(key = pheno2, value = cor)

# Set the 'pheno1' column as row names
row.names(mf_df_cor) <- mf_df_cor$pheno1
# Remove the 'pheno1' column
mf_df_cor$pheno1 <- NULL

###################################################
# get p value dataframe

male_df_p <- male_df %>% select(c("pheno1", "pheno2", "p_value"))
female_df_p <- female_df %>% select(c("pheno1", "pheno2", "p_value"))
female_df_p <- female_df_p %>% rename(pheno1 = pheno2, pheno2 = pheno1)
mf_df_p <- rbind(male_df_p, female_df_p)

phenos <- unique(mf_df_p$pheno1)
for (i in phenos) {
  new_row <- data.frame(
    pheno1 = i,
    pheno2 = i,
    p_value = 0
  )
  mf_df_p <- rbind(mf_df_p, new_row)
}

mf_df_p <- mf_df_p %>%
  spread(key = pheno2, value = p_value)

# Set the 'pheno1' column as row names
row.names(mf_df_p) <- mf_df_p$pheno1
# Remove the 'pheno1' column
mf_df_p$pheno1 <- NULL

# cluster
hclust_out <- hclust(dist(mf_df_cor), method = "complete")
cluster_order <- order.dendrogram(as.dendrogram(hclust_out))

# Reorder the male and female correlation matrices to match the female order
cor_mat_ordered_r <- mf_df_cor[cluster_order, cluster_order]
cor_mat_ordered_p <- mf_df_p[cluster_order, cluster_order]

# Convert matrices to dataframes
cor_df <- as.data.frame(cor_mat_ordered_r)
pvalue_df <- as.data.frame(cor_mat_ordered_p)

# Add variable names
cor_df$pheno <- rownames(cor_df)
pvalue_df$pheno <- rownames(pvalue_df)

# Reshape dataframes from wide to long format
cor_df <- reshape2::melt(cor_df, id.vars = "pheno") %>% rename(cor = value)
pvalue_df <- reshape2::melt(pvalue_df, id.vars = "pheno") %>% rename(p = value)

# Merge correlation and p-value dataframes
df <- merge(cor_df, pvalue_df, by = c("variable", "pheno"))

# remove NA
df[is.na(df)] <- 0

# Rename columns
colnames(df) <- c("pheno1", "pheno2", "correlation", "p_value")

# Order phenotypes
df$pheno1 <- factor(df$pheno1, levels = colnames(cor_mat_ordered_r))
df$pheno2 <- factor(df$pheno2, levels = colnames(cor_mat_ordered_r))

# Bonferroni-corrected p-value threshold
n <- dim(df)[1] - 12
alpha <- 0.05
bonferroni_threshold <- alpha / n

# Add column indicating whether Bonferroni-corrected p-value is significant
df$significant_bonferroni <- df$p_value < bonferroni_threshold
  

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = pheno1, y = pheno2)) +
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
  ggtitle("Genetic correlation of male and female")

ggsave(filename = "out_fig/genetic_correlation_male_and_female.pdf", plot = p,
       width = 7, height = 6, units = "in")


####################################################
#########           male vs female         #########
####################################################

# read data
fvsm_df <- read.csv("key_results/genetic_correlation/rg_20230909/rg_hip_female_vs_male.csv")

###################################################
# get cor dataframe

fvsm_df_cor <- fvsm_df %>% select(c("pheno_f", "pheno_m", "cor"))

fvsm_df_cor <- fvsm_df_cor %>%
  spread(key = pheno_m, value = cor)

# Set the 'pheno_f' column as row names
row.names(fvsm_df_cor) <- fvsm_df_cor$pheno_f
# Remove the 'pheno_f' column
fvsm_df_cor$pheno_f <- NULL

###################################################
# get p value dataframe

fvsm_df_p <- fvsm_df %>% select(c("pheno_f", "pheno_m", "p_value"))

fvsm_df_p <- fvsm_df_p %>%
  spread(key = pheno_m, value = p_value)

# Set the 'pheno_f' column as row names
row.names(fvsm_df_p) <- fvsm_df_p$pheno_f
# Remove the 'pheno_f' column
fvsm_df_p$pheno_f <- NULL

# cluster
hclust_out <- hclust(dist(fvsm_df_cor), method = "complete")
cluster_order <- order.dendrogram(as.dendrogram(hclust_out))

# Reorder the male and female correlation matrices to match the female order
cor_mat_ordered_r <- fvsm_df_cor[cluster_order, cluster_order]
cor_mat_ordered_p <- fvsm_df_p[cluster_order, cluster_order]

# Convert matrices to dataframes
cor_df <- as.data.frame(cor_mat_ordered_r)
pvalue_df <- as.data.frame(cor_mat_ordered_p)

# Add variable names
cor_df$pheno <- rownames(cor_df)
pvalue_df$pheno <- rownames(pvalue_df)

# Reshape dataframes from wide to long format
cor_df <- reshape2::melt(cor_df, id.vars = "pheno") %>% rename(cor = value)
pvalue_df <- reshape2::melt(pvalue_df, id.vars = "pheno") %>% rename(p = value)

# Merge correlation and p-value dataframes
df <- merge(cor_df, pvalue_df, by = c("variable", "pheno"))

# remove NA
df[is.na(df)] <- 0

# Rename columns
colnames(df) <- c("pheno1", "pheno2", "correlation", "p_value")

# Order phenotypes
df$pheno1 <- factor(df$pheno1, levels = colnames(cor_mat_ordered_r))
df$pheno2 <- factor(df$pheno2, levels = colnames(cor_mat_ordered_r))

# Bonferroni-corrected p-value threshold
n <- dim(df)[1]
alpha <- 0.05
bonferroni_threshold <- alpha / n

# Add column indicating whether Bonferroni-corrected p-value is significant
df$significant_bonferroni <- df$p_value < bonferroni_threshold

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = pheno1, y = pheno2)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 6), fill = df$correlation), shape = 21) +
  # geom_text(aes(label=round(correlation,2), size = 2.5), show.legend = FALSE) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(-1, -0.5, 0, 0.5, 1), 
                       labels=c(-1, -0.5, 0, 0.5, 1), 
                       limits = c(-1, 1), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, significant_bonferroni),
             aes(x = pheno1, y = pheno2),
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Male", y = "Female", fill = "Correlation", size = "-log10(p-value)") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  ggtitle("Genetic correlation of male vs female")

ggsave(filename = "out_fig/genetic_correlation_male_vs_female.pdf", plot = p,
       width = 8, height = 7, units = "in")

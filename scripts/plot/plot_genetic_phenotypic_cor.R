library(ggplot2)
library(dplyr)
library(patchwork)

select_pheno <- c("pelvic_height", "pelvic_width", "pelvic_inlet_width", "subpubic_angle",
                  "oblique_pelvic_inlet_length", "ear_left2ear_right", "arm_devide_torso",
                  "shoulder_width", "iliac_isthmus_breadth", "acetabular_diameter")

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


#----------

# male vs female genetic correlation
setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/key_results")

df <- read.csv("genetic_correlation/male_vs_female_genetic_cor.csv")

n <- nrow(df)

# Bonferroni-corrected p-value threshold
alpha <- 0.05
bonferroni_threshold <- alpha / n

# Add column indicating whether Bonferroni-corrected p-value is significant
df$significant_bonferroni <- df$p_value < bonferroni_threshold

color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = pheno1, y = pheno2)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 8), fill = cor), shape = 21) +
  scale_fill_gradientn(colors = color_palette, limits = c(-1, 1)) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, significant_bonferroni), color = "black", shape = 8, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Female genetic correlation", y = "Male genetic correlation", fill = "Correlation", size = "-log10(p-value)")

ggsave(filename = "../out_fig/male_vs_female_genetic_cor.pdf", plot = p, width = 8, height = 7, units = "in")

print(p)

########################################################################################
# male or female self not correlated, but when combine them together they are correlated
########################################################################################

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")
# regress together
male_female_df <- read.csv("key_results/hip_select_pheno_cm_residual_regress_on_female_and_male.csv")
p1 <- ggplot(male_female_df %>% filter(sex == 0), aes(x = pelvic_inlet_width, y = ear_left2ear_right)) +
  geom_point(color = "#E15759") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_minimal()
p2<- ggplot(male_female_df %>% filter(sex == 1), aes(x = pelvic_inlet_width, y = ear_left2ear_right)) +
  geom_point(color = "#1F78B4") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_minimal()
p3 <- ggplot(male_female_df, aes(x = pelvic_inlet_width, y = ear_left2ear_right)) +
  geom_point(aes(color = as.factor(sex))) +
  geom_smooth(method = "lm", se = FALSE, color = "black", aes(group = 1)) +  # Regression line for all people
  scale_color_manual(values = c("1" = "#1F78B4", "0" = "#E15759"), 
                     labels = c("1" = "Male", "0" = "Female")) +
  labs(color = "Sex") +
  theme_minimal()
pdf("out_fig/added_effect.pdf",
    width = 15, 
    height = 5)
p1 + p2 + p3
dev.off()


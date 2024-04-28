rm(list= ls())
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsci)
library(gridExtra)
library(grid)
library(dplyr)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(stringr)

##########################################
#             Phenotypic
##########################################

domi_fid <- read.csv('/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/UKB_xray_image_info/fids/initial_patient_flt/fid1707.csv')
domi_fid <- domi_fid %>% 
  rename("handedness" = "X1707.0.0") %>% 
  select(c("eid", "handedness"))

# Replace the numbers with the corresponding text
domi_fid$handedness[domi_fid$handedness == 1] <- "right"
domi_fid$handedness[domi_fid$handedness == 2] <- "left"
domi_fid$handedness[domi_fid$handedness == 3] <- "both"

# Convert the column to a factor
domi_fid$handedness <- as.factor(domi_fid$handedness)

# Remove rows where handedness is not "right", "left", or "both"
domi_fid <- domi_fid[domi_fid$handedness %in% c("right", "left"), ] # , "both"

pheno <- read.csv('/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/key_results/hip_pheno_23_cm.csv')

pheno <- pheno %>% 
  mutate("lr_ratio_iliac_isthmus_breadth" = pheno$sciatic_notch_left2inferior_iliac_spine_right / pheno$sciatic_notch_right2inferior_iliac_spine_left) %>% 
  mutate("lr_ratio_acetabular_diameter" = pheno$iliopubic_eminence_left2acetabular_inferior_right / pheno$iliopubic_eminence_right2acetabular_inferior_left) %>% 
  # mutate("lr_ratio_acetabular_inclination" = pheno$acetabular_inclination_left / pheno$acetabular_inclination_right) %>% 
  select(c("eid", "sex", "age_imaging_visit", "standing_height_imaging_visit", "lr_ratio_iliac_isthmus_breadth", "lr_ratio_acetabular_diameter")) # , "lr_ratio_acetabular_inclination"

df <- merge(pheno, domi_fid, by = 'eid')

df$z_iliac <- scale(df$lr_ratio_iliac_isthmus_breadth)
df$z_acetabular_diameter <- scale(df$lr_ratio_acetabular_diameter)
# df$z_acetabular_inclination <- scale(df$lr_ratio_acetabular_inclination)

# Remove rows where z-score is > 4 or < -4
df <- df[!(df$z_iliac > 4 | df$z_iliac < -4), ]
df <- df[!(df$z_acetabular_diameter > 4 | df$z_acetabular_diameter < -4), ]
# df <- df[!(df$z_acetabular_inclination > 4 | df$z_acetabular_inclination < -4), ]

# Remove z-score columns
df$z_iliac <- NULL
df$z_acetabular_diameter <- NULL
# df$z_acetabular_inclination <- NULL


### PLOT ###
custom_colors <- c(left = "#fc8d62", right = "#8da0cb") # "both" = 'gray', 

my_comparisons <- list( c("right", "left")) # , c("right", "both"), c("left", "both") 

# Create boxplot for lr_ratio_iliac_isthmus_breadth
p1 <- ggviolin(df, x = "handedness", y = "lr_ratio_iliac_isthmus_breadth", 
                fill = "handedness", width = 0.7) +
  geom_boxplot(width = 0.25, fill = "white", alpha = 0.9, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "t.test", 
                     ref.group = NULL, hide.ns = TRUE) +
  geom_hline(yintercept = 1, linetype="dashed", color = "black", size = 0.5) +
  scale_fill_manual(values = custom_colors) +
  ylab("Left right ratio of iliac isthmus breadth") +
  xlab("") 

# Create boxplot for lr_ratio_acetabular_diameter
p2 <- ggviolin(df, x = "handedness", y = "lr_ratio_acetabular_diameter", 
               fill = "handedness", width = 0.7) +
  geom_boxplot(width = 0.25, fill = "white", alpha = 0.9, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "t.test", 
                     ref.group = NULL, hide.ns = TRUE) +
  geom_hline(yintercept = 1, linetype="dashed", color = "black", size = 0.5) +
  scale_fill_manual(values = custom_colors) +
  ylab("Left right ratio of acetabular diameter") +
  xlab("Handedness")+
  theme(legend.position = "none")

# Create boxplot for lr_ratio_acetabular_inclination
# p3 <- ggviolin(df, x = "handedness", y = "lr_ratio_acetabular_inclination", 
#                 fill = "handedness", width = 0.7) +
#   geom_boxplot(width = 0.25, fill = "white", alpha = 0.9, outlier.shape = NA) +
#   stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "t.test", 
#                      ref.group = NULL, hide.ns = TRUE) +
#   geom_hline(yintercept = 1, linetype="dashed", color = "black", size = 0.5) +
#   scale_fill_manual(values = custom_colors) +
#   ylab("Left right ratio of acetabular inclination") +
#   xlab("") +
#   theme(legend.position = "none")
# Display the plots
pdf("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/out_fig/relationship_handedness_lr.pdf", width=7, height=5)
p1 + p2 # + p3
dev.off()

##########################################
#         Genetic correlation
##########################################

pelvic_df <- read.csv("key_results/genetic_correlation/rg_20230830/rg_hip_handedness.csv")
handedness <- c("equally_handed_rg", "left_handed_rg", "right_handed_rg")
df <- pelvic_df %>% filter(pheno2 %in% handedness)

df$pheno2[str_detect(df$pheno2, "equally")] <- "Both"
df$pheno2[str_detect(df$pheno2, "left")] <- "Left"
df$pheno2[str_detect(df$pheno2, "right")] <- "Right"

set2_palette <- brewer.pal(n = 8, name = "Set2")

p <- ggplot(df, aes(x = cor, y = pheno1, col = pheno2)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(y = pheno1, xmin = cor - se, xmax = cor + se), 
                 height = 0, position = position_dodge(width = 0.5), size = 1) +
  theme_bw() +
  theme(legend.position = c(1, 1), legend.justification = c(1.1, 1.1), 
        legend.title = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  +
  scale_color_manual(
    values = c("gray", set2_palette[2:3]), 
    labels = c("Both", "Left", "Right")
  ) + 
  labs(x = "Gentic correlation", y = "")

ggsave("out_fig/rg_btw_hip_lr_ratio_and_handedness.pdf", p, height = 6, width = 5)



##########################################
#       Association by regression
##########################################

domi_fid <- read.csv('/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/UKB_xray_image_info/fids/initial_patient_flt/fid1707.csv')
domi_fid <- domi_fid %>% 
  rename("handedness" = "X1707.0.0") %>% 
  select(c("eid", "handedness"))

# Replace the numbers with the corresponding text
# domi_fid$handedness[domi_fid$handedness == 1] <- "right"
domi_fid$handedness[domi_fid$handedness == 2] <- 0
# domi_fid$handedness[domi_fid$handedness == 3] <- "both"

# Convert the column to a factor
# domi_fid$handedness <- as.factor(domi_fid$handedness)

# Remove rows where handedness is not "right", "left", or "both"
# domi_fid <- domi_fid[domi_fid$handedness %in% c("right", "left", "both"), ]

pheno <- read.csv('/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/key_results/hip_pheno_23_cm.csv')

pheno <- pheno %>% 
  mutate("lr_ratio_iliac_isthmus_breadth" = pheno$sciatic_notch_left2inferior_iliac_spine_left / pheno$sciatic_notch_right2inferior_iliac_spine_right) %>% 
  mutate("lr_ratio_acetabular_diameter" = pheno$iliopubic_eminence_left2acetabular_inferior_left / pheno$iliopubic_eminence_right2acetabular_inferior_right) %>% 
  mutate("lr_ratio_acetabular_inclination" = pheno$acetabular_inclination_left / pheno$acetabular_inclination_right) %>% 
  select(c("eid", "sex", "age_imaging_visit", "standing_height_imaging_visit", "lr_ratio_iliac_isthmus_breadth", "lr_ratio_acetabular_diameter", "lr_ratio_acetabular_inclination"))

df <- merge(pheno, domi_fid, by = 'eid')

df$z_iliac <- scale(df$lr_ratio_iliac_isthmus_breadth)
df$z_acetabular_diameter <- scale(df$lr_ratio_acetabular_diameter)
# df$z_acetabular_inclination <- scale(df$lr_ratio_acetabular_inclination)

# Remove rows where z-score is > 4 or < -4
df <- df[!(df$z_iliac > 4 | df$z_iliac < -4), ]
df <- df[!(df$z_acetabular_diameter > 4 | df$z_acetabular_diameter < -4), ]
# df <- df[!(df$z_acetabular_inclination > 4 | df$z_acetabular_inclination < -4), ]

# Remove z-score columns
df$z_iliac <- NULL
df$z_acetabular_diameter <- NULL
# df$z_acetabular_inclination <- NULL

# only keep left and right handedness
df <- df %>% filter(handedness %in% c(1, 0))

# regression
iliac_isthmus_breadth <- lm(lr_ratio_iliac_isthmus_breadth ~ sex + age_imaging_visit + handedness, data = df)
summary(iliac_isthmus_breadth)
### output
# Call:
#   lm(formula = lr_ratio_iliac_isthmus_breadth ~ sex + age + handedness, 
#      data = df)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.279978 -0.046782 -0.001463  0.045225  0.286376 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  9.110e-01  2.868e-03 317.691  < 2e-16 ***
#   sex          3.987e-03  7.196e-04   5.540 3.04e-08 ***
#   age          4.814e-04  4.816e-05   9.995  < 2e-16 ***
#   handedness  -2.401e-03  1.217e-03  -1.972   0.0486 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.06938 on 37591 degrees of freedom
# Multiple R-squared:  0.003876,	Adjusted R-squared:  0.003796 
# F-statistic: 48.75 on 3 and 37591 DF,  p-value: < 2.2e-16


acetabular_diameter <- lm(lr_ratio_acetabular_diameter ~ sex + age_imaging_visit + handedness, data = df)
summary(acetabular_diameter)
### output
# Call:
#   lm(formula = lr_ratio_acetabular_diameter ~ sex + age + handedness, 
#      data = df)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.176221 -0.027799 -0.000268  0.027421  0.173622 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 1.009e+00  1.719e-03 586.745  < 2e-16 ***
#   sex         2.073e-03  4.313e-04   4.805 1.55e-06 ***
#   age         1.085e-04  2.886e-05   3.761  0.00017 ***
#   handedness  3.263e-03  7.296e-04   4.472 7.77e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04159 on 37591 degrees of freedom
# Multiple R-squared:  0.001596,	Adjusted R-squared:  0.001516 
# F-statistic: 20.03 on 3 and 37591 DF,  p-value: 5.793e-13




























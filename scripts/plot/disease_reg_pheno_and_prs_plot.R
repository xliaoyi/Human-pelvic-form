rm(list = ls())
library(ggplot2)
library(reshape2)
library(patchwork)
library(tidyverse)
setwd("/Users/louis/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")
# outline

# 1. OA & walking phenotype
# 2. OA & walking PRS
# 3. pregnant phenotype & PRS

select_pheno <- c("pelvic_height", "pelvic_width", "pelvic_inlet_width", "subpubic_angle",
                  "oblique_pelvic_inlet_length", "iliac_isthmus_breadth", "acetabular_diameter", 
                  "leg_divide_torso") #  "bi_acetabular_width" "pelvic_inlet_area", 

common <- c('M54', 'M16', 'M17', 'M23',
            'back_pain', 'hip_pain', 'knee_pain', 'walking_pace')
select_common <- c('M16', 'M17', 'M23',
                   'hip_pain', 'knee_pain', 'walking_pace')
female_specific <- c('delivery_emergency_caesarean_section',
                     'birth_weight_first_child', 'n_weeks_gestation')

# rename phenotypes
rename_map_pheno <- list(
  "acetabular_diameter" = "Acetabular diameter",
  "arm_devide_torso" = "Arm:torso",
  "ear_left2ear_right" = "Head width",
  "iliac_isthmus_breadth" = "Iliac isthmus breadth",
  "oblique_pelvic_inlet_length" = "Oblique pelvic inlet length",
  "pelvic_height" = "Pelvic height",
  "pelvic_inlet_width" = "Pelvic inlet width",
  "pelvic_width" = "Pelvic width",
  "shoulder_width" = "Shoulder width",
  "subpubic_angle" = "Subpubic angle",
  "pelvic_inlet_area" = "Pelvic inlet area",
  "head_divide_inlet_width" = "Head:pelvic inlet width",
  "head_divide_oblique_inlet_length" = "Head:pelvic oblique inlet length",
  "head_area_divide_pelvic_inlet_area" = "Head area:pelvic inlet area",
  "shoulder_divide_inlet_width" = "Shoulder:pelvic inlet width",
  "shoulder_divide_oblique_inlet_length" = "Shoulder:pelvic oblique inlet length",
  "shoulder_area_divide_pelvic_inlet_area" = "Shoulder area:pelvic inlet area",
  "leg_length" = "Leg length",
  "leg_divide_torso" = "Leg:torso",
  "bi_acetabular_width" = "Biacetabular width",
  "BW_devide_pelvic_inlet_width" = "Baby birth weight:pelvic inlet width",
  "BW_devide_oblique_pelvic_inlet_length" = "Baby birth weight:oblique pelvic inlet length",
  "BW_devide_pelvic_inlet_area" = "Baby birth weight:pelvic inlet area"
)
rename_map_outcome <- list(
  "M54" = "Dorsalgia (M54)",
  "back_pain" = "Back pain",
  "M16" = "Hip OA (M16)",
  "hip_pain" = "Hip pain",
  "M17" = "Knee OA (M17)",
  "knee_pain" = "Knee pain",
  "M23" = "Internal derangement of knee (M23)",
  "walking_pace" = "Walking pace",
  "delivery_spontaneous_vertex" = "Spontaneous vertex",
  "delivery_emergency_caesarean_section" = "Emergency caesarean section",
  "birth_weight_first_child" = "Baby birth weight",
  "n_weeks_gestation" = "Number of gestation weeks",
  "N81" = "Genital prolapse (N81)",
  "keep_allna_incontinence" = "Incontinence"
)

newcommon <- c("Dorsalgia (M54)", "Back pain", "Hip OA (M16)", "Hip pain",
               "Knee OA (M17)", "Knee pain", "Internal derangement of knee (M23)", "Walking pace")
# newcommon <- c("Hip OA (M16)", "Hip pain",
#                "Knee OA (M17)", "Knee pain", "Internal derangement of knee (M23)", "Walking pace")
newpheno <- c("Leg:torso",
              "Acetabular diameter", "Arm:torso", "Biacetabular width", "Iliac isthmus breadth", 
              "Pelvic height", "Pelvic width", "Head width", "Shoulder width", "Subpubic angle",
              "Pelvic inlet width", "Oblique pelvic inlet length","Pelvic inlet area", 
              "Head:pelvic inlet width", "Head:pelvic oblique inlet length", "Head area:pelvic inlet area",
              "Shoulder:pelvic inlet width", "Shoulder:pelvic oblique inlet length", "Shoulder area:pelvic inlet area")
newfemale <- c("Emergency caesarean section", "Baby birth weight", "Number of gestation weeks", "Genital prolapse (N81)", "Incontinence")

pheno2remove <- c("Head:pelvic inlet width", "Head:pelvic oblique inlet length", "Head area:pelvic inlet area",
                  "Shoulder:pelvic inlet width", "Shoulder:pelvic oblique inlet length", "Shoulder area:pelvic inlet area", 
                  "Arm:torso", "Head width", "Shoulder width")
##################################
##### OA & walking phenotype ##### 
##################################

##################################
# female

df <- read.csv("key_results/disease_phenotypes_association_full_female_20231027.csv")
df <- df %>% filter(outcome %in% c(common)) %>% filter(pheno %in% select_pheno)

# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno),
    outcome = ifelse(outcome %in% names(rename_map_outcome), rename_map_outcome[outcome], outcome)
  )

df$outcome <- factor(df$outcome, levels = c(newcommon))
df$pheno <- factor(df$pheno, levels = c(newpheno))

df <- df %>% filter(!pheno %in% pheno2remove)

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))

# remove back OA back pain phenotypes
df <- df %>% filter(outcome %in% newcommon)

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p1 <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 3), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.8,0.9,1,1.1,1.2), 
                       labels=c(0.8,0.9,1,1.1,1.2), 
                       limits = c(0.8, 1.2), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Female") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = "out_fig/disease_regression_female_phenotype_231027.pdf", plot = p1,
       width = 5.6, height = 5.6, units = "in")


##################################
# male
df <- read.csv("key_results/disease_phenotypes_association_full_male_20231027.csv")
df <- df %>% filter(outcome %in% c(common)) %>% filter(pheno %in% select_pheno)

# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno),
    outcome = ifelse(outcome %in% names(rename_map_outcome), rename_map_outcome[outcome], outcome)
  )

df$outcome <- factor(df$outcome, levels = c(newcommon))
df$pheno <- factor(df$pheno, levels = c(newpheno))

df <- df %>% filter(!pheno %in% pheno2remove)

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))

# remove back OA back pain phenotypes
df <- df %>% filter(outcome %in% newcommon)

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p2 <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 3), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.8,0.9,1,1.1,1.2), 
                       labels=c(0.8,0.9,1,1.1,1.2), 
                       limits = c(0.8, 1.2), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Male") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggsave(filename = "out_fig/disease_regression_male_phenotype_231027.pdf", plot = p2,
       width = 5.6, height = 5.6, units = "in")





##################################
# both
df <- read.csv("key_results/disease_phenotypes_association_full_both_20240112.csv")
df <- df %>% filter(outcome %in% c(common)) %>% filter(pheno %in% select_pheno)

# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno),
    outcome = ifelse(outcome %in% names(rename_map_outcome), rename_map_outcome[outcome], outcome)
  )

df$outcome <- factor(df$outcome, levels = c(newcommon))
df$pheno <- factor(df$pheno, levels = c(newpheno))

df <- df %>% filter(!pheno %in% pheno2remove)

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))

# remove back OA back pain phenotypes
df <- df %>% filter(outcome %in% newcommon)

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 3), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.8,0.9,1,1.1,1.2), 
                       labels=c(0.8,0.9,1,1.1,1.2), 
                       limits = c(0.8, 1.2), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Both") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = "out_fig/disease_regression_both_phenotype_240204.pdf", plot = p,
       width = 5.6, height = 5.0, units = "in")


##################################
######## OA & walking PRS ######## 
##################################

##################################
# female

# df <- read.csv("key_results/disease_phenotypes_association_400k_female_prs_240111.csv")
df <- read.csv("key_results/disease_phenotypes_association_400k_female_prs_240222.csv")
df <- df %>% filter(outcome %in% c(common)) %>% filter(pheno %in% select_pheno)

# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno),
    outcome = ifelse(outcome %in% names(rename_map_outcome), rename_map_outcome[outcome], outcome)
  )

df$outcome <- factor(df$outcome, levels = c(newcommon))
df$pheno <- factor(df$pheno, levels = c(newpheno))

df <- df %>% filter(!pheno %in% pheno2remove)

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))

# remove back OA back pain phenotypes
df <- df %>% filter(outcome %in% newcommon)

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p1 <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 3), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.8,0.9,1,1.1,1.2), 
                       labels=c(0.8,0.9,1,1.1,1.2), 
                       limits = c(0.8, 1.2), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Female") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = "out_fig/disease_regression_female_prs_240111.pdf", plot = p1,
       width = 5.6, height = 5.6, units = "in")


##################################
# male
df <- read.csv("key_results/disease_phenotypes_association_400k_male_prs_231027.csv")
df <- df %>% filter(outcome %in% c(common)) %>% filter(pheno %in% select_pheno)

# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno),
    outcome = ifelse(outcome %in% names(rename_map_outcome), rename_map_outcome[outcome], outcome)
  )

df$outcome <- factor(df$outcome, levels = c(newcommon))
df$pheno <- factor(df$pheno, levels = c(newpheno))

df <- df %>% filter(!pheno %in% pheno2remove)

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))

# remove back OA back pain phenotypes
df <- df %>% filter(outcome %in% newcommon)

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p2 <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 3), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.8,0.9,1,1.1,1.2), 
                       labels=c(0.8,0.9,1,1.1,1.2), 
                       limits = c(0.8, 1.2), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Male") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = "out_fig/disease_regression_male_prs_231027.pdf", plot = p2,
       width = 5.6, height = 5.6, units = "in")





##################################
# both
# df <- read.csv("key_results/disease_phenotypes_association_400k_both_prs_231031.csv")
df <- read.csv("key_results/disease_phenotypes_association_400k_both_prs_240112.csv")
df <- df %>% filter(outcome %in% c(common)) %>% filter(pheno %in% select_pheno)

# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno),
    outcome = ifelse(outcome %in% names(rename_map_outcome), rename_map_outcome[outcome], outcome)
  )

df$outcome <- factor(df$outcome, levels = c(newcommon))
df$pheno <- factor(df$pheno, levels = c(newpheno))

df <- df %>% filter(!pheno %in% pheno2remove)

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))

# remove back OA back pain phenotypes
df <- df %>% filter(outcome %in% newcommon)

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 3), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.8,0.9,1,1.1,1.2), 
                       labels=c(0.8,0.9,1,1.1,1.2), 
                       limits = c(0.8, 1.2), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Both") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = "out_fig/disease_regression_both_prs_240204.pdf", plot = p,
       width = 5.6, height = 5.0, units = "in")



##################################
#### Pregnant phenotype & PRS #### 
##################################

##################################
# since number of idividual in imaging data set is too small, so I didn't perform phenotype association 
newpheno <- c("Acetabular diameter", "Iliac isthmus breadth", #, "Biacetabular width"
              "Pelvic height", "Pelvic width", "Subpubic angle",
              "Pelvic inlet width", "Oblique pelvic inlet length", #"Pelvic inlet area", 
              "Baby birth weight:pelvic inlet width", "Baby birth weight:oblique pelvic inlet length") # , "Baby birth weight:pelvic inlet area"

# female_specific <- c('delivery_emergency_caesarean_section', 'delivery_spontaneous_vertex',
#                      'birth_weight_first_child', 'n_weeks_gestation')


newfemale <- c("Emergency caesarean section", "Baby birth weight", "Number of gestation weeks") # , "Spontaneous vertex"

# PRS

# df <- read.csv("key_results/disease_phenotypes_association_400k_female_prs_231028.csv")
# df <- read.csv("key_results/disease_phenotypes_association_400k_female_prs_240112.csv")
df <- read.csv("key_results/disease_phenotypes_association_400k_female_prs_240222.csv")

# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno),
    outcome = ifelse(outcome %in% names(rename_map_outcome), rename_map_outcome[outcome], outcome)
  )

# df <- df %>% filter(outcome %in% c(newcommon, newfemale)) %>% filter(pheno %in% newpheno) # %>% filter(!pheno %in% pheno2remove)
df <- df %>% filter(outcome %in% c(newfemale)) %>% filter(pheno %in% newpheno)

df$outcome <- factor(df$outcome, levels = c(newfemale))
df$pheno <- factor(df$pheno, levels = c(newpheno))

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))

# keep pregnant phenotypes
# df <- df %>% filter(outcome %in% newfemale)

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p1 <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 3), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.8,0.9,1,1.1,1.2), 
                       labels=c(0.8,0.9,1,1.1,1.2), 
                       limits = c(0.8, 1.2), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Female PRS association") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = "out_fig/disease_regression_female_pregnant_prs_240204.pdf", plot = p1,
       width = 5.15, height = 5, units = "in")


##################################
########  Pelvic floor PRS ####### 
##################################
# df <- read.csv("key_results/disease_phenotypes_association_400k_female_prs_231028.csv")
# df <- read.csv("key_results/disease_phenotypes_association_400k_female_prs_240112.csv")
df <- read.csv("key_results/disease_phenotypes_association_400k_female_prs_240222.csv")
df <- df %>% filter(outcome %in% c("N81", "keep_allna_incontinence"))

# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno),
    outcome = ifelse(outcome %in% names(rename_map_outcome), rename_map_outcome[outcome], outcome)
  )

pheno2remove <- c("Baby birth weight:pelvic inlet width", "Baby birth weight:oblique pelvic inlet length", "Baby birth weight:pelvic inlet area", "Biacetabular width", "Pelvic inlet area")
df <- df %>% filter(!pheno %in% pheno2remove)

newfemale <- c("Genital prolapse (N81)", "Incontinence") # , "Spontaneous vertex"

df$pheno <- factor(df$pheno, levels = c(newpheno))
df$outcome <- factor(df$outcome, levels = c(newfemale))

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))
# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p1 <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 3), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.8,0.9,1,1.1,1.2), 
                       labels=c(0.8,0.9,1,1.1,1.2), 
                       limits = c(0.8, 1.2), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Female PRS association") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "out_fig/pelvic_floor_female_prs_240204.pdf", plot = p1,
       width = 3.9, height = 4.0, units = "in")


##################################
###### female skeletal PRS ####### 
##################################
df_skeletal <- read.csv("key_results/disease_phenotypes_association_400k_female_prs_240316_skeletal_blood_chem.csv")
df_pelvis <- read.csv("key_results/disease_phenotypes_association_400k_female_prs_240316_pelvis_blood_chem.csv")
df <- bind_rows(df_skeletal, df_pelvis)
write.csv(df,
          "key_results/disease_phenotypes_association_400k_female_prs_240316_skeletal_pelvis_blood_chem.csv",
          row.names = FALSE)


# Apply renaming to the df dataframe
# df <- df %>%
#   mutate(
#     pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno),
#     outcome = ifelse(outcome %in% names(rename_map_outcome), rename_map_outcome[outcome], outcome)
#   )
# 
# pheno2remove <- c("Baby birth weight:pelvic inlet width", "Baby birth weight:oblique pelvic inlet length", "Baby birth weight:pelvic inlet area", "Biacetabular width", "Pelvic inlet area")
# df <- df %>% filter(!pheno %in% pheno2remove)
# 
# newfemale <- c("Genital prolapse (N81)", "Incontinence") # , "Spontaneous vertex"
# 
# df$pheno <- factor(df$pheno, levels = c(newpheno))
newfemale <- c("pelvic_inlet_width", "oblique_pelvic_inlet_length","subpubic_angle",
               "pelvic_height","pelvic_width", "iliac_isthmus_breadth" ,"acetabular_diameter",
               "trochanter_distance", "shoulder_width","head_diameter","humerus","tibia",
               "femur","forearm" ,"torso_length","standing_height_imaging_visit")
df$pheno <- factor(df$pheno, levels = rev(newfemale))

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))
# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p1 <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 3), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.8,0.9,1,1.1,1.2), 
                       labels=c(0.8,0.9,1,1.1,1.2), 
                       limits = c(0.8, 1.2), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Female PRS association") +
  theme(plot.title = element_text(hjust = 0.5)) #, panel.grid.major = element_blank(), panel.grid.minor = element_blank()
ggsave(filename = "out_fig/skeletal_pelvis_blood_chem_female_prs_240316.pdf", plot = p1,
       width = 6, height = 9, units = "in")


##################################
####### male skeletal PRS ######## 
##################################
df <- read.csv("key_results/disease_phenotypes_association_400k_male_prs_240322_skeletal_pelvis_blood_chem.csv")


# Apply renaming to the df dataframe
# df <- df %>%
#   mutate(
#     pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno),
#     outcome = ifelse(outcome %in% names(rename_map_outcome), rename_map_outcome[outcome], outcome)
#   )
# 
# pheno2remove <- c("Baby birth weight:pelvic inlet width", "Baby birth weight:oblique pelvic inlet length", "Baby birth weight:pelvic inlet area", "Biacetabular width", "Pelvic inlet area")
# df <- df %>% filter(!pheno %in% pheno2remove)
# 
# newfemale <- c("Genital prolapse (N81)", "Incontinence") # , "Spontaneous vertex"
# 
# df$pheno <- factor(df$pheno, levels = c(newpheno))
newmale <- c("pelvic_inlet_width", "oblique_pelvic_inlet_length","subpubic_angle",
               "pelvic_height","pelvic_width", "iliac_isthmus_breadth" ,"acetabular_diameter",
               "trochanter_distance", "shoulder_width","head_diameter","humerus","tibia",
               "femur","forearm" ,"torso_length","standing_height_imaging_visit")
df$pheno <- factor(df$pheno, levels = rev(newmale))

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))
# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p1 <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 3), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.8,0.9,1,1.1,1.2), 
                       labels=c(0.8,0.9,1,1.1,1.2), 
                       limits = c(0.8, 1.2), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Male PRS association") +
  theme(plot.title = element_text(hjust = 0.5)) #, panel.grid.major = element_blank(), panel.grid.minor = element_blank()
ggsave(filename = "out_fig/skeletal_pelvis_blood_chem_male_prs_240322.pdf", plot = p1,
       width = 5, height = 9, units = "in")

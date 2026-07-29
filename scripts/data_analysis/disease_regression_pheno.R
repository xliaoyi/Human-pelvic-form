rm(list = ls())

library(dplyr)

setwd("/Users/louis/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape")
# df <- read.csv("key_results/phenotype_with_disease.csv")
# df_female <- read.csv("key_results/phenotype_with_disease_female.csv")
# df_male <- read.csv("key_results/phenotype_with_disease_male.csv")

# outcome phenotypes
fid_info <- read.csv("../UKB_xray_image_info/fids/fid_disease/fid_info_400k_20240222.csv")



##################################
############# Female ############# 
##################################
# residual
df_female <- read.csv("UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_female_cm_20240104.txt", sep = " ")
df_female <- subset(df_female, select = -c(FID))
df_female <- df_female %>% rename(eid = 'IID')
df_female <- merge(df_female, fid_info, by = 'eid')

# logistic regression
# outcome_phenos <- c('M54', 'M16', 'M17', 'M23', 'had_menopause', 'mom_has_child', 'taken_oral_contraceptive_pill', 'has_stillbirths', 'had_same_sex_intercourse', 
#                     'N72', 'N73', 'N80', 'N81', 'N83', 'N84', 'N85', 'N86', 'N87', 'N88', 'N89', 'N90', 'N92', 'N93', 'N94',
#                     'N95', 'N97', 'O03', 'O04', 'O20', 'O32', 'O24', 'O36', 'O42', 'O47', 'O48', 'O63', 'O68', 'O70', 'O72', 
#                     'no_pain', 'back_pain', 'neck_or_shoulder_pain', 'headache', 'knee_pain', 'stomach_or_abdominal_pain', 'hip_pain', 'facial_pain', 'pain_all_over_the_body', 
#                     'fracture_spine', 'fracture_other_bones', 'fracture_wrist', 'fracture_arm', 'fracture_leg', 'fracture_ankle', 'delivery_elective_caesarean_section', 
#                     'delivery_spontaneous_vertex', 'delivery_other_than_specified', 'delivery_emergency_caesarean_section', 'delivery_ventouse', 'delivery_onset_caesarean_section', 
#                     'delivery_onset_spontaneous', 'delivery_onset_medical_induction', 'male_baby', 'female_baby', 'has_caesarean', 'pregnancy_termination', 'walking_pace')
outcome_phenos <- c('M54', 'M16', 'M17', 'M23',
                    'back_pain', 'knee_pain', 'hip_pain',
                    'delivery_spontaneous_vertex', 'delivery_emergency_caesarean_section','walking_pace')

# measured_phenos <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
#                      'iliac_flare_ratio', 'oblique_pelvic_inlet_length', 'iliac_flare_angle',
#                      'subpubic_angle', 'ear_left2ear_right', 'iliac_isthmus_breadth',
#                      'acetabular_diameter', 'acetabular_inclination', 'trochanter_left2trochanter_right',
#                      'shoulder_width', 'arm_devide_torso', 
#                      'head_divide_inlet_width', 'shoulder_divide_inlet_width',
#                      'head_divide_oblique_inlet_length',
#                      'shoulder_divide_oblique_inlet_length', 'pelvic_inlet_area',
#                      'head_area_divide_pelvic_inlet_area',
#                      'shoulder_area_divide_pelvic_inlet_area')
measured_phenos <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
                     'oblique_pelvic_inlet_length',
                     'subpubic_angle', 'iliac_isthmus_breadth',
                     'acetabular_diameter')

# covar <- c("age", "bmi", "standing_height", "E11", "F32", "F33", "I25", "pelvis_bmd",
#            "femur_bmd_left", "femur_bmd_right", "spine_bmd", 'smoking', 
#            "father_death_age", "mother_death_age", "blood_pressure", "sleep_duration", "household_income", "n_live_birth")

covar <- c("age", "sex", "standing_height", "weight", "household_income")

pheno_scores <- df_female[measured_phenos]
pheno_z_scores <- scale(pheno_scores)
pheno_z_scores <- as.data.frame(pheno_z_scores)
df_female[measured_phenos] <- pheno_z_scores[measured_phenos]

results_female <- data.frame()

for (outcome in outcome_phenos) {
  print(outcome)
  tmp_df <- na.omit(df_female[c(outcome, measured_phenos, covar)])
  print(dim(tmp_df))
  for (pheno in measured_phenos) {
    
    # for ratio pheno, shouldn't adjusting for height
    if (grepl("devide|divide", pheno)) {
      covar <- c("age", "weight", "household_income")
    }  else {
      covar <- c("age", "standing_height", "weight", "household_income")
    }
    
    # Run logistic regression
    formula <- as.formula(paste(outcome, " ~ ", pheno, 
                                " + ", paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "binomial")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))[pheno, "Estimate"]
    p_value <- coef(summary(res))[pheno, "Pr(>|z|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results_female data frame
    results_female <- rbind(results_female, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
}

# linear regression
# outcome_phenos = c('patient_birth_weight', 'hip_circumference', 'periods_started_age', 'n_live_birth', 'birth_weight_first_child', 
#                    'first_live_birth_age', 'last_live_birth_age', 'oral_contraceptive_pill_started_age', 'oral_contraceptive_pill_stopped_age', 
#                    'periods_stopped_age', 'menstrual_cycle_length', 'n_spontaneous_miscarriages', 'n_pregnancy_terminations', 
#                    'primiparous_age', 'igf_1', 'testosterone', 'shbg', 'oestradiol', 'n_same_sex_sexual_partner', 
#                    'n_sexual_partner', 'walking_pace', 'n_pregnancy', 'n_weeks_gestation', 'baby_birth_weight', 'postnatal_stay_duration')

outcome_phenos = c('birth_weight_first_child', 'n_weeks_gestation', "n_live_birth")

for (outcome in outcome_phenos) {
  print(outcome)
  tmp_df <- na.omit(df_female[c(outcome, measured_phenos, covar)])
  print(dim(tmp_df))
  for (pheno in measured_phenos) {
    
    # for ratio pheno, shouldn't adjusting for height
    if (grepl("devide|divide", pheno)) {
      covar <- c("age", "weight", "household_income")
    }  else {
      covar <- c("age", "standing_height", "weight", "household_income")
    }
    
    # Run linear regression
    formula <- as.formula(paste(outcome, " ~ ", pheno, 
                                " + ", paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "gaussian")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))[pheno, "Estimate"]
    p_value <- coef(summary(res))[pheno, "Pr(>|t|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results_female data frame
    results_female <- rbind(results_female, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
}

print(results_female)

# save
write.csv(results_female, file = "key_results/disease_phenotypes_association_full_female_20231027.csv", row.names = FALSE)

################################
############# Male #############
################################

df_male <- read.csv("UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_male_cm_231025.txt", sep = " ")
df_male <- subset(df_male, select = -c(FID))
df_male <- df_male %>% rename(eid = 'IID')
df_male <- merge(df_male, fid_info, by = 'eid')

# logistic regression
# outcome_phenos <- c('M54', 'M16', 'M17', 'M23', 'dad_has_child', 'had_same_sex_intercourse', 'no_pain', 
#                     'back_pain', 'neck_or_shoulder_pain', 'headache', 'knee_pain', 'stomach_or_abdominal_pain', 
#                     'hip_pain', 'facial_pain', 'pain_all_over_the_body', 'fracture_spine', 'fracture_other_bones', 
#                     'fracture_wrist', 'fracture_arm', 'fracture_leg', 'fracture_ankle', 'walking_pace')
outcome_phenos <- c('M54', 'M16', 'M17', 'M23',
                    'back_pain', 'knee_pain', 'hip_pain',
                    'walking_pace')

# measured_phenos <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
#                      'iliac_flare_ratio', 'oblique_pelvic_inlet_length', 'iliac_flare_angle',
#                      'subpubic_angle', 'ear_left2ear_right', 'iliac_isthmus_breadth',
#                      'acetabular_diameter', 'acetabular_inclination', 'trochanter_left2trochanter_right',
#                      'shoulder_width', 'arm_devide_torso', 
#                      'head_divide_inlet_width', 'shoulder_divide_inlet_width',
#                      'head_divide_oblique_inlet_length',
#                      'shoulder_divide_oblique_inlet_length', 'pelvic_inlet_area',
#                      'head_area_divide_pelvic_inlet_area',
#                      'shoulder_area_divide_pelvic_inlet_area')
measured_phenos <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
                     'oblique_pelvic_inlet_length',
                     'subpubic_angle', 'ear_left2ear_right', 'iliac_isthmus_breadth',
                     'acetabular_diameter', 'trochanter_left2trochanter_right',
                     'shoulder_width', 'arm_devide_torso', 
                     'head_divide_inlet_width', 'shoulder_divide_inlet_width',
                     'head_divide_oblique_inlet_length',
                     'shoulder_divide_oblique_inlet_length', 'pelvic_inlet_area',
                     'head_area_divide_pelvic_inlet_area',
                     'shoulder_area_divide_pelvic_inlet_area',
                     'arm_divide_leg', 'leg_divide_torso', 'bi_acetabular_width')


# covar <- c("age", "bmi", "standing_height", "E11", "F32", "F33", "I25", "pelvis_bmd",
#            "femur_bmd_left", "femur_bmd_right", "spine_bmd", 'smoking', 
#            "father_death_age", "mother_death_age", "blood_pressure", "sleep_duration", "household_income")

covar <- c("age", "standing_height", "weight", "household_income")

pheno_scores <- df_male[measured_phenos]
pheno_z_scores <- scale(pheno_scores)
pheno_z_scores <- as.data.frame(pheno_z_scores)
df_male[measured_phenos] <- pheno_z_scores[measured_phenos]

results_male <- data.frame()

for (outcome in outcome_phenos) {
  print(outcome)
  tmp_df <- na.omit(df_male[c(outcome, measured_phenos, covar)])
  for (pheno in measured_phenos) {
    
    # for ratio pheno, shouldn't adjusting for height
    if (grepl("devide|divide", pheno)) {
      covar <- c("age", "weight", "household_income")
    }  else {
      covar <- c("age", "standing_height", "weight", "household_income")
    }
    
    # Run logistic regression
    formula <- as.formula(paste(outcome, " ~ ", pheno, 
                                " + ", paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "binomial")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))[pheno, "Estimate"]
    p_value <- coef(summary(res))[pheno, "Pr(>|z|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results_male data frame
    results_male <- rbind(results_male, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
}

# linear regression
# outcome_phenos = c('patient_birth_weight', 'hip_circumference','n_child_fathered', 
#                    'igf_1', 'testosterone', 'shbg', 'oestradiol', 'n_same_sex_sexual_partner', 'n_sexual_partner', 'walking_pace')
# 
# for (outcome in outcome_phenos) {
#   print(outcome)
#   tmp_df <- na.omit(df_male[c(outcome, measured_phenos, covar)])
#   for (pheno in measured_phenos) {
#     
#     # Run linear regression
#     formula <- as.formula(paste(outcome, " ~ ", pheno, 
#                                 " + age + bmi + standing_height + 
#                                 E11 + F32 + F33 + I25 + pelvis_bmd +
#                                 femur_bmd_left + femur_bmd_right + spine_bmd +
#                                 smoking + 
#                                 father_death_age + mother_death_age + 
#                                 blood_pressure + sleep_duration + household_income"))
#     res <- glm(formula, data = tmp_df, family = "gaussian")
#     
#     # Save beta coefficient and p-value
#     beta <- coef(summary(res))[pheno, "Estimate"]
#     p_value <- coef(summary(res))[pheno, "Pr(>|t|)"]
#     neglog10p <- -log10(p_value)
#     odds_ratio <- exp(beta)
#     
#     # Append to results_male data frame
#     results_male <- rbind(results_male, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
#   }
# }

print(results_male)

# save
write.csv(results_male, file = "key_results/disease_phenotypes_association_full_male_20231027.csv", row.names = FALSE)


##################################
############# Both ############### 
##################################
leg_torso <- read.csv("key_results/hip_select_pheno_cm_eid_flt_z_flt.csv") %>% select(c("eid", 'leg_divide_torso'))

df <- read.csv("UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_both_cm_20240104.txt", sep = " ")
df <- subset(df, select = -c(FID))
df <- df %>% rename(eid = 'IID')
df <- merge(df, fid_info, by = 'eid')

df <- merge(df, leg_torso, by = 'eid')

# logistic regression
# outcome_phenos <- c('M54', 'M16', 'M17', 'M23', 'had_same_sex_intercourse', 'no_pain', 'back_pain', 'neck_or_shoulder_pain', 'headache', 
#                     'knee_pain', 'stomach_or_abdominal_pain', 'hip_pain', 'facial_pain', 'pain_all_over_the_body', 'fracture_spine', 
#                     'fracture_other_bones', 'fracture_wrist', 'fracture_arm', 'fracture_leg', 'fracture_ankle', 'walking_pace')
outcome_phenos <- c('M54', 'M16', 'M17', 'M23',
                    'back_pain', 'knee_pain', 'hip_pain',
                    'walking_pace')

# measured_phenos <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
#                      'iliac_flare_ratio', 'oblique_pelvic_inlet_length', 'iliac_flare_angle',
#                      'subpubic_angle', 'ear_left2ear_right', 'iliac_isthmus_breadth',
#                      'acetabular_diameter', 'acetabular_inclination', 'trochanter_left2trochanter_right',
#                      'shoulder_width', 'arm_devide_torso', 
#                      'head_divide_inlet_width', 'shoulder_divide_inlet_width',
#                      'head_divide_oblique_inlet_length',
#                      'shoulder_divide_oblique_inlet_length', 'pelvic_inlet_area',
#                      'head_area_divide_pelvic_inlet_area',
#                      'shoulder_area_divide_pelvic_inlet_area')
measured_phenos <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
                     'oblique_pelvic_inlet_length', 'pelvic_inlet_area',
                     'subpubic_angle', 'iliac_isthmus_breadth',
                     'acetabular_diameter', 'leg_divide_torso')

# covar <- c("age", "sex", "bmi", "standing_height", "weight", "E11", "F32", "F33", "I25", "pelvis_bmd",
#            "femur_bmd_left", "femur_bmd_right", "spine_bmd", 'smoking',
#            "father_death_age", "mother_death_age", "blood_pressure", "sleep_duration", "household_income")
covar <- c("age_imaging_visit", "sex", "weight_imaging_visit", "standing_height_imaging_visit",
           "household_income", "E11", "F32", "F33", "I25", 'smoking', "sleep_duration")
# covar <- c("age", "sex", "standing_height", "weight", "household_income")

pheno_scores <- df[measured_phenos]
pheno_z_scores <- scale(pheno_scores)
pheno_z_scores <- as.data.frame(pheno_z_scores)
df[measured_phenos] <- pheno_z_scores[measured_phenos]

results <- data.frame()

for (outcome in outcome_phenos) {
  print(outcome)
  for (pheno in measured_phenos) {
    
    # # for ratio pheno, shouldn't adjusting for height
    if (grepl("leg_over_torso", pheno)) {
      # covar <- c("age", "sex", "weight", "household_income")
      covar <- c("age_imaging_visit", "sex", "weight_imaging_visit", # "standing_height_imaging_visit",
                 "household_income", "E11", "F32", "F33", "I25", 'smoking', "sleep_duration")
      tmp_df <- na.omit(df[c(outcome, measured_phenos, covar)])
      print(dim(tmp_df))
    }  else {
      # covar <- c("age", "sex", "standing_height", "weight", "household_income")
      covar <- c("age_imaging_visit", "sex", "weight_imaging_visit", "standing_height_imaging_visit",
                 "household_income", "E11", "F32", "F33", "I25", 'smoking', "sleep_duration")
      tmp_df <- na.omit(df[c(outcome, measured_phenos, covar)])
      print(dim(tmp_df))
    }
    
    # Run logistic regression
    formula <- as.formula(paste(outcome, " ~ ", pheno, 
                                " + ", paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "binomial")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))[pheno, "Estimate"]
    p_value <- coef(summary(res))[pheno, "Pr(>|z|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results data frame
    results <- rbind(results, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
}

# linear regression
# outcome_phenos = c('patient_birth_weight', 'hip_circumference', 'igf_1', 'testosterone', 'shbg', 'oestradiol',
#                    'n_same_sex_sexual_partner', 'n_sexual_partner',  'walking_pace')
# outcome_phenos = c('acceleration')
# 
# for (outcome in outcome_phenos) {
#   print(outcome)
#   tmp_df <- na.omit(df[c(outcome, measured_phenos, covar)])
#   for (pheno in measured_phenos) {
# 
#     # Run linear regression
#     formula <- as.formula(paste(outcome, " ~ ", pheno, 
#                                 " + ", paste(covar, collapse = " + ")))
#     res <- glm(formula, data = tmp_df, family = "gaussian")
# 
#     # Save beta coefficient and p-value
#     beta <- coef(summary(res))[pheno, "Estimate"]
#     p_value <- coef(summary(res))[pheno, "Pr(>|t|)"]
#     neglog10p <- -log10(p_value)
#     odds_ratio <- exp(beta)
# 
#     # Append to results data frame
#     results <- rbind(results, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
#   }
# }

print(results)

# save
# write.csv(results, file = "key_results/disease_phenotypes_association_full_both_20231027.csv", row.names = FALSE)
write.csv(results, file = "key_results/disease_phenotypes_association_full_both_20240112.csv", row.names = FALSE)

##########################################################################################
###                                   Plot                                             ###
##########################################################################################

library(ggplot2)
library(reshape2)
library(patchwork)

select_pheno <- c("pelvic_height", "pelvic_width", "pelvic_inlet_width", "subpubic_angle",
                  "oblique_pelvic_inlet_length", "ear_left2ear_right", "arm_devide_torso",
                  "shoulder_width", "iliac_isthmus_breadth", "acetabular_diameter", "pelvic_inlet_area", 
                  "head_divide_inlet_width", "shoulder_divide_inlet_width",
                  "head_divide_oblique_inlet_length", "shoulder_divide_oblique_inlet_length", 
                  "head_area_divide_pelvic_inlet_area", "shoulder_area_divide_pelvic_inlet_area")

common <- c('M54', 'M16', 'M17', 'M23',
            'back_pain', 'knee_pain', 'hip_pain', 'walking_pace')
female_specific <- c('delivery_spontaneous_vertex', 'delivery_emergency_caesarean_section',
                     'birth_weight_first_child', 'n_weeks_gestation')

##################################
############# Female ############# 
##################################

df <- read.csv("key_results/disease_phenotypes_association_full_female_20231008.csv")
df <- df %>% filter(outcome %in% c(common, female_specific)) %>% filter(pheno %in% select_pheno)
df$outcome <- factor(df$outcome, levels = c(common, female_specific))

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 6), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.7,0.8,0.9,1,1.1,1.2,1.3), 
                       labels=c(0.7,0.8,0.9,1,1.1,1.2,1.3), 
                       limits = c(0.7, 1.3), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Female")

ggsave(filename = "out_fig/disease_regression_female.pdf", plot = p,
       width = 14, height = 8, units = "in")

# Print the plot
print(p)

##################################
############# Male ############### 
##################################

df <- read.csv("key_results/disease_phenotypes_association_full_male_20231008.csv")
df <- df %>% filter(outcome %in% c(common, male_specific)) %>% filter(pheno %in% select_pheno)
df$outcome <- factor(df$outcome, levels = c(common, male_specific))

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 6), fill = df$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.7,0.8,0.9,1,1.1,1.2,1.3), 
                       labels=c(0.7,0.8,0.9,1,1.1,1.2,1.3), 
                       limits = c(0.7, 1.3), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Male")

ggsave(filename = "out_fig/disease_regression_male.pdf", plot = p,
       width = 9.73, height = 8, units = "in")

# Print the plot
print(p)


##################################
############# Both ############### 
##################################

df <- read.csv("key_results/disease_phenotypes_association_full_both_20231008.csv")
df <- df %>% filter(outcome %in% c(common)) %>% filter(pheno %in% select_pheno)
df$outcome <- factor(df$outcome, levels = c(common))

df$FDR <- p.adjust(df$p_value, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(df$p_value), 6), fill = df$odds_ratio), shape = 21) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.7,0.8,0.9,1,1.1,1.2,1.3), 
                       labels=c(0.7,0.8,0.9,1,1.1,1.2,1.3), 
                       limits = c(0.7, 1.3), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(df, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Odds Ratio", size = "-log10(p-value)") +
  ggtitle("Both male and female")

ggsave(filename = "out_fig/disease_regression_both.pdf", plot = p,
       width = 9.5, height = 8, units = "in", limitsize = FALSE)

# Print the plot
print(p)

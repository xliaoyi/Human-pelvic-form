rm(list = ls())
library(dplyr)
setwd("/Users/louis/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/PRS_OUTPUT")

PCs <- read.delim("../fid22009.csv", sep = ',')
colnames(PCs) <- c("eid", paste0("PC",1:40)) 
pcs <- PCs[,c("eid", paste0("PC", 1:20))]

fid.info <- read.delim("../../../UKB_xray_image_info/fids/fid_disease/fid_info_400k_20240222.csv", sep = ',')

# merge pcs and fid.info (outcome phenotypes)
pcs_fid <- merge(pcs, fid.info, by.x = "eid", by.y = "eid")

male_train_eids_old <- read.csv("/Users/louis/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_male_cm_230713.txt", sep = " ")
male_train_eids_new <- read.csv("/Users/louis/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/GWAS_INPUT_DATA/male_skeletal_gwas_cm_240318.txt", sep = " ")
# female_train_eids <- read.csv("/Users/louis/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_female_cm_230713.txt", sep = " ")
female_train_eids <- read.csv("/Users/louis/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/GWAS_INPUT_DATA/female_skeletal_gwas_cm_240310.txt", sep = " ")

both_train_eids <- read.csv("/Users/louis/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_cm_230713.txt", sep = " ")
male_train_eids_old <- as.vector(male_train_eids_old[,"FID"])
male_train_eids_new <- as.vector(male_train_eids_new[,"FID"])
male_train_eids <- unique(c(male_train_eids_old, male_train_eids_new))
female_train_eids <- as.vector(female_train_eids[,"FID"])
both_train_eids <- as.vector(both_train_eids[,"FID"])
# train_eids <- c(male_train_eids, female_train_eids)
# train_eids <- read.csv("../GWAS_INPUT_DATA/hip_select_pheno_gwas_both_cm_231109.txt", sep = " ")
# train_eids <- as.vector(train_eids[,"FID"])

##################################
############# Female ############# 
##################################

# union R32 and R15
pcs_fid <- pcs_fid %>% mutate(R32_and_R15 = ifelse(R32 == 1 | R15 == 1, 1, 0))

# logistic regression
# outcome_phenos <- c('R15', 'R32', 'N80', 'N81')
outcome_phenos <- c('M54', 'M16', 'M17', 'M23',
                    'back_pain', 'knee_pain', 'hip_pain',
                    'delivery_spontaneous_vertex', 'delivery_emergency_caesarean_section','walking_pace', 
                    'N81', 'keep_allna_incontinence') # , 'remove_allna_incontinence'

# measured_phenos <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
#                      'oblique_pelvic_inlet_length',
#                      'subpubic_angle', 'ear_left2ear_right', 'iliac_isthmus_breadth',
#                      'acetabular_diameter', 'trochanter_left2trochanter_right',
#                      'shoulder_width', 'arm_devide_torso',
#                      'head_divide_inlet_width', 'shoulder_divide_inlet_width',
#                      'head_divide_oblique_inlet_length',
#                      'shoulder_divide_oblique_inlet_length', 'pelvic_inlet_area',
#                      'head_area_divide_pelvic_inlet_area',
#                      'shoulder_area_divide_pelvic_inlet_area',
#                      'arm_divide_leg', 'leg_divide_torso', 'bi_acetabular_width',
#                      'BW_devide_pelvic_inlet_width', 'BW_devide_oblique_pelvic_inlet_length', 'BW_devide_pelvic_inlet_area')

measured_phenos <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
                     'oblique_pelvic_inlet_length',
                     'subpubic_angle', 'iliac_isthmus_breadth',
                     'acetabular_diameter') # , 'BW_devide_pelvic_inlet_width', 'BW_devide_oblique_pelvic_inlet_length', 'BW_devide_pelvic_inlet_area'
                                            # ,'pelvic_inlet_area', 'bi_acetabular_width', 'BW_devide_pelvic_inlet_width', 'BW_devide_oblique_pelvic_inlet_length', 'BW_devide_pelvic_inlet_area'

covar <- c("SCORESUM", "weight", "household_income", "n_live_birth","E11", "F32", "F33", "I25", 'smoking', "sleep_duration") #, "blood_pressure" 
# covar <- c("SCORESUM", "weight", "household_income", "n_live_birth")

results_female <- data.frame()

for (outcome in outcome_phenos) {
  print(outcome)
  for (pheno in measured_phenos) {
    prs_df <- read.table(paste0("hip_female_400k_240109/prs_output_", # hip_female_400k_230913
                                pheno, 
                                "_pst_eff_a1_b0.5_phiauto_all.profile"), header = T)[, c("FID", "SCORESUM")]
    # filter imaged patients
    prs_df <- prs_df %>% filter(!FID %in% female_train_eids)
    
    # merge prs, outcome phenotypes, pcs
    tmp_df <- merge(prs_df, pcs_fid, by.x = "FID", by.y = "eid")
    tmp_df <- tmp_df %>% filter(sex == 0)
    tmp_df <- na.omit(tmp_df[c(outcome, covar)])
    
    pheno_scores <- tmp_df["SCORESUM"]
    pheno_z_scores <- scale(pheno_scores)
    pheno_z_scores <- as.data.frame(pheno_z_scores)
    tmp_df["SCORESUM"] <- pheno_z_scores["SCORESUM"]
    
    # Run logistic regression
    formula <- as.formula(paste(outcome, " ~ ",  paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "binomial")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))["SCORESUM", "Estimate"]
    p_value <- coef(summary(res))["SCORESUM", "Pr(>|z|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results_female data frame
    results_female <- rbind(results_female, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
  print(dim(tmp_df))
}

# linear regression
outcome_phenos = c('birth_weight_first_child', 'n_weeks_gestation', "n_live_birth")

# outcome_phenos = c('first_live_birth_age', 'last_live_birth_age', 'n_live_birth', 'n_pregnancy_terminations', 
#                    'n_sexual_partner', 'n_spontaneous_miscarriages', 'n_stillbirths')

for (outcome in outcome_phenos) {
  print(outcome)
  for (pheno in measured_phenos) {
    prs_df <- read.table(paste0("hip_female_400k_240109/prs_output_", 
                                pheno, 
                                "_pst_eff_a1_b0.5_phiauto_all.profile"), header = T)[, c("FID", "SCORESUM")]
    # filter imaged patients
    prs_df <- prs_df %>% filter(!FID %in% female_train_eids)
    
    # merge prs, outcome phenotypes, pcs
    tmp_df <- merge(prs_df, pcs_fid, by.x = "FID", by.y = "eid")
    tmp_df <- tmp_df %>% filter(sex == 0)
    tmp_df <- na.omit(tmp_df[c(outcome, covar)])
    
    pheno_scores <- tmp_df["SCORESUM"]
    pheno_z_scores <- scale(pheno_scores)
    pheno_z_scores <- as.data.frame(pheno_z_scores)
    tmp_df["SCORESUM"] <- pheno_z_scores["SCORESUM"]
    
    # Run linear regression
    formula <- as.formula(paste(outcome, " ~ ",  paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "gaussian")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))["SCORESUM", "Estimate"]
    p_value <- coef(summary(res))["SCORESUM", "Pr(>|t|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results_female data frame
    results_female <- rbind(results_female, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
  print(dim(tmp_df))
}

print(results_female)

# save
write.csv(results_female, 
          file = "../../key_results/disease_phenotypes_association_400k_female_prs_240222.csv", 
          row.names = FALSE)

##################################
####### Female blood chem ######## 
##################################

# skeletal traits

results_female <- data.frame()
measured_phenos <- c('standing_height_imaging_visit','head_diameter','trochanter_distance','shoulder_width',
                     'torso_length','humerus','femur','forearm','tibia') 
covar <- c("SCORESUM", "weight", "household_income", "n_live_birth","E11", "F32", "F33", "I25", 'smoking', "sleep_duration" ) 

outcome_phenos = c('testosterone', 'igf_1', 'shbg', 'oestradiol',
                   'periods_started_age', 'periods_stopped_age')

for (outcome in outcome_phenos) {
  print(outcome)
  for (pheno in measured_phenos) {
    prs_df <- read.table(paste0("hip_female_400k_240313/prs_output_", 
                                pheno, 
                                "_pst_eff_a1_b0.5_phiauto_all.profile"), header = T)[, c("FID", "SCORESUM")]
    # filter imaged patients
    prs_df <- prs_df %>% filter(!FID %in% female_train_eids)
    
    # merge prs, outcome phenotypes, pcs
    tmp_df <- merge(prs_df, pcs_fid, by.x = "FID", by.y = "eid")
    tmp_df <- tmp_df %>% filter(sex == 0)
    tmp_df <- na.omit(tmp_df[c(outcome, covar)])
    
    pheno_scores <- tmp_df["SCORESUM"]
    pheno_z_scores <- scale(pheno_scores)
    pheno_z_scores <- as.data.frame(pheno_z_scores)
    tmp_df["SCORESUM"] <- pheno_z_scores["SCORESUM"]
    
    # Run linear regression
    formula <- as.formula(paste(outcome, " ~ ",  paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "gaussian")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))["SCORESUM", "Estimate"]
    p_value <- coef(summary(res))["SCORESUM", "Pr(>|t|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results_female data frame
    results_female <- rbind(results_female, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
  print(dim(tmp_df))
}

print(results_female)

# save
write.csv(results_female, 
          file = "../../key_results/disease_phenotypes_association_400k_female_prs_240316_skeletal_blood_chem.csv", 
          row.names = FALSE)

# pelvic traits

pcs_fid$duration.able.give.birth <- pcs_fid$periods_stopped_age - pcs_fid$periods_started_age
pcs_fid$mean.period <- (pcs_fid$periods_stopped_age + pcs_fid$periods_started_age)/2

results_female <- data.frame()
measured_phenos <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
                     'oblique_pelvic_inlet_length',
                     'subpubic_angle', 'iliac_isthmus_breadth',
                     'acetabular_diameter') 

covar <- c("SCORESUM", "weight", "household_income", "n_live_birth","E11", "F32", "F33", "I25", 'smoking', "sleep_duration" ) 

# outcome_phenos = c('testosterone', 'igf_1', 'shbg', 'oestradiol',
#                    'periods_started_age', 'periods_stopped_age')
outcome_phenos = c('mean.period')

for (outcome in outcome_phenos) {
  print(outcome)
  for (pheno in measured_phenos) {
    prs_df <- read.table(paste0("hip_female_400k_240109/prs_output_", 
                                pheno, 
                                "_pst_eff_a1_b0.5_phiauto_all.profile"), header = T)[, c("FID", "SCORESUM")]
    # filter imaged patients
    prs_df <- prs_df %>% filter(!FID %in% female_train_eids)
    
    # merge prs, outcome phenotypes, pcs
    tmp_df <- merge(prs_df, pcs_fid, by.x = "FID", by.y = "eid")
    tmp_df <- tmp_df %>% filter(sex == 0)
    tmp_df <- na.omit(tmp_df[c(outcome, covar)])
    
    pheno_scores <- tmp_df["SCORESUM"]
    pheno_z_scores <- scale(pheno_scores)
    pheno_z_scores <- as.data.frame(pheno_z_scores)
    tmp_df["SCORESUM"] <- pheno_z_scores["SCORESUM"]
    
    # Run linear regression
    formula <- as.formula(paste(outcome, " ~ ",  paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "gaussian")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))["SCORESUM", "Estimate"]
    p_value <- coef(summary(res))["SCORESUM", "Pr(>|t|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results_female data frame
    results_female <- rbind(results_female, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
  print(dim(tmp_df))
}

print(results_female)

# save
write.csv(results_female, 
          file = "../../key_results/disease_phenotypes_association_400k_female_prs_240316_pelvis_blood_chem.csv", 
          row.names = FALSE)
##################################
############# Male ############# 
##################################

# logistic regression
# outcome_phenos <- c('M54', 'M16', 'M17', 'M23', 'dad_has_child', 'had_same_sex_intercourse', 'no_pain', 
#                     'back_pain', 'neck_or_shoulder_pain', 'headache', 'knee_pain', 'stomach_or_abdominal_pain', 
#                     'hip_pain', 'facial_pain', 'pain_all_over_the_body', 'fracture_spine', 'fracture_other_bones', 
#                     'fracture_wrist', 'fracture_arm', 'fracture_leg', 'fracture_ankle', 'walking_pace')
outcome_phenos <- c('M54', 'M16', 'M17', 'M23',
                    'back_pain', 'knee_pain', 'hip_pain',
                    'walking_pace', 'erectile_dysfunction')

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

# covar <- c("SCORESUM", "weight", "household_income", "E11", "F32", "F33", "I25", 'smoking', "sleep_duration")
covar <- c("SCORESUM", "weight", "household_income", "E11", "F32", "F33", "I25", 'smoking', "sleep_duration" ) 

results_male <- data.frame()

for (outcome in outcome_phenos) {
  print(outcome)
  for (pheno in measured_phenos) {
    prs_df <- read.table(paste0("hip_male_400k_230913/prs_output_", 
                                pheno, 
                                "_pst_eff_a1_b0.5_phiauto_all.profile"), header = T)[, c("FID", "SCORESUM")]
    # filter imaged patients
    prs_df <- prs_df %>% filter(!FID %in% male_train_eids)
    
    # merge prs, outcome phenotypes, pcs
    tmp_df <- merge(prs_df, pcs_fid, by.x = "FID", by.y = "eid")
    tmp_df <- tmp_df %>% filter(sex == 1)
    tmp_df <- na.omit(tmp_df[c(outcome, covar)])
    
    pheno_scores <- tmp_df["SCORESUM"]
    pheno_z_scores <- scale(pheno_scores)
    pheno_z_scores <- as.data.frame(pheno_z_scores)
    tmp_df["SCORESUM"] <- pheno_z_scores["SCORESUM"]
    
    # Run logistic regression
    formula <- as.formula(paste(outcome, " ~ ",  paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "binomial")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))["SCORESUM", "Estimate"]
    p_value <- coef(summary(res))["SCORESUM", "Pr(>|z|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results_male data frame
    results_male <- rbind(results_male, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
  print(dim(tmp_df))
}

# linear regression
# outcome_phenos = c('weight', 'patient_birth_weight', 'hip_circumference', 'n_child_fathered',
#                    'igf_1', 'testosterone', 'n_same_sex_sexual_partner', 'n_sexual_partner', 
#                    'shbg', 'oestradiol')
# 
# for (outcome in outcome_phenos) {
#   print(outcome)
#   for (pheno in measured_phenos) {
#     prs_df <- read.table(paste0("hip_male_400k_230913/prs_output_", 
#                                 pheno, 
#                                 "_pst_eff_a1_b0.5_phiauto_all.profile"), header = T)[, c("FID", "SCORESUM")]
#     # filter imaged patients
#     prs_df <- prs_df %>% filter(!FID %in% male_train_eids)
#     
#     # merge prs, outcome phenotypes, pcs
#     tmp_df <- merge(prs_df, pcs_fid, by.x = "FID", by.y = "eid")
#     tmp_df <- tmp_df %>% filter(sex == 1)
#     tmp_df <- na.omit(tmp_df[c(outcome, covar)])
#     
#     pheno_scores <- tmp_df["SCORESUM"]
#     pheno_z_scores <- scale(pheno_scores)
#     pheno_z_scores <- as.data.frame(pheno_z_scores)
#     tmp_df["SCORESUM"] <- pheno_z_scores["SCORESUM"]
#     
#     # Run linear regression
#     formula <- as.formula(paste(outcome, " ~ ",  paste(covar, collapse = " + ")))
#     res <- glm(formula, data = tmp_df, family = "gaussian")
#     
#     # Save beta coefficient and p-value
#     beta <- coef(summary(res))["SCORESUM", "Estimate"]
#     p_value <- coef(summary(res))["SCORESUM", "Pr(>|t|)"]
#     neglog10p <- -log10(p_value)
#     odds_ratio <- exp(beta)
#     
#     # Append to results_male data frame
#     results_male <- rbind(results_male, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
#   }
#   print(dim(tmp_df))
# }

print(results_male)

# save
write.csv(results_male, 
          file = "../../key_results/disease_phenotypes_association_400k_male_prs_240222.csv", 
          row.names = FALSE)


##################################
######## Male blood chem #########
##################################

# skeletal traits

results_male <- data.frame()
measured_phenos <- c('standing_height_imaging_visit','head_diameter',
                     'trochanter_distance','shoulder_width',
                     'torso_length','humerus','femur','forearm','tibia',
                     'pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
                     'oblique_pelvic_inlet_length',
                     'subpubic_angle', 'iliac_isthmus_breadth',
                     'acetabular_diameter') 
covar <- c("SCORESUM", "weight", "household_income","E11", "F32", "F33", "I25", 'smoking', "sleep_duration" ) 

outcome_phenos = c('testosterone', 'igf_1', 'shbg', 'oestradiol')

for (outcome in outcome_phenos) {
  print(outcome)
  for (pheno in measured_phenos) {
    prs_df <- read.table(paste0("hip_male_400k_240321/prs_output_", 
                                pheno, 
                                "_pst_eff_a1_b0.5_phiauto_all.profile"), header = T)[, c("FID", "SCORESUM")]
    # filter imaged patients
    prs_df <- prs_df %>% filter(!FID %in% male_train_eids)
    
    # merge prs, outcome phenotypes, pcs
    tmp_df <- merge(prs_df, pcs_fid, by.x = "FID", by.y = "eid")
    tmp_df <- tmp_df %>% filter(sex == 0)
    tmp_df <- na.omit(tmp_df[c(outcome, covar)])
    
    pheno_scores <- tmp_df["SCORESUM"]
    pheno_z_scores <- scale(pheno_scores)
    pheno_z_scores <- as.data.frame(pheno_z_scores)
    tmp_df["SCORESUM"] <- pheno_z_scores["SCORESUM"]
    
    # Run linear regression
    formula <- as.formula(paste(outcome, " ~ ",  paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "gaussian")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))["SCORESUM", "Estimate"]
    p_value <- coef(summary(res))["SCORESUM", "Pr(>|t|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results_female data frame
    results_male <- rbind(results_male, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
  print(dim(tmp_df))
}

print(results_male)

# save
write.csv(results_male, 
          file = "../../key_results/disease_phenotypes_association_400k_male_prs_240322_skeletal_pelvis_blood_chem.csv", 
          row.names = FALSE)



##################################
############# Both ############# 
##################################

# logistic regression
# outcome_phenos <- c('M54', 'M16', 'M17', 'M23', 'had_same_sex_intercourse', 'no_pain', 'back_pain', 'neck_or_shoulder_pain', 'headache', 
#                     'knee_pain', 'stomach_or_abdominal_pain', 'hip_pain', 'facial_pain', 'pain_all_over_the_body', 'fracture_spine', 
#                     'fracture_other_bones', 'fracture_wrist', 'fracture_arm', 'fracture_leg', 'fracture_ankle', 'walking_pace')
outcome_phenos <- c('M54', 'M16', 'M17', 'M23',
                    'back_pain', 'knee_pain', 'hip_pain',
                    'walking_pace')

measured_phenos <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
                     'oblique_pelvic_inlet_length',
                     'subpubic_angle', 'iliac_isthmus_breadth',
                     'acetabular_diameter', 'pelvic_inlet_area',
                     'leg_divide_torso')
# measured_phenos <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width',
#                      'oblique_pelvic_inlet_length',
#                      'subpubic_angle', 'ear_left2ear_right', 'iliac_isthmus_breadth',
#                      'acetabular_diameter', 'trochanter_left2trochanter_right',
#                      'shoulder_width', 'arm_devide_torso', 
#                      'head_divide_inlet_width', 'shoulder_divide_inlet_width',
#                      'head_divide_oblique_inlet_length',
#                      'shoulder_divide_oblique_inlet_length', 'pelvic_inlet_area',
#                      'head_area_divide_pelvic_inlet_area',
#                      'shoulder_area_divide_pelvic_inlet_area',
#                      'arm_divide_leg', 'leg_divide_torso', 'bi_acetabular_width')

# covar <- c("SCORESUM", "weight", "household_income", "E11", "F32", "F33", "I25", 'smoking', "sleep_duration")
covar <- c("SCORESUM", "weight", "household_income", "E11", "F32", "F33", "I25", 'smoking', "sleep_duration" ) 


results_both <- data.frame()

for (outcome in outcome_phenos) {
  print(outcome)
  for (pheno in measured_phenos) {
    prs_df <- read.table(paste0("hip_both_400k_230913/prs_output_", 
                                pheno, 
                                "_pst_eff_a1_b0.5_phiauto_all.profile"), header = T)[, c("FID", "SCORESUM")]
    # filter imaged patients
    prs_df <- prs_df %>% filter(!FID %in% both_train_eids)
    
    # merge prs, outcome phenotypes, pcs
    tmp_df <- merge(prs_df, pcs_fid, by.x = "FID", by.y = "eid")
    
    tmp_df <- na.omit(tmp_df[c(outcome, covar)])
    
    pheno_scores <- tmp_df["SCORESUM"]
    pheno_z_scores <- scale(pheno_scores)
    pheno_z_scores <- as.data.frame(pheno_z_scores)
    tmp_df["SCORESUM"] <- pheno_z_scores["SCORESUM"]
    
    # Run logistic regression
    formula <- as.formula(paste(outcome, " ~ ",  paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "binomial")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))["SCORESUM", "Estimate"]
    p_value <- coef(summary(res))["SCORESUM", "Pr(>|z|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results_both data frame
    results_both <- rbind(results_both, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
  print(dim(tmp_df))
}

# linear regression
# outcome_phenos = c('patient_birth_weight', 'hip_circumference', 'igf_1', 'testosterone', 'shbg', 'oestradiol',
#                    'n_same_sex_sexual_partner', 'n_sexual_partner')
# outcome_phenos = c('acceleration')
# 
# for (outcome in outcome_phenos) {
#   print(outcome)
#   for (pheno in measured_phenos) {
#     prs_df <- read.table(paste0("hip_both_400k_230913/prs_output_",
#                                 pheno,
#                                 "_pst_eff_a1_b0.5_phiauto_all.profile"), header = T)[, c("FID", "SCORESUM")]
#     # filter imaged patients
#     prs_df <- prs_df %>% filter(!FID %in% both_train_eids)
# 
#     # merge prs, outcome phenotypes, pcs
#     tmp_df <- merge(prs_df, pcs_fid, by.x = "FID", by.y = "eid")
# 
#     tmp_df <- na.omit(tmp_df[c(outcome, covar)])
# 
#     pheno_scores <- tmp_df["SCORESUM"]
#     pheno_z_scores <- scale(pheno_scores)
#     pheno_z_scores <- as.data.frame(pheno_z_scores)
#     tmp_df["SCORESUM"] <- pheno_z_scores["SCORESUM"]
# 
#     # Run linear regression
#     formula <- as.formula(paste(outcome, " ~ ",  paste(covar, collapse = " + ")))
#     res <- glm(formula, data = tmp_df, family = "gaussian")
# 
#     # Save beta coefficient and p-value
#     beta <- coef(summary(res))["SCORESUM", "Estimate"]
#     p_value <- coef(summary(res))["SCORESUM", "Pr(>|t|)"]
#     neglog10p <- -log10(p_value)
#     odds_ratio <- exp(beta)
# 
#     # Append to results_both data frame
#     results_both <- rbind(results_both, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
#   }
#   print(dim(tmp_df))
# }

print(results_both)

# save
write.csv(results_both, 
          file = "../../key_results/disease_phenotypes_association_400k_both_prs_240222.csv", 
          row.names = FALSE)




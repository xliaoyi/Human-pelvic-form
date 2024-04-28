setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")

library(dplyr)

length_pheno <- c('pelvic_height',
                  'pelvic_width',
                  'pelvic_inlet_width',
                  'oblique_pelvic_inlet_length',
                  'ear_left2ear_right',
                  'trochanter_left2trochanter_right',
                  'iliac_isthmus_breadth',
                  'acetabular_diameter',
                  'shoulder_width')

standing_height <- read.csv('key_results/hip_select_pheno_cm_eid_flt.csv') %>% select(c("eid", "sex", "standing_height"))

##################################
############## Both ##############
##################################

both_df <- read.csv('UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_cm_230908.txt', sep = " ") %>% select(-FID) %>% rename(eid = "IID")
both_df <- merge(standing_height, both_df, by = 'eid')

both_df_residual <- both_df

for (i in length_pheno) {
  lm.fit <- lm(both_df[,i] ~ both_df$standing_height)
  both_df_residual[,i] <- residuals(lm.fit)
}

write.csv(both_df_residual, "key_results/hip_select_pheno_cm_residual.csv", row.names = FALSE)

# test for ratio of height
# t <- both_df
# for (i in length_pheno) {
#   t[,i] <- both_df[,i] / both_df$standing_height
# }
##################################
############## Male ##############
##################################

male_df <- read.csv('UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_male_cm_230908.txt', sep = " ") %>% select(-FID) %>% rename(eid = "IID")
male_df <- merge(standing_height, male_df, by = 'eid', )

male_df_residual <- male_df

for (i in length_pheno) {
  lm.fit <- lm(male_df[,i] ~ male_df$standing_height)
  male_df_residual[,i] <- residuals(lm.fit)
}

write.csv(male_df_residual, "key_results/hip_select_pheno_cm_male_residual.csv", row.names = FALSE)

##################################
############# Female ############# 
##################################
female_df <- read.csv('UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_female_cm_230908.txt', sep = " ") %>% select(-FID) %>% rename(eid = "IID")
female_df <- merge(standing_height, female_df, by = 'eid')

female_df_residual <- female_df

for (i in length_pheno) {
  lm.fit <- lm(female_df[,i] ~ female_df$standing_height)
  female_df_residual[,i] <- residuals(lm.fit)
}

write.csv(female_df_residual, "key_results/hip_select_pheno_cm_female_residual.csv", row.names = FALSE)


##################################
# Regress on both male and female #
##################################
male_df <- read.csv('UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_male_cm_230908.txt', sep = " ") %>% select(-FID) %>% rename(eid = "IID")
male_df <- merge(standing_height, male_df, by = 'eid')

female_df <- read.csv('UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_female_cm_230908.txt', sep = " ") %>% select(-FID) %>% rename(eid = "IID")
female_df <- merge(standing_height, female_df, by = 'eid')

df <- rbind(male_df, female_df)
df_residual <- df

for (i in length_pheno) {
  lm.fit <- lm(df[,i] ~ df$standing_height)
  df_residual[,i] <- residuals(lm.fit)
}

write.csv(df_residual, "key_results/hip_select_pheno_cm_residual_regress_on_female_and_male.csv", row.names = FALSE)


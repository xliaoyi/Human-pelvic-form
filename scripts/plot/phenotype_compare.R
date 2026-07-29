rm(list = ls())

library(ggplot2)
library(ggpubr)
library(tidyverse)

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape")

df <- read.csv("key_results/hip_select_pheno_both_flt_melt.csv")

df <- df %>% mutate(sex = ifelse(sex == 1, "Male", "Female"))

pdf("out_fig/phenotypes_compare_between_male_female.pdf",
    width = 12, 
    height = 10)

ggboxplot(df, x = "sex", y = "value",
                color = "sex", palette =c("steelblue", "indianred")) +
  facet_wrap(~phenotype, scales = "free") +
  labs(title="Comparison of Different Sexes on Different Phenotypes",
       x="", y="Phenotype value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(method = "t.test")

dev.off()



################# use residual ################# 
rm(list = ls())

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(reshape2)
library(RColorBrewer)

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape")

male_eid <- read.csv('key_results/hip_select_pheno_cm_eid_flt_z_flt_male.csv') %>% select(c("eid", 'sex', "standing_height"))
female_eid <- read.csv('key_results/hip_select_pheno_cm_eid_flt_z_flt_female.csv') %>% select(c("eid", 'sex', "standing_height"))
eid <- rbind(male_eid, female_eid)

male_df <- read.csv("UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_male_cm_230908.txt", sep = " ")
female_df <- read.csv("UKB_Imaging_Genetics/GWAS_INPUT_DATA/hip_select_pheno_gwas_female_cm_230908.txt", sep = " ")
df <- rbind(male_df, female_df) %>% rename(eid = "IID") %>% select(-FID)

df <- merge(eid, df, by = 'eid') %>% select(-eid)

length_pheno <- c('pelvic_height',
                  'pelvic_width',
                  'pelvic_inlet_width',
                  'oblique_pelvic_inlet_length',
                  'ear_left2ear_right',
                  'trochanter_left2trochanter_right',
                  'shoulder_width',
                  'iliac_isthmus_breadth',
                  'acetabular_diameter')

df_residual <- df

for (i in length_pheno) {
  lm.fit <- lm(df[,i] ~ df$standing_height)
  df_residual[,i] <- residuals(lm.fit)
}

df_residual <- df_residual %>% mutate(sex = ifelse(sex == 1, "Male", "Female"))

# melt
df_melted <- melt(df_residual, id.vars = "sex", variable.name = "pheno", value.name = "value")
# 
# measured_phenos <- c('pelvic_height', 'ear_left2ear_right', 'shoulder_width', 'arm_devide_torso', 'pelvic_width', 
#                      'pelvic_inlet_width', 'oblique_pelvic_inlet_length', 'subpubic_angle', 'iliac_isthmus_breadth', 'acetabular_diameter')
measured_phenos <- c('pelvic_height', 'pelvic_width', 
                     'pelvic_inlet_width', 'oblique_pelvic_inlet_length', 'subpubic_angle', 'iliac_isthmus_breadth', 'acetabular_diameter')

df_melted <- df_melted %>% filter(pheno %in% measured_phenos)
# write.csv(df_melted, "key_results/phenotype_residual_gender_compare.csv")

pdf("out_fig/phenotypes_compare_between_male_female_residuals.pdf",
    width = 10, 
    height = 10)


ggboxplot(df_melted, x = "sex", y = "value", 
          fill = "sex", palette =c("#E15759", "#1F78B4"),
          outlier.shape=4, outlier.alpha = 0.1) +
  facet_wrap(~pheno, scales = "free", nrow = 2) +
  labs(x="", y="Phenotype (residual/degree/ratio)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major.x = element_blank()) +
  guides(fill = guide_legend(title = NULL)) # +
  # stat_compare_means(method = "t.test")

dev.off()

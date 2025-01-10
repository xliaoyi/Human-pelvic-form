rm(list = ls())

library(ggplot2)
library(tidyverse)

setwd("/Users/louis/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")

fid.info <- read.delim("../UKB_xray_image_info/fids/fid_disease/fid_info_400k_20240222.csv", 
                       sep = ',')
pheno <- read.csv("key_results/hip_select_pheno_cm_eid_flt_z_flt_female.csv")

skeletal.pheno <- c('standing_height_imaging_visit', 'head_diameter', 
                    'trochanter_distance', 'shoulder_width',
                    'torso_length', 'humerus', 'femur', 'forearm', 'tibia')
pelvic.pheno <- c('pelvic_height', 'pelvic_width', 'pelvic_inlet_width', 
                  'oblique_pelvic_inlet_length', 'subpubic_angle', 
                  'iliac_isthmus_breadth', 'acetabular_diameter')
outcome.pheno <- c('testosterone', 'igf_1', 'shbg', 'oestradiol',
                   'periods_started_age', 'periods_stopped_age')
covar <- c("age_imaging_visit", "weight", "household_income", "n_live_birth") # standing_height_imaging_visit
all.pheno <- c("eid", skeletal.pheno, pelvic.pheno, outcome.pheno, covar)

###
### merge data and select subset and scale
###

pheno = pheno[,c("eid", skeletal.pheno, pelvic.pheno)]
fid.info = fid.info[,c("eid",outcome.pheno, covar)]

df = merge(pheno, fid.info, by = 'eid')

measured.pheno <- c(skeletal.pheno, pelvic.pheno)

df[measured.pheno] = scale(df[measured.pheno])

# calc delta_t = age at imaging - period stop age
df$year.post.period.stop = df$age_imaging_visit - df$periods_stopped_age
df$duration.able.give.birth = df$periods_stopped_age - df$periods_started_age

###
### regression
###

results <- data.frame()

for (outcome in c(outcome.pheno, "year.post.period.stop", "duration.able.give.birth")) {
  print(outcome)
  for (pheno in measured.pheno) {
    
    # for height pheno, shouldn't adjust for height
    if (pheno == "standing_height_imaging_visit") {
      covar <- c("weight", "household_income", "n_live_birth")
    }  else {
      covar <- c("weight", "household_income", "n_live_birth", 
                 "standing_height_imaging_visit")
    }
    
    tmp_df <- na.omit(df[c(outcome, pheno, covar)])
    print(dim(tmp_df))
    
    # Run linear regression
    formula <- as.formula(paste(outcome, " ~ ", pheno, 
                                " + ", paste(covar, collapse = " + ")))
    res <- glm(formula, data = tmp_df, family = "gaussian")
    
    # Save beta coefficient and p-value
    beta <- coef(summary(res))[pheno, "Estimate"]
    p_value <- coef(summary(res))[pheno, "Pr(>|t|)"]
    neglog10p <- -log10(p_value)
    odds_ratio <- exp(beta)
    
    # Append to results data frame
    results <- rbind(results, data.frame(outcome, pheno, beta, p_value, neglog10p, odds_ratio))
  }
}



###
### plot
###

results$FDR <- p.adjust(results$p_value, method = "fdr")
results <- results %>% mutate("SIGNIF" = ifelse((FDR < 0.05), TRUE, FALSE))

# order phenotypes

newfemale <- c("pelvic_inlet_width", "oblique_pelvic_inlet_length","subpubic_angle",
               "pelvic_height","pelvic_width", "iliac_isthmus_breadth" ,"acetabular_diameter",
               "trochanter_distance", "shoulder_width","head_diameter","humerus","tibia",
               "femur","forearm" ,"torso_length","standing_height_imaging_visit")

results$pheno <- factor(results$pheno, levels = rev(newfemale))

color_palette <- colorRampPalette(c("steelblue", 
                                    "steelblue1", 
                                    "white", 
                                    "darkorange", 
                                    "indianred"))(100)
p1 <- ggplot(results, aes(x = outcome, y = pheno)) +
  geom_point(aes(size = pmin(-log10(results$p_value), 3), 
                 fill = results$odds_ratio), shape = 22) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(0.8,0.9,1,1.1,1.2), 
                       labels=c(0.8,0.9,1,1.1,1.2), 
                       limits = c(0.8, 1.2), 
                       oob = scales::squish) +
  scale_size_continuous(range = c(0, 8)) +
  geom_point(data = subset(results, SIGNIF), 
             aes(x = outcome, y = pheno), 
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "exp(beta)", size = "-log10(p-value)") +
  ggtitle("Female phenotypic association") +
  theme(plot.title = element_text(hjust = 0.5)) # panel.grid.major = element_blank(), panel.grid.minor = element_blank()

ggsave(filename = "out_fig/disease_regression_female_phenotype_20319.pdf", 
       plot = p1,
       width = 6, height = 9, units = "in")

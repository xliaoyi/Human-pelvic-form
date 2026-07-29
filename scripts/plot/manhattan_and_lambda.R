rm(list = ls())
setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/LMM_OUTPUT/bolt_lmm_20230817_cm_both")
library(tidyverse)

all_pheno_sumstats <- data.frame()

files <- c("acetabular_diameter_lmm_combo.tab", "iliac_isthmus_breadth_lmm_combo.tab",
           # "arm_devide_torso_lmm_combo.tab", "shoulder_width_lmm_combo.tab", "ear_left2ear_right_lmm_combo.tab",
           "subpubic_angle_lmm_combo.tab", "oblique_pelvic_inlet_length_lmm_combo.tab", "pelvic_inlet_width_lmm_combo.tab",
           "pelvic_width_lmm_combo.tab", "pelvic_height_lmm_combo.tab")

# Loop through each file and rbind to the master dataframe
for (file in files) {
  print(paste0("start reading ", file))
  temp_df <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  temp_df <- temp_df[, !(names(temp_df) %in% c("CHISQ_BOLT_LMM", "P_BOLT_LMM"))]
  print(paste0("finished reading ", file))
  all_pheno_sumstats <- rbind(all_pheno_sumstats, temp_df)
  
  # Remove duplicates based on SNP column and keep the one with the lowest P_BOLT_LMM_INF
  print(paste0("start filtering ", file))
  all_pheno_sumstats <- all_pheno_sumstats %>%
    group_by(SNP) %>%
    arrange(P_BOLT_LMM_INF) %>%
    slice(1) %>%
    ungroup()
  print(paste0("finished filtering ", file))
  print(dim(all_pheno_sumstats))
}

write.csv(all_pheno_sumstats, file = "all_pheno_sumstats.txt")


# make manhattan plot
library(qqman)

df <- read.csv("all_pheno_sumstats.txt", sep = ",", header = TRUE)

png(filename="../../../out_fig/manhattan.png", height = 10000, width = 35000, res=1200)
par(mar=c(5.1, 6.1, 4.1, 2.1))
manhattan(df,
          chr = "CHR",
          bp = "BP",
          snp = "SNP",
          p = "P_BOLT_LMM_INF",
          col = c("#4E79A7", "#F28E2B"),
          suggestiveline = FALSE,
          genomewideline = -log10(5e-08))
dev.off()


# calculate lambda
# Initialize a data frame to store results
results_df <- data.frame(file = character(0), lambda = numeric(0))

for (file in files) {
  print(paste0("start reading ", file))
  
  # Read the data
  temp_df <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  print(paste0("finished reading ", file))
  
  # Convert P-values to chi-squared statistics
  # temp_df$chisq <- qchisq(1 - temp_df$P_BOLT_LMM_INF, df=1)
  
  # Calculate lambda
  lambda <- median(temp_df$CHISQ_BOLT_LMM_INF) / qchisq(0.5, df=1)
  
  print(paste0(file, " - ", lambda))
  # Append the results to the results data frame
  results_df <- rbind(results_df, data.frame(file = file, lambda = lambda))
}

# Print the results
print(results_df)

write.csv(results_df, "phenotype_lambdas.csv")
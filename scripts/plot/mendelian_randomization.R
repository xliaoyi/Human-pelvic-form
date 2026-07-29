rm(list = ls())

setwd("/Users/louis/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/LMM_OUTPUT/bolt_lmm_20230817_cm_female")

library(TwoSampleMR)
library(MRInstruments)
library(tidyverse)

all_pheno_sumstats <- data.frame()
files <- c("pelvic_inlet_area_lmm_combo.tab", "pelvic_inlet_width_lmm_combo.tab", "oblique_pelvic_inlet_length_lmm_combo.tab")
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
write.csv(all_pheno_sumstats, file = "all_inlet_sumstats.txt")

df <- read.csv("acetabular_diameter_lmm_combo.tab", sep = '\t')
# use snp clumped by plink
snp_list = read.table("../../PLINK_CLUMP_OUTPUT/lmm_female_clump_20240319_gwas/acetabular_diameter.clumped.ranges", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)$SNP
df <- df %>% filter(SNP %in% snp_list)

df$CHRPOS <- paste0(df$CHR, ":", df$BP)
write.csv(df, file = "/Users/louis/Downloads/acetabular_diameter_flt.txt")

out <- read.csv("/Users/louis/Downloads/reprogen_ANM_201K_170621.txt", sep = "\t")
out$CHRPOS <- paste0(out$CHR, ":", out$POS)
write.csv(out, file = "/Users/louis/Downloads/menarche_age.txt")

# t <- merge(df, out, by = 'CHRPOS')

exp <- read_exposure_data(
  # filename = "UKB_Imaging_Genetics/LMM_OUTPUT/bolt_lmm_20230817_cm_female/pelvic_inlet_area_lmm_combo.tab",
  # filename = "/Users/louis/Downloads/pelvic_inlet_width_flt.txt",
  # filename = "/Users/louis/Downloads/subpubic_angle_flt.txt",
  filename = "/Users/louis/Downloads/oblique_pelvic_inlet_length_flt.txt",
  # filename = "/Users/louis/Downloads/pelvic_width_flt.txt",
  # filename = "/Users/louis/Downloads/pelvic_height_flt.txt",
  # filename = "/Users/louis/Downloads/acetabular_diameter_flt.txt",
  sep = ",",
  snp_col = "CHRPOS",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_BOLT_LMM_INF",
  units_col = "",
  gene_col = "",
  samplesize_col = ""
)

exp_instruments <- subset(exp, pval.exposure < 1e-6)
# clumped_data <- clump_data(exp_instruments)
# instrument_snps <- clumped_data$SNP


df <- read.table("../../DOWNLOAD_GWAS_SUMSTATS/Maternal_BW_European_meta.NG2019.txt", sep = " ", header = T)

out <- df %>% filter(SNP %in% instrument_snps)
write.csv(out, file = "../../Maternal_BW_European_meta.NG2019.txt.test")

out <- read_outcome_data(
  filename = "/Users/louis/Downloads/menarche_age.txt",
  sep = ",",
  snp_col = "CHRPOS",
  beta_col = "Effect",
  se_col = "SE",
  effect_allele_col = "Effect_Allele",
  other_allele_col = "Other_Allele",
  eaf_col = "",
  pval_col = "Pval",
  units_col = "",
  gene_col = "",
  samplesize_col = "N"
)

# out.flt <- out %>% filter(SNP %in% instrument_snps)
out.flt <- out %>% filter(SNP %in% exp_instruments$SNP)

dat <- harmonise_data(
  exposure_dat = exp, 
  outcome_dat = out.flt
)

res <- mr(dat)
mr_scatter_plot(res, dat)

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]] + theme_bw()



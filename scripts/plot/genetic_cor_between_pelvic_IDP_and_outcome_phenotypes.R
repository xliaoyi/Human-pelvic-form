rm(list = ls())

library(ggplot2)
library(dplyr)
library(pheatmap)
library(Hmisc)

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")

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
  "bi_acetabular_width" = "Biacetabular width",
  "BW_devide_pelvic_inlet_width" = "Baby birth weight:pelvic inlet width",
  "BW_devide_oblique_pelvic_inlet_length" = "Baby birth weight:oblique pelvic inlet length",
  "BW_devide_pelvic_inlet_area" = "Baby birth weight:pelvic inlet area"
)
rename_map_outcome <- list(
  "androgen_phenotypes.female.CBAT" = "Female CBAT",
  "androgen_phenotypes.female.SHBG" = "Female SHBG",
  "androgen_phenotypes.female.Testosterone" = "Female testosterone",
  "androgen_phenotypes.male.CBAT" = "Male CBAT",
  "androgen_phenotypes.male.SHBG" = "Male SHBG",
  "androgen_phenotypes.male.Testosterone" = "Male testosterone",
  "EGG_Maternal_GWAS_META_gestational_duration" = "Gestational duration (maternal)",
  "EGG_Maternal_GWAS_META_postterm_delivery" = "Postterm delivery (maternal)",
  "EGG_Maternal_GWAS_META_preterm_delivery" = "Preterm delivery (maternal)",
  "Fetal_BW_European_meta.NG2019" = "Birth weight (European Fetal)",
  "Fetal_BW_transethnic_meta.NG2019" = "Birth weight (transethnic Fetal)",
  "Fetal_gest_duration_NComms2019" = "Gestational duration (fetal)",
  "Fetal_postterm_birth_NComms2019" = "Postterm delivery (fetal)",
  "Fetal_preterm_birth_NComms2019" = "Preterm delivery (fetal)",
  "Maternal_BW_European_meta.NG2019" = "Birth weight (European maternal)",
  "Maternal_BW_transethnic_meta.NG2019" = "Birth weight (transethnic maternal)",
  "finngen_R9_E4_PCOS" = "Polycystic ovarian syndrome",
  "finngen_R9_N14_ENDOMETRIOSIS" = "Endometriosis",
  "finngen_R9_N14_FEMALEGENINF" = "Inflammatory diseases of female pelvic organs",
  "finngen_R9_O15_DELIV_CAESAR" = "Delivery by caesarean section",
  "finngen_R9_O15_HAEMORRH_EARLY_PREG" = "Haemorrhage in early pregnancy",
  "finngen_R9_O15_LABOUR_FETAL_STRESS" = "Labour and delivery complicated by fetal stress",
  "finngen_R9_O15_LABOUR_LONG" = "Long labour",
  "finngen_R9_O15_LABOUR_MALPOS" = "Obstructed labour due to malposition and malpresentation of fetus",
  "finngen_R9_O15_LABOUR_PELVIC_ABNORM" = "Obstructed labour due to maternal pelvic abnormality",
  "finngen_R9_O15_MEMBR_PREMAT_RUPT" = "Premature rupture of membranes",
  "finngen_R9_O15_PLAC_DISORD" = "Placental disorders",
  "finngen_R9_O15_POLYHYDR" = "Polyhydramnios",
  "finngen_R9_O15_POOR_FETGRO" = "Poor fetal growth",
  "finngen_R9_O15_PREECLAMPS" = "Pre-eclampsia",
  "finngen_R9_O15_PREG_ABORT" = "Pregnancy with abortive outcome",
  "finngen_R9_O15_PREG_ECTOP" = "Ectopic pregnancy",
  "finngen_R9_O15_PREG_PROLONGED" = "Prolonged pregnancy",
  "finngen_R9_O15_PRETERM" = "Preterm labour and delivery",
  "neale_Oestradiol_pmol_L" = "Oestradiol",
  "neale_usual_walking_pace" = "Walking pace"
)

# newoutcome <- c("Female CBAT", "Female SHBG", "Female testosterone", "Male CBAT", "Male SHBG", "Male testosterone",
#                 "Gestational duration (maternal)", "Postterm delivery (maternal)",  "Preterm delivery (maternal)",
#                 "Birth weight (European Fetal)", "Birth weight (transethnic Fetal)",
#                 "Gestational duration (fetal)", "Postterm delivery (fetal)", "Preterm delivery (fetal)",
#                 "Birth weight (European maternal)", "Birth weight (transethnic maternal)",
#                 "Polycystic ovarian syndrome", "Endometriosis", "Inflammatory diseases of female pelvic organs",
#                 "Delivery by caesarean section", "Haemorrhage in early pregnancy", "Labour and delivery complicated by fetal stress",
#                 "Long labour", "Obstructed labour due to malposition and malpresentation of fetus", "Obstructed labour due to maternal pelvic abnormality",
#                 "Premature rupture of membranes", "Placental disorders", "Polyhydramnios", "Poor fetal growth", "Pre-eclampsia", "Pregnancy with abortive outcome",
#                 "Ectopic pregnancy", "Prolonged pregnancy", "Preterm labour and delivery", "Oestradiol", "Walking pace")
newoutcome <- c("Obstructed labour due to malposition and malpresentation of fetus", "Obstructed labour due to maternal pelvic abnormality")
newpheno <- c("Acetabular diameter", "Iliac isthmus breadth", # , "Biacetabular width"
              "Pelvic height", "Pelvic width", "Subpubic angle",
              "Pelvic inlet width", "Oblique pelvic inlet length", # , "Pelvic inlet area"
              "Baby birth weight:pelvic inlet width", "Baby birth weight:oblique pelvic inlet length") # , "Baby birth weight:pelvic inlet area"

### female

df <- read.csv("key_results/genetic_correlation/rg_20230909/rg_hip_pelvic_vs_outcome_female.csv")

# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno_pelvic = ifelse(pheno_pelvic %in% names(rename_map_pheno), rename_map_pheno[pheno_pelvic], pheno_pelvic),
    pheno_outcome = ifelse(pheno_outcome %in% names(rename_map_outcome), rename_map_outcome[pheno_outcome], pheno_outcome)
  )

# filteration
df <- df %>% filter(pheno_pelvic %in% newpheno) %>% filter(pheno_outcome %in% newoutcome)

df$pheno_pelvic <- factor(df$pheno_pelvic, levels = c(newpheno))
df$pheno_outcome <- factor(df$pheno_outcome, levels = c(newoutcome))

# Add column indicating whether Bonferroni-corrected p-value is significant
df$FDR <- p.adjust(df$p_value, method = "fdr") < 0.05

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = pheno_outcome, y = pheno_pelvic)) +
  geom_point(data = df, 
             aes(size = pmin(neglog10p, 3), fill = cor), shape = 22) +
  # scale_size_continuous(range = c(3, 7)) +
  scale_size_continuous(range = c(0, 8), breaks = c(1, 2, 3), labels = c(1, 2, 3)) +
  # geom_text(data = subset(df, neglog10p >= -log10(0.5)), 
  #           aes(label = paste0(round(cor, 2), "\n", round(neglog10p, 2)), size = 0.5), 
  #           show.legend = FALSE, vjust = 1.6) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(-1, -0.5, 0, 0.5, 1), 
                       labels=c(-1, -0.5, 0, 0.5, 1), 
                       limits = c(-1, 1), 
                       oob = scales::squish) +
  geom_point(data = subset(df, FDR),
             aes(x = pheno_outcome, y = pheno_pelvic),
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Correlation", size = "-log10(p-value)") +
  ggtitle("Genetic correlation") +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = "out_fig/rg_pelvic_pheno_and_outcome_pheno_female_240204.pdf", plot = p,
       width = 4.9, height = 7.1, units = "in")


### male


df <- read.csv("key_results/genetic_correlation/rg_20230909/rg_hip_pelvic_vs_outcome_male.csv")

# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno_pelvic = ifelse(pheno_pelvic %in% names(rename_map_pheno), rename_map_pheno[pheno_pelvic], pheno_pelvic),
    pheno_outcome = ifelse(pheno_outcome %in% names(rename_map_outcome), rename_map_outcome[pheno_outcome], pheno_outcome)
  )

# filteration
df <- df %>% filter(pheno_pelvic %in% newpheno) %>% filter(pheno_outcome %in% newoutcome)

df$pheno_pelvic <- factor(df$pheno_pelvic, levels = c(newpheno))
df$pheno_outcome <- factor(df$pheno_outcome, levels = c(newoutcome))

# Add column indicating whether Bonferroni-corrected p-value is significant
df$FDR <- p.adjust(df$p_value, method = "fdr") < 0.05

# keep p < 0.05
df$p_value <- ifelse(df$neglog10p < -log10(0.5), 0, df$neglog10p)

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = pheno_outcome, y = pheno_pelvic)) +
  geom_point(data = subset(df, neglog10p >= -log10(0.5)), 
             aes(size = pmin(neglog10p, 6), fill = cor), shape = 21) +
  scale_size_continuous(range = c(0, 8)) +
  # geom_text(data = subset(df, neglog10p >= -log10(0.5)), 
  #           aes(label = paste0(round(cor, 2), "\n", round(neglog10p, 2)), size = 0.5), 
  #           show.legend = FALSE, vjust = 1.6) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(-1, -0.5, 0, 0.5, 1), 
                       labels=c(-1, -0.5, 0, 0.5, 1), 
                       limits = c(-1, 1), 
                       oob = scales::squish) +
  geom_point(data = subset(df, FDR),
             aes(x = pheno_outcome, y = pheno_pelvic),
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Correlation", size = "-log10(p-value)") +
  ggtitle("Genetic correlation between male pelvic phenotypes and outcome phenotypes")

ggsave(filename = "out_fig/rg_pelvic_pheno_and_outcome_pheno_male.pdf", plot = p,
       width = 14, height = 9, units = "in")


### both


df <- read.csv("key_results/genetic_correlation/rg_20230909/rg_hip_pelvic_vs_outcome_both.csv")

# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno_pelvic = ifelse(pheno_pelvic %in% names(rename_map_pheno), rename_map_pheno[pheno_pelvic], pheno_pelvic),
    pheno_outcome = ifelse(pheno_outcome %in% names(rename_map_outcome), rename_map_outcome[pheno_outcome], pheno_outcome)
  )

# filteration
df <- df %>% filter(pheno_pelvic %in% newpheno) %>% filter(pheno_outcome %in% newoutcome)

df$pheno_pelvic <- factor(df$pheno_pelvic, levels = c(newpheno))
df$pheno_outcome <- factor(df$pheno_outcome, levels = c(newoutcome))

# Add column indicating whether Bonferroni-corrected p-value is significant
df$FDR <- p.adjust(df$p_value, method = "fdr") < 0.05

# keep p < 0.05
df$p_value <- ifelse(df$neglog10p < -log10(0.5), 0, df$neglog10p)

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = pheno_outcome, y = pheno_pelvic)) +
  geom_point(data = subset(df, neglog10p >= -log10(0.5)), 
             aes(size = pmin(neglog10p, 6), fill = cor), shape = 21) +
  scale_size_continuous(range = c(0, 8)) +
  # geom_text(data = subset(df, neglog10p >= -log10(0.5)), 
  #           aes(label = paste0(round(cor, 2), "\n", round(neglog10p, 2)), size = 0.5), 
  #           show.legend = FALSE, vjust = 1.6) +
  scale_fill_gradientn(colors = color_palette, 
                       breaks=c(-1, -0.5, 0, 0.5, 1), 
                       labels=c(-1, -0.5, 0, 0.5, 1), 
                       limits = c(-1, 1), 
                       oob = scales::squish) +
  geom_point(data = subset(df, FDR),
             aes(x = pheno_outcome, y = pheno_pelvic),
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Correlation", size = "-log10(p-value)") +
  ggtitle("Genetic correlation between both sex pelvic phenotypes and outcome phenotypes")

ggsave(filename = "out_fig/rg_pelvic_pheno_and_outcome_pheno_both.pdf", plot = p,
       width = 14, height = 9, units = "in")


### difference between male and female
df <- read.csv("key_results/genetic_correlation/rg_20230909/rg_hip_pelvic_vs_outcome_male_and_female.csv")
df$z_difference <- (df$female_cor - df$male_cor) / sqrt(df$female_se^2 + df$male_se^2)
df$p_value_difference <- 2 * (1 - pnorm(abs(df$z_difference)))

# keep significant result
df_sorted <- df[order(df$p_value_difference), ]
df.flt <- df_sorted %>% filter(df_sorted$p_value_difference < 0.05)
df.flt_sorted <- df.flt[order(df.flt$female_cor, decreasing = TRUE), ]

write.csv(df.flt_sorted, "key_results/genetic_correlation/rg_20230909/male_vs_female.csv", row.names = FALSE)


rm(list=ls())
setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape")

library(ggplot2)
library(dplyr)

gcta_h2g_df <- read.csv("key_results/h2g/h2g_all_20230918.csv")
gcta_h2g_df <- gcta_h2g_df %>% rename(gcta.h2g = var, gcta.h2g.se = se)

ldsc_h2g_df <- read.csv("key_results/h2g/h2g_ldsc_all_20230909.csv")
ldsc_h2g_df <- ldsc_h2g_df %>% rename(ldsc.h2g = h2g, ldsc.h2g.se = se)

# # load in heritability and genetic correlation estimates
# df_mvsf <- read.csv("key_results/genetic_correlation/rg_20230824/rg_hip_female_vs_male.csv")
# df_mvsf <- df_mvsf %>% filter(pheno_f == pheno_m)
# 
# # get both df
# both_df <- read.csv("key_results/genetic_correlation/rg_20230824/rg_hip_both.csv")
# 
# both_df_1 <- both_df %>% dplyr::select(c("pheno1", "pheno1_h2g", "pheno1_h2g_se")) %>%
#   rename(pheno = pheno1, pheno_h2g = pheno1_h2g, pheno_h2g_se = pheno1_h2g_se)
# both_df_2 <- both_df %>% dplyr::select(c("pheno2", "pheno2_h2g", "pheno2_h2g_se")) %>%
#   rename(pheno = pheno2, pheno_h2g = pheno2_h2g, pheno_h2g_se = pheno2_h2g_se)
# both_df_melt <- rbind(both_df_1, both_df_2) %>% distinct() %>% mutate(sex = "both")
# 
# # get male and female df
# male_df <- df_mvsf %>%
#   dplyr::select(c("pheno_m", "pheno_m_h2g", "pheno_m_h2g_se")) %>%
#   rename(pheno = pheno_m, pheno_h2g = pheno_m_h2g, pheno_h2g_se = pheno_m_h2g_se) %>%
#   mutate(sex = "male")
# female_df <- df_mvsf %>%
#   dplyr::select(c("pheno_f", "pheno_f_h2g", "pheno_f_h2g_se")) %>%
#   rename(pheno = pheno_f, pheno_h2g = pheno_f_h2g, pheno_h2g_se = pheno_f_h2g_se) %>%
#   mutate(sex = "female")
# 
# # merge male and female df
# mf_df <- rbind(male_df, female_df) %>% distinct()
# 
# # merge three dfs
# df_final <- rbind(mf_df, both_df_melt)
# 
# ldsc_h2g_df <- df_final %>% rename(ldsc.h2g = pheno_h2g, ldsc.h2g.se = pheno_h2g_se)

# merge two h2g tables
df <- merge(gcta_h2g_df, ldsc_h2g_df, by = c("pheno", "sex"))

# filter phenotypes
select_pheno <- c("pelvic_height", "pelvic_width", "pelvic_inlet_width", "subpubic_angle",
                  "oblique_pelvic_inlet_length", "iliac_isthmus_breadth", "acetabular_diameter")
df <- df %>% filter(pheno %in% select_pheno)


# rename phenotypes
rename_map_pheno <- list(
  "acetabular_diameter" = "Acetabular diameter",
  "ear_left2ear_right" = "Head width",
  "iliac_isthmus_breadth" = "Iliac isthmus breadth",
  "oblique_pelvic_inlet_length" = "Oblique pelvic inlet length",
  "pelvic_height" = "Pelvic height",
  "pelvic_inlet_width" = "Pelvic inlet width",
  "pelvic_width" = "Pelvic width",
  "shoulder_width" = "Shoulder width",
  "subpubic_angle" = "Subpubic angle",
  "pelvic_inlet_area" = "Pelvic inlet area",
  "bi_acetabular_width" = "Biacetabular width"
)
# Apply renaming to the df dataframe
df <- df %>%
  mutate(
    pheno = ifelse(pheno %in% names(rename_map_pheno), rename_map_pheno[pheno], pheno)
  )
df$pheno <- factor(df$pheno, levels = unique(df$pheno))

custom_colors <- c("Pelvic height" = "#76B7B2", "Pelvic width" = "#B07AA1", "Pelvic inlet width" = "#4E79A7",
                   "Oblique pelvic inlet length" = "#FF9DA7", "Iliac isthmus breadth" = "#59A14F", 
                   "Acetabular diameter" = "#E15759", "Subpubic angle" = "#EDC948")
# Create the scatter plot with error bars
p <- ggplot(df %>% filter(sex == "both"), aes(x = gcta.h2g, y = ldsc.h2g)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", size = 0.5) +
  # geom_abline(slope = 1, intercept = 0, size = 0.5, ) + 
  geom_point(aes(color = pheno), size = 3) +
  geom_errorbar(aes(ymin = ldsc.h2g - ldsc.h2g.se, ymax = ldsc.h2g + ldsc.h2g.se, color = pheno, width = 0), size = 1) +
  geom_errorbarh(aes(xmin = gcta.h2g - gcta.h2g.se, xmax = gcta.h2g + gcta.h2g.se, color = pheno, height = 0), size = 1) +
  scale_color_manual(values = custom_colors) +
  xlab(" Heritability estimated in GCTA") +
  ylab("Heritability estimated in LDSC") +
  # xlim(c(0.2, 0.54)) +
  # ylim(c(0.24, 0.6)) +
  theme_minimal() +
  theme(legend.position = c(1, 1), legend.justification = c(1.8, 1), 
        legend.title = element_blank()) 
ggsave(filename = "out_fig/h2g_gcta_vs_ldsc.pdf", plot = p,
       width = 4, height = 4, units = "in")


### from Arbel
rm(list = ls())
library("dplyr")
library("tidyr")
library("ggplot2")
library("grid")
library("gridExtra")
library("ggpubr")
library("reshape2")


setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape")

# h2g from LDSC
#
# # load in heritability and genetic correlation estimates
# df_mvsf <- read.csv("key_results/genetic_correlation/rg_20230701/rg_hip_female_vs_male.csv")
# df_mvsf <- df_mvsf %>% filter(pheno_f == pheno_m)
# 
# # get both df
# both_df <- read.csv("key_results/genetic_correlation/rg_20230701/rg_hip_both.csv")
# 
# both_df_1 <- both_df %>% select(c("pheno1", "pheno1_h2g", "pheno1_h2g_se")) %>% 
#   rename(pheno = pheno1, pheno_h2g = pheno1_h2g, pheno_h2g_se = pheno1_h2g_se)
# both_df_2 <- both_df %>% select(c("pheno2", "pheno2_h2g", "pheno2_h2g_se")) %>% 
#   rename(pheno = pheno2, pheno_h2g = pheno2_h2g, pheno_h2g_se = pheno2_h2g_se)
# both_df_melt <- rbind(both_df_1, both_df_2) %>% distinct() %>% mutate(sex = "both")
# 
# # get male and female df
# male_df <- df_mvsf %>% 
#   select(c("pheno_m", "pheno_m_h2g", "pheno_m_h2g_se")) %>% 
#   rename(pheno = pheno_m, pheno_h2g = pheno_m_h2g, pheno_h2g_se = pheno_m_h2g_se) %>% 
#   mutate(sex = "male")
# female_df <- df_mvsf %>% 
#   select(c("pheno_f", "pheno_f_h2g", "pheno_f_h2g_se")) %>% 
#   rename(pheno = pheno_f, pheno_h2g = pheno_f_h2g, pheno_h2g_se = pheno_f_h2g_se) %>% 
#   mutate(sex = "female")
# 
# # merge male and female df
# mf_df <- rbind(male_df, female_df) %>% distinct()
# 
# # merge three dfs
# df_final <- rbind(mf_df, both_df_melt)

# h2g from GCTA
h2g_df <- read.csv("key_results/h2g/h2g_all_20230918.csv")
h2g_df <- h2g_df %>% rename(pheno_h2g = var, pheno_h2g_se = se)

# add genetic correlation

# load in heritability and genetic correlation estimates
df_mvsf <- read.csv("key_results/genetic_correlation/rg_20230909/rg_hip_female_vs_male.csv")
df_mvsf <- df_mvsf %>% filter(pheno_f == pheno_m)
cor_df <- df_mvsf[, c("pheno_f", "cor", "se")]
df <- merge(h2g_df, cor_df, by.x = "pheno", by.y = "pheno_f")


# create dataframe of both sex heritabilities for later calculation
both_h2 <- rep(df[df$sex == "both", "pheno_h2g"], each = 3)

# calculate relative heritability and format df
df <- df %>%
  mutate(relative_h2 = pheno_h2g / both_h2) %>%
  mutate(relative_h2_se = ((pheno_h2g + pheno_h2g_se) / both_h2) - relative_h2 ) %>%
  mutate(sex = factor(sex, levels = c("female", "both", "male"))) %>%
  arrange(cor, pheno, sex) %>%
  mutate(pheno = factor(pheno, levels = unique(pheno)))

# write to table
write.table(df, file = "key_results/relative_h2.txt", quote=FALSE, sep="\t", row.names=FALSE) 


# create PDF
pdf(file = "out_fig/r2_by_h2.pdf", width = 9, height = 10)

# plot relative heritabilities for each phenotype
# p1 <- ggplot(df[df$sex != 'both',], aes(x = relative_h2, y = pheno, col = sex)) +
#   geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
#   geom_point(size = 3, position = position_dodge(width = 0.5), color = "black") +
#   geom_errorbarh(aes(y = pheno, xmin = relative_h2 - relative_h2_se, xmax = relative_h2 + relative_h2_se), 
#                  height = 0, position = position_dodge(width = 0.5), size = 1) +
#   scale_x_continuous(breaks = c(0.9,1,1.1,1.2,1.3), limits=c(0.85, 1.35), expand=c(0,0)) +
#   theme_bw() +
#   theme(legend.position = c(1, 1), legend.justification = c(5.9, 1.1), 
#         legend.title = element_blank(), axis.title.y = element_blank(),
#         axis.text.y = element_text(hjust = 0)) +
#   scale_color_manual(labels = c("Female", "Male"), values = c("#E15759", "#1F78B4")) +
#   labs(x = "Heritability Relative to Heritability of Both-sex Sample") 

p1 <- ggplot(df[df$sex != 'both',], aes(x = relative_h2, y = pheno, col = sex, group = sex)) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(y = pheno, xmin = relative_h2 - relative_h2_se, xmax = relative_h2 + relative_h2_se), 
                 height = 0, position = position_dodge(width = 0.7), size = 3) +
  scale_x_continuous(breaks = c(0.9,1,1.1,1.2,1.3), limits=c(0.85, 1.35), expand=c(0,0)) +
  theme_bw() +
  theme(legend.position = c(1, 1), legend.justification = c(5.9, 1.1), 
        legend.title = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = 0)) +
  scale_color_manual(labels = c("Female", "Male"), values = c("#E15759", "#1F78B4")) +
  geom_point(size = 2.5, position = position_dodge(width = 0.7), color = 'black') +
  labs(x = "Heritability Relative to Heritability of Both-sex Sample")

# get correlation and adjust positions for plot

measured_phenos <- c('pelvic_height', 'pelvic_width', 
                     'pelvic_inlet_width', 'oblique_pelvic_inlet_length', 'subpubic_angle', 'iliac_isthmus_breadth', 'acetabular_diameter')
select_pheno_df <- df %>% filter(pheno %in% measured_phenos)
corr_df <- select_pheno_df %>%
  select(c(1, 5)) %>%
  distinct() %>%
  mutate(pheno_point = seq(-0.022, 1.222, by = 1.244 / (length(pheno) - 1)))

# correlations
p2 <- ggplot(corr_df, aes(x = 0, y = cor)) +
  geom_segment(aes(x = 0, y = -0.05, xend = 0, yend = 1.25),
               arrow = arrow(length = unit(0.3, "cm"), end = "both"), size = 1, color = "black") +
  geom_segment(aes(x = 0, y = cor, xend = 1, yend = pheno_point), alpha = 0.3) +
  geom_point(shape = 1, size = 3) +
  geom_point(aes(x = 1, y = pheno_point), size = 0.5) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2), limits = c(-0.05, 1.25), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.1, 1.01)) +
  theme(plot.margin = margin(10,0,19,5), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_blank(), axis.title.y = element_text(size = 17), axis.text.y = element_text(size = 15)) +
  labs(x = "", y = "Genetic correlation between males and females")
# ggsave("out_fig/genetic_correlation.pdf", p2, height = 10, width = 1.5, units = "in")

# place plots next to each other
lay <- rbind( c(1,2,2,2,2,2,2))
p <- grid.arrange(p2, p1, ncol = 2, layout_matrix = lay,
                  top = textGrob("", gp = gpar(fontsize = 16)))

p
dev.off()

## STATS
# relative heritability diff by correlation
stat_df <- df[df$sex != "both", c(1,3,6,8,9)]
f <- stat_df[stat_df$Sex == "female", ]
m <- stat_df[stat_df$Sex == "male", ]
stat_df <- data.frame(Code = f$Code, Correlation = f$Correlation, h2_diff = abs(f$relative_h2 - m$relative_h2))
# heritability difference
model <- cor.test(stat_df$Correlation, stat_df$h2_diff)
print(model)




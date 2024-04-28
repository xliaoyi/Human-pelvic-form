rm(list = ls())

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")

library(tidyverse)
library(patchwork)
library(ggthemes)


############################################# 
### for two visits correlation plot
############################################# 



visit_1_df <- read.csv("key_results/visit_1_df.csv")
visit_2_df <- read.csv("key_results/visit_2_df.csv")

length.pheno <- c("pelvic_height", 
                  "pelvic_width",
                  "pelvic_inlet_width",
                  "oblique_pelvic_inlet_length",
                  "iliac_isthmus_breadth",
                  # "bi_acetabular_width",
                  "acetabular_diameter")

visit_1_df <- visit_1_df %>% 
  select(c("eid", length.pheno)) %>% 
  pivot_longer(cols = length.pheno,
               names_to = "Phenotypes",
               values_to = "Length_1")
visit_2_df <- visit_2_df %>% 
  select(c("eid", length.pheno)) %>% 
  pivot_longer(cols = length.pheno,
               names_to = "Phenotypes",
               values_to = "Length_2")
two_visit_df <- merge(visit_1_df, visit_2_df, by = c("eid", "Phenotypes"))
# write.csv(two_visit_df, "key_results/two_visit_df_avg_left_right.csv")


cor_value <- cor(two_visit_df$Length_1, 
                 two_visit_df$Length_2, 
                 method = "pearson")


custom_colors <- c("pelvic_height" = "#76B7B2", "pelvic_width" = "#B07AA1", "pelvic_inlet_width" = "#4E79A7",
                   "oblique_pelvic_inlet_length" = "#FF9DA7", "iliac_isthmus_breadth" = "#59A14F", 
                   "acetabular_diameter" = "#E15759") # "bi_acetabular_width" = "#F28E2B", 

p <- ggplot(two_visit_df, 
            aes(x=Length_1, 
                y=Length_2,
                color = Phenotypes)) +
  geom_point(alpha = 1, size = 3, shape = 1) +
  xlim(c(1, 38)) + ylim(c(1, 38)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color = 'black', size = 0.3) +
  labs(title="",
       x="First visit measurements (cm)",
       y="Second visit measurements (cm)") +
  annotate("text", x = Inf, y = -Inf, 
           label = sprintf("italic(R^2) == %.3f", cor_value**2), 
           hjust = 1.2, vjust = -2, 
           size = 6, color = "black", parse = TRUE) +
  # xlim(c(2.5, 8)) + ylim(c(2.5, 8)) +
  scale_color_manual(values = custom_colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_classic(base_size = 16) +
  guides(color = guide_legend(override.aes=list(shape = 15, size = 5)))

pdf("out_fig/two_visits_correlation.pdf",
    width = 7.5,
    height = 5)
p
dev.off()

############################################# 
### for left and right correlation plot
############################################# 

df <- read.csv("key_results/hip_pheno_23_cm.csv")

lr_cor_df <- df %>% dplyr::select(c('image_id',
                                    'sciatic_notch_left2inferior_iliac_spine_left',
                                    'sciatic_notch_right2inferior_iliac_spine_right',
                                    'iliopubic_eminence_left2acetabular_inferior_left',
                                    'iliopubic_eminence_right2acetabular_inferior_right',
                                    'acetabular_inclination_left',
                                    'acetabular_inclination_right'))
# Compute z-scores for each column excluding 'image_id'
z_scores <- as.data.frame(lapply(lr_cor_df[-1], function(x) scale(x)))

# Determine rows to keep based on z-scores
rows_to_keep <- apply(z_scores, 1, function(row) all(abs(row) <= 4))

# Filter rows based on z-scores
filtered_df <- lr_cor_df[rows_to_keep, ]


cor_value <- cor(filtered_df$sciatic_notch_left2inferior_iliac_spine_left,
                 filtered_df$sciatic_notch_right2inferior_iliac_spine_right,
                 method = "pearson")
p1 <- ggplot(filtered_df,
             aes(x=sciatic_notch_left2inferior_iliac_spine_left,
                 y=sciatic_notch_right2inferior_iliac_spine_right)) +
  geom_point(color="#59A14F", alpha = 1, size = 3, shape = 1) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  labs(title="Iliac isthmus breadth",
       x="Right measurements (cm)",
       y="Left measurements (cm)") +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("italic(r) == %.3f", cor_value),
           hjust = 1, vjust = -1.2, 
           size = 6, color = "black", parse = TRUE) +
  xlim(c(2.5, 6.7)) + ylim(c(2.5, 6.7)) +
  theme_classic(base_size = 16)

# Acetabular diameter

cor_value <- cor(filtered_df$iliopubic_eminence_left2acetabular_inferior_left,
                 filtered_df$iliopubic_eminence_right2acetabular_inferior_right,
                 method = "pearson")

p2 <- ggplot(filtered_df,
             aes(x=iliopubic_eminence_left2acetabular_inferior_left,
                 y=iliopubic_eminence_right2acetabular_inferior_right)) +
  geom_point(color="#E15759", alpha = 1, size = 3, shape = 1) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  labs(title="Acetabular diameter",
       x="Right measurements (cm)",
       y="") +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("italic(r) == %.3f", cor_value),
           hjust = 1, vjust = -1.2, 
           size = 6, color = "black", parse = TRUE) +
  xlim(c(4, 8)) + ylim(c(4, 8)) +
  theme_classic(base_size = 16)

# Acetabular inclination

# cor_value <- cor(filtered_df$acetabular_inclination_left,
#                  filtered_df$acetabular_inclination_right,
#                  method = "pearson")
# 
# p3 <- ggplot(filtered_df,
#              aes(x=acetabular_inclination_left,
#                  y=acetabular_inclination_right)) +
#   geom_point(color="#8CD17D", alpha = 0.3, size = 3) +
#   geom_abline(intercept=0, slope=1, linetype="dashed") +
#   labs(title="Acetabular inclination",
#        x="Left measurements (degree)",
#        y="Right measurements (degree)") +
#   annotate("text", x = Inf, y = -Inf,
#            label = sprintf("italic(r) == %.3f", cor_value),
#            hjust = 1.2, vjust = -2,
#            size = 5, color = "black", parse = TRUE) +
#   xlim(c(39, 81)) + ylim(c(39, 81)) +
#   theme_classic(base_size = 12)




############################################# 
### for two visits discrepancy correlation plot
############################################# 

visit_1_df <- read.csv("key_results/visit_1_df_w_lr.csv")
visit_2_df <- read.csv("key_results/visit_2_df_w_lr.csv")

visit_1_df <- visit_1_df %>% 
  mutate(iliac_isthmus_breadth_diff_1 = abs(sciatic_notch_left2inferior_iliac_spine_left - sciatic_notch_right2inferior_iliac_spine_right)) %>% 
  mutate(acetabular_diameter_diff_1 = abs(iliopubic_eminence_left2acetabular_inferior_left - iliopubic_eminence_right2acetabular_inferior_right)) %>% 
  mutate(acetabular_inclination_diff_1 = abs(acetabular_inclination_left - acetabular_inclination_right)) %>% 
  mutate(iliac_isthmus_breadth_mean = (sciatic_notch_left2inferior_iliac_spine_left + sciatic_notch_right2inferior_iliac_spine_right)/2) %>% 
  mutate(acetabular_diameter_mean = (iliopubic_eminence_left2acetabular_inferior_left + iliopubic_eminence_right2acetabular_inferior_right)/2) %>% 
  mutate(acetabular_inclination_mean = (acetabular_inclination_left + acetabular_inclination_right)/2) %>% 
  mutate(iliac_isthmus_breadth_proportion = iliac_isthmus_breadth_diff_1/iliac_isthmus_breadth_mean) %>% 
  mutate(acetabular_diameter_proportion = acetabular_diameter_diff_1/acetabular_diameter_mean) %>% 
  select(c("eid", "iliac_isthmus_breadth_diff_1", "acetabular_diameter_diff_1", "acetabular_inclination_diff_1",
           "iliac_isthmus_breadth_proportion", "acetabular_diameter_proportion"))

visit_2_df <- visit_2_df %>% 
  mutate(iliac_isthmus_breadth_diff_2 = abs(sciatic_notch_left2inferior_iliac_spine_left - sciatic_notch_right2inferior_iliac_spine_right)) %>% 
  mutate(acetabular_diameter_diff_2 = abs(iliopubic_eminence_left2acetabular_inferior_left - iliopubic_eminence_right2acetabular_inferior_right)) %>% 
  mutate(acetabular_inclination_diff_2 = abs(acetabular_inclination_left - acetabular_inclination_right)) %>% 
  mutate(iliac_isthmus_breadth_mean = (sciatic_notch_left2inferior_iliac_spine_left + sciatic_notch_right2inferior_iliac_spine_right)/2) %>% 
  mutate(acetabular_diameter_mean = (iliopubic_eminence_left2acetabular_inferior_left + iliopubic_eminence_right2acetabular_inferior_right)/2) %>% 
  mutate(acetabular_inclination_mean = (acetabular_inclination_left + acetabular_inclination_right)/2) %>% 
  mutate(iliac_isthmus_breadth_proportion = iliac_isthmus_breadth_diff_2/iliac_isthmus_breadth_mean) %>% 
  mutate(acetabular_diameter_proportion = acetabular_diameter_diff_2/acetabular_diameter_mean) %>% 
  select(c("eid", "iliac_isthmus_breadth_diff_2", "acetabular_diameter_diff_2", "acetabular_inclination_diff_2",
           "iliac_isthmus_breadth_proportion", "acetabular_diameter_proportion"))

visit_df <- merge(visit_1_df, visit_2_df, by = 'eid')

visit_df <- visit_df %>% 
  mutate(iliac_isthmus_breadth_mean = (sciatic_notch_left2inferior_iliac_spine_left + sciatic_notch_right2inferior_iliac_spine_right)/2) %>% 
  mutate(iliac_isthmus_breadth_mean = (sciatic_notch_left2inferior_iliac_spine_left + sciatic_notch_right2inferior_iliac_spine_right)/2) %>% 
  mutate(iliac_isthmus_breadth_mean = (sciatic_notch_left2inferior_iliac_spine_left + sciatic_notch_right2inferior_iliac_spine_right)/2)


# Iliac isthmus breadth

cor_value <- cor(visit_df$iliac_isthmus_breadth_diff_1, 
                 visit_df$iliac_isthmus_breadth_diff_2,
                 method = "pearson")

p3 <- ggplot(visit_df, 
             aes(x=iliac_isthmus_breadth_diff_1, 
                 y=iliac_isthmus_breadth_diff_2)) +
  geom_point(color="#59A14F", alpha = 1, size = 3, shape = 1) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  labs(title="Iliac isthmus breadth",
       x="First visit discrepancy (cm)",
       y="Second visit discrepancy (cm)") +
  annotate("text", x = Inf, y = -Inf, 
           label = sprintf("italic(r) == %.3f", cor_value), 
           hjust = 1, vjust = -1.2, 
           size = 6, color = "black", parse = TRUE) +
  xlim(c(-1.5, 1)) + ylim(c(-1.5, 1)) +
  theme_classic(base_size = 16)

# Acetabular diameter

cor_value <- cor(visit_df$acetabular_diameter_diff_1, 
                 visit_df$acetabular_diameter_diff_2,
                 method = "pearson")

p4 <- ggplot(visit_df, 
             aes(x=acetabular_diameter_diff_1, 
                 y=acetabular_diameter_diff_2)) +
  geom_point(color="#E15759", alpha = 1, size = 3, shape = 1) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  labs(title="Acetabular diameter",
       x="First visit discrepancy (cm)",
       y="") +
  annotate("text", x = Inf, y = -Inf, 
           label = sprintf("italic(r) == %.3f", cor_value), 
           hjust = 1, vjust = -1.2, 
           size = 6, color = "black", parse = TRUE) +
  xlim(c(-1, 1)) + ylim(c(-1, 1)) +
  theme_classic(base_size = 16)

# Acetabular inclination

# cor_value <- cor(visit_df$acetabular_inclination_diff_1, 
#                  visit_df$acetabular_inclination_diff_2,
#                  method = "pearson")

# p3 <- ggplot(visit_df, 
#              aes(x=acetabular_inclination_diff_1, 
#                  y=acetabular_inclination_diff_2)) +
#   geom_point(color="#8CD17D", alpha = 0.3, size = 3) +
#   geom_abline(intercept=0, slope=1, linetype="dashed") +
#   theme_classic() +
#   labs(title="Acetabular inclination",
#        x="1st visit discrepancy (degree)",
#        y="2nd visit discrepancy (degree)") +
#   annotate("text", x = Inf, y = -Inf, 
#            label = sprintf("italic(r) == %.3f", cor_value), 
#            hjust = 1.2, vjust = -2, 
#            size = 5, color = "black", parse = TRUE) 

pdf("out_fig/left_right_two_visits_discrepancy_correlation.pdf",
    width = 10, 
    height = 5)
p1 | p2 | p3 | p4
dev.off()
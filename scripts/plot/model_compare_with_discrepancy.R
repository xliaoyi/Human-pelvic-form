rm(list = ls())

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")

library(ggplot2)
library(patchwork)

me_arm_ratio_flt <- read.csv("key_results/me_arm_left_right_ratio_flt.csv", row.names = 1)
ek_arm_ratio_flt <- read.csv("key_results/ek_arm_left_right_ratio_flt.csv", row.names = 1)

# Calculate r and p-value for ek dataset
ek_cor_test <- cor.test(ek_arm_ratio_flt$arm_ratio_visit_1, ek_arm_ratio_flt$arm_ratio_visit_2)
ek_r2 <- ek_cor_test$estimate^2
ek_p <- ek_cor_test$p.value

# Calculate r and p-value for me dataset
me_cor_test <- cor.test(me_arm_ratio_flt$arm_ratio_visit_1, me_arm_ratio_flt$arm_ratio_visit_2)
me_r2 <- me_cor_test$estimate^2
me_p <- me_cor_test$p.value


ek <- ggplot(ek_arm_ratio_flt, aes(x = arm_ratio_visit_1, y = arm_ratio_visit_2)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(title = "Kun et al., 2023",
       x = "Left right arm ratio at first visit",
       y = "Left right arm ratio at second visit") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = Inf, y = -Inf, label = sprintf("italic(r)^2 == %.3f", ek_r2), hjust = 1.78, vjust = -4, parse = TRUE, size=26) +
  annotate("text", x = Inf, y = -Inf, label = sprintf("italic(P) == %.3e", ek_p), hjust = 1.2, vjust = -2.5, parse = TRUE, size=26) +
  xlim(c(0.94, 1.06)) + ylim(c(0.94, 1.06)) +
  theme_bw()

me <- ggplot(me_arm_ratio_flt, aes(x = arm_ratio_visit_1, y = arm_ratio_visit_2)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(title = "Optimized model in this study",
       x = "Left right arm ratio at first visit",
       y = "") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = Inf, y = -Inf, label = sprintf("italic(r)^2 == %.3f", me_r2), hjust = 1.95, vjust = -4, parse = TRUE, size=16) +
  annotate("text", x = Inf, y = -Inf, label = sprintf("italic(P) == %.3e", me_p), hjust = 1.2, vjust = -2.5, parse = TRUE, size=26) +
  xlim(c(0.94, 1.06)) + ylim(c(0.94, 1.06)) +
  theme_bw()

pdf("out_fig/model_compare_with_ratio.pdf",
    width = 10, 
    height = 5)
ek + me
dev.off()
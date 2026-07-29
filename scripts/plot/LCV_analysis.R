rm(list = ls())

library(dplyr)
library(ggplot2)

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/LCV_OUTPUT/")

df <- read.csv("LCV_female_output_results.csv")

# filteration
df <- df %>% mutate(p_value = 10^(log10_p)) %>% mutate(neglog10p = -log10_p)

# Add column indicating whether Bonferroni-corrected p-value is significant
df$FDR <- p.adjust(df$p_value, method = "fdr") < 0.05

# keep p < 0.05
df$p_value <- ifelse(df$neglog10p < -log10(0.5), 0, df$neglog10p)

# Define the color palette
color_palette <- colorRampPalette(c("steelblue", "steelblue1", "white", "darkorange", "indianred"))(100)

p <- ggplot(df, aes(x = trait2, y = trait1)) +
  geom_point(data = subset(df, neglog10p >= -log10(0.5)), 
             aes(size = pmin(neglog10p, 6), fill = gcp), shape = 21) +
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
             aes(x = trait2, y = trait1),
             color = "black", shape = 8, size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "", y = "", fill = "Correlation", size = "-log10(p-value)") +
  ggtitle("LCV between female pelvic phenotypes and outcome phenotypes")

ggsave(filename = "../../out_fig/LCV_female.pdf", plot = p,
       width = 15, height = 9, units = "in")

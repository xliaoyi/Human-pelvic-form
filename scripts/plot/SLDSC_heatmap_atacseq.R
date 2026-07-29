rm(list=ls())

library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(reshape2)

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/")

select_pheno <- c("pelvic_height", "pelvic_width", "pelvic_inlet_width", "subpubic_angle",
                  "oblique_pelvic_inlet_length", "ear_left2ear_right", "arm_devide_torso",
                  "shoulder_width", "iliac_isthmus_breadth", "acetabular_diameter", "pelvic_inlet_area", 
                  "head_divide_inlet_width", "shoulder_divide_inlet_width",
                  "head_divide_oblique_inlet_length", "shoulder_divide_oblique_inlet_length", 
                  "head_area_divide_pelvic_inlet_area", "shoulder_area_divide_pelvic_inlet_area")
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
  "shoulder_area_divide_pelvic_inlet_area" = "Shoulder area:pelvic inlet area"
)

######### for both
df <- read.csv("key_results/SLDSC_results_from_sa_paper_both.csv")

df <- df %>% filter(Phenotype %in% select_pheno)

df$fdr <- p.adjust(df$Enrichment_p, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((fdr < 0.05), "*", ""))

# Convert data to wide format suitable for the heatmap
df_wide <- dcast(df, Annotation ~ Phenotype, value.var = "Enrichment")

# rename phenotypes
names(df_wide) <- ifelse(names(df_wide) %in% names(rename_map_pheno),
                         rename_map_pheno[names(df_wide)],
                         names(df_wide))

df_wide <- merge(df_wide, unique(df %>% select(c("Annotation", "Group"))), by = "Annotation")
mat <- as.data.frame(df_wide[,-1])
rownames(mat) <- df_wide$Annotation

# Prepare the annotation data
df_wide_annotation <- dcast(df, Annotation ~ Phenotype, value.var = "SIGNIF")
annotation_mat <- as.matrix(df_wide_annotation[,-1])
rownames(annotation_mat) <- df_wide_annotation$Annotation

# PLOT
col_fun = colorRamp2(c(-5, 1, 50), c("lightblue", "white", "darkorange"))
hm <- Heatmap(mat %>% select(-Group), name = "mat", col = col_fun,
              heatmap_legend_param = list(title = "Enrichment"),
              row_split = factor(mat$Group, levels = c("Human", "Human_Brain_Filtered", "Human_Gain", "Mouse", "Mouse_Brain_Filtered", "Mouse_Gain", "Human_Mouse_Intersect", "Human_Mouse_Shared")), 
              row_title = NULL,
              row_gap = unit(5, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) { 
                grid.text(sprintf(annotation_mat[i, j]), x, y, just=c("center","center"), gp = gpar(fontsize = 15, col = "black", fontface = "bold"))
              }, cluster_columns = FALSE, cluster_rows = FALSE,
              cluster_row_slices = FALSE,
              border = TRUE,
              rect_gp = gpar(col = "white", lwd = 2),
              show_heatmap_legend = FALSE
)
lgd = Legend(col_fun = col_fun, direction = "horizontal", title = 'Enrichment')

pdf(file = "out_fig/SLDSC_with_ATAC-seq_from_sa_paper_both.pdf", width = 14, height = 14)
draw(hm, column_title = "SLDSC for ATAC-Seq Both Male and Female", padding = unit(c(22, 20, 2, 20), "mm"))
draw(lgd, x = unit(0.85, "npc"), y = unit(0.17, "npc"))
dev.off()



######### for female
df <- read.csv("key_results/SLDSC_results_from_sa_paper_female.csv")

df <- df %>% filter(Phenotype %in% select_pheno)

df$fdr <- p.adjust(df$Enrichment_p, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((fdr < 0.05), "*", ""))

# Convert data to wide format suitable for the heatmap
df_wide <- dcast(df, Annotation ~ Phenotype, value.var = "Enrichment")

# rename phenotypes
names(df_wide) <- ifelse(names(df_wide) %in% names(rename_map_pheno),
                         rename_map_pheno[names(df_wide)],
                         names(df_wide))

df_wide <- merge(df_wide, unique(df %>% select(c("Annotation", "Group"))), by = "Annotation")
mat <- as.data.frame(df_wide[,-1])
rownames(mat) <- df_wide$Annotation

# Prepare the annotation data
df_wide_annotation <- dcast(df, Annotation ~ Phenotype, value.var = "SIGNIF")
annotation_mat <- as.matrix(df_wide_annotation[,-1])
rownames(annotation_mat) <- df_wide_annotation$Annotation

# PLOT
col_fun = colorRamp2(c(-5, 1, 50), c("lightblue", "white", "darkorange"))
hm <- Heatmap(mat %>% select(-Group), name = "mat", col = col_fun,
              heatmap_legend_param = list(title = "Enrichment"),
              row_split = factor(mat$Group, levels = c("Human", "Human_Brain_Filtered", "Human_Gain", "Mouse", "Mouse_Brain_Filtered", "Mouse_Gain", "Human_Mouse_Intersect", "Human_Mouse_Shared")), 
              row_title = NULL,
              row_gap = unit(5, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf(annotation_mat[i, j]), x, y, just=c("center","center"), gp = gpar(fontsize = 15, col = "black", fontface = "bold"))
              }, cluster_columns = FALSE, cluster_rows = FALSE,
              cluster_row_slices = FALSE,
              border = TRUE,
              rect_gp = gpar(col = "white", lwd = 2),
              show_heatmap_legend = FALSE
)
lgd = Legend(col_fun = col_fun, direction = "horizontal", title = 'Enrichment')

pdf(file = "out_fig/SLDSC_with_ATAC-seq_from_sa_paper_female.pdf", width = 14, height = 14)
draw(hm, column_title = "SLDSC for ATAC-Seq Female", padding = unit(c(22, 20, 2, 20), "mm"))
draw(lgd, x = unit(0.85, "npc"), y = unit(0.17, "npc"))
dev.off()



######### for male
df <- read.csv("key_results/SLDSC_results_from_sa_paper_male.csv")

df <- df %>% filter(Phenotype %in% select_pheno)

df$fdr <- p.adjust(df$Enrichment_p, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((fdr < 0.05), "*", ""))

# Convert data to wide format suitable for the heatmap
df_wide <- dcast(df, Annotation ~ Phenotype, value.var = "Enrichment")


# rename phenotypes
names(df_wide) <- ifelse(names(df_wide) %in% names(rename_map_pheno),
                         rename_map_pheno[names(df_wide)],
                         names(df_wide))

df_wide <- merge(df_wide, unique(df %>% select(c("Annotation", "Group"))), by = "Annotation")
mat <- as.data.frame(df_wide[,-1])
rownames(mat) <- df_wide$Annotation

# Prepare the annotation data
df_wide_annotation <- dcast(df, Annotation ~ Phenotype, value.var = "SIGNIF")
annotation_mat <- as.matrix(df_wide_annotation[,-1])
rownames(annotation_mat) <- df_wide_annotation$Annotation

# PLOT
col_fun = colorRamp2(c(-5, 1, 50), c("lightblue", "white", "darkorange"))
hm <- Heatmap(mat %>% select(-Group), name = "mat", col = col_fun,
              heatmap_legend_param = list(title = "Enrichment"),
              row_split = factor(mat$Group, levels = c("Human", "Human_Brain_Filtered", "Human_Gain", "Mouse", "Mouse_Brain_Filtered", "Mouse_Gain", "Human_Mouse_Intersect", "Human_Mouse_Shared")), 
              row_title = NULL,
              row_gap = unit(5, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf(annotation_mat[i, j]), x, y, just=c("center","center"), gp = gpar(fontsize = 15, col = "black", fontface = "bold"))
              }, cluster_columns = FALSE, cluster_rows = FALSE,
              cluster_row_slices = FALSE,
              border = TRUE,
              rect_gp = gpar(col = "white", lwd = 2),
              show_heatmap_legend = FALSE
)
lgd = Legend(col_fun = col_fun, direction = "horizontal", title = 'Enrichment')

pdf(file = "out_fig/SLDSC_with_ATAC-seq_from_sa_paper_male.pdf", width = 14, height = 14)
draw(hm, column_title = "SLDSC for ATAC-Seq Male", padding = unit(c(22, 20, 2, 20), "mm"))
draw(lgd, x = unit(0.85, "npc"), y = unit(0.17, "npc"))
dev.off()
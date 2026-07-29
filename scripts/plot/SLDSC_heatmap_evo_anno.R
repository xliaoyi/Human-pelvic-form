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
df <- read.csv("key_results/SLDSC_results_evo_anno_both.csv")
df$Annotation <- factor(df$Annotation, levels = c("GSE63648_7pcw_ac_me2_Hu_gain",
                                                  "GSE63648_8_5pcw_ac_me2_Hu_gain",
                                                  "GSE63648_12Fpcw_ac_me2_Hu_gain",
                                                  "GSE63648_12Opcw_ac_me2_Hu_gain",
                                                  "enhancers_human_macaque_prefontal_cortex_hg19",
                                                  "promoters_human_macaque_prefontal_cortex_hg19",
                                                  "nchaes_merged_hg19",
                                                  "selsweeps_233putativeregions_peyregne",
                                                  "depleted_neanderthal_vernot",
                                                  "depleted_neanderthal_denisovan_vernot"))
df <- df %>% filter(Phenotype %in% select_pheno)
df <- df %>% filter(!Annotation %in% "nchaes_merged_hg19")

df$fdr <- p.adjust(df$Enrichment_p, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((fdr < 0.05), "*", ""))

# Convert data to wide format suitable for the heatmap
df_wide <- dcast(df, Annotation ~ Phenotype, value.var = "Enrichment")

# rename phenotypes
names(df_wide) <- ifelse(names(df_wide) %in% names(rename_map_pheno),
                         rename_map_pheno[names(df_wide)],
                         names(df_wide))

mat <- as.data.frame(df_wide[,-1])
rownames(mat) <- df_wide$Annotation
mat$Avg <- rowMeans(mat)


# Prepare the annotation data
df_wide_annotation <- dcast(df, Annotation ~ Phenotype, value.var = "SIGNIF")
annotation_mat <- as.matrix(df_wide_annotation[,-1])
rownames(annotation_mat) <- df_wide_annotation$Annotation

# PLOT
col_fun = colorRamp2(c(-4, -1, 1, 11), c("steelblue", "steelblue1", "white", "darkorange"))
row_ha = rowAnnotation(avg=anno_barplot(mat$Avg), width = unit(4, "cm"))
hm <- Heatmap(mat %>% select(-Avg), name = "mat", col = col_fun,
              heatmap_legend_param = list(title = "Enrichment"),
              # row_split = factor(mat$Group, levels = c("Human", "Human_Brain_Filtered", "Human_Gain", "Mouse", "Mouse_Brain_Filtered", "Mouse_Gain", "Human_Mouse_Intersect", "Human_Mouse_Shared")), 
              # row_title = NULL,
              # row_gap = unit(5, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) { 
                grid.text(sprintf(annotation_mat[i, j]), x, y, just=c("center","center"), gp = gpar(fontsize = 15, col = "black", fontface = "bold"))
              }, cluster_columns = FALSE, cluster_rows = FALSE,
              cluster_row_slices = FALSE,
              border = TRUE,
              rect_gp = gpar(col = "white", lwd = 2),
              show_heatmap_legend = FALSE,
              column_names_side = "bottom", row_names_side = "left",
              right_annotation = row_ha
)
lgd = Legend(col_fun = col_fun, direction = "horizontal", title = 'Enrichment', at = c(-4, 1, 6, 11), legend_width = unit(1.5, "inches"))

pdf(file = "out_fig/SLDSC_with_evo_anno_both.pdf", width = 15, height = 10)
draw(hm, padding = unit(c(20, 50, 20, 20), "mm"))
draw(lgd, x = unit(0.2, "npc"), y = unit(0.25, "npc"))
dev.off()


######### for female 
df <- read.csv("key_results/SLDSC_results_evo_anno_female.csv")
df$Annotation <- factor(df$Annotation, levels = c("GSE63648_7pcw_ac_me2_Hu_gain",
                                                  "GSE63648_8_5pcw_ac_me2_Hu_gain",
                                                  "GSE63648_12Fpcw_ac_me2_Hu_gain",
                                                  "GSE63648_12Opcw_ac_me2_Hu_gain",
                                                  "enhancers_human_macaque_prefontal_cortex_hg19",
                                                  "promoters_human_macaque_prefontal_cortex_hg19",
                                                  "nchaes_merged_hg19",
                                                  "selsweeps_233putativeregions_peyregne",
                                                  "depleted_neanderthal_vernot",
                                                  "depleted_neanderthal_denisovan_vernot"))
df <- df %>% filter(Phenotype %in% select_pheno)
df <- df %>% filter(!Annotation %in% "nchaes_merged_hg19")

df$fdr <- p.adjust(df$Enrichment_p, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((fdr < 0.05), "*", ""))

# Convert data to wide format suitable for the heatmap
df_wide <- dcast(df, Annotation ~ Phenotype, value.var = "Enrichment")
mat <- as.data.frame(df_wide[,-1])
rownames(mat) <- df_wide$Annotation

mat$Avg <- rowMeans(mat)


# Prepare the annotation data
df_wide_annotation <- dcast(df, Annotation ~ Phenotype, value.var = "SIGNIF")
annotation_mat <- as.matrix(df_wide_annotation[,-1])
rownames(annotation_mat) <- df_wide_annotation$Annotation

# PLOT
col_fun = colorRamp2(c(-4, -1, 1, 11), c("steelblue", "steelblue1", "white", "darkorange"))
row_ha = rowAnnotation(avg=anno_barplot(mat$Avg), width = unit(4, "cm"))
hm <- Heatmap(mat %>% select(-Avg), name = "mat", col = col_fun,
              heatmap_legend_param = list(title = "Enrichment"),
              # row_split = factor(mat$Group, levels = c("Human", "Human_Brain_Filtered", "Human_Gain", "Mouse", "Mouse_Brain_Filtered", "Mouse_Gain", "Human_Mouse_Intersect", "Human_Mouse_Shared")), 
              # row_title = NULL,
              # row_gap = unit(5, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) { 
                grid.text(sprintf(annotation_mat[i, j]), x, y, just=c("center","center"), gp = gpar(fontsize = 15, col = "black", fontface = "bold"))
              }, cluster_columns = FALSE, cluster_rows = FALSE,
              cluster_row_slices = FALSE,
              border = TRUE,
              rect_gp = gpar(col = "white", lwd = 2),
              show_heatmap_legend = FALSE,
              column_names_side = "bottom", row_names_side = "left",
              right_annotation = row_ha
)
lgd = Legend(col_fun = col_fun, direction = "horizontal", title = 'Enrichment', at = c(-4, 1, 6, 11), legend_width = unit(1.5, "inches"))

pdf(file = "out_fig/SLDSC_with_evo_anno_female.pdf", width = 15, height = 10)
draw(hm, padding = unit(c(20, 50, 20, 20), "mm"))
draw(lgd, x = unit(0.2, "npc"), y = unit(0.25, "npc"))
dev.off()

######### for male 
df <- read.csv("key_results/SLDSC_results_evo_anno_male.csv")
df$Annotation <- factor(df$Annotation, levels = c("GSE63648_7pcw_ac_me2_Hu_gain",
                                                  "GSE63648_8_5pcw_ac_me2_Hu_gain",
                                                  "GSE63648_12Fpcw_ac_me2_Hu_gain",
                                                  "GSE63648_12Opcw_ac_me2_Hu_gain",
                                                  "enhancers_human_macaque_prefontal_cortex_hg19",
                                                  "promoters_human_macaque_prefontal_cortex_hg19",
                                                  "nchaes_merged_hg19",
                                                  "selsweeps_233putativeregions_peyregne",
                                                  "depleted_neanderthal_vernot",
                                                  "depleted_neanderthal_denisovan_vernot"))
df <- df %>% filter(Phenotype %in% select_pheno)
df <- df %>% filter(!Annotation %in% "nchaes_merged_hg19")

df$fdr <- p.adjust(df$Enrichment_p, method = "fdr")
df <- df %>% mutate("SIGNIF" = ifelse((fdr < 0.05), "*", ""))

# Convert data to wide format suitable for the heatmap
df_wide <- dcast(df, Annotation ~ Phenotype, value.var = "Enrichment")
mat <- as.data.frame(df_wide[,-1])
rownames(mat) <- df_wide$Annotation

mat$Avg <- rowMeans(mat)


# Prepare the annotation data
df_wide_annotation <- dcast(df, Annotation ~ Phenotype, value.var = "SIGNIF")
annotation_mat <- as.matrix(df_wide_annotation[,-1])
rownames(annotation_mat) <- df_wide_annotation$Annotation

# PLOT
col_fun = colorRamp2(c(-4, -1, 1, 11), c("steelblue", "steelblue1", "white", "darkorange"))
row_ha = rowAnnotation(avg=anno_barplot(mat$Avg), width = unit(4, "cm"))
hm <- Heatmap(mat %>% select(-Avg), name = "mat", col = col_fun,
              heatmap_legend_param = list(title = "Enrichment"),
              # row_split = factor(mat$Group, levels = c("Human", "Human_Brain_Filtered", "Human_Gain", "Mouse", "Mouse_Brain_Filtered", "Mouse_Gain", "Human_Mouse_Intersect", "Human_Mouse_Shared")), 
              # row_title = NULL,
              # row_gap = unit(5, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) { 
                grid.text(sprintf(annotation_mat[i, j]), x, y, just=c("center","center"), gp = gpar(fontsize = 15, col = "black", fontface = "bold"))
              }, cluster_columns = FALSE, cluster_rows = FALSE,
              cluster_row_slices = FALSE,
              border = TRUE,
              rect_gp = gpar(col = "white", lwd = 2),
              show_heatmap_legend = FALSE,
              column_names_side = "bottom", row_names_side = "left",
              right_annotation = row_ha
)
lgd = Legend(col_fun = col_fun, direction = "horizontal", title = 'Enrichment', at = c(-4, 1, 6, 11), legend_width = unit(1.5, "inches"))

pdf(file = "out_fig/SLDSC_with_evo_anno_male.pdf", width = 15, height = 10)
draw(hm, padding = unit(c(20, 50, 20, 20), "mm"))
draw(lgd, x = unit(0.2, "npc"), y = unit(0.25, "npc"))
dev.off()


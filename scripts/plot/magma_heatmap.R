rm(list= ls())

library(ggplot2)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(viridis)

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
#######################################################
# MAGMA with RNA-Seq data from science advance paper
####################################################### 

select_pheno <- c("pelvic_height", "pelvic_width", "pelvic_inlet_width", "subpubic_angle",
                  "oblique_pelvic_inlet_length", "iliac_isthmus_breadth", "acetabular_diameter")

# for different subelements
# for both sex
setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/MAGMA_OUTPUT/20230730_hip_both/diff_subelements/")
file_list <- list.files(pattern = "\\.gsa.out$")

all_pheno_magma <- NULL

for (file in file_list) {
  basename <- sub("2023-07-30_", "", file)
  basename <- sub("_geneset_analysis.gsa.out", "", basename)
  
  if (basename %in% select_pheno) {
    
    df <- read.csv(file, sep = '', skip = 4, header = T)
    df <- df %>% 
      mutate("PHENO" = basename) %>% 
      mutate("NEGLOG10P" = -log10(P))
    
    all_pheno_magma <- rbind(all_pheno_magma, df)
    
  }
}

# Apply renaming to the df dataframe
# all_pheno_magma <- all_pheno_magma %>%
#   mutate(
#     PHENO = ifelse(PHENO %in% names(rename_map_pheno), rename_map_pheno[PHENO], PHENO)
#   )

all_pheno_magma$FDR <- p.adjust(all_pheno_magma$P, method = "fdr")
all_pheno_magma <- all_pheno_magma %>% mutate("SIGNIF" = ifelse((FDR < 0.05), "*", ""))
write.csv(all_pheno_magma, file = "/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/key_results/magma_subelements.csv")
# Convert data to wide format suitable for the heatmap
df_wide <- dcast(all_pheno_magma, PHENO ~ VARIABLE, value.var = "NEGLOG10P")
df_wide$PHENO <- factor(df_wide$PHENO, levels = unique(all_pheno_magma$PHENO))

# Convert df_wide into a matrix, keeping the PHENO column as row names
mat <- data.matrix(df_wide[,-1])
rownames(mat) <- df_wide$PHENO

# Prepare the annotation data
df_wide_annotation <- dcast(all_pheno_magma, PHENO ~ VARIABLE, value.var = "SIGNIF")
annotation_mat <- as.matrix(df_wide_annotation[,-1])
rownames(annotation_mat) <- df_wide_annotation$PHENO

# Plot
col_fun = colorRamp2(c(min(mat), median(mat), max(mat)), c("white", "darkolivegreen3", "darkgreen"))
hm <- Heatmap(mat, name = "mat", col = col_fun,
              heatmap_legend_param = list(title = "-log10(q-value)"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf(annotation_mat[i, j]), x, y, just=c("center","center"), gp = gpar(fontsize = 15, col = "white", fontface = "bold"))
              }, cluster_columns = TRUE,
              show_heatmap_legend = FALSE)
lgd = Legend(col_fun = col_fun, direction = "horizontal", title = '-log10(q-value)')
pdf(file = "../../../../out_fig/MAGMA_pelvic_subelements_RNA_seq_diff_subelements_diff_subelements_240213.pdf", width = 6, height = 7)
draw(hm, column_title = "MAGMA with pelvic subelements RNA-Seq")
draw(lgd, x = unit(0.78, "npc"), y = unit(0.1, "npc"))
dev.off()


# for different subelements at different time points
# for both sex
setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/MAGMA_OUTPUT/20230730_hip_both/diff_timepoints/")
file_list <- list.files(pattern = "\\.gsa.out$")

all_pheno_magma <- NULL

for (file in file_list) {
  basename <- sub("2023-07-30_", "", file)
  basename <- sub("_geneset_analysis.gsa.out", "", basename)
  
  if (basename %in% select_pheno) {
    
    df <- read.csv(file, sep = '', skip = 4, header = T)
    df <- df %>% 
      mutate("PHENO" = basename) %>% 
      mutate("NEGLOG10P" = -log10(P))
    
    all_pheno_magma <- rbind(all_pheno_magma, df)
    
  }
}
# Apply renaming to the df dataframe
# all_pheno_magma <- all_pheno_magma %>%
#   mutate(
#     PHENO = ifelse(PHENO %in% names(rename_map_pheno), rename_map_pheno[PHENO], PHENO)
#   )
all_pheno_magma$FDR <- p.adjust(all_pheno_magma$P, method = "fdr")
all_pheno_magma <- all_pheno_magma %>% mutate("SIGNIF" = ifelse((FDR < 0.05), "*", ""))
write.csv(all_pheno_magma, file = "/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/key_results/magma_subelements_diff_time.csv")
# Convert data to wide format suitable for the heatmap
df_wide <- dcast(all_pheno_magma, PHENO ~ VARIABLE, value.var = "NEGLOG10P")
df_wide$PHENO <- factor(df_wide$PHENO, levels = unique(all_pheno_magma$PHENO))

# Convert df_wide into a matrix, keeping the PHENO column as row names
mat <- data.matrix(df_wide[,-1])
rownames(mat) <- df_wide$PHENO

# Prepare the annotation data
df_wide_annotation <- dcast(all_pheno_magma, PHENO ~ VARIABLE, value.var = "SIGNIF")
annotation_mat <- as.matrix(df_wide_annotation[,-1])
rownames(annotation_mat) <- df_wide_annotation$PHENO

# Plot
col_fun = colorRamp2(c(min(mat), median(mat), max(mat)), c("white", "darkolivegreen3", "darkgreen"))
hm <- Heatmap(mat, name = "mat", col = col_fun,
              heatmap_legend_param = list(title = "-log10(q-value)"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf(annotation_mat[i, j]), x, y, gp = gpar(fontsize = 15))
              }, cluster_columns = FALSE,
              show_heatmap_legend = FALSE)
lgd = Legend(col_fun = col_fun, direction = "horizontal", title = '-log10(q-value)')
pdf(file = "../../../../out_fig/MAGMA_pelvic_subelements_RNA_seq_diff_subelements_diff_timepoints_240213.pdf", width = 14, height = 8)
draw(hm, column_title = "MAGMA with pelvic subelements RNA-Seq")
draw(lgd, x = unit(0.88, "npc"), y = unit(0.1, "npc"))
dev.off()
#######################################################
#                  MAGMA with GTEx v8
#######################################################

# for both sex
setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics/FUMA_OUTPUT/2023-07-26_hip_both")
file_list <- list.files()

all_pheno_magma <- NULL

for (file in file_list) {
  path <- paste0(file, "/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out")
  
  df <- read.csv(path , sep = '', skip = 5, header = T)
  df <- df %>% 
    mutate("PHENO" = file) %>% 
    mutate("NEGLOG10P" = -log10(P))
  
  all_pheno_magma <- rbind(all_pheno_magma, df)
}

all_pheno_magma <- all_pheno_magma %>% filter(PHENO %in% select_pheno)
# Apply renaming to the df dataframe
# all_pheno_magma <- all_pheno_magma %>%
#   mutate(
#     PHENO = ifelse(PHENO %in% names(rename_map_pheno), rename_map_pheno[PHENO], PHENO)
#   )
all_pheno_magma$FDR <- p.adjust(all_pheno_magma$P, method = "fdr")
all_pheno_magma <- all_pheno_magma %>% mutate("SIGNIF" = ifelse((FDR < 0.05), "*", ""))
write.csv(all_pheno_magma, file = "/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/key_results/magma_gtex_v8.csv")
# Convert data to wide format suitable for the heatmap
df_wide <- dcast(all_pheno_magma, PHENO ~ VARIABLE, value.var = "NEGLOG10P")
df_wide$PHENO <- factor(df_wide$PHENO, levels = unique(all_pheno_magma$PHENO))

# Convert df_wide into a matrix, keeping the PHENO column as row names
mat <- data.matrix(df_wide[,-1])
rownames(mat) <- df_wide$PHENO

# Prepare the annotation data
df_wide_annotation <- dcast(all_pheno_magma, PHENO ~ VARIABLE, value.var = "SIGNIF")
annotation_mat <- as.matrix(df_wide_annotation[,-1])
rownames(annotation_mat) <- df_wide_annotation$PHENO

# Plot
col_fun = colorRamp2(c(min(mat), median(mat), max(mat)), c("white", "plum3", "purple3"))

hm <- Heatmap(mat, name = "mat", col = col_fun, 
              # heatmap_legend_param = list(title = "FDR-adjusted\n-log10(p)"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf(annotation_mat[i, j]), x, y, just = c("center", "top"), gp = gpar(fontsize = 15))
              },
              show_heatmap_legend = FALSE)
lgd = Legend(col_fun = col_fun, direction = "horizontal", title = '-log10(q-value)')
pdf(file = "../../../out_fig/MAGMA_GTEx_v8_240213.pdf", width = 16, height = 8)
draw(hm, column_title = "MAGMA with GTEx v8", padding = unit(c(10, 2, 2, 2), "mm"))
draw(lgd, x = unit(0.9, "npc"), y = unit(0.25, "npc"))
dev.off()

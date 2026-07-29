rm(list=ls())
library(biomaRt)
library(dplyr)

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/from_sci_adv")

# Set up the mart object
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh="GRCh37")

# Get a list of all files in the directory
files <- list.files(path = "RNA-Seq_data", pattern = "*.genes.results.txt", full.names = TRUE)

all_df <- NULL

# Loop through all files
for (file.name in files) {
  # Read in the data
  sa_data <- read.csv(file.name, sep = "\t")
  
  # Extract the gene_ids
  gene_ids <- as.vector(sa_data$gene_id)
  
  # Map gene_ids to ENSG IDs
  mapped_id <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), 
                     filters = 'external_gene_name', 
                     values = gene_ids, 
                     mart = ensembl)
  
  # Merge the mapping with the original data
  df <- merge(mapped_id, sa_data, by.x = 'external_gene_name', by.y = 'gene_id')
  
  # get tissue name
  split_string <- strsplit(file.name, "/")[[1]]
  final_split <- strsplit(split_string[2], "\\.")[[1]]
  tissue.name <- final_split[1]
  
  # Filter and select columns
  df <- df %>% dplyr::select(c("ensembl_gene_id", "TPM")) %>%
    rename_with(~ "gene_id", "ensembl_gene_id") %>%
    rename_with(~ tissue.name, "TPM") 
    # filter(get(tissue.name) >= 1)

  # Merge with the all_df dataframe
  if (is.null(all_df)) {
    all_df <- df
  } else {
    all_df <- merge(all_df, df, by = 'gene_id')
  }
  
  print(dim(all_df))
}
write.table(all_df, file = "RNA-Seq_data/all_tissues_TPM.txt", sep = "\t", row.names = FALSE)

######################################################################
#     average of different time points of different subelements
######################################################################
all_df <- read.table("RNA-Seq_data/all_tissues_TPM.txt", sep = " ", header = T)

# remove GSM5058379_Ilium59, the author mentioned this data is not accurate

# "Seven biological replicates of each pelvic subelement were sequenced, but 
# one acetabulum sample failed in sequencing because of low concentration, 
# and one ilium sample was removed in downstream processing because of highly 
# outlying values. Therefore, N = 6 for acetabulum and ilium."
cols <- colnames(all_df)
cols <- cols[!cols %in% c("gene_id", "GSM5058379_Ilium59")]

# Acetabulum
remain.cols <- cols[!cols %in% c("GSM5058367_Acetabulum53A")]
all_df <- mutate(all_df, Acetabulum53 = all_df[, c("GSM5058367_Acetabulum53A")] - rowMeans(all_df[, remain.cols]))
remain.cols <- cols[!cols %in% c("GSM5058368_Acetabulum54A", "GSM5058369_Acetabulum54B")]
all_df <- mutate(all_df, Acetabulum54 = rowMeans(all_df[, c("GSM5058368_Acetabulum54A", "GSM5058369_Acetabulum54B")]) - rowMeans(all_df[, remain.cols]))
remain.cols <- cols[!cols %in% c("GSM5058370_Acetabulum57A", "GSM5058371_Acetabulum57B")]
all_df <- mutate(all_df, Acetabulum57 = rowMeans(all_df[, c("GSM5058370_Acetabulum57A", "GSM5058371_Acetabulum57B")]) - rowMeans(all_df[, remain.cols]))
remain.cols <- cols[!cols %in% c("GSM5058372_Acetabulum59")]
all_df <- mutate(all_df, Acetabulum59 = all_df[, c("GSM5058372_Acetabulum59")] - rowMeans(all_df[, remain.cols]))

# Ilium
remain.cols <- cols[!cols %in% c("GSM5058373_Ilium53A")]
all_df <- mutate(all_df, Ilium53 = all_df[, c("GSM5058367_Acetabulum53A")] - rowMeans(all_df[, remain.cols]))
remain.cols <- cols[!cols %in% c("GSM5058374_Ilium54A", "GSM5058375_Ilium54B")]
all_df <- mutate(all_df, Ilium54 = rowMeans(all_df[, c("GSM5058374_Ilium54A", "GSM5058375_Ilium54B")]) - rowMeans(all_df[, remain.cols]))
remain.cols <- cols[!cols %in% c("GSM5058376_Ilium57A", "GSM5058377_Ilium57B", "GSM5058378_Ilium57C")]
all_df <- mutate(all_df, Ilium57 = rowMeans(all_df[, c("GSM5058376_Ilium57A", "GSM5058377_Ilium57B", "GSM5058378_Ilium57C")]) - rowMeans(all_df[, remain.cols]))

# Ischium
remain.cols <- cols[!cols %in% c("GSM5058380_Ischium53A")]
all_df <- mutate(all_df, Ischium53 = all_df[, c("GSM5058380_Ischium53A")] - rowMeans(all_df[, remain.cols]))
remain.cols <- cols[!cols %in% c("GSM5058381_Ischium54A", "GSM5058382_Ischium54B")]
all_df <- mutate(all_df, Ischium54 = rowMeans(all_df[, c("GSM5058381_Ischium54A", "GSM5058382_Ischium54B")]) - rowMeans(all_df[, remain.cols]))
remain.cols <- cols[!cols %in% c("GSM5058383_Ischium57A", "GSM5058384_Ischium57B", "GSM5058385_Ischium57C")]
all_df <- mutate(all_df, Ischium57 = rowMeans(all_df[, c("GSM5058383_Ischium57A", "GSM5058384_Ischium57B", "GSM5058385_Ischium57C")]) - rowMeans(all_df[, remain.cols]))
remain.cols <- cols[!cols %in% c("GSM5058386_Ischium59")]
all_df <- mutate(all_df, Ischium59 = all_df[, c("GSM5058386_Ischium59")] - rowMeans(all_df[, remain.cols]))

# Pubis
remain.cols <- cols[!cols %in% c("GSM5058387_Pubis53A")]
all_df <- mutate(all_df, Pubis53 = all_df[, c("GSM5058387_Pubis53A")] - rowMeans(all_df[, remain.cols]))
remain.cols <- cols[!cols %in% c("GSM5058388_Pubis54A", "GSM5058389_Pubis54B")]
all_df <- mutate(all_df, Pubis54 = rowMeans(all_df[, c("GSM5058388_Pubis54A", "GSM5058389_Pubis54B")]) - rowMeans(all_df[, remain.cols]))
remain.cols <- cols[!cols %in% c("GSM5058390_Pubis57A", "GSM5058391_Pubis57B", "GSM5058392_Pubis57C")]
all_df <- mutate(all_df, Pubis57 = rowMeans(all_df[, c("GSM5058390_Pubis57A", "GSM5058391_Pubis57B", "GSM5058392_Pubis57C")]) - rowMeans(all_df[, remain.cols]))
remain.cols <- cols[!cols %in% c("GSM5058393_Pubis59")]
all_df <- mutate(all_df, Pubis59 = all_df[, c("GSM5058393_Pubis59")] - rowMeans(all_df[, remain.cols]))

# remove columns startswith "GSM"
all_df <- all_df[, !grepl("^GSM", names(all_df))]

write.table(all_df, file = "RNA-Seq_data/all_tissues_TPM_diff_time_points.txt", sep = " ", row.names = FALSE)


######################################################################
#                 average of different subelements
######################################################################
all_df <- read.table("RNA-Seq_data/all_tissues_TPM.txt", sep = " ", header = T)

# Define the vector of subelements
subelements <- c("Acetabulum", "Ilium", "Ischium", "Pubis")

new_all_df <- all_df
# For each subelement
for (subelement in subelements) {
  
  # Find columns containing the subelement
  cols_subelement <- grep(subelement, names(all_df), value = TRUE)
  
  # Find the remaining columns excluding "gene_id"
  cols_remain <- setdiff(names(all_df), c("gene_id", cols_subelement))
  
  # Subtract mean of remaining columns
  new_all_df[subelement] <- rowMeans(all_df[cols_subelement]) - rowMeans(all_df[cols_remain])
}

# remove columns startswith "GSM"
new_all_df <- new_all_df[, !grepl("^GSM", names(new_all_df))]

write.table(new_all_df, file = "RNA-Seq_data/all_tissues_TPM_diff_subelements.txt", sep = " ", row.names = FALSE)


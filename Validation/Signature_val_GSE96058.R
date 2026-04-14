#> In this script we obtain and preprocess the metadata and the counts data

library(GEOquery)
library(tidyverse)

# 1.- Download data -----------------------------------------------------------


# 1.1 Download the supplementary file (The actual expression matrix)
# getGEOSuppFiles("GSE96058", baseDir = "D:/GSE96058")

# 1.1.2 Read the specific expression file (SCAN-B typically provides a large .txt or .csv)

raw <- read.csv("D:/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz", 
                          row.names = 1, check.names = FALSE)
counts_data <- raw
# 1.2 Download metadata

gse <- getGEO("GSE96058", GSEMatrix = TRUE)

# 1.2.2 Asign to object

pheno <- pData(gse[[1]])



# 2.- Preprocess metadata -------------------------------------------------

pheno$characteristics_ch1.3 <- gsub("\\D", "", pheno$characteristics_ch1.3)


metadata <-
  pheno %>%
  mutate(
    tissue = source_name_ch1,
    age = as.numeric(`age at diagnosis:ch1`),
    tumor_size = as.numeric(`tumor size:ch1`),
    lymph_group = `lymph node group:ch1`,
    lymph_status =  `lymph node status:ch1`,
    er_status = as.numeric(`er status:ch1`),
    pgr_status = as.numeric(`pgr status:ch1`),
    her2_status = as.numeric(`her2 status:ch1`),
    ki67_status = as.numeric(`ki67 status:ch1`),
    nhg = as.factor(`nhg:ch1`),
    er_pred_mgc = as.numeric(`er prediction mgc:ch1`),
    # This are predictions made by RNA if mgc its Molecular Gene Classifier which is older than SCN which is Single Sample Classifier (SCAN-B) and if SGC its single gene classifier
    er_pred_sgc = as.numeric(`er prediction sgc:ch1`),
    pgr_pred_mfc = as.numeric(`pgr prediction mgc:ch1`),
    pgr_pred_sgc = as.numeric(`pgr prediction sgc:ch1`),
    her2_pred_mfc = as.numeric(`her2 prediction mgc:ch1`),
    her2_pred_sgc = as.numeric(`her2 prediction sgc:ch1`),
    ki67_pred_mfc = as.numeric(`ki67 prediction mgc:ch1`),
    ki67_pred_sgc = as.numeric(`ki67 prediction sgc:ch1`),
    nhg_pred_mgc = as.numeric(`nhg prediction mgc:ch1`),
    pam50 = as.factor(`pam50 subtype:ch1`),
    os_months = as.numeric(`overall survival days:ch1`),
    os_status = as.numeric(`overall survival event:ch1`),
    endocrine_tx = as.numeric(`endocrine treated:ch1`),
    chemo_tx = as.numeric(`chemo treated:ch1`)
    
    
  ) %>%
  select(
    -c(
      source_name_ch1,
      characteristics_ch1.2,
      characteristics_ch1.3,
      characteristics_ch1.4,
      characteristics_ch1.5,
      characteristics_ch1.6,
      characteristics_ch1.7,
      characteristics_ch1.8,
      characteristics_ch1.9,
      characteristics_ch1.10,
      characteristics_ch1.11,
      characteristics_ch1.12,
      characteristics_ch1.13,
      characteristics_ch1.14,
      characteristics_ch1.15,
      characteristics_ch1.16,
      characteristics_ch1.17,
      characteristics_ch1.18,
      characteristics_ch1.19,
      characteristics_ch1.20,
      characteristics_ch1.21,
      characteristics_ch1.22,
      characteristics_ch1.23,
      characteristics_ch1.24,
      # Each one of the characteristics_ch1. corresponds to its equivalent in the next lines and both correspond in orther to its characteristic in mutate
      `age at diagnosis:ch1`,
      `tumor size:ch1`,
      `lymph node group:ch1`,
      `lymph node status:ch1`,
      `er status:ch1`,
      `pgr status:ch1`,
      `her2 status:ch1`,
      `ki67 status:ch1`,
      `nhg:ch1`,
      `er prediction mgc:ch1`,
      `er prediction sgc:ch1`,
      `pgr prediction mgc:ch1`,
      `pgr prediction sgc:ch1`,
      `her2 prediction mgc:ch1`,
      `her2 prediction sgc:ch1`,
      `ki67 prediction mgc:ch1`,
      `ki67 prediction sgc:ch1`,
      `nhg prediction mgc:ch1`,
      `pam50 subtype:ch1`,
      `overall survival days:ch1`,
      `overall survival event:ch1`,
      `endocrine treated:ch1`,
      `chemo treated:ch1`
      
    )
  )

counts_data[1:5,1:5]

  
metadata_her2_low <- 
  metadata %>% 
  filter(her2_status == 0)

# 3.- Preprocess data -----------------------------------------------------

# 3.1 Match it with metadata

# 3.1.1 Identify patients in both sets

common_samples <- intersect(colnames(counts_data), metadata$title)

common_samples_her2_low <- intersect(colnames(counts_data), metadata_her2_low$title)

# 3.1.2 Keep the patients in counts data that also have metadata

counts_data <- counts_data[, common_samples]


counts_her2_low <- counts_data[, common_samples_her2_low]

#> //Dictionary//################################################
#> 
#> counts_data <-  Gene counts of the patients with available metadata
#> 
#> metadata_her2_low <- Processed metadata

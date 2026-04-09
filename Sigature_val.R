# Script to perform normalization, annotation and batch correction for PCA and cluster analysis

library(oligo)
library(GEOquery)
library(tidyverse)
library(limma)
library(hgu133plus2.db)
library(gridExtra)
library(readr)



# 1.- Load data and metadata ----------------------------------------------

# 1.1 Get supplementary files
# The use of supp files (.cel) is for future analysis of multiple data bases

# getGEOSuppFiles("GSE2034", baseDir = "D:/tcga")

# 1.1.2 Untar files
#untar("D:/tcga/GSE2034/GSE2034_RAW.tar", exdir = "D:/tcga/GSE2034/")

# 1.1.3 Listing .cel files

cel_files <- list.celfiles("D:/tcga/GSE2034/", full.names = TRUE, listGzipped = TRUE)

# 1.1.4 Reading in cel files

pre_raw_data <- read.celfiles(cel_files)

raw_data <- pre_raw_data

# 1.2 Metadata

pre_metadata <- read.delim("C:/R/METABRIC/Validation/GSE2034/Metadatos/GSE2034_metadata", sep = ",")
metadata <- pre_metadata

#> pData of raw_data initially only contains the ID for the counts and an index
#> meanwhile metadata contains the full metadata but the IDs differ from the count data IDs

# Object that specifies which GSE is in use

gse_obj <- "GSE2034"



# 3.- Preprocessing metadata --------------------------------------------------

#3.1 Clean metadata

metadata <- metadata %>%
  mutate(
    SURVIVAL = ifelse(relapse..1.True. == 1, 1, 0),
    SURVIVAL = as.numeric(SURVIVAL),
    SURVIVAL_MON = time.to.relapse.or.last.follow.up..months.,
    SURVIVAL_MON = as.numeric(SURVIVAL_MON),
    id = GEO.asscession.number
    )

# 3.2 Object with the names of each file

id <- sampleNames(raw_data)
sample_names <- gsub("\\..*", "", id)
# 3.2.2 Add the object to the phenotipic data

pData(raw_data)$id <- sample_names

# 3.4 Join metadata 

pData(raw_data) <- 
  pData(raw_data) %>% 
  rownames_to_column("file_name") %>% 
  full_join(metadata, by = "id", keep = FALSE) %>% 
  mutate(comp_file_name = file_name) %>% 
  column_to_rownames("file_name") %>% 
  mutate(file_name = comp_file_name,
         comp_file_name = NULL)
  
metadata <- pData(raw_data)


# 4.- Preprocess data -----------------------------------------------------


# Normalize
norm_data <- rma(raw_data)

# Boxplot after normalization
boxplot(exprs(norm_data), 
        las = 2, 
        main = paste0("RMA normalized - ", gse_obj))

# 4.3 Expression matrix

expr_matrix <- exprs(norm_data)

# 5.- Probe ID to symbol --------------------------------------------------

# 5.1 Get mapping to change probe IDs to gene names
library(hgu133a.db)
probe_gene <- AnnotationDbi::select(
  hgu133a.db,
  keys = rownames(expr_matrix),
  columns = "SYMBOL",
  keytype = "PROBEID"
)

# 5.2 Delete NA symbols and duplicate symbols and add to rownames

# 5.2.1 Transform matrix to a tidy data frame and add Symbols

gene_expres_matrix <- 
  expr_matrix %>%
  as.data.frame() %>%
  rownames_to_column("PROBEID") %>%
  inner_join(probe_gene, by = "PROBEID") %>%  # Join with your mapping object
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%   # Remove NAs and empty symbols
  
  # 5.2.2 Calculate variance for each probe across all samples
  
  mutate(variance = apply(dplyr::select(., -PROBEID, -SYMBOL), 1, var)) %>%
  
  # 5.2.3 Keep only the probe with the highest variance per Gene Symbol
  
  group_by(SYMBOL) %>%
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% 
  ungroup() %>%
  
  # 5.2.4 Remove columns, add symbol to rownames and reformat to matrix
  
  dplyr::select(-PROBEID, -variance) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()


# 1. Transpose so genes are rows for scaling
# 2. scale() works on columns, so we scale, then transpose back 
# so patients are rows for the prediction
gse2034_scaled <- t(scale(t(gene_expres_matrix))) 

# Convert to a data frame so you can add clinical columns
final_df <- as.data.frame(t(gse2034_scaled)) 

metadata <- metadata %>% 
  filter(ER.Status == "ER+")

final_df <- final_df[rownames(final_df) %in% metadata$file_name,] %>% 
  rownames_to_column("file_name")


library(survival)
# Join with the metadata you cleaned earlier
final_df <- final_df %>%
  left_join(metadata, by = "file_name") %>%
  mutate(
    # Create the Survival Object inside the dataframe
    surv_obj = Surv(time = SURVIVAL_MON, event = SURVIVAL)
  ) %>% 
  column_to_rownames("file_name")



# Create a vector of the genes your model expects
# required_genes <- c("GBP5", "CDCA5", "FGD3") # Add any others that fail
# 
# # Add them to the dataframe as 0 (the mean of scaled data)
# for(gen in required_genes) {
#   if(!(gen %in% colnames(final_df))) {
#     final_df[[gen]] <- 0
#   }
# }


# This is the "testing" phase
gse2034_results <- predict(final_fit, new_data = final_df) %>%
  bind_cols(final_df)


# This is where you get your "P-value" for the test
validation_test <- coxph(Surv(SURVIVAL_MON, SURVIVAL) ~ .pred_time, data = gse2034_results)

summary(validation_test)


library(survminer)

# Create risk groups based on the median of your predictions
gse2034_results <- gse2034_results %>%
  mutate(risk_group = ifelse(.pred_time < median(.pred_time), "High Risk", "Low Risk"))

# Fit the KM curve
km_fit <- survfit(Surv(SURVIVAL_MON, SURVIVAL) ~ risk_group, data = gse2034_results)

# Plot
ggsurvplot(km_fit, 
           data = gse2034_results, 
           pval = TRUE, 
           risk.table = TRUE,
           title = "Validation in GSE2034 (Untreated Cohort)",
           palette = c("#E41A1C", "#377EB8"))


gse2034_results <- gse2034_results %>%
  mutate(pred_z = scale(.pred_time))

# Run Cox again
summary(coxph(Surv(SURVIVAL_MON, SURVIVAL) ~ pred_z, data = gse2034_results))
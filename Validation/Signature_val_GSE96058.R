#> In this script we obtain and preprocess the metadata.gse96058 and the counts data

library(GEOquery)
library(tidyverse)

# 1.- Download data -----------------------------------------------------------


# 1.1 Download the supplementary file (The actual expression matrix)
# getGEOSuppFiles("GSE96058", baseDir = "D:/GSE96058")

# 1.1.2 Read the specific expression file (SCAN-B typically provides a large .txt or .csv)

raw <- read.csv("D:/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz", 
                row.names = 1, check.names = FALSE)
counts_data.gse96058 <- raw
# 1.2 Download metadata.gse96058

gse <- getGEO("GSE96058", GSEMatrix = TRUE)

# 1.2.2 Asign to object

pheno <- pData(gse[[1]])



# 2.- Preprocess metadata.gse96058 -------------------------------------------------

pheno$characteristics_ch1.3 <- gsub("\\D", "", pheno$characteristics_ch1.3)


metadata.gse96058 <-
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
    os_months = as.numeric(`overall survival days:ch1`) / 30.4166667,
    os_status = as.numeric(`overall survival event:ch1`),
    endocrine_tx = as.numeric(`endocrine treated:ch1`),
    chemo_tx = as.numeric(`chemo treated:ch1`)
    
    
  ) %>%
  dplyr::select(
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

counts_data.gse96058[1:5,1:5]


metadata.gse96058_er_pos <-
  metadata.gse96058 %>% 
  filter(er_status == 1) %>% 
  rownames_to_column("id") %>% 
  mutate(EVENT = os_status,
         EVENT_MON = os_months,
         id = NULL) 

# 3.- Preprocess data -----------------------------------------------------

# 3.1 Match it with metadata.gse96058

# 3.1.1 Identify patients in both sets

common_samples <- intersect(colnames(counts_data.gse96058), metadata.gse96058_er_pos$title)

# 3.1.2 Keep the patients in counts data that also have metadata.gse96058

counts_data.gse96058_erpos <- counts_data.gse96058[, common_samples]

counts_data.gse96058_erpos <- scale(t(counts_data.gse96058_erpos))




#> //Dictionary//################################################
#> 
#> counts_data.gse96058 <-  Gene counts of the patients with available metadata.gse96058
#> 
#> metadata.gse96058_her2_low <- Processed metadata.gse96058


# Find genes present in both data sets
common_genes_meta.gse96058 <- intersect(late_death.genes, colnames(counts_data.gse96058_erpos))


# Check how many are lost
print(paste("Original:", length(late_death.genes), "Common:", length(common_genes_meta.gse96058)))

counts_data.gse96058_erpos <- counts_data.gse96058_erpos[ , common_genes_meta.gse96058]

final_df <- 
  counts_data.gse96058_erpos %>% 
  as.data.frame() %>% 
  rownames_to_column("title") %>% 
  left_join(metadata.gse96058_er_pos, by = "title") %>% 
  mutate(surv_obj = Surv(time = EVENT_MON, event = EVENT, type = "right")) %>% 
  dplyr::select(all_of(late_death.genes),
                EVENT,
                EVENT_MON,
                title,
                surv_obj) %>% 
  column_to_rownames("title")

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# 6.- Validation ----------------------------------------------------------


# This is the predict phase

gse96058_results <- predict(final_fit, new_data = final_df, type = "linear_pred") %>%
  bind_cols(final_df)


# To get the p value

validation_test <- coxph(Surv(EVENT_MON, EVENT) ~ .pred_linear_pred, data = gse96058_results)

summary(validation_test)


library(survminer)

# Create risk groups based on the median of the predictions

gse96058_results <- gse96058_results %>%
  mutate(risk_group = as.factor(ifelse(.pred_linear_pred < median(.pred_linear_pred), "High Risk", "Low Risk")))


# Fit the KM curve

km_fit <- survfit(Surv(EVENT_MON, EVENT) ~ risk_group, data = gse96058_results)

# Plot

ggsurvplot(km_fit, 
           data = gse96058_results, 
           pval = TRUE, 
           risk.table = TRUE,
           title = "Validation in GSE2034 (Untreated Cohort)",
           palette = c("#E41A1C", "#377EB8"))


gse96058_results <- gse96058_results %>%
  mutate(pred_z = scale(.pred_linear_pred))

# Run Cox again
gse96058_results$risk_group <- relevel(gse96058_results$risk_group, ref = "Low Risk")

summary_gse96058 <- summary(coxph(Surv(EVENT_MON, EVENT) ~ risk_group, data = gse96058_results))

# Calculate the actual Concordance Index
c_index_results.gse96058 <- concordance(Surv(EVENT_MON, EVENT) ~ .pred_linear_pred, 
                               data = gse96058_results)



library(timeROC)

# Area under the curve per time

res_auc <- timeROC(T = gse96058_results$EVENT_MON,
                   delta = gse96058_results$EVENT,
                   marker = -gse96058_results$pred_z,
                   cause = 1, # The event code
                   times = c(36, 60, 80), # 3, 5, and 10 years
                   iid = TRUE)

# View the AUC values

res_auc_gse96058 <- res_auc$AUC %>% 
  as.data.frame()


# Cox multivariado con clinica

late_genes.patients_gse96058.cox <- 
  final_df %>% 
  as.data.frame() %>% 
  rownames_to_column("title") %>% 
  left_join(metadata.gse96058_er_pos, by = "title", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y")) %>% 
  column_to_rownames("title") %>% 
  mutate(HER2 = her2_pred_sgc, 
         LYMPH = lymph_group,
         PAM50 = pam50,
         AGE = as.numeric(age),
         KI67 = ki67_pred_sgc,
         SCORE = gse96058_results$.pred_linear_pred * -1
  ) %>% 
  dplyr::select(all_of(late_death.genes),
                surv_obj,
                LYMPH,
                PAM50,
                AGE,
                HER2,
                KI67,
                SCORE
  ) %>% 
  na.omit()

independent_prog.gse96058 <- coxph(surv_obj ~ PAM50 + KI67 + HER2 + AGE + LYMPH + SCORE, 
                                   data = late_genes.patients_gse96058.cox) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

text_validation.gse96058 <- paste("Al dividir en grupos de alto y bajo riesgo, esta firma en GSE96058 consiguio un HR de ",
      round(summary_gse96058$coefficients[2], 3),
      " (IC 95% de ",
      round(summary_gse96058$conf.int[3], 3),
      " - ",
      round(summary_gse96058$conf.int[4], 3),
      ", pval de ",
      summary_gse96058$coefficients[5],
      "), C score de ",
      round(c_index_results.gse96058$concordance, 2),
      ", área bajo la curva a los 3 años de ",
      round(res_auc_gse96058[1,1], 3),
      ", a los 5 años de ",
      round(res_auc_gse96058[2,1], 3),
      ", y a los 6 años de",
      round(res_auc_gse96058[3,1], 3)
)

cat(text_signature, text_validation.gse96058, sep = ". ")





lymph_rows <- independent_prog.gse96058[grepl("LYMPH", independent_prog.gse96058$term),]

best_lymph <- lymph_rows[which.min(lymph_rows$p.value),]

significance.gse96058 <- if((independent_prog.gse96058[independent_prog.gse96058$term == "SCORE",]$p.value < 0.05) == TRUE){
  "se mantuvo como un predictor de la supervivencia independiente significativo "
}else{
  "no se mantuvo como un predictor de la supervivencia independiente significativo " 
}

text_independent_prog.gse96058 <- paste0(
  "Con esta base de datos, al realizar un cox multivariado junto a edad, ganglios linfaticos, fenotipos de PAM50, y KI67, la firma ",
  significance.gse96058,
  "obteniendo un HR de ",
  round(independent_prog.gse96058$estimate[independent_prog.gse96058$term == "SCORE"], 3),
  " (IC 95% de ",
  round(independent_prog.gse96058$conf.low[independent_prog.gse96058$term == "SCORE"], 2),
  " - ",
  round(independent_prog.gse96058$conf.high[independent_prog.gse96058$term == "SCORE"], 2),
  " pvalue de ",
  independent_prog.gse96058$p.value[independent_prog.gse96058$term == "SCORE"],
  ")",
  if((independent_prog.gse96058[independent_prog.gse96058$term == "SCORE",]$p.value < best_lymph$p.value) == TRUE){
    paste0(" superando a los ganglios linfaticos como predictor (HR de ",
           round(best_lymph$estimate, 3),
           " pval de ",
           best_lymph$p.value,
           ")")
  }else{
    paste0(" sin lograr superar a los ganglios linfaticos como predictor (HR de ",
           round(best_lymph$estimate, 3),
           " pval de ",
           best_lymph$p.value,
           ")")
  },
  if((independent_prog.gse96058[independent_prog.gse96058$term == "SCORE",]$p.value < independent_prog.gse96058[independent_prog.gse96058$term == "KI67",]$p.value) == TRUE & (independent_prog.gse96058[independent_prog.gse96058$term == "SCORE",]$p.value < best_lymph$p.value) == TRUE){
    paste0(" e igualmente superando a KI67 como predictor (HR de ",
           round(independent_prog.gse96058$estimate[independent_prog.gse96058$term == "KI67"], 3),
           " pval de ",
           independent_prog.gse96058$p.value[independent_prog.gse96058$term == "KI67"],
           ")")
  }else if((independent_prog.gse96058[independent_prog.gse96058$term == "SCORE",]$p.value < independent_prog.gse96058[independent_prog.gse96058$term == "KI67",]$p.value) == FALSE & (independent_prog.gse96058[independent_prog.gse96058$term == "SCORE",]$p.value < best_lymph$p.value) == TRUE){
    paste0(" sin lograr superar a KI67 como predictor (HR de ",
           round(independent_prog.gse96058$estimate[independent_prog.gse96058$term == "KI67"], 3),
           " pval de ",
           independent_prog.gse96058$p.value[independent_prog.gse96058$term == "KI67"],
           ").")
  }
  else if((independent_prog.gse96058[independent_prog.gse96058$term == "SCORE",]$p.value < independent_prog.gse96058[independent_prog.gse96058$term == "KI67",]$p.value) == TRUE & (independent_prog.gse96058[independent_prog.gse96058$term == "SCORE",]$p.value < best_lymph$p.value) == FALSE){
  paste0(" pero si superando a KI67 como predictor (HR de ",
         round(independent_prog.gse96058$estimate[independent_prog.gse96058$term == "KI67"], 3),
         " pval de ",
         independent_prog.gse96058$p.value[independent_prog.gse96058$term == "KI67"],
         ")")
  }
)

text_gse96058 <- paste(text_validation.gse96058, text_independent_prog.gse96058, sep = ". ")

cat(text_metabric, text_gse96058, sep = "\n")

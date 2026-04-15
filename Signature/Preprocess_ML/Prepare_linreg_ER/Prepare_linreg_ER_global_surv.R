library(survival)


# In this file we prepare the data for the linear regression model using onl ER+ patients
# and with survival parameters

label <- "Para predecir supervivencia en pacientes ER positivo de la base de METABRIC "

# 1.- Preparing metadata --------------------------------------------------

er_patients <- alive_brca.death %>% 
  filter(ER_IHC == "Positve")

lm_metadata <- er_patients

# 1.2 List of genes to use (check dictionary below to understand the different variables that are used)

late_death.genes <-  common_genes_meta.gse96058 # boruta_signature #significant_genes$term #rownames(res_sig) common_genes_meta.gse96058


# 1.3 Object with all the patients ER+ and expression of only the genes of interest
rownames(counts_data) <- make.names(rownames(counts_data))

late_genes.patients <- counts_data[late_death.genes, er_patients$PATIENT_ID]

# 1.4 Transpose only since scaling is done in the lin regression recipe

late_genes.patients <- t(late_genes.patients) 

# 1.5.1 Check that the patients are in the same order

all(rownames(late_genes.patients) == er_patients$PATIENT_ID)

# 1.5.2 Add a column of EVENT as a binary term for it to be the outcome and the months of survival

late_genes.patients <- 
  late_genes.patients %>% 
  as.data.frame() %>% 
  rownames_to_column("PATIENT_ID") %>% 
  left_join(er_patients, by = "PATIENT_ID") %>% 
  column_to_rownames("PATIENT_ID") %>% 
  mutate(EVENT = as.numeric(SURVIVAL_STAT),
         EVENT_MON = as.numeric(OS_MONTHS)
  ) %>%  # Turn to factor for machine learning
  dplyr::select(all_of(late_death.genes),
         EVENT_MON,
         EVENT) %>% 
  filter(EVENT_MON > 0) %>% # Eliminate those with 0 survival months
  drop_na() %>% 
    mutate(surv_obj = Surv(
      time  = EVENT_MON,
      event = EVENT,
      type  = "right"))




# /Dictionary/ ##########################
#>  VARIABLES FOR 1.2 late_death.genes
#>  
#> rownames(res_sig) <- List of genes from differential expression
#> significant_genes$term <- Contains the genes determined by the cox analysis as significant so as to be used as signature input in ML models


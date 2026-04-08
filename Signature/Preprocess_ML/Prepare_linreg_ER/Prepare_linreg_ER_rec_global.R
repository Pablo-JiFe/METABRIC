# In this file we prepare the data for the linear regression model using only ER+ patients
# and with recurrence parameters


# 1.- Preparing metadata --------------------------------------------------

er_patients_recu <- metadata %>% 
  filter(ER_IHC == "Positve")

# 1.2 List of genes to use (check dictionary below to understand the different variables that are used)

late_death.genes <-  significant_genes$term #rownames(res_sig)
significant_genes$term %>% 
  cat(sep = ", ")


# 1.3 Object with all the patients ER + and expression of only the genes of interest

late_genes.patients <- counts_data[late_death.genes, er_patients_recu$PATIENT_ID]

# 1.4  Scaling is done in the linear regression recipe


late_genes.patients <- t(late_genes.patients) 

# 1.5.1 Check that the patients are in the same order


all(rownames(late_genes.patients) == er_patients_recu$PATIENT_ID)

# 1.5.2 Add a column called SURVIVAL and SURVIVAL_MON to create the surv object
# NOTE that this file is recurrence, it is still stored in survival so as to not have to change the 
# linear regression file

late_genes.patients <- 
  late_genes.patients %>% 
  as.data.frame() %>% 
  mutate(SURVIVAL = er_patients_recu$RECURR_STAT, 
         SURVIVAL = as.numeric(SURVIVAL),
         SURVIVAL_MON = er_patients_recu$RFS_MONTHS,
         SURVIVAL_MON = as.numeric(SURVIVAL_MON)
  ) %>%  # Turn to factor for machine learning
  filter(SURVIVAL_MON > 0) %>% 
  drop_na()  

library(survival)

late_genes.patients$surv_obj <- Surv(
  time  = late_genes.patients$SURVIVAL_MON,
  event = late_genes.patients$SURVIVAL,
  type  = "right"
)

# /Dictionary/ ##########################
#>  VARIABLES FOR 1.2 late_death.genes
#>  
#> rownames(res_sig) <- List of genes from differential expression
#> significant_genes$term <- Contains the genes determined by the cox analysis as significant so as to be used as signature input in ML models

# In this file we prepare the data for the linear regression model using onl ER+ patients
# and with survival parameters


# 1.- Preparing metadata --------------------------------------------------

er_patients <- alive_brca.death %>% 
  filter(ER_IHC == "Positve")

# 1.2 List of genes to use (check dictionary below to understand the different variables that are used)

late_death.genes <-  significant_genes$term #rownames(res_sig)
significant_genes$term %>% 
  cat(sep = ", ")


# 1.3 Object with all the patients ER+ and expression of only the genes of interest

late_genes.patients <- counts_data[late_death.genes, er_patients$PATIENT_ID]

# 1.4 Transpose only since scaling is done in the lin regression recipe

late_genes.patients <- t(late_genes.patients) 

# 1.5.1 Check that the patients are in the same order

all(rownames(late_genes.patients) == er_patients$PATIENT_ID)

# 1.5.2 Add a column of SURVIVAL as a binary term for it to be the outcome and the months of survival

late_genes.patients <- 
  late_genes.patients %>% 
  as.data.frame() %>% 
  mutate(SURVIVAL = er_patients$SURVIVAL_STAT, 
         SURVIVAL = as.numeric(SURVIVAL),
         SURVIVAL_MON = er_patients$OS_MONTHS,
         SURVIVAL_MON = as.numeric(SURVIVAL_MON)
  ) %>%  # Turn to factor for machine learning
  filter(SURVIVAL_MON > 0) %>% # Eliminate those with 0 survival months
  drop_na() 

library(survival)

# 2.- Create surv object for the different ML models as outcome

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


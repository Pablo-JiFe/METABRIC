# In this file we prepare the data for the linear regression model global

# 1.- Preparing metadata --------------------------------------------------

# 1.2 List of genes to use (check dictionary below to understand the different variables that are used)

late_death.genes <-  significant_genes$term #rownames(res_sig)


# 1.3 Object with the patients that are alive or died because of breast cancer and expression of only the genes of interest

late_genes.patients <- counts_data[late_death.genes, alive_brca.death$PATIENT_ID]

# 1.4 Scaling is done in the linear regresison recipe

late_genes.patients <- t(late_genes.patients) 

# 1.5.1 Check that the patients are in the same order

all(rownames(late_genes.patients) == alive_brca.death$PATIENT_ID)

# 1.5.2 Add a column of SURVIVAL as a binary term for it to be the outcome and the months of survival

late_genes.patients <- 
  late_genes.patients %>% 
  as.data.frame() %>% 
  mutate(SURVIVAL = alive_brca.death$SURVIVAL_STAT, 
         SURVIVAL = as.numeric(SURVIVAL),
         SURVIVAL_MON = alive_brca.death$OS_MONTHS,
         SURVIVAL_MON = as.numeric(SURVIVAL_MON)
  ) %>%  # Turn to factor for machine learning
  filter(SURVIVAL_MON > 0) %>% 
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


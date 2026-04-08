

library(survival)
library(survminer)
library(broom)
library(survcomp)

# In this file we continue with feature selectiuon using cox survival analysis


# 1.- Preparing metadata --------------------------------------------------


# 1.1 Object with live patients and patients diseased by breast cancer

alive_brca.death <- metadata %>% 
  filter(VITAL_STATUS != "Died of Other Causes")

# 1.2 List of genes to use (check dictionary below to understand the different variables that are used)

late_death.genes <-  # rownames(res_sig) # coef_tbl$term


# 1.3 Object with the patients that are either alive or died from breast cancer and expression of only the genes of interest

late_genes.patients <- counts_data[late_death.genes, alive_brca.death$PATIENT_ID]

late_genes.patients <- 
  late_genes.patients %>% 
  drop_na()

# 1.4 Scale so that they are comparable

late_genes.patients <- scale(t(late_genes.patients)) 

# 1.5.1 Check that the patients are in the same order

all(rownames(late_genes.patients) == alive_brca.death$PATIENT_ID)

# 1.5.2 Add a column of SURVIVAL as a binary term for it to be the outcome

late_genes.patients <- 
  late_genes.patients %>% 
  as.data.frame() %>% 
  mutate(SURVIVAL = alive_brca.death$SURVIVAL_STAT, 
         SURVIVAL = as.numeric(SURVIVAL),
         SURVIVAL_MON = alive_brca.death$OS_MONTHS, 
         SURVIVAL_MON = as.numeric(SURVIVAL_MON)
  ) 



# 2.- Survival ------------------------------------------------------------


# 2.1 Survival object

surv <- Surv(late_genes.patients$SURVIVAL_MON, late_genes.patients$SURVIVAL)


# 2.3 Cox survival

fit <- coxph(surv ~ . - SURVIVAL_MON - SURVIVAL, data = late_genes.patients)

# 2.4 Risk scores

risk_scores <- predict(fit, type = "risk")

# 2.5 C score

c_index_val <- concordance.index(x = risk_scores, 
                                 surv.time = alive_brca.death$OS_MONTHS, 
                                 surv.event = alive_brca.death$SURVIVAL_STAT, 
                                 method = "noether")

print(c_index_val$c.index)


# 3.- Preparing for input into machine learning models --------------------

# 3.1 Convert the model object to a tidy data frame
tidy_results <- tidy(fit)

# 3.2 Filter for significant genes
significant_genes <- tidy_results %>%
  filter(p.value < 0.05) %>%
  arrange(p.value) # Optional: sort by significance

# 3.3 View results
print(significant_genes)
significant_genes$term %>% 
  cat(sep = ", ")

# /Dictionary/ ##########################
#>  VARIABLES FOR 1.2 late_death.genes
#>  
#> rownames(res_sig) <- List of genes from differential expression
#> coef_tbl$term <-  List of genes obained from penalized linear regression
#> 
#> 
#> VARIABLES FOR MACHINE LEARNING MODELS
#> significant_genes$term <- Contains the genes determined by the cox analysis as significant so as to be used as signature input in ML models

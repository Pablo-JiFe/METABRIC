# In this file we prepare the data for the logistical regression model


# 1.- Preparing metadata --------------------------------------------------

er_patients_recu <- metadata %>% 
  filter(ER_IHC == "Positve")

# 1.2 List of genes described in the differential expression as being of prognosis for late death

late_death.genes <-  significant_genes$term #[res_sig$logFC > 0]#rownames(res_sig)
significant_genes$term %>% 
  cat(sep = ", ")


# 1.3 Object with all the patients and expression of only the genes of interest

late_genes.patients <- counts_data[late_death.genes, er_patients_recu$PATIENT_ID]

# 1.4 Scale so that they are comparable

#late_genes.patients <- t(late_genes.patients)

late_genes.patients <- t(late_genes.patients) 

# 1.5.1 Check that the patients are in the same order


all(rownames(late_genes.patients) == er_patients_recu$PATIENT_ID)

# 1.5.2 Add a column of SURVIVAL as a binary term for it to be the outcome

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


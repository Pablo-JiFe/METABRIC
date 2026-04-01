# In this file we prepare the data for the logistical regression model


# 1.- Preparing metadata --------------------------------------------------


# 1.1 Object with live patients and patients diseased by breast canbcer

alive_brca.death <- metadata[metadata$VITAL_STATUS != "Died of Other Causes" ,]

# 1.2 List of genes described in the differential expression as being of prognosis for late death

late_death.genes <- rownames(res_sig)[res_sig$logFC > 0]


# 1.3 Object with all the patients and expression of only the genes of interest

late_genes.patients <- counts_data[late_death.genes, alive_brca.death$PATIENT_ID]
  
# 1.4 Scale so that they are comparable

#late_genes.patients <- t(late_genes.patients)

late_genes.patients <- scale(t(late_genes.patients))

# 1.5.1 Check that the patients are in the same order

alive_brca.death <- alive_brca.death[match(rownames(late_genes.patients),alive_brca.death$PATIENT_ID),]

all(rownames(late_genes.patients) == alive_brca.death$PATIENT_ID)

# 1.5.2 Add a column of SURVIVAL as a binary term for it to be the outcome

late_genes.patients <- 
  late_genes.patients %>% 
  as.data.frame() %>% 
  mutate(SURVIVAL = alive_brca.death$SURVIVAL_STAT, 
         SURVIVAL = as.factor(SURVIVAL)) # Turn to factor for machine learning

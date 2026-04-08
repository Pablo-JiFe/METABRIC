library(limma)
# 1.- Dividing by lymph nodes ------------------------------
library(pheatmap)

# 1.- Dividing by median of time until death ------------------------------

# 1.1 Dividing time until death form median

# 1.1.2 Obtain median of overall survival in months

median_disease.brca <- median(metadata_diseased.brca$OS_MONTHS)

# 1.2 Generate column corresponding to time until death if it was earlier or later than the median

metadata_diseased.brca <- metadata_diseased.brca %>% 
  mutate(BOOLEAN_TEMP = ifelse(OS_MONTHS <= median_disease.brca,
                               yes = "Early",
                               no = "Late"))


# 1.1 Generate column corresponding to lymph nodes, those that have 0 in one group and those with more than 0 in another

col_data <- metadata %>%
  mutate(LYMPH = ifelse(LYMPH_NODES_EXAMINED_POSITIVE == 0, 0, 1),
         LYMPH = as.factor(LYMPH)) %>%
  drop_na() %>%
  select(PATIENT_ID, LYMPH) %>% # Create only the object to use for Limma
  column_to_rownames("PATIENT_ID")


# 2.- Differential expression -----------------------------------------------


# 2.2 Data counts of the patients that had lymph node information in the metadata

count_data <- counts_data[colnames(counts_data) %in% rownames(col_data)]

# 2.2.2 Making sure they are in the same order

count_data <- count_data[match(rownames(col_data), colnames(count_data))]

all(colnames(count_data) == rownames(col_data))


# 2.3 Generate limma object

library(limma)

# 2.4 Design based on object separating on lymph nodes

design <- model.matrix(~ 0 + LYMPH, data = col_data)

# 2.4.2 Asign make.names objects as colnames

colnames(design) <- make.names(colnames(design)) 

# 2.5 Fit

fit <- lmFit(count_data, design)

# 2.5.2 Contrast matrix comparing lymph 0 to > 0 lymph

contrast.matrix <- makeContrasts(LYMPH0 - LYMPH1,
                                 levels = design)

# 2.5.3 Fit based on contrasts

fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

topTable(fit)

# 2.6 Results

res <- topTable(fit, coef = 1, number = Inf)

# 2.6.2 Results that correspond to a signfiicant p value and log fold change

res_sig <- res %>%
  filter(adj.P.Val < 0.01 & abs(logFC) > 0.2) # 0.1
res_sig 


library(pheatmap)

# 1.- Dividing by median of time until death ------------------------------

HAY QUE ELIMINAR NA!!!!

# 1.1 Dividing time until death form median

# 1.1.2 Obtain median of overall surbival in months

median_disease.brca <- median(metadata_diseased.brca$OS_MONTHS)

# 1.2 Generate column corresponding to time until death if it was earlier or later than the median

metadata_diseased.brca <- metadata_diseased.brca %>% 
  mutate(BOOLEAN_TEMP = ifelse(OS_MONTHS <= median_disease.brca,
                               yes = "Early",
                               no = "Late"))





# 2.- Differential expression -----------------------------------------------

#> 2.1 Data frame where row names correspond to the names of the patients and adittionaly contains the column
#> previously established as being correspondant to early or late death with respect to the median

col_data <- data.frame(metadata_diseased.brca$BOOLEAN_TEMP, row.names = metadata_diseased.brca$PATIENT_ID)

# 2.2 Data counts of the patients that are either alive or died because of breast cancer

count_data <- counts_data[colnames(counts_data) %in% metadata_diseased.brca$PATIENT_ID]

# 2.2.2 Making shure they are in the same order

count_data <- count_data[match(rownames(col_data), colnames(count_data))]

all(colnames(count_data) == rownames(col_data))

# 2.2.3 Changing the name

col_data <-
  col_data %>% 
  mutate(BOOLEAN_TEMP = metadata_diseased.brca.BOOLEAN_TEMP,
         .keep = "unused")


# 2.3 Generate limma object


library(limma)

# 2.4 Design based on object separating on alive and diseased

design <- model.matrix(~ 0 + BOOLEAN_TEMP, data = col_data)

# 2.4.2 Asign make.names objects as colnames

colnames(design) <- make.names(colnames(design)) 

# 2.5 Fit

fit <- lmFit(count_data, design)

# 2.5.2 Contrast matrix comparing late vs early

contrast.matrix <- makeContrasts(BOOLEAN_TEMPLate - BOOLEAN_TEMPEarly,
                                 levels = design)

# 2.5.3 Fit based on contrasts

fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

topTable(fit)

# 2.6 Results

res <- topTable(fit, coef = 1, number = Inf)

# 2.6.2 Results that correspond to a signfiicant p value and log fold change

res_sig <- res %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)


# 2.7 QQ plot

qqt(
  fit$t,
  df = fit$df.prior + fit$df.residual,
  pch = 16,
  cex = 0.2
)
abline(-0.1, 2.8, col = "red", lwd = 2)

# 2.8 HEATMAP

#2.8.1 Vector with the names of significant genes
genes_sig <- rownames(res_sig)

#2.8.2 Expression submatrix with only significant genes
expr_sig <- count_data[rownames(count_data) %in% genes_sig, ]
expr_sig <- expr_sig[genes_sig, ]

#2.8.3 Create column annotations
annotation_col <- data.frame(Group = col_data[colnames(expr_sig), "BOOLEAN_TEMP"])
rownames(annotation_col) <- colnames(expr_sig)

#2.8.4 Plot heatmap
pheatmap(
  expr_sig,
  scale = "row",
  color = colorRampPalette(c("darkblue", "lightblue", "#FFF7FB", "#F9C5D5", "deeppink", "darkviolet"))(100),
  annotation_col = annotation_col,
  annotation_colors = list(
    Group = c(Early = "#F15BB5",Late = "purple")),
  show_rownames = TRUE,
  show_colnames = FALSE,
  clustering_method = "complete",
  main = "Differentially Expressed Genes: Early vs Late"
)



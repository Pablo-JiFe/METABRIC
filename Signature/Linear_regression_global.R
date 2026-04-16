
late_genes.patients.cox <- late_genes.patients[rownames(test_data),]

late_genes.patients.cox <- 
  late_genes.patients.cox %>% 
  as.data.frame() %>% 
  rownames_to_column("PATIENT_ID") %>% 
  left_join(lm_metadata, by = "PATIENT_ID", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y")) %>% 
  column_to_rownames("PATIENT_ID") %>% 
  mutate(HER2 = HER2_SNP6, 
         LYMPH = LYMPH_NODES_EXAMINED_POSITIVE,
         AGE = as.numeric(AGE_AT_DIAGNOSIS ),
         MENO = INFERRED_MENOPAUSAL_STATE    ,
         SCORE = test_pred$.pred_linear_pred,
         HORMONE = HORMONE_THERAPY ,
         CHEMO = CHEMOTHERAPY,
         SURGERY = BREAST_SURGERY
  ) %>% 
  dplyr::select(all_of(late_death.genes),
                surv_obj,
                AGE,
                LYMPH,
                HER2,
                MENO,
                SCORE,
                HORMONE,
                CHEMO,
                SURGERY) %>%  
  na.omit()

independent_prog <- coxph(surv_obj ~ HORMONE + CHEMO + SURGERY + MENO + HER2 + AGE + LYMPH + SCORE, 
                          data = late_genes.patients.cox) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE)


#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################


input_genes <- paste0(late_death.genes, collapse = ", ")

genes_terms <- paste0(coef_tbl$term, collapse = ", ")

input <- paste0(
  "Input: ",
  input_genes,
  ". ")

text_signature <- paste0(
  label,
  "se evaluó una firma con regresión lineal penalizada con un alpha de ",
  best_params$mixture,
  " un lambda de ",
  best_params$penalty,
  " y se obtuvo la firma con "
  ,
  length(coef_tbl$term),
  " genes ",
  " (",
  genes_terms,
  ")",
  " que al estratificar en alto y bajo riesgo se obtuvo un HR de ",
  round(summary_cox$coefficients[2], 3),
  " (IC 95% de ",
  round(summary_cox$conf.int[3], 3),
  " - ",
  round(summary_cox$conf.int[4], 3),
  ", pval de ",
  summary_cox$coefficients[5],
  "), C score de ",
  round(concordancia$concordance, 2),
  ", área bajo la curva a los 3 años de ",
  round(auc[1, 1], 3),
  ", a los 5 años de ",
  round(auc[2, 1], 3),
  ", y a los 10 años de ",
  round(auc[3, 1], 3)
)

cat(
  input ,
text_signature,
sep = "\n"
)




significance <- if ((independent_prog[independent_prog$term == "SCORE", ]$p.value < 0.05) == TRUE) {
  "se mantuvo como un predictor de la supervivencia independiente significativo "
} else{
  "no se mantuvo como un predictor de supervivencia independiente significativo "
}

text_independent_prog <- paste0(
  "Al realizar un cox multivariado junto a edad, tratamiento, ganglios linfáticos  y menopausia, la firma ",
  significance,
  "obteniendo un HR de ",
  round(independent_prog$estimate[independent_prog$term == "SCORE"], 3),
  " (IC 95% de ",
  round(independent_prog$conf.low[independent_prog$term == "SCORE"], 2),
  " - ",
  round(independent_prog$conf.high[independent_prog$term == "SCORE"], 2),
  " pvalue de ",
  independent_prog$p.value[independent_prog$term == "SCORE"],
  ")",
  if ((independent_prog[independent_prog$term == "SCORE", ]$p.value < independent_prog[independent_prog$term == "LYMPH", ]$p.value) == TRUE) {
    paste0(
      " superando a los ganglios linfáticos  como predictor (HR de ",
      round(independent_prog$estimate[independent_prog$term == "LYMPH"], 3),
      " pval de ",
      independent_prog$p.value[independent_prog$term == "LYMPH"],
      ")"
    )
  } else{
    paste0(
      " sin lograr superar a los ganglios linfaticos como predictor (HR de ",
      round(independent_prog$estimate[independent_prog$term == "LYMPH"], 3),
      " pval de ",
      independent_prog$p.value[independent_prog$term == "LYMPH"],
      ")."
    )
  }
)

text_metabric <- paste(text_signature, text_independent_prog, sep = ". ")

library(tidymodels)
library(readr)       # for importing data
library(vip)
library(censored)

# Initial object late_genes.patients created in preprocessing
# Its made so that the modifications have to be done on preprocessing so even if the
# thing to be studied is recurrence the object will stay as EVENT and EVENT_MON


# 1.- Split ---------------------------------------------------------------


set.seed(123)

# 1.1 Split

split <- initial_split(
  late_genes.patients,
  prop = 0.8,
  strata = EVENT   # maintain proportions
)

# > split
# <Training/Testing/Total>
#   <1186/298/1484>

# 1.2 Asign to objects for training and testing

train_data <- training(split)
test_data  <- testing(split)


# 2.- Preparing recipe and model ------------------------------------------

# 2.1 Recipe

lr_rec <- recipe(surv_obj ~ ., data = train_data) %>% # Survival object created in preprocessing for linear regressions
  update_role(EVENT_MON, EVENT, new_role = "non_predictor") %>%
  step_dummy(all_nominal_predictors()) %>% 
  step_zv(all_predictors()) %>% # Eliminates variables with a single value
  step_nzv(all_predictors()) %>% # Eliminates higly sparsed variables
  step_impute_mean(all_predictors()) %>% # Imputes NAs to mean of those variables
  step_normalize(all_predictors())  # Normalize all



# 2.2 Model

lr_mod <- proportional_hazards(
  penalty = tune(),    # lambda establishes the severith of the penalty
  mixture = tune()     # alpha establishes the type, 1 being lasso, 0 being ridge, and 0.5 being elasticnet
) %>%
  set_engine("glmnet") # Engine that permits penalizing by elasticnet, ridge, and lasso

# 2.3 Workflow

lr_wf <- workflow() %>%
  add_model(lr_mod) %>%
  add_recipe(lr_rec)


# 3.- Selecting best penalizing values ------------------------------------


set.seed(123)

# 3.1 Folds for evaluating with resamples

folds <- vfold_cv(
  train_data,
  v = 10,
  strata = EVENT
)

# 3.2 Grid for range of penalizing

grid <- grid_regular(
  penalty(range = c(-4, 1)),   
  mixture(range = c(0, 1)),
  levels = 10
)

# 3.3 Running the different penalization methods

res_ml <- tune_grid(
  lr_wf,
  resamples = folds,
  grid = grid,
  metrics = metric_set(concordance_survival), # Evaluates the different penalizing methods by c score
  control = control_grid(save_pred = TRUE)
)

# 3.3.1 Observe metrics

collect_metrics(res_ml)

# 3.3.2 Object with best parameters for penalizing

best_params <- select_best(res_ml, metric = "concordance_survival")



# 4.- Actual training -----------------------------------------------------

# 4.1 Final workflow with the best parameter selected previously

final_wf <- finalize_workflow(lr_wf, best_params)

# 4.2 Final fit with training data

final_fit <- fit(final_wf, data = train_data)

# 4.2.2 Observing genes that are maintained after penalziation

library(broom)

coef_tbl <- tidy(final_fit) %>%
  filter(estimate != 0) %>%
  arrange(desc(abs(estimate)))

cat(coef_tbl$term, sep = ", ")


# 5.- Testing -------------------------------------------------------------

# 5.1 Predictions on test data

test_pred <- predict(final_fit, new_data = test_data, type = "linear_pred")

# 5.2 Creating groups for kaplan mayer curves

# 5.2.1 Creating column on test data with its prediction

test_data$risk_score <- test_pred$.pred_linear_pred


# 5.2.2 Dividing the groups by median so as to establish a high and low risk and create a column

test_data$risk_group <- as.factor(ifelse(
  test_data$risk_score < median(test_data$risk_score),
  "High",
  "Low"
))

test_data$risk_group <- relevel(test_data$risk_group, ref = "Low")

# 5.3 Creating kapan meier curve

library(survival)
library(survminer)

# 5.3.1 Based on the survival object compare the risk_groups created previously

fit_km <- survfit(surv_obj ~ risk_group, data = test_data)

# 5.3.2 Plot

ggsurvplot(fit_km, data = test_data, title = "Recurrence ER+  chemotherapy", ylab = "Recurrence probability")

# 5.4 Tables

survdiff(surv_obj ~ risk_group, data = test_data)

# 5.4.1 Cox analysis table with HR

summary_cox <- summary(coxph(surv_obj ~ risk_group, data = test_data))

# 6.- C score

# 6.1 Combine the test and training data to compare concordance

eval_df <- dplyr::bind_cols(
  test_data,
  test_pred
)

# 6.2 Compare concordance

concordancia <- concordance(
  surv_obj ~ .pred_linear_pred,
  data = eval_df
)



# 7.- AUC per time --------------------------------------------------------


library(timeROC)

# 7.1 Creating object on testing data

time_roc <- timeROC(
  T = test_data$EVENT_MON,
  delta = test_data$EVENT,
  marker = -test_pred$.pred_linear_pred,
cause = 1,
  times = c(36, 60, 72, 120),  # 3y, 5y, 10y
  iid = TRUE
)

auc <- time_roc$AUC %>% 
  as.data.frame()


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

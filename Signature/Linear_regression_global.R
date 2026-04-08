library(tidymodels)
library(readr)       # for importing data
library(vip)
library(censored)

# Initial object late_genes.patients created in preprocessing
# Its made so that the modifications have to be done on preprocessing so even if the
# thing to be studied is recurrence the object will stay as SURVIVAL and SURVIVAL_MON


# 1.- Split ---------------------------------------------------------------


set.seed(123)

# 1.1 Split

split <- initial_split(
  late_genes.patients,
  prop = 0.8,
  strata = SURVIVAL   # maintain proportions
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
  step_rm(SURVIVAL, SURVIVAL_MON) %>% # Eliminate the columns of the outcomes since the outcome is the surv object
  step_zv(all_predictors()) %>% # Eliminates variables with a single value
  step_nzv(all_predictors()) %>% # Eliminates higly sparsed variables
  step_impute_mean(all_predictors()) %>% # Imputes NAs to mean of those variables
  step_normalize(all_predictors()) # Normalize all


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
  v = 5,
  strata = SURVIVAL
)

# 3.2 Grid for range of penalizing

grid <- grid_regular(
  penalty(range = c(-4, 0)),   
  mixture(range = c(0, 1)),
  levels = 5
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

test_data$risk_group <- ifelse(
  test_data$risk_score > median(test_data$risk_score),
  "High",
  "Low"
)

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

coxph(surv_obj ~ risk_group, data = test_data)

# 6.- C score

# 6.1 Combine the test and training data to compare concordance

eval_df <- dplyr::bind_cols(
  test_data,
  test_pred
)

# 6.2 Compare concordance

concordance(
  surv_obj ~ .pred_linear_pred,
  data = eval_df
)



# 7.- AUC per time --------------------------------------------------------


library(timeROC)

# 7.1 Creating object on testing data

time_roc <- timeROC(
  T = test_data$SURVIVAL_MON,
  delta = test_data$SURVIVAL,
  marker = -test_pred$.pred_linear_pred,
  cause = 1,
  times = c(36, 60, 120),  # 3y, 5y, 10y
  iid = TRUE
)

time_roc$AUC 

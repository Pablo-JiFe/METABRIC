library(tidymodels)
library(readr)       # for importing data
library(vip)


# Data split

set.seed(123)

data_split <- initial_split(late_genes.patients, strata = SURVIVAL)

train_data <- training(data_split)
test_data  <- testing(data_split)

# Folds split

folds <- vfold_cv(train_data, v = 5, strata = SURVIVAL)


# Model specification

lr_mod <- logistic_reg(
  penalty = tune(),
  mixture = 0
) %>%
  set_engine("glmnet")

# Recipe

lr_recipe <- recipe(SURVIVAL ~ ., data = train_data)

# Workflow

lr_workflow <- workflow() %>%
  add_model(lr_mod) %>%
  add_recipe(lr_recipe)

# Grid for tuning

lr_reg_grid <- grid_space_filling(
  penalty(),
  size = 30
)


# Hyperparameter tuning

lr_res <- tune_grid(
  lr_workflow,
  resamples = folds,
  grid = lr_reg_grid,
  metrics = metric_set(roc_auc),
  control = control_grid(save_pred = TRUE)  
)


# Plot tuning results

lr_plot <- lr_res %>%
  collect_metrics() %>%
  ggplot(aes(x = penalty, y = mean)) +
  geom_point() +
  geom_line() +
  ylab("Area under the ROC Curve") +
  scale_x_log10(labels = scales::label_number())

lr_plot


# Best models


top_models <- lr_res %>%
  show_best(metric = "roc_auc", n = 25) %>%
  arrange(penalty)

top_models


# Select best hyperparameter


lr_best <- lr_res %>%
  select_best(metric = "roc_auc")

lr_best


# ROC curve from validation results


lr_auc <- lr_res %>%
  collect_predictions(parameters = lr_best) %>%
  roc_curve(SURVIVAL, .pred_1) %>%   
  mutate(model = "Logistic Regression")

autoplot(lr_auc)


# Final model 

set.seed(345)

final_lr_workflow <- lr_workflow %>%
  finalize_workflow(lr_best)

last_lr_fit <- final_lr_workflow %>%
  last_fit(data_split)

lr_final_fit <- final_lr_workflow %>%
  fit(data = train_data)

# Final metrics


last_lr_fit %>%
  collect_metrics()


# Variable importance

last_lr_fit %>%
  extract_fit_parsnip() %>%
  vip(num_features = 22)

last_lr_fit %>%
  extract_fit_parsnip() %>% 
  tidy()


# Final ROC curve on test set


last_lr_fit %>%
  collect_predictions() %>%
  roc_curve(SURVIVAL, .pred_1) %>%
  autoplot()




last_lr_fit %>%
  collect_predictions() %>% 
  roc_auc(truth = SURVIVAL, .pred_0)







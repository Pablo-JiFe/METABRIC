# Fix the random numbers by setting the seed 

# This enables the analysis to be reproducible when random numbers are used 

set.seed(222)

# Put 3/4 of the data into the training set 

data_split <- initial_split(late_genes.patients, strata = SURVIVAL)

# Create data frames for the two sets:

late_genes.other <- training(data_split)
late_genes.test  <- testing(data_split)


set.seed(234)
val_set <- initial_validation_split(late_genes.other, 
                            strata = SURVIVAL)

rf_val_resample <- validation_set(val_set)


rf_mod <- 
  rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
  set_engine("ranger") %>% 
  set_mode("classification")


rf_recipe <- 
  recipe(SURVIVAL ~ ., data = late_genes.other)

rf_workflow <- 
  workflow() %>% 
  add_model(rf_mod) %>% 
  add_recipe(rf_recipe)


set.seed(345)
rf_res <- 
  rf_workflow %>% 
  tune_grid(rf_val_resample,
            grid = 25,
            control = control_grid(save_pred = TRUE),
            metrics = metric_set(roc_auc))

rf_res %>% 
  show_best(metric = "roc_auc")
autoplot(rf_res)

rf_best <- 
  rf_res %>% 
  select_best(metric = "roc_auc")
rf_best


rf_res %>% 
  collect_predictions()


rf_auc <- 
  rf_res %>% 
  collect_predictions(parameters = rf_best) %>% 
  roc_curve(SURVIVAL, .pred_1) %>% 
  mutate(model = "Random Forest")

bind_rows(rf_auc, lr_auc) %>% 
  ggplot(aes(x = 1 - specificity, y = sensitivity, col = model)) + 
  geom_path(lwd = 1.5, alpha = 0.8) +
  geom_abline(lty = 3) + 
  coord_equal() + 
  scale_color_viridis_d(option = "plasma", end = .6)

rf_res %>% 
  collect_predictions(parameters = rf_best) %>% 
  roc_auc(truth = SURVIVAL, .pred_0)

library(tidymodels)
library(vip)

late_genes.patients <- late_genes.patients %>% 
  mutate(SURVIVAL = as.factor(SURVIVAL))

# Fix the random numbers by setting the seed 

# This enables the analysis to be reproducible when random numbers are used 

set.seed(222)

# Put 3/4 of the data into the training set 

data_split <- initial_split(late_genes.patients, strata = SURVIVAL)

# Create data frames for the two sets:

train_data <- training(data_split)
test_data  <- testing(data_split)

# Create recipe that defines the model formula
late_genes.rec <- 
  recipe(SURVIVAL ~ ., data = train_data) # Survival is the outcome and all variables are predictors

# Set model (glm for logistic regression)
late_genes.mod <- 
  logistic_reg() %>% 
  set_engine("glm")

# Workflow with the recipe and model
late_genes.wf <- 
  workflow() %>% 
  add_model(late_genes.mod) %>% 
  add_recipe(late_genes.rec)

# Fit logistic regression model with training data
late_genes.fit <- 
  late_genes.wf %>% 
  fit(data = train_data)

# Extract fitted model
late_genes.fit %>% 
  extract_fit_parsnip() %>% 
  tidy()

# Predictions on the test set
late_genes.aug <- 
  augment(late_genes.fit, test_data)

# ROC 
late_genes.aug %>% 
  roc_curve(truth = SURVIVAL, .pred_0) %>% 
  autoplot()

lr_auc <- late_genes.aug %>%
  roc_curve(SURVIVAL, .pred_0) %>%   
  mutate(model = "Logistic Regression")

# Plot
autoplot(lr_auc)

# AUC
late_genes.aug %>% 
  roc_auc(truth = SURVIVAL, .pred_0)

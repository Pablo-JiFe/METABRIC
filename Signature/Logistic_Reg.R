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

#

late_genes.rec <- 
  recipe(SURVIVAL ~ ., data = train_data)

late_genes.mod <- 
  logistic_reg() %>% 
  set_engine("glm")


late_genes.wf <- 
  workflow() %>% 
  add_model(late_genes.mod) %>% 
  add_recipe(late_genes.rec) 

late_genes.fit <- 
  late_genes.wf %>% 
  fit(data = train_data)

late_genes.fit %>% 
  extract_fit_parsnip() %>% 
  tidy()


late_genes.aug <- 
  augment(late_genes.fit, test_data)

late_genes.aug %>% 
  roc_curve(truth = SURVIVAL, .pred_0) %>% 
  autoplot()


lr_auc <- late_genes.aug %>%
  roc_curve(SURVIVAL, .pred_0) %>%   
  mutate(model = "Logistic Regression")

autoplot(lr_auc)


late_genes.aug %>% 
  roc_auc(truth = SURVIVAL, .pred_0)

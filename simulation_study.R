print("Working Directory:")
print(getwd())
print("Sys.Date():")
print(Sys.Date())

# Set the library path using the full file path
.libPaths("/users/4/neher015/R/x86_64-pc-linux-gnu-library/4.3")

# Ensure library path has been set
print(.libPaths())

# Load libraries
if (!("pacman" %in% installed.packages()[,"Package"])) {
  install.packages("pacman", lib = "/users/4/neher015/R/x86_64-pc-linux-gnu-library/4.3")
}

pacman::p_load(tidyverse, Rcpp, RcppArmadillo,
               parallel, gridExtra, rstan,
               scales, coda, MASS, lme4)

# Source scripts
source("src/utility_functions.R")
source("src/generate_data_v2.R")

# First read in the arguments listed at the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 1) {
  slurm_array_task_id <- args[1]
  print("class(slurm_array_task_id:")
  class(slurm_array_task_id)
  slurm_array_task_id <- as.numeric(args[1])
  } else {
    print("commandArgs have not been supplied.")
    # Supply default values
    slurm_array_task_id <- 16
  }

output_dir_name <- file.path("simulation_study_results",
                             paste0(Sys.Date(), "_simulation_study"))

print("output_dir_name:")
print(output_dir_name)

# Check if the directory exists, and if not, create it
if (!dir.exists(output_dir_name)) {
  dir.create(output_dir_name)
}

print("Slurm Array Task ID:")
print(slurm_array_task_id)

print("Main Output Directory:")
print(output_dir_name)

# Save the R environment details to a file in the output directory
session_info_file <- file.path(output_dir_name, "session_info.txt")
sessionInfo() %>%
  capture.output() %>%
  writeLines(session_info_file)

# Define fixed parameters
S <- 20 # Number of datasets to simulate/ scenario
n_iter <- 10000
n_burnin <- 5000
n_sample <- n_iter - n_burnin

# Define simulation scenarios
scenarios <- expand.grid(
  sigma2_ksi_true = c(0, .5, 1),
  sigma2_theta_true = c(0, .5, 1),
  covars = c(F), 
  n_views = 4,
  r = 4, # 6,
  features_per_view = 500, # 150, # c(150, 1000),
  method = c("BIP", "BIPmixed", "RandMVLearn", "PCA2step", "Integrative2step", "CooperativeLearning")
)

# Filter to retain only the scenarios relevant to the manuscript
scenarios <- scenarios %>%
  filter(
    (sigma2_theta_true == 0 & sigma2_ksi_true == 0) |       # Scenario 1
      (sigma2_theta_true == 1 & sigma2_ksi_true == .5) |  # Scenario 2
      (sigma2_theta_true == .5 & sigma2_ksi_true == 1)    # Scenario 3
  )

n_scenarios <- scenarios %>%
  group_by(method) %>%
  summarise(n = n()) %>%
  pull(n) %>%
  unique()

if (length(n_scenarios) !=1 ) {
  stop("n_scenarios can only be length 1.")
}

n_methods <- scenarios %>%
  pull(method) %>%
  unique() %>%
  length()

scenarios <- scenarios %>%
  mutate(scenario_id = rep(1:n_scenarios, times = n_methods))

# Repeat the scenarios S times and add the train_seed column
expanded_scenarios <- do.call(rbind, replicate(S, scenarios, simplify = FALSE))
expanded_scenarios$train_seed <- rep(1:S, each = nrow(scenarios))

# Create subdirectory for the slurm_array_task_id
subdir_names <- file.path(output_dir_name, paste0(Sys.Date(), "_task_", 1:nrow(expanded_scenarios)))
expanded_scenarios$subdir_name <- subdir_names

print("Number of tasks for this simulation study:")
print(nrow(expanded_scenarios))

# Save the scenarios
scenarios_file <- file.path(output_dir_name, "scenarios.csv")
expanded_scenarios %>% write.csv(scenarios_file)

# Extract simulation scenario for current job to global environment
subdir_name <- subdir_names[slurm_array_task_id]

print("subdir_name:")
print(subdir_name)

dir.create(subdir_name, showWarnings = FALSE, recursive = TRUE)
simulation_scenario <- expanded_scenarios[slurm_array_task_id, ]
train_seed <- simulation_scenario$train_seed
test_seed <- train_seed + S
covars_flag <- simulation_scenario$covars
sigma2_ksi_true <- simulation_scenario$sigma2_ksi_true
sigma2_theta_true <- simulation_scenario$sigma2_theta_true # Site-level variance fixed across sites
r <- simulation_scenario$r
n_views <- simulation_scenario$n_views
features_per_view <- simulation_scenario$features_per_view
method <- simulation_scenario$method

# Generate data
train_test_sets <- generate_train_test_sets(
  scenario = 1, setting = 1, 
  train_seed = train_seed,
  test_seed = test_seed,
  sigma2_ksi_true = sigma2_ksi_true, 
  sigma2_theta_true = sigma2_theta_true,
  development = T
)

train_set <- train_test_sets$train_set
test_set <- train_test_sets$test_set

if (method == "BIP" | method == "BIPmixed") {
  
  source("src/BIP.R") # Includes BIPmixed implementation
  source("src/BIPpredict.R")
  
  # # Reshape train_set for BIP/ BIPmixed
  # reshaped_data <- reshape_for_BIP(train_set, n_views, covars_flag)
  # trainList <- reshaped_data$trainList
  # IndicVar <- reshaped_data$IndicVar
  trainList <- list(train_set$X1, train_set$X2, train_set$X3, train_set$X4, train_set$Y)
  IndicVar <- c(rep(0, n_views),1)
  
  # # Reshape test_set for BIP_predict
  # test_data_reshaped <- reshape_for_BIPpredict(test_set, n_views, covars_flag)
  # testList <- test_data_reshaped$testList
  testList <- list(test_set$X1, test_set$X2, test_set$X3, test_set$X4)
  
  if (method == "BIP") {
    
    BIP_result <- BIP(dataList = trainList, IndicVar = IndicVar, Method = "BIP",
                      nbrcomp = r, sample = n_sample, burnin = n_burnin)
    
    y_preds <- BIPpredict(dataListNew = testList, Result = BIP_result, meth = "BMA")$ypredict

  } else if (method == "BIPmixed") {
    
    BIP_result <- BIP(dataList = trainList, IndicVar = IndicVar, Method = "BIPmixed",
                           nbrcomp = r, sample = n_sample, burnin = n_burnin,
                           Z_family = train_set$Z_family, Z_site = train_set$Z_site)
    
    y_preds <- BIPpredict(dataListNew = testList, Result=BIP_result, meth="BMA", 
                                   Z_site = test_set$Z_site)$ypredict
    
    # # Check convergence
    # convergence_check_result <- check_BIPmixed_coverage(BIP_result, train_set, n_iter, n_burnin, n_views)
    # write.csv(convergence_check_result$variance_param_coverage, file = file.path(subdir_name, paste0("variance_param_coverage_", slurm_array_task_id, ".csv")), row.names = FALSE)
    
  }
  
  # Get variable selection performance
  variable_selection_performance <- check_global_variable_selection(BIP_result, IndicVar, train_set, method_name = method)
  write.csv(variable_selection_performance, file = file.path(subdir_name, paste0("variable_selection_performance_", slurm_array_task_id, ".csv")), row.names = FALSE)
  
  if (train_seed == 1) {
    # Save the Model Fit
    saveRDS(BIP_result, file = file.path(subdir_name, "BIP_result.rds"))
  }
  
} 

if (method == "PCA2step") {
  
  # Early Fusion
  library(stats)
  library(lme4)
  
  # Concatenate X matrices (assuming they are the multi-view data)
  train_df <- cbind(train_set$X1, train_set$X2, train_set$X3, train_set$X4) %>%
    as.data.frame()
  
  # Perform PCA on the concatenated matrices
  pca_train <- prcomp(train_df, scale. = TRUE)
  
  W <- pca_train$rotation[, 1:r]
  
  # Use the first r principal components
  train_pca <- as.data.frame(pca_train$x[, 1:r])
  
  # Add random effects to the PCA-transformed data
  train_pca$Y <- train_set$Y
  train_pca$Z_site <- apply(train_set$Z_site, 1, function(r) which(r == 1)) %>%
    as.factor()
  train_pca$Z_family <- apply(train_set$Z_family, 1, function(r) which(r == 1)) %>%
    as.factor()
  
  # Fit the linear mixed model with family effects nested within site effects
  # Create a formula string by pasting together the PC terms
  pc_terms <- paste("PC", 1:r, sep = "", collapse = " + ")
  formula_string <- paste("Y ~", pc_terms, "+ (1 | Z_site) + (1 | Z_site:Z_family)")
  # Convert the string to a formula
  lmm_formula <- as.formula(formula_string)
  # Fit the model
  lmm_model <- lmer(lmm_formula, data = train_pca)
  
  # Now for the test set
  test_df <- cbind(test_set$X1, test_set$X2, test_set$X3, test_set$X4) %>% 
    as.data.frame()
  
  # Calculate test PCs
  test_pca <- as.matrix( test_df ) %*% W %>%
    as.data.frame()
  
  # Add random effects to the PCA-transformed test data
  test_pca$Z_site <- apply(test_set$Z_site, 1, function(r) which(r == 1)) %>%
    as.factor()
  test_pca$Z_family <- apply(test_set$Z_family, 1, function(r) which(r == 1)) %>%
    as.factor()
  
  # Predict using the fitted mixed model
  y_preds <- predict(lmm_model, newdata = test_pca, allow.new.levels = TRUE)
  
}

# if (method == "MOFA") {
#   
#   library(data.table)
#   reticulate::use_python("~/miniconda3/bin/python", required = TRUE)
#   library(MOFA2)
#   x_list <- train_set[c("X1", "X2", "X3", "X4")]  # Replace with your view names
#   xt_list <- lapply(x_list, function(x) t(x)) # MOFA requires transposed inputs
#   MOFAobject <- create_mofa(xt_list)
#   # plot_data_overview(MOFAobject)
#   data_opts <- get_default_data_options(MOFAobject)
#   model_opts <- get_default_model_options(MOFAobject)
#   train_opts <- get_default_training_options(MOFAobject)
#   MOFAobject <- prepare_mofa(
#     object = MOFAobject,
#     data_options = data_opts,
#     model_options = model_opts,
#     training_options = train_opts
#   )
#   output_dir_name = file.path(subdir_name,"model.hdf5")
#   MOFAobject.trained <- run_mofa(MOFAobject, output_dir_name, use_basilisk = T)
#   
# }

if (method == "Integrative2step") {
  
  library(r.jive)
  library(glmnet) # Load glmnet for LASSO
  library(lme4)
  views <- train_set[c("X1", "X2", "X3", "X4")] %>%
    lapply(function(x) t(x))
  jive_fit <- jive(views, scale = TRUE)
  joint_structure <- jive_fit$joint
  # Prepare the matrix of predictors and response
  X <- cbind(joint_structure[[1]], joint_structure[[2]], 
                    joint_structure[[3]], joint_structure[[4]]) %>%
    apply(2, scale) %>% t()
  Y <- matrix(train_set$Y, ncol = 1)
  # Apply LASSO with 10-fold cross-validation
  lasso_cv <- cv.glmnet(X, Y, alpha = 1, family = "gaussian", nfolds = 10)
  # Extract optimal lambda
  optimal_lambda <- lasso_cv$lambda.min
  # Fit the LASSO model with optimal lambda
  lasso_fit <- glmnet(X, Y, alpha = 1, lambda = optimal_lambda)
  # Identify selected features (excluding intercept)
  selected_features <- which(coef(lasso_fit)[-1] != 0)
  # Get names of selected features (excluding intercept)
  selected_feature_names <- paste0("V", selected_features)
  # Create a formula using selected features
  selected_terms <- paste(selected_feature_names, collapse = " + ")
  formula_string <- paste("Y ~", selected_terms, "+ (1 | Z_site) + (1 | Z_site:Z_family)")
  # Convert to formula and fit the model
  lmm_formula <- as.formula(formula_string)
  # Add random effects to the joint structure data
  train_df <- as.data.frame(X)[, selected_features]
  train_df$Y <- train_set$Y
  train_df$Z_site <- apply(train_set$Z_site, 1, function(r) which(r == 1)) %>%
    as.factor()
  train_df$Z_family <- apply(train_set$Z_family, 1, function(r) which(r == 1)) %>%
    as.factor()
  
  # Fit the linear mixed model
  lmm_model <- lmer(lmm_formula, data = train_df)
  
  # Now for the test set
  views <- test_set[c("X1", "X2", "X3", "X4")] %>%
    lapply(function(x) t(x))
  jive_fit <- jive(views, scale = TRUE)
  joint_structure <- jive_fit$joint
  # Prepare the matrix of predictors and response
  X <- cbind(joint_structure[[1]], joint_structure[[2]], 
             joint_structure[[3]], joint_structure[[4]]) %>%
    apply(2, scale) %>% t()
  
  # Add random effects to the joint structure data
  test_df <- as.data.frame(X)[, selected_features]
  test_df$Z_site <- apply(test_set$Z_site, 1, function(r) which(r == 1)) %>%
    as.factor()
  test_df$Z_family <- apply(test_set$Z_family, 1, function(r) which(r == 1)) %>%
    as.factor()
  
  # Predict using the fitted mixed model
  y_preds <- predict(lmm_model, newdata = test_df, allow.new.levels = TRUE)
  
  # Evaluate feature selection
  coefficients <- coef(lasso_fit)[-1] # Exclude the intercept
  p_m <- views %>% sapply(ncol)
  view_index <- rep(1:length(views), p_m)
  eta_hat <- rep(0, sum(p_m))
  eta_hat[selected_features] <- 1
  feature_selection_results <- data.frame(view = view_index,
                                          eta_hat = eta_hat)
  
  # Extract the truly important features
  eta_true_global <- train_set$TrueVar1 
  
  false_pos_rate <- c()
  false_neg_rate <- c()
  f1_measure <- c()
  view_label <- c()
  for (view in 1:n_views) {
    # True positives, false positives, true negatives, false negatives
    eta_hat_view <- feature_selection_results %>% filter(view == view) %>% pull(eta_hat)
    true_positive <- sum(eta_hat_view == 1 & eta_true_global == 1)
    false_positive <- sum(eta_hat_view == 1 & eta_true_global == 0)
    true_negative <- sum(eta_hat_view == 0 & eta_true_global == 0)
    false_negative <- sum(eta_hat_view == 0 & eta_true_global == 1)
    # Calculate False Positive Rate (FPR), False Negative Rate (FNR), and F1 measure
    fpr <- false_positive / (false_positive + true_negative)
    fnr <- false_negative / (false_negative + true_positive)
    precision <- true_positive / (true_positive + false_positive)
    recall <- true_positive / (true_positive + false_negative)
    f1 <- ( 2 * (precision * recall) / (precision + recall) ) * 100 # Report as a percentage
    # Append the calculated values to the respective lists
    false_pos_rate <- c(false_pos_rate, fpr)
    false_neg_rate <- c(false_neg_rate, fnr)
    f1_measure <- c(f1_measure, f1)
    view_label <- c(view_label, view)
  }
  
  variable_selection_performance <- data.frame(
    FalsePosRate = false_pos_rate,
    FalseNegRate = false_neg_rate,
    F1measure = f1_measure,
    AUC = NA,
    Method = factor(method, levels = levels(scenarios$method)),
    View = as.integer(view_label)
  )
  
  write.csv(variable_selection_performance, file = file.path(subdir_name, paste0("variable_selection_performance_", slurm_array_task_id, ".csv")), row.names = FALSE)
  
  
}

if (method == "CooperativeLearning") {
  
  library(multiview)
  x_list <- list(X1 = train_set$X1, X2 = train_set$X2, X3 = train_set$X3, X4 = train_set$X4)
  y <- train_set$Y
  cv_fit <- cv.multiview(x_list = x_list, y = y, rho = .5, family = gaussian())
  # plot(cv_fit)
  
  test_x_list <- list(X1 = test_set$X1, X2 = test_set$X2, X3 = test_set$X3, X4 = test_set$X4)
  y_preds <- predict(cv_fit, newx = test_x_list, s = "lambda.1se") %>%
    as.vector()
  print("y_preds:")
  print(y_preds)
  
  coefficients <- coef(cv_fit, s = "lambda.1se")[-1] # Exclude the intercept
  p_m <- x_list %>% sapply(ncol)
  view_index <- rep(1:length(x_list), p_m)
  feature_selection_results <- data.frame(view = view_index,
                                           coef = coefficients) %>%
    mutate(eta_hat = ifelse(coef==0, 0, 1))
  
  # Extract the truly important features
  eta_true_global <- train_set$TrueVar1 
  
  false_pos_rate <- c()
  false_neg_rate <- c()
  f1_measure <- c()
  view_label <- c()
  for (view in 1:n_views) {
    # True positives, false positives, true negatives, false negatives
    eta_hat_view <- feature_selection_results %>% filter(view == view) %>% pull(eta_hat)
    true_positive <- sum(eta_hat_view == 1 & eta_true_global == 1)
    false_positive <- sum(eta_hat_view == 1 & eta_true_global == 0)
    true_negative <- sum(eta_hat_view == 0 & eta_true_global == 0)
    false_negative <- sum(eta_hat_view == 0 & eta_true_global == 1)
    # Calculate False Positive Rate (FPR), False Negative Rate (FNR), and F1 measure
    fpr <- false_positive / (false_positive + true_negative)
    fnr <- false_negative / (false_negative + true_positive)
    precision <- true_positive / (true_positive + false_positive)
    recall <- true_positive / (true_positive + false_negative)
    f1 <- ( 2 * (precision * recall) / (precision + recall) ) * 100 # Report as a percentage
    # Append the calculated values to the respective lists
    false_pos_rate <- c(false_pos_rate, fpr)
    false_neg_rate <- c(false_neg_rate, fnr)
    f1_measure <- c(f1_measure, f1)
    view_label <- c(view_label, view)
  }
  
  variable_selection_performance <- data.frame(
    FalsePosRate = false_pos_rate,
    FalseNegRate = false_neg_rate,
    F1measure = f1_measure,
    AUC = NA,
    Method = factor(method, levels = levels(scenarios$method)),
    View = as.integer(view_label)
  )
  
  write.csv(variable_selection_performance, file = file.path(subdir_name, paste0("variable_selection_performance_", slurm_array_task_id, ".csv")), row.names = FALSE)
  
  }

if (method == "RandMVLearn") {
  
  # Center our outcomes to ensure comparability between RandMVLearn and BIP(mixed)
  Y_train_mean <- mean(train_set$Y)
  train_set$Y <- train_set$Y - Y_train_mean
  test_set$Y <- test_set$Y - Y_train_mean
  
  reticulate::use_python("~/miniconda3/bin/python", required = TRUE)
  library(RandMVLearn)
  torch <- reticulate::import("torch") # Import Python's torch library
  
  # Training data excluding Z_site and Z_family
  Xdata_train <- list(
    train_set$X1, 
    train_set$X2, 
    train_set$X3,
    train_set$X4
  )
  
  Y_train <- train_set$Y  # Outcome matrix
  
  # Testing data excluding Z_site and Z_family
  Xdata_test <- list(
    test_set$X1, 
    test_set$X2,
    test_set$X3,
    test_set$X4
  )
  
  Y_test <- test_set$Y  # Outcome matrix
  
  # Convert Xdata and Y to torch tensors for training
  Xdata_train_tensors <- lapply(Xdata_train, function(x) torch$from_numpy(x))
  Y_train_tensor <- torch$from_numpy(Y_train)
  
  # Convert Xdata and Y to torch tensors for testing
  Xdata_test_tensors <- lapply(Xdata_test, function(x) torch$from_numpy(x))
  Y_test_tensor <- torch$from_numpy(Y_test)
  
  # Training the model using RandMVLearnR
  RandMVLearn_train <- RandMVLearnR(myseed = as.integer(test_seed+1), # Ensures it's diff to any seeds used in data generation
                                    Xdata = Xdata_train_tensors,
                                    Y = Y_train_tensor,
                                    ncomponents = as.integer(r),
                                    outcometype = 'continuous',
                                    standardize_X = TRUE)
  
  # Predicting using RandMVPredict
  
  RandMVLearn_prediction_results <- RandMVPredict(Ytest = Y_test_tensor, Ytrain = Y_train_tensor, 
                                                  Xtest = Xdata_test_tensors, Xtrain = Xdata_train_tensors, 
                                                  myEstimates = RandMVLearn_train,
                                                  outcometype = 'continuous', standardize_X = TRUE)
  
  # Convert the tensor back to a NumPy array
  y_preds_np <- RandMVLearn_prediction_results$predictedEstimates$numpy()
  
  # Convert the NumPy array back to an R matrix
  y_preds <- as.matrix(y_preds_np)
  
  # We also extract the test MSE from
  RandMVLearn_MSE <- RandMVLearn_prediction_results$TestError$numpy() %>%
    as.numeric()
  
  # Get features selected by RandMVLearn
  # Note, the method returns a list of indicators: 
  # 1col matrices of length p_m where 1 implies the feature is selected
  features_selection_results <- RandMVLearn_train$Var_selection
  
  # Feature selection predictions
  features_selection_probs <- lapply(RandMVLearn_train$gamma, function(g) {
    g_np <- g$numpy()
    g <- as.vector(g_np)
    return(g)
  })
  
  # Extract the truly important features
  eta_true_global <- train_set$TrueVar1 
  
  # Assuming eta_true_global contains the ground truth (0/1) for all features
  # And features_selection_probs is a list with length equal to the number of views
  auc_results <- sapply(features_selection_probs, function(probabilities) {
    # Compute AUC for the current view
    roc_curve <- pROC::roc(eta_true_global, probabilities)
    pROC::auc(roc_curve)
  })
  
  # Get variable selection performance
  
  # Calculate the true selected features globally (eta_true_global)
  # Initialize empty lists to store performance metrics
  false_pos_rate <- c()
  false_neg_rate <- c()
  f1_measure <- c()
  view_label <- c()
  for (view in 1:n_views) {
    # True positives, false positives, true negatives, false negatives
    eta_hat_view <- as.vector(features_selection_results[[view]])
    true_positive <- sum(eta_hat_view == 1 & eta_true_global == 1)
    false_positive <- sum(eta_hat_view == 1 & eta_true_global == 0)
    true_negative <- sum(eta_hat_view == 0 & eta_true_global == 0)
    false_negative <- sum(eta_hat_view == 0 & eta_true_global == 1)
    # Calculate False Positive Rate (FPR), False Negative Rate (FNR), and F1 measure
    fpr <- false_positive / (false_positive + true_negative)
    fnr <- false_negative / (false_negative + true_positive)
    precision <- true_positive / (true_positive + false_positive)
    recall <- true_positive / (true_positive + false_negative)
    f1 <- ( 2 * (precision * recall) / (precision + recall) ) * 100 # Report as a percentage
    # Append the calculated values to the respective lists
    false_pos_rate <- c(false_pos_rate, fpr)
    false_neg_rate <- c(false_neg_rate, fnr)
    f1_measure <- c(f1_measure, f1)
    view_label <- c(view_label, view)
  }
  
  variable_selection_performance <- data.frame(
    FalsePosRate = false_pos_rate,
    FalseNegRate = false_neg_rate,
    F1measure = f1_measure,
    AUC = auc_results,
    Method = factor(method, levels = levels(scenarios$method)),
    View = as.integer(view_label)
  )
  
  write.csv(variable_selection_performance, file = file.path(subdir_name, paste0("variable_selection_performance_", slurm_array_task_id, ".csv")), row.names = FALSE)
  
} 

# Estimate prediction performance metrics
family_ids <- apply(test_set$Z_family, 1, function(row) which(row==1))

prediction_data <- data.frame(
  True_Y = test_set$Y,
  Predicted_Y = y_preds,
  Method = method, length(test_set$Y),
  Family = c(family_ids)
)

print("head(prediction_data):")
print(head(prediction_data))

# Prediction performance overall
prediction_performance_by_method <- prediction_data %>% 
  group_by(Method) %>% 
  summarize(calculate_prediction_metrics(True_Y, Predicted_Y))

if (method=="RandMVLearn") {
  prediction_performance_by_method$MSE <- RandMVLearn_MSE
}

# Prediction performance stratified by cluster
prediction_performance_by_family <- prediction_data %>% 
  group_by(Method, Family) %>% 
  summarize(calculate_prediction_metrics(True_Y, Predicted_Y))

if (method=="RandMVLearn") {
  prediction_performance_by_family$MSE <- NA
}

# Write-out prediction results
write.csv(prediction_data, file = file.path(subdir_name, paste0("prediction_data_", slurm_array_task_id, ".csv")), row.names = FALSE)
write.csv(prediction_performance_by_method, file = file.path(subdir_name, paste0("prediction_performance_by_method_", slurm_array_task_id, ".csv")), row.names = FALSE)
write.csv(prediction_performance_by_family, file = file.path(subdir_name, paste0("prediction_performance_by_family_", slurm_array_task_id, ".csv")), row.names = FALSE)

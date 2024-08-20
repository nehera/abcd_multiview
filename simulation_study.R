# First read in the arguments listed at the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  print("Not enough arguments supplied.")
  # Supply default values
  slurm_array_task_id <- 1
  slurm_job_id <- 999999
} else {
  slurm_array_task_id <- args[1]
  slurm_job_id <- args[2]
}

output_dir_name <- file.path("simulation_study_results",
                             paste0(Sys.Date(), "_simulation_study"))

# Check if the directory exists, and if not, create it
if (!dir.exists(output_dir_name)) {
  dir.create(output_dir_name)
}

print("Slurm Array Task ID:")
print(slurm_array_task_id)
print("Slurm Job ID:")
print(slurm_job_id)
print("Main Output Directory:")
print(output_dir_name)

# Create subdirectory for the slurm_array_task_id and slurm_job_id
subdir_name <- file.path(output_dir_name, paste0(Sys.Date(), "_task_", slurm_array_task_id, "_job_", slurm_job_id))
dir.create(subdir_name, showWarnings = FALSE, recursive = TRUE)

# Load libraries
if (("pacman" %in% installed.packages()[,"Package"]) == FALSE) { install.packages("pacman") }
pacman::p_load(tidyverse, Rcpp, RcppArmadillo,
               parallel, gridExtra, rstan,
               scales, coda, MASS, lme4)

if (("BIPnet" %in% installed.packages()[,"Package"]) == FALSE) {
  pacman::p_load(devtools)
  devtools::install_github('chekouo/BIPnet')
}

# Source scripts
source("src/simulation_study_data_generator.R")
source("src/BIP.R") # Includes BIPmixed implementation
source("src/BIPpredict.R")

# Save the R environment details to a file in the output directory
session_info_file <- file.path(output_dir_name, "session_info.txt")
sessionInfo() %>%
  capture.output() %>%
  writeLines(session_info_file)

# Define fixed parameters
S <- 100 # Number of datasets to simulate/ scenario
n_iter <- 10000
n_burnin <- floor(n_iter*0.5)
n_sample <- n_iter - n_burnin

# Define simulation scenarios
scenarios <- expand.grid(
  sigma2_ksi_true = c(0, 0.75, 1.5),
  sigma2_theta_true = c(0, 0.75, 1.5),
  covars = c(F), 
  n_views = 4, 
  r = 6,
  features_per_view = 150
) 

scenarios <- scenarios %>%
  mutate(scenario_id = 1:nrow(scenarios))

# Repeat the scenarios S times and add the train_seed column
expanded_scenarios <- do.call(rbind, replicate(S, scenarios, simplify = FALSE))
expanded_scenarios$train_seed <- rep(1:S, each = nrow(scenarios))

print("Number of tasks for this simulation study:")
print(nrow(expanded_scenarios))

# Save the scenarios
scenarios_file <- file.path(output_dir_name, "scenarios.csv")
expanded_scenarios %>% write.csv(scenarios_file)

# Extract simulation scenario for current job to global environment 
simulation_scenario <- expanded_scenarios[slurm_array_task_id, ]
train_seed <- simulation_scenario$train_seed
test_seed <- train_seed + S
covars_flag <- simulation_scenario$covars
sigma2_ksi_true <- simulation_scenario$sigma2_ksi_true
N_sites <- 30 # Assume n sites fixed across train and test sets
sigma2_theta_true <- rep(simulation_scenario$sigma2_theta_true, N_sites)
r <- simulation_scenario$r
n_views <- simulation_scenario$n_views
features_per_view <- simulation_scenario$features_per_view

# Generate train data
train_set <- simulation_study_data_generator(seed = train_seed,
                                             sigma2_ksi_true,
                                             sigma2_theta_true,
                                             covars = covars_flag,
                                             n_views,
                                             features_per_view,
                                             r,
                                             train_set = T)

# Reshape train data
trainList <- list(train_set$Y)
for (m in 1:n_views) {
  trainList[[m+1]] <- train_set$X[[m]]
}
IndicVar <- c(1, rep(0, n_views))
if (covars_flag == T) {
  covar_index <- length(trainList) + 1
  trainList[[covar_index]] <- train_set$W
  IndicVar[covar_index] <- 2
}

# Fit models
BIP_start_time <- Sys.time()
BIP_result <- BIP(dataList = trainList, IndicVar = IndicVar, Method = "BIP",
                  nbrcomp = r, sample = n_sample, burnin = n_burnin)
BIP_end_time <- Sys.time()
print("BIP required:")
print(BIP_end_time - BIP_start_time)

BIPmixed_start_time <- Sys.time()
BIPmixed_result <- BIP(dataList = trainList, IndicVar = IndicVar, Method = "BIPmixed",
                    nbrcomp = r, sample = n_sample, burnin = n_burnin,
                    Z_family = train_set$Z_family, Z_site = train_set$Z_site)
BIPmixed_end_time <- Sys.time()
print("BIPmixed required")
print(BIPmixed_end_time - BIPmixed_start_time)

# Check BIPmixed convergence

# Combine samples into a 3D array for rstan::monitor
combine_samples_nested <- function(samples_list, n_iter, n_chains) {
  n_params <- ncol(samples_list[[1]]$mu_samples) +
    ncol(samples_list[[1]]$beta_samples) +  # FIGURE OUT WHAT TO DO ABOUT THIS IF THERE'S NO FIXED EFFECTS
    ncol(samples_list[[1]]$ksi_samples) +
    ncol(samples_list[[1]]$theta_samples) +
    ncol(samples_list[[1]]$sigma2_ksi_samples) +
    ncol(samples_list[[1]]$sigma2_samples) +
    ncol(samples_list[[1]]$sigma2_theta_samples)  # beta, ksi, theta, sigma2_ksi, sigma2, sigma2_theta
  
  # Create an empty array
  combined_array <- array(NA, dim = c(n_iter, n_chains, n_params))
  
  # Fill in the array
  for (chain in 1:n_chains) {
    samples <- samples_list[[chain]]
    combined_array[, chain, 1] <- samples$mu
    combined_array[, chain, 2:(1 + ncol(samples$beta_samples))] <- samples$beta_samples
    combined_array[, chain, (1 + ncol(samples$beta_samples) + 1):(1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples))] <- samples$ksi_samples
    combined_array[, chain, (1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + 1):(1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples))] <- samples$theta_samples
    combined_array[, chain, 1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + 1] <- samples$sigma2_ksi_samples
    combined_array[, chain, (1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + 2):(1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + 1 + ncol(samples$sigma2_theta_samples))] <- samples$sigma2_theta_samples
    combined_array[, chain, (2 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + ncol(samples$sigma2_theta_samples) + 1)] <- samples$sigma2_samples
  }
  
  return(combined_array)
}
samples_list <- list(BIPmixed_result)
n_chains <- length(samples_list)
combined_samples <- combine_samples_nested(samples_list, n_iter, n_chains)

# Parameter names
param_names <- c("mu",
                 paste0("beta_", 1:ncol(samples_list[[1]]$beta_samples)),
                 paste0("ksi_", 1:ncol(samples_list[[1]]$ksi_samples)),
                 paste0("theta_", 1:ncol(samples_list[[1]]$theta_samples)),
                 "sigma2_ksi",
                 paste0("sigma2_theta_", 1:ncol(samples_list[[1]]$sigma2_theta_samples)),
                 "sigma2")

# Assign parameter names
dimnames(combined_samples) <- list(
  iterations = NULL,
  chains = NULL,
  parameters = param_names
)

# Use rstan::monitor
mcmc_summary <- monitor(combined_samples)

# Extract summary statistics
mean_values <- round(mcmc_summary$mean, 4)
median_values <- round(mcmc_summary$`50%`, 4)
lower_bounds <- round(mcmc_summary$`2.5%`, 4)
upper_bounds <- round(mcmc_summary$`97.5%`, 4)

# Combine the true values into a single vector, ordered according to the MCMC parameters
mu_true <- train_set$mu
ksi_true <- train_set$ksi
theta_true <- train_set$theta
beta_true <- train_set$beta
U_true <- train_set$U
alpha_true <- train_set$alpha
sigma2_ksi_true <- train_set$nu2$sigma2_ksi
sigma2_theta_true <- train_set$nu2$sigma2_theta
sigma2_true <- train_set$sigma2
true_values <- c(mu_true, beta_true, ksi_true, theta_true, sigma2_ksi_true, sigma2_theta_true, sigma2_true)

# Initial values
# Extract initial values
init_values <- samples_list[[1]]$initial_values
mu_init <- init_values$mu_init
beta_init <- init_values$beta_init
ksi_init <- init_values$ksi_init
theta_init <- init_values$theta_init
sigma2_ksi_init <- init_values$sigma2_ksi_init
sigma2_theta_init <- init_values$sigma2_theta_init
sigma2_init <- init_values$sigma2_init
initial_values <- c(mu_init, beta_init, ksi_init, theta_init,
                    sigma2_ksi_init, sigma2_theta_init, sigma2_init)

# Initialize a dataframe to hold the comparisons
comparison <- data.frame(
  param_name = param_names,
  lower = lower_bounds,
  mean = mean_values,
  median = median_values,
  upper = upper_bounds,
  true_value = round(true_values, 4),
  initial = round(initial_values, 4)
)

# Add a logical vector to see if true values are within the credible intervals
comparison$within_credible_interval <- with(comparison, true_value >= lower & true_value <= upper)

# Filter for fixed effect parameters
fixed_effect_comparison <- comparison %>% filter(grepl("^mu|^beta_", param_name))

# Filter for random intercept parameters (those starting with "ksi_")
random_intercept_comparison <- comparison %>% filter(grepl("^ksi_|^theta_", param_name))

# % of random intercept parameters with correct sign
random_intercept_comparison$correct_sign <- sign(random_intercept_comparison$mean) == sign(random_intercept_comparison$true_value)

# Filter the comparison dataframe for variance parameters
variance_comparison <- comparison %>% filter(grepl("^sigma2$|^sigma2_ksi$|^sigma2_theta_", param_name))

cat("% of all parameters within credible interval: ", mean(comparison$within_credible_interval), "\n")
cat("% of random intercept parameters within credible interval: ", mean(random_intercept_comparison$within_credible_interval), "\n")
cat("% of random intercept parameters with correct sign: ", mean(random_intercept_comparison$correct_sign), "\n")
cat("% of variance parameters within credible interval: ", mean(variance_comparison$within_credible_interval), "\n")

# Store variance parameter coverage, calculating the proportion of site variances where the truth is covered in CI
variance_param_coverage <- data.frame(
  param_name = variance_comparison$param_name,
  within_credible_interval = variance_comparison$within_credible_interval
) %>%
  mutate(param_type = case_when(
    param_name == "sigma2_ksi" ~ "sigma2_ksi",
    grepl("^sigma2_theta_", param_name) ~ "sigma2_theta",
    param_name == "sigma2" ~ "sigma2"
  )) %>%
  group_by(param_type) %>%
  summarise(
    total = n(),
    within_ci = sum(within_credible_interval),
    proportion_within_ci = within_ci / total
  )

print(variance_param_coverage)

# Create a combined data frame in one step
variance_df <- data.frame(
  sigma2 = as.vector(BIPmixed_result$sigma2_samples),
  sigma2_ksi = as.vector(BIPmixed_result$sigma2_ksi_samples)
) %>%
  mutate(iter = 1:n()) %>%
  gather(key = "parameter", value = "value", -iter)

# Plot the variance trace plot
variance_trace_plot <- ggplot(variance_df, aes(x = iter, y = value, color = parameter)) +
  geom_line(alpha = 0.7) +
  geom_vline(xintercept = n_burnin, linetype = "dashed", color = "black") + 
  labs(x = "Iteration", y = "Value", color = "Parameter") +
  theme_minimal()

# Display the plot
print(variance_trace_plot)

# Generate test data
test_set <- simulation_study_data_generator(seed = test_seed,
                                            sigma2_ksi_true,
                                            sigma2_theta_true,
                                            covars = covars_flag,
                                            n_views,
                                            features_per_view,
                                            r,
                                            train_set = F)

# Reshape test data
testList <- test_set$X
IndicVar <- rep(0, n_views)
if (covars_flag == T) {
  covar_index <- length(testList) + 1
  testList[[covar_index]] <- test_set$W
}

# Get predictions
y_preds_BIP <- BIPpredict(dataListNew = testList, Result=BIP_result, meth="BMA")$ypredict

y_preds_BIPmixed <- BIPpredict(dataListNew = testList, Result=BIPmixed_result, meth="BMA", 
                               Z_site = test_set$Z_site, Z_family = test_set$Z_family)$ypredict

## -- Estimate performance metrics

# Estimate variable selection performance

gamma_true <- test_set$gamma
outcome_index <- which(IndicVar == 1)
BIP_OutCompoSelPerformance <- BIPnet::ComputeVarCriteria(BIP_result$VarSelMean[[outcome_index]], gamma_true, thres = 0.5) %>%
  as.data.frame() %>%
  mutate(Method="BIP")
BIPmixed_OutCompoSelPerformance <- BIPnet::ComputeVarCriteria(BIPmixed_result$VarSelMean[[outcome_index]], gamma_true, thres = 0.5) %>%
  as.data.frame() %>%
  mutate(Method="BIPmixed")
component_selection_performance <- rbind(
  BIP_OutCompoSelPerformance,
  BIPmixed_OutCompoSelPerformance
)

# Let's assess variable selection globally
eta_true_global <- colSums(test_set$Eta)
omics_index <- which(IndicVar == 0)
BIP_VarSelGlobalPerformance_list <- lapply(BIP_result$VarSelMeanGlobal[omics_index], function(pred.prob) {
  BIPnet::ComputeVarCriteria(pred.prob, eta_true_global, thres = 0.5) %>%
    as.data.frame() %>%
    mutate(Method = "BIP")
})
BIPmixed_VarSelGlobalPerformance_list <- lapply(BIPmixed_result$VarSelMeanGlobal[omics_index], function(pred.prob) {
  BIPnet::ComputeVarCriteria(pred.prob, eta_true_global, thres = 0.5) %>%
    as.data.frame() %>%
    mutate(Method = "BIPmixed")
})

# Function to combine list of data frames and add index column
combine_dataframes_with_index <- function(df_list, index_col_name = "View") {
  # Add index to each data frame
  for (i in seq_along(df_list)) {
    df_list[[i]][[index_col_name]] <- i
  }
  # Combine all data frames into one
  combined_df <- do.call(rbind, df_list)
  return(combined_df)
}

variable_selection_performance <- rbind(
  combine_dataframes_with_index(BIP_VarSelGlobalPerformance_list),
  combine_dataframes_with_index(BIPmixed_VarSelGlobalPerformance_list)
)

# Estimate prediction performance metrics

# Define a function to calculate MSE, bias squared, prediction variance
calculate_prediction_metrics <- function(Y, y_preds) {
  mse <- mean((Y - y_preds)^2)
  bias2 <- (mean(y_preds) - mean(Y))^2
  var_pred <- var(y_preds)
  data.frame(MSE = mse, Bias2 = bias2, Variance = var_pred)
}

# We reshape data for prediction analysis
family_ids <- apply(test_set$Z_family, 1, function(row) which(row==1))
prediction_data <- data.frame(
  True_Y = c(test_set$Y, test_set$Y),
  Predicted_Y = c(y_preds_BIP, y_preds_BIPmixed),
  Method = rep(c("BIP", "BIPmixed"), each = length(test_set$Y)),
  Family = c(family_ids, family_ids)
)

# Prediction performance overall
prediction_performance_by_method <- prediction_data %>% 
  group_by(Method) %>% 
  summarize(calculate_prediction_metrics(True_Y, Predicted_Y))

# Prediction performance stratified by cluster
prediction_performance_by_family <- prediction_data %>% 
  group_by(Method, Family) %>% 
  summarize(calculate_prediction_metrics(True_Y, Predicted_Y))

# Create the scatterplot with the 1:1 reference line
prediction_scatter_plot <-ggplot(prediction_data, aes(x = True_Y, y = Predicted_Y, color = Method)) +
  geom_point(alpha = 0.8) +  # Scatter plot with points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 1:1 line
  labs(x = "True Y",
       y = "Predicted Y",
       color = "Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(prediction_scatter_plot)

# Prediction summary statistics overall
prediction_data %>% group_by(Method) %>% summarise(cor(True_Y, Predicted_Y))
prediction_data %>% group_by(Method) %>% summarise(mean(Predicted_Y))

# Prediction distribution relative to truth overall
pred_distribution_data <- data.frame(
  Y = c(test_set$Y, y_preds_BIP, y_preds_BIPmixed),
  Y_type = rep(c("True_Y", "BIP_Preds", "BIPmixed_Preds"), each = length(test_set$Y))
)

pred_distribution_plot <- ggplot(pred_distribution_data, aes(x = Y, color = Y_type)) + 
  geom_density() +
  geom_vline(xintercept = 1)

print(pred_distribution_plot)

# Create density plot for MSE stratified by family
mse_plot <- ggplot(prediction_performance_by_family, aes(x = MSE, fill = Method, color = Method)) +
  geom_density(alpha = 0.5) +
  labs(x = "Mean Squared Error (MSE)",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend from this plot

# Create density plot for Mean Bias stratified by family
bias_plot <- ggplot(prediction_performance_by_family, aes(x = Bias2, fill = Method, color = Method)) +
  geom_density(alpha = 0.5) +
  labs(x = "Bias2",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend from this plot

# Create density plot for Variance
variance_plot <- ggplot(prediction_performance_by_family, aes(x = Variance, fill = Method, color = Method)) +
  geom_density(alpha = 0.5) +
  labs(x = "Variance in Predictions",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")  # Include legend only in this plot

# Arrange the plots in a grid
prediction_family_stratified_plot <- grid.arrange(mse_plot, bias_plot, variance_plot, ncol = 1)

print(prediction_family_stratified_plot)

## -- Write-out results

# Save each data frame as a .csv file with the slurm_array_task_id in the file name
write.csv(variance_param_coverage, file = file.path(subdir_name, paste0("variance_param_coverage_", slurm_array_task_id, ".csv")), row.names = FALSE)
write.csv(component_selection_performance, file = file.path(subdir_name, paste0("component_selection_performance_", slurm_array_task_id, ".csv")), row.names = FALSE)
write.csv(variable_selection_performance, file = file.path(subdir_name, paste0("variable_selection_performance_", slurm_array_task_id, ".csv")), row.names = FALSE)
write.csv(prediction_performance_by_method, file = file.path(subdir_name, paste0("prediction_performance_by_method_", slurm_array_task_id, ".csv")), row.names = FALSE)
write.csv(prediction_performance_by_family, file = file.path(subdir_name, paste0("prediction_performance_by_family_", slurm_array_task_id, ".csv")), row.names = FALSE)

# Save the plots as .png files with the slurm_array_task_id in the file name
png(file = file.path(subdir_name, paste0("prediction_scatter_plot_", slurm_array_task_id, ".png")))
print(prediction_scatter_plot)
dev.off()

png(file = file.path(subdir_name, paste0("prediction_family_stratified_plot_", slurm_array_task_id, ".png")))
print(prediction_family_stratified_plot)
dev.off()

png(file = file.path(subdir_name, paste0("pred_distribution_plot_", slurm_array_task_id, ".png")))
print(pred_distribution_plot)
dev.off()

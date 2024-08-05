# First read in the arguments listed at the command line
args <- commandArgs(trailingOnly=TRUE) # args is now a list of character vectors 

if(length(args) == 0){ 
  print("No arguments supplied.") 
  # Supply default value
  slurm_id <- 1 
} else{ 
  slurm_id = args[1] 
} 
cat("Slurm ID:", slurm_id)

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
source("simulation_study_data_generator")
source("BIPpredict.R")

# Set-up output directory
current_date <- Sys.Date()
# Create a directory name with the current date
dir_name <- paste0("simulation_output_", current_date)
dir.create(dir_name)
# Save the R environment details to a file in the new directory
session_info_file <- file.path(dir_name, "session_info.txt")
sessionInfo() %>%
  capture.output() %>%
  writeLines(session_info_file)

# Define fixed parameters
S <- 100 # Number of datasets to simulate/ scenario
N_sites <- 30
n_families_per_site_train <- 100
n_families_per_site_test <- 25
n_individs_per_family <- 2
r <- 4 # Number of components in latent factor
n_chains <- 2
n_iter <- 10000
n_burnin <- floor(n_iter*0.5)
n_sample <- n_iter - n_burnin

# Define simulation scenarios
scenarios <- expand.grid(
  high_signal_to_noise = c(F,T),
  high_kurtosis = c(F,T),
  covars = c(F,T)
)

# Extract simulation scenario for current job to global environment 
simulation_scenario <- scenarios[slurm_id, ]
high_signal_to_noise_flag <- simulation_scenario$high_signal_to_noise
high_kurtosis_flag <- simulation_scenario$high_kurtosis
covars_flag <- simulation_scenario$covars

# Generate data
train_set <- simulation_study_data_generator(seed = slurm_id,
                                             high_signal_to_noise = high_signal_to_noise_flag,
                                             high_kurtosis = high_kurtosis_flag,
                                             covars = covars_flag,
                                             train_set = T)

test_set <- simulation_study_data_generator(seed = slurm_id,
                                            high_signal_to_noise = high_signal_to_noise_flag,
                                            high_kurtosis = high_kurtosis_flag,
                                            covars = covars_flag,
                                            train_set = F)

# Reshape data
trainList <- list(train_set$Y, train_set$X[[1]], train_set$X[[2]], train_set$X[[3]], train_set$X[[4]])
IndicVar <- c(0, 1, 1, 1, 1)
if (covars_flag == T) {
  trainList[[6]] <- train_set$W
  IndicVar[6] <- 2
}

# Fit models
BIP_result <- BIPnet::BIP(dataList = trainList, IndicVar = IndicVar, Method = "BIP",
                          # TODO remove n_sample workaround by using our version of BIP
                          nbrcomp = r, sample = n_sample+1, burnin = n_burnin)

BIPmixed_result <- BIP(dataList = dataList, IndicVar = IndicVar, Method = "BIPmixed",
                    nbrcomp = r, sample = n_sample, burnin = n_burnin,
                    Z_family = Z_family, Z_site = Z_site, Z_family_to_site = Z_family_to_site,
                    # TODO remove prior specification
                    mu_prior_var = priors$mu_prior_var, mu_beta = priors$mu_beta,
                    beta_prior_var = priors$beta_prior_var,
                    sigma_ksi_prior = c(priors$sigma_ksi_prior_a, priors$sigma_ksi_prior_b),
                    sigma_theta_prior = c(priors$sigma_theta_prior_a, priors$sigma_theta_prior_b),
                    sigma_prior = c(priors$sigma_prior_a, priors$sigma_prior_b))

# Get predictions
# TODO Add method to result upstream
BIP_result$Method <- 0 
y_preds_BIP <- BIPpredict(dataListNew=list(simulation_results_new$X[[1]], simulation_results_new$W), 
                          Result=BIP_result, meth="BMA")$ypredict

y_preds_BIPmixed <- BIPpredict(dataListNew=list(simulation_results_new$X[[1]]), 
                               Result=BIPmixed_result, meth="BMA", 
                               Wnew = simulation_results_new$W, 
                               Z_family_to_site = simulation_results_new$Z_family_to_site, 
                               Z_family = simulation_results_new$Z_family)$ypredict

# Estimate performance metrics
  # Estimate variable selection performance
  # Estimate prediction performance metrics

# Write-out results
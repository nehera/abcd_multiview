dev <- FALSE
if (dev == FALSE) {
  # Get CLI arguments
  args = commandArgs(trailingOnly=TRUE)
  # test if there is at least one argument: if not, return an error
  if (length(args)==0) {
    stop("At least one argument must be supplied in CLI", call.=FALSE)
  } else if (length(args)==1) {
    # default output file
    args[2] = "out.txt"
  }
  job <- as.numeric(args[1])
} else {
  job <- 1
}

# Simulate data
source("code/simulate_simple_data.R")
simulation_results <- simulate_iid_data(seed=job)
# Get BIP results
data_list <- list(simulation_results$X_list[[1]], 
                  simulation_results$X_list[[2]], 
                  simulation_results$Y)
# Scale data
data_list <- lapply(data_list, scale)
indic_var <- c(0, 0, 1) 
method <- "BIP" # method without grouping information to start
group_list <- NULL # default when grouping information not included
bip_0 <- BIPnet::BIP(dataList = data_list, IndicVar = indic_var, Method = method)
# Store simulated data and results
today_date <- Sys.Date()
fname_simulation_results <- paste0("data/", today_date, "_job_", job, "_simulation_results.rds")
saveRDS(simulation_results, fname_simulation_results)
fname_data_list <- paste0("data/", today_date, "_job_", job, "_simulation_data_list.rds")
saveRDS(data_list, fname_data_list)
fname_bip_0 <- paste0("data/", today_date, "_job_", job, "_simulation_BIP_results.rds")
saveRDS(bip_0, fname_bip_0)

# Load packages
library(tidyverse)
# Set parameters
r <- 4
n_obs <- 200
p_m <- 10
prob_component_selection <- 0.5
prob_feature_selection <- 0.5
n_sample <- 5000
n_burnin <- 1000
n_iterations <- n_sample + n_burnin
# Load data
m <- 1 # Start with the 1st view
# data_list <- readRDS("data/2023-09-08_simulation_data_list.rds")
x <- data_list[[m]]

# Set covariates
sigma2 <- rep(1, p_m)
tau2 <- matrix(1, nrow = r, ncol = p_m)
# simulation_results <- readRDS("data/2023-09-08_simulation_results.rds")
U <- simulation_results$U # U is a fixed set of covariates

# Define sub-routines
initialize_gamma <- function(r=4, prior_component_selection=0.5) {
  gamma <- rbinom(n = r, size = 1, prob = prior_component_selection)
  return(gamma)
}

initialize_eta <- function(gamma, r=4, p_m=10, prior_feature_selection=0.5) {
  eta <- matrix(nrow = r, ncol = p_m)
  for (l in 1:r) {
    if (gamma[l] == 1) {
      eta[l, ] <- rbinom(n = p_m, size = 1, prior_feature_selection)
    } else {
      eta[l, ] <- rep(0, p_m)
    }
  }
  return(eta)
}

## -- calculate P_lj
calculate_mvnorm_var_j <- function(eta_j, sigma2_j, tau2_j, U, n_obs) {
  n_components_active_in <- sum(eta_j)
  if (n_components_active_in > 0) {
    components_active_in <- which(eta_j==1)
    Sigma2_j <- U[, components_active_in, drop = FALSE] %*% 
      diag(n_components_active_in) %*% 
      t(U[, components_active_in, drop = FALSE]) + # TODO use Woodbury matrix formula for inversion
      diag(n_obs) 
  } else {
    Sigma2_j <- diag(n_obs)
  }
  mvnorm_var <- sigma2_j * Sigma2_j
  return(mvnorm_var)
}

calculate_log_dmvnorm_j <- function(eta_j, sigma2_j, tau2_j, U, n_obs, x_j) {
  mvnorm_var_j <- calculate_mvnorm_var_j(eta_j, sigma2_j, tau2_j, U, n_obs)
  log_dmvnorm_j <- mvtnorm::dmvnorm(x = x_j, mean = rep(0, n_obs), sigma = mvnorm_var_j, log = TRUE)
  return(log_dmvnorm_j)
}

calculate_log_G_j <- function(gamma, eta_j, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection) {
  log_dmvnorm_j <- calculate_log_dmvnorm_j(eta_j, sigma2_j, tau2_j, U, n_obs, x_j)
  active_gamma <- which(gamma == 1)
  n_active_gamma <- length(active_gamma) # TODO account for all inactive
  n_active_eta_given_gamma_1 <- eta_j[active_gamma] %>% sum()
  log_G_j <- log_dmvnorm_j + 
    n_active_eta_given_gamma_1 * log(prob_feature_selection) + 
    (n_active_gamma - n_active_eta_given_gamma_1) * log(1 - prob_feature_selection)
  return(log_G_j) 
}

calculate_log_PQ_lj <- function(l, gamma, eta_j, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection) {
  gamma_1 <- replace(gamma, l, 1)
  eta_1 <- replace(eta_j, l, 1)
  eta_0 <- replace(eta_j, l, 0)
  log_G_j_1 <- calculate_log_G_j(gamma_1, eta_1, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  log_G_j_0 <- calculate_log_G_j(gamma_1, eta_0, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  # Use logsumexp trick
  max_arg <- max(log_G_j_1, log_G_j_0) 
  a <- log_G_j_1 - max_arg
  b <- log_G_j_0 - max_arg
  # TODO understand if the log(exp(...)) is ok.
  log_P_lj <- a - log(exp(a)+exp(b)) 
  log_Q_lj <- b - log(exp(a)+exp(b))
  return(c(log_P_lj, log_Q_lj))
}

calculate_eta_lj_threshold <- function(l, gamma, eta_j, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection) {
  gamma_1 <- replace(gamma, l, 1)
  eta_1 <- replace(eta_j, l, 1)
  eta_0 <- replace(eta_j, l, 0)
  log_G_j_1 <- calculate_log_G_j(gamma_1, eta_1, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  log_G_j_0 <- calculate_log_G_j(gamma_1, eta_0, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  return(log_G_j_1 - log_G_j_0)
}

## -- Calculate log_target_density_l 
log_target_density_l <- function(l, gamma, eta, sigma2, tau2, U, n_obs, p_m, 
                                 x, prob_component_selection, prob_feature_selection) {
  sum_log_target_density_lj <- 0
  for (j in 1:p_m) {
    log_dmvnorm_j <- calculate_log_dmvnorm_j(eta[, j], sigma2[j], tau2[,j], U, n_obs, x[,j])
    sum_log_target_density_lj <- sum_log_target_density_lj + log_dmvnorm_j + 
      eta[l,j]*log(prob_feature_selection) + (1-eta[l,j])*log(1-prob_feature_selection)
  }
  log_target_density_l <- gamma[l]*log(prob_component_selection) + 
    (1-gamma[l])*log(1-prob_component_selection) + sum_log_target_density_lj
  return(log_target_density_l)
}

## -- Calculate log_proposal_l
log_proposal_l_density <- function(l, gamma_prime, eta_prime, gamma, eta,
                                   sigma2, tau2, U, n_obs, p_m, x, prob_feature_selection) {
  if (gamma[l]==1 & gamma_prime[l]==0 & sum(eta_prime[l,])==0) {
    return(0)
  } else if (gamma[l]==0 & gamma_prime[l]==1) {
    log_proposal_l_density <- 0
    for (j in 1:p_m) {
      PQ_lj <- calculate_log_PQ_lj(l, gamma_prime, eta_prime[, j], sigma2[j], tau2[, j], U, n_obs, x[, j], prob_feature_selection)
      log_proposal_l_density <- log_proposal_l_density + eta_prime[l,j] * PQ_lj[1] + (1-eta_prime[l,j]) * PQ_lj[2]
    }
    return(log_proposal_l_density)
  } else {
    stop("Error in log_proposal_l_density evaluation.")
  }
}

# Define MCMC
# seed=1
sample_gamma_eta <- function(job, seed, r, n_obs, p_m, 
                             prob_component_selection, prob_feature_selection, 
                             n_iterations, x, sigma2, tau2, U, verbose=FALSE) {
  cat("Job", job, "\n")
  # Set seed for starting point
  set.seed(seed)
  # Set parameters' starting point
  gamma_chain <- matrix(nrow = r, ncol = n_iterations)
  eta_chain <- array(dim = c(r, p_m, n_iterations))
  gamma <- initialize_gamma()
  eta <- initialize_eta(gamma)
  
  if (verbose==TRUE) {
    initial_conditions <- list(gamma=gamma, eta=eta, sigma2=sigma2, tau2=tau2, U=U)
    # r x n_iteration proposals and accept/ rejects occur, so last two dimensions represent this
    log_proposal_forward_chain <- matrix(nrow = r, ncol = n_iterations)
    log_proposal_backward_chain <- matrix(nrow = r, ncol = n_iterations)
    log_target_chain <- matrix(nrow = r, ncol = n_iterations)
    log_target_new_chain <- matrix(nrow = r, ncol = n_iterations)
    log_acceptance_ratio_chain <- matrix(nrow = r, ncol = n_iterations)
    gamma_new_chain <- array(dim = c(r, r, n_iterations))
    eta_new_chain <- array(dim = c(r, p_m, r, n_iterations))
  }
  
  # Sample for n_iterations
  for (iter in 1:n_iterations) {
    cat("iter", iter, "\n")
    # Sample across components
    for (l in 1:r) {
      # Propose values
      gamma_new <- replace(gamma, l, 1-gamma[l])
      eta_new <- eta
      if (gamma_new[l] == 0) {
        eta_new[l, ] <- rep(0, p_m)
      } else {
        eta_new[l, ] <- rbinom(n = p_m, size = 1, prob = prob_feature_selection)
      }
      # Calculate log acceptance ratio
      log_target <- log_target_density_l(l, gamma, eta, sigma2, tau2, U, n_obs, p_m, x, prob_component_selection, prob_feature_selection)
      log_target_new <- log_target_density_l(l, gamma_new, eta_new, sigma2, tau2, U, n_obs, p_m, x, prob_component_selection, prob_feature_selection)
      log_proposal_forward <- log_proposal_l_density(l, gamma_new, eta_new, gamma, eta, sigma2, tau2, U, n_obs, p_m, x, prob_feature_selection)
      log_proposal_backward <- log_proposal_l_density(l, gamma, eta, gamma_new, eta_new, sigma2, tau2, U, n_obs, p_m, x, prob_feature_selection)
      log_acceptance_ratio <- log_target_new + log_proposal_backward - log_target - log_proposal_forward
      
      if (verbose==TRUE) {
        # Note, these items are available in the stored images
        print("gamma:")
        print(gamma)
        print("gamma new:")
        print(gamma_new)
        print("log_acceptance_ratio:")
        print(log_acceptance_ratio)
        log_proposal_backward_chain[l, iter] <- log_proposal_backward
        log_proposal_forward_chain[l, iter] <- log_proposal_forward
        log_target_chain[l,iter] <- log_target_new
        log_target_new_chain[l,iter] <- log_target_new
        log_acceptance_ratio_chain[l,iter] <- log_acceptance_ratio
        gamma_new_chain[,l,iter] <- gamma_new
        eta_new_chain[,,l,iter] <- eta_new
      }
      
      if (log(runif(1)) < log_acceptance_ratio) {
        # Accept proposed gamma and eta
        gamma_chain[,iter] <- gamma <- gamma_new
        eta_chain[,,iter] <- eta <- eta_new
      } else {
        # Reject proposed gamma and eta
        gamma_chain[,iter] <- gamma
        eta_chain[,,iter] <- eta
      }
      if (gamma[l]==1) {
        # Gibbs sample to mix feature activation indicators
        # TODO remove for loops. Note, calculating the thresholds should be step 1 and accept/ reject step 2.
        eta_lj_thresholds <- rep(NA, p_m) 
        for (j in 1:p_m) {
          eta_lj_thresholds[j] <- calculate_eta_lj_threshold(l, gamma, eta[, j], sigma2[j], tau2[, j], U, n_obs, x[, j], prob_feature_selection)
        }
        for (j in 1:p_m) {
          u <- runif(1)
          if (log(u/(1-u)) < eta_lj_thresholds[j]) {
            # Turn on feature
            eta[l,j] <- 1
          } else {
            # Turn off feature
            eta[l,j] <- 0
          }
        }
        # Store Gibbs feature activation mixing result
        eta_chain[,,iter] <- eta
      }
    }
    # Store intermediates
    if (iter %% 500 == 0) {
      # Save progress
      fileConn <- file(paste0("data/", Sys.Date(), "_job_", job, "_seed_", seed, "_progress.txt"))
      writeLines(paste("Job", job, "seed", seed, "is on iter", iter), fileConn)
      close(fileConn)
      # Save output
      output_name <- paste0("data/", Sys.Date(), "_job_", job, "_seed_", seed, "_gamma_chain.rds")
      saveRDS(gamma_chain, file = output_name)
      output_name <- paste0("data/", Sys.Date(), "_job_", job, "_seed_", seed, "_eta_chain.rds")
      saveRDS(eta_chain, file = output_name)
      save.image(file = paste0("data/", Sys.Date(), "_job_", job, "_seed_", seed, "_image.RData"))
    }
  }
  return(list(gamma_chain=gamma_chain, eta_chain=eta_chain))
}

job_results <- list()
# TODO parallelize chain starting points
for (s in 1:3) {
  job_results[[s]] <- sample_gamma_eta(job, seed=s, r, n_obs, p_m, 
                                      prob_component_selection, prob_feature_selection,
                                      n_iterations, x, sigma2, tau2, U, verbose=TRUE)
}

# Save output
output_name <- paste0("data/", Sys.Date(), "_job_", job, "_results.rds")
saveRDS(job_results, file = output_name)
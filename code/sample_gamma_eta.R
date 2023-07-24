# Set-up environment
library(tidyverse)

dev <- TRUE # Flag for indicating active development
verbose <- TRUE 

if (dev == TRUE) { 
  seed <- 1 } else {
    # Get CLI arguments
    args = commandArgs(trailingOnly=TRUE)
    # test if there is at least one argument: if not, return an error
    if (length(args)==0) {
      stop("At least one argument must be supplied (input file).n", call.=FALSE)
    } else if (length(args)==1) {
      # default output file
      args[2] = "out.txt"
    }
    seed <- as.numeric(args[1])
  }

set.seed(seed)

# Load simulated data
simulation_results <- readRDS("data/2023-07-03_simulation_results.RDS")
data_list <- readRDS("data/2023-07-03_simulation_data_list.rds")

# Start with the 1st view
m <- 1 
x <- data_list[[m]]
n <- sapply(data_list, nrow)[m]
p_m <- sapply(data_list, ncol)[m]

# Specify parameters
r <- 4 
prior_component_selection <- 0.5
prior_variable_selection <- 0.05
n_sample <- 5000
n_burnin <- 1000
n_iterations <- n_sample + n_burnin

# Initialize data structures
gamma <- rbinom(n = r, size = 1, prob = prior_component_selection)
eta <- matrix(nrow = r, ncol = p_m)
for (l in 1:r) {
  if (gamma[l] == 1) {
    eta[l, ] <- rbinom(n = p_m, size = 1, prior_variable_selection)
  } else {
    eta[l, ] <- rep(0, p_m)
  }
}
gamma_chain <- matrix(nrow = r, ncol = n_iterations+1)
gamma_chain[, 1] <- gamma
eta_chain <- array(dim = c(r, p_m, n_iterations+1))
eta_chain[,,1] <- eta

sigma2 <- rep(1, p_m)
tau2 <- matrix(1, nrow = r, ncol = p_m)
U <- simulation_results$U # U is a fixed set of covariates

if (verbose==TRUE) {
  # r x n_iteration proposals and accept/ rejects occur, so last two dimensions represent this
  log_target_prime_chain <- matrix(nrow = r, ncol = n_iterations)
  log_target_chain <- matrix(nrow = r, ncol = n_iterations)
  log_acceptance_ratio_chain <- matrix(nrow = r, ncol = n_iterations)
  acceptance_indicator_chain <- matrix(nrow = r, ncol = n_iterations)
  gamma_prime_chain <- array(dim = c(r, r, n_iterations))
  eta_prime_chain <- array(dim = c(r, p_m, r, n_iterations))
}

# Specify functions
get_mvnorm_var <- function(j, gamma, U, sigma2, n) {
  Sigma2_j <- U[, gamma] %*% diag(tau2[gamma,j]) %*% t(U[, gamma]) + diag(n) 
  mvnorm_var <- sigma2[j] * Sigma2_j
  return(mvnorm_var)
}

get_logG <- function(gamma, x_j, mvnorm_var, 
                     prior_variable_selection, n) {
  log_dmvnorm_j <- mvtnorm::dmvnorm(x = x_j, mean = rep(0, n), 
                                    sigma = mvnorm_var, log = TRUE)
  n_components_active <- sum(gamma)
  # TODO implement fix RE probability of features being on
  log_p <- rep(prior_variable_selection, n_components_active) %>% 
    log() %>% sum()
  return(log_dmvnorm_j + log_p)
}

get_P_lj <- function(log_G0, log_G1) {
  x_0 <- max(log_G0, log_G1)
  x <- log_G1 - x_0
  y <- log_G0 - x_0
  return(exp(x) / (exp(x) + exp(y)))
}

log_target_density <- function(gamma, x_list, mvnorm_var_list,
                               prior_component_selection, prior_variable_selection,
                               r, n) {
  log_prod_G <- mapply(get_logG, gamma, 
                       x_list, mvnorm_var_list, 
                       prior_variable_selection, n) %>% sum()
  # TODO update calculation of log_prod_prior_component_selection because gamma is bernoulli! 
  log_prod_prior_component_selection <- prior_component_selection %>% 
    rep(r) %>% log() %>% sum()
  return(log_prod_G + log_prod_prior_component_selection)
}

# Reshape input
x_list <- as.list(as.data.frame(x))

# Begin MCMC
if (dev == TRUE) {
  iter <- 1
  l <- 1
  j <- 1
}

for (iter in 1:n_iterations) {
  
  if (verbose==TRUE) { print(iter) }
  
  # Initialize candidate values
  gamma_prime <- gamma
  eta_prime <- eta
  
  for (l in 1:r) {
    
    if (verbose==TRUE) { 
      print("component:")
      print(l) 
      print("gamma:")
      print(gamma)
      }
    
    mvnorm_var_list <- lapply(1:p_m, get_mvnorm_var, gamma, U, sigma2, n)

    if (gamma[l] == 1) {
      # Propose deactivation
      gamma_prime[l] <- 0
      eta_prime[l,] <- rep(0, p_m)
    } else if (gamma[l] == 0) {
      # Propose activation
      gamma_prime[l] <- 1
      for (j in 1:p_m) {
        if (verbose==TRUE) { print(j) }
        eta_1 <- eta_prime[, j]
        eta_1[l] <- 1
        eta_0 <- eta_prime[, j]
        eta_0[l] <- 0
        log_G1 <- get_logG(gamma_prime, x[, j], mvnorm_var_list[[j]], 
                           prior_variable_selection, n)
        # TODO understand if logG calculation should be truly independent of eta
        log_G0 <- get_logG(gamma_prime, x[, j], mvnorm_var_list[[j]], 
                           prior_variable_selection, n)
        P_lj <- get_P_lj(log_G0, log_G1)
        if (verbose==TRUE) { print(P_lj) }
        eta_prime[l,j] <- rbinom(1, 1, prob = P_lj)
      }
    }
    
    mvnorm_var_list_prime <- lapply(1:p_m, get_mvnorm_var, gamma_prime, U, sigma2, n)
    
    # TODO verify log_target calculations are correct
    log_target_prime <- log_target_density(gamma_prime, x_list, 
                                           mvnorm_var_list_prime, 
                                           prior_component_selection, 
                                           prior_variable_selection, r, n)
    log_target_previous <- log_target_density(gamma, x_list, 
                                              mvnorm_var_list, 
                                              prior_component_selection, 
                                              prior_variable_selection, r, n)
    
    log_acceptance_ratio <- log_target_prime - log_target_previous
    
    if (verbose==TRUE) { 
      print("gamma_prime:")
      print(gamma_prime)
      print("log_target_prime:")
      print(log_target_prime)
      print("log_target_previous:")
      print(log_target_previous)
      print("log_acceptance_ratio:")
      print(log_acceptance_ratio)
      log_target_prime_chain[l,iter] <- log_target_prime
      log_target_chain[l,iter] <- log_target_previous
      log_acceptance_ratio_chain[l,iter] <- log_acceptance_ratio
      gamma_prime_chain[,l,iter] <- gamma_prime
      eta_prime_chain[,,l,iter] <- eta_prime
    }
    
    if (log(runif(1)) < log_acceptance_ratio) {
      # Accept
      gamma_chain[,iter] <- gamma_prime
      eta_chain[,,iter] <- eta_prime
      mvnorm_var_list <- mvnorm_var_list_prime
      if (verbose == TRUE) {
        acceptance_indicator_chain[l,iter] <- 1
      }
    } else {
      # Reject
      gamma_chain[,iter] <- gamma_chain[,iter-1]
      eta_chain[,,iter] <- eta_chain[,,iter-1]
      if (verbose == TRUE) {
        acceptance_indicator_chain[l,iter] <- 0
      }
    }
    
  }
  
  if (iter %% 500 == 0) {
    # Save progress
    fileConn <- file(paste0("data/", Sys.Date(), "_progress_seed_", seed, ".txt"))
    writeLines(paste("Job seed", seed, "is on iter", iter), fileConn)
    close(fileConn)
    # Save output
    output_name <- paste0("data/", Sys.Date(), "_gamma_chain_seed_", seed, ".rds")
    saveRDS(gamma_chain, file = output_name)
    output_name <- paste0("data/", Sys.Date(), "_eta_chain_seed_", seed, ".rds")
    saveRDS(eta_chain, file = output_name)
    if (verbose==TRUE) {
      output_name <- paste0("data/", Sys.Date(), "_log_target_list_", seed, ".rds")
      saveRDS(list(log_target_prime_chain, log_target_chain, log_acceptance_ratio), 
              file = output_name)
      saveRDS(acceptance_indicator_chain, file = "data/acceptance_indicator_chain.rds")
      saveRDS(gamma_prime_chain, file = "data/gamma_prime_chain.rds")
      saveRDS(eta_prime_chain, file = "data/eta_prime_chain.rds")
      save.image()
    }
  }
  
}

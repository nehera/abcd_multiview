library(tidyverse)

set.seed(1)
verbose <- TRUE

data_list <- readRDS("data/2023-07-03_simulation_data_list.rds")

# TODO start with 1st view
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
gamma_previous <- rbinom(n = r, size = 1, prob = prior_component_selection)
eta_previous <- matrix(nrow = r, ncol = p_m)
for (l in 1:r) {
  if (gamma_previous[l] == 1) {
    eta_previous[l, ] <- rbinom(n = p_m, size = 1, prior_variable_selection)
  } else {
    eta_previous[l, ] <- rep(0, p_m)
  }
}
gamma_chain <- matrix(nrow = r, ncol = n_iterations)
gamma_chain[, 1] <- gamma_previous
eta_chain <- array(dim = c(r, p_m, n_iterations))
eta_chain[,,1] <- eta_previous

sigma2 <- rep(1, p_m)
tau2 <- matrix(1, nrow = r, ncol = p_m)
U <- matrix(rnorm(n*r), nrow = n, ncol = r)

# Specify functions
get_mvnorm_var <- function(j, gamma, U, sigma2, n) {
  Sigma2_j <- U[, gamma] %*% diag(tau2[gamma,j]) %*% t(U[, gamma]) + diag(n) # TODO store as array for acceptance/ rejection
  mvnorm_var <- sigma2[j] * Sigma2_j
  return(mvnorm_var)
}

get_logG <- function(gamma, x_j, mvnorm_var, 
                     prior_variable_selection, n) {
  log_dmvnorm_j <- mvtnorm::dmvnorm(x = x_j, mean = rep(0, n), sigma = mvnorm_var, log = TRUE)
  n_components_active <- sum(gamma)
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
  log_prod_prior_component_selection <- prior_component_selection %>% 
    rep(r) %>% log() %>% sum()
  return(log_prod_G + log_prod_prior_component_selection)
}

# Reshape input
x_list <- as.list(as.data.frame(x))

# Begin MCMC
for (iter in 2:n_iterations) {
  
  if (verbose==TRUE) { print(iter) }
  
  # Initialize candidate values
  gamma_prime <- gamma_previous
  eta_prime <- eta_previous
  
  for (l in 1:r) {
    
    if (verbose==TRUE) { 
      print(l) 
      print("gamma_previous:")
      print(gamma_previous)
      print("eta_previous:")
      print(eta_previous)
      }
    
    mvnorm_var_list <- lapply(j:p_m, get_mvnorm_var, gamma_previous, U, sigma2, n)

    if (gamma[l] == 1) {
      # Propose deactivation
      gamma_prime[l] <- 0
      eta_prime[l, ] <- rep(0, p_m)
    } else if (gamma[l] == 0) {
      # Propose activation
      gamma_prime[l] <- 1
      for (j in 1:p_m) {
        # TODO propose based on P_lj
        # TODO start with 1st feature
        j <- 1
        if (verbose==TRUE) { print(j) }
        gamma_1 <- gamma_prime # TODO understand if fine, but should be since set to 1 pre-feature loop?
        eta_1 <- eta_prime[, j]
        eta_1[l] <- 1
        eta_0 <- eta_prime[, j]
        eta_0[l] <- 0
        log_G1 <- get_logG(gamma_1, eta_1, x[, j], mvnorm_var_list[[j]], 
                           prior_variable_selection, n, r)
        log_G0 <- get_logG(gamma_1, eta_0, x[, j], mvnorm_var_list[[j]], 
                           prior_variable_selection, n, r)
        P_lj <- get_P_lj(log_G0, log_G1)
        if (verbose==TRUE) { print(P_lj) }
        eta_prime[l,j] <- rbinom(1, 1, prob = P_lj)
      }
    }
    
    mvnorm_var_list_prime <- lapply(j:p_m, get_mvnorm_var, gamma_prime, U, sigma2, n)
    
    # TODO verify log_target calculations are correct
    log_target_prime <- log_target_density(gamma_prime, x_list, 
                                           mvnorm_var_list_prime, 
                                           prior_component_selection, 
                                           prior_variable_selection, r, n)
    log_target_previous <- log_target_density(gamma_previous, x_list, 
                                              mvnorm_var_list_previous, 
                                              prior_component_selection, 
                                              prior_variable_selection, r, n)
    
    log_acceptance_ratio <- log_target_prime - log_target_previous
    
    if (verbose==TRUE) { 
      print("gamma_prime:")
      print(gamma_prime)
      print("eta_prime:")
      print(eta_prime)
      print("log_target_prime:")
      print(log_target_prime)
      print("log_target_previous:")
      print(log_target_previous)
      print("log_acceptance_ratio:")
      print(log_acceptance_ratio)
    }
    
    if (log(runif(1)) < log_acceptance_ratio) {
      # Accept
      gamma_chain[,iter] <- gamma_prime
      eta_chain[,,iter] <- eta_prime
      mvnorm_var_list_previous <- mvnorm_var_list_prime
    } else {
      # Reject
      gamma_chain[,iter] <- gamma_chain[,iter-1]
      eta_chain[,,iter] <- eta_chain[,,iter-1]
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
  }
  
}

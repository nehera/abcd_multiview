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

sigma2 <- rep(1, p_m)
tau2 <- matrix(1, nrow = r, ncol = p_m)
U <- matrix(rnorm(n*r), nrow = n, ncol = r)

# Specify functions
get_logG <- function(gamma_j, eta_j, x_j, mvnorm_var, 
                     prior_variable_selection, n, r) {
  log_dmvnorm_j <- mvtnorm::dmvnorm(x = x_j, mean = rep(0, n), sigma = mvnorm_var, log = TRUE)
  n_features_active <- sum(eta_j)
  log_prod <- n_features_active * log(prior_variable_selection) + 
    (r-n_features_active) * log(1-prior_variable_selection)
  return(log_dmvnorm_j + log_prod)
}

get_P_lj <- function(log_G0, log_G1) {
  x_0 <- max(log_G0, log_G1)
  x <- log_G1 - x_0
  y <- log_G0 - x_0
  return(exp(x) / (exp(x) + exp(y)))
}

#### YOU ARE HERE
log_target_density <- function(gamma_prime, eta_prime, U_iter, Sigma2_m, Tau2_m, 
                               X_m, prior_component_selection, prior_variable_selection) {
  r <- ncol(U_iter)
  p_m <- length(Sigma2_m)
  log_prod_prior_component_selection <- prior_component_selection %>% 
    rep(r) %>% log() %>% sum()
  logG <- numeric(p_m)
  for (j in 1:p_m) { # TODO: Remove for loop
    logG[j] <- get_logG(j, gamma_prime, eta_prime[,j], U_iter, 
                        Sigma2_m, Tau2_m, X_m, prior_variable_selection)
  }
  return(sum(logG))
}

# Begin mcmc step 1
gamma_prime <- gamma
eta_prime <- eta
for (l in 1:r) {
  # TODO start with the 1st component
  l <- 1
  if (verbose==TRUE) { print(l) }
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
      Sigma2_j <- U[, gamma] %*% diag(tau2[gamma,j]) %*% t(U[, gamma]) + diag(n)
      mvnorm_var <- sigma2[j] * Sigma2_j
      gamma_1 <- gamma_prime # TODO understand if fine, but should be since set to 1 pre-feature loop?
      eta_1 <- eta_prime[, j]
      eta_1[l] <- 1
      eta_0 <- eta_prime[, j]
      eta_0[l] <- 0
      log_G1 <- get_logG(gamma_1, eta_1, x[, j], mvnorm_var, 
               prior_variable_selection, n, r)
      log_G0 <- get_logG(gamma_1, eta_0, x[, j], mvnorm_var, 
                         prior_variable_selection, n, r)
      P_lj <- get_P_lj(log_G0, log_G1)
      if (verbose==TRUE) { print(P_lj) }
      eta_prime[l,j] <- rbinom(1, 1, prob = P_lj)
    }
  }
  
  # TODO accept/ reject
  
  
}

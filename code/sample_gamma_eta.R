# Set-up environment
library(tidyverse)

dev <- FALSE # Flag for indicating active development
verbose <- TRUE 

if (dev == TRUE) { 
  seed <- 1 
  } else {
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
simulation_results <- readRDS("data/2023-07-03_simulation_results.rds")
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
if (dev==TRUE) { n_iterations <- 3 }

# Initialize data structures
gamma_chain <- matrix(nrow = r, ncol = n_iterations)
eta_chain <- array(dim = c(r, p_m, n_iterations))

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
U <- simulation_results$U # U is a fixed set of covariates

log_target_chain <- numeric(n_iterations)

# Specify functions
get_mvnorm_var <- function(j, gamma, U, sigma2, n) {
  Sigma2_j <- U[, gamma, drop = FALSE] %*% diag(tau2[gamma,j]) %*% t(U[, gamma, drop = FALSE]) + diag(n) 
  mvnorm_var <- sigma2[j] * Sigma2_j
  return(mvnorm_var)
}

get_log_dmvnorm_j <- function(x_j, mvnorm_var_j, n) {
  log_dmvnorm_j <- mvtnorm::dmvnorm(x = x_j, mean = rep(0, n), 
                          sigma = mvnorm_var_j, log = TRUE)
  return(log_dmvnorm_j)
}

get_log_dmvnorm_vector <- function(gamma, x, U, sigma2, p_m, n) {
  x_list <- as.list(as.data.frame(x))
  # TODO breakout from function to reduce redundant calculations
  mvnorm_var_list <- lapply(1:p_m, get_mvnorm_var, gamma, U, sigma2, n) 
  log_dmvnorm_vector <- mapply(get_log_dmvnorm_j, x_list, mvnorm_var_list, n)
  return(unname(log_dmvnorm_vector))
}

get_log_G_j <- function(log_dmvnorm_j, eta_j, prior_variable_selection, r) {
  n_components_feature_active_in <- sum(eta_j)
  log_p_eta_j <- n_components_feature_active_in * log(prior_variable_selection) + 
    (r - n_components_feature_active_in) * log(1 - prior_variable_selection)
  return(log_dmvnorm_j + log_p_eta_j)
}

get_log_G <- function(eta, log_dmvnorm_vector, prior_variable_selection, r, p_m) {
  n_active_features <- sum(eta)
  log_G <- sum(log_dmvnorm_vector) +
    log(prior_variable_selection) * n_active_features +
    log(1 - prior_variable_selection) * (r * p_m - n_active_features) # TODO verify whether this line is appropriate
  return(log_G)
}

get_P_lj <- function(log_G0, log_G1) {
  x_0 <- max(log_G0, log_G1)
  x <- log_G1 - x_0
  y <- log_G0 - x_0
  return(exp(x) / (exp(x) + exp(y)))
}

get_P_l <- function(l, log_dmvnorm_vector, eta, r, p_m) {
  P_l <- rep(NA, p_m)
  for (j in 1:p_m) {
    eta_1 <- eta[, j]
    eta_1[l] <- 1
    eta_0 <- eta[, j]
    eta_0[l] <- 0
    log_G1 <- get_log_G_j(log_dmvnorm_vector[j], eta_1, 
                          prior_variable_selection, r)
    log_G0 <- get_log_G_j(log_dmvnorm_vector[j], eta_0, 
                          prior_variable_selection, r)
    P_l[j] <- get_P_lj(log_G0, log_G1)
  }
  return(P_l)
}

# TODO should I be calculating this just for the lth component worth of data? 
log_target_density <- function(gamma, eta, log_dmvnorm_vector, 
                               prior_component_selection, 
                               prior_variable_selection, r, p_m) {
  log_G <- get_log_G(eta, log_dmvnorm_vector, prior_variable_selection, r, p_m)
  n_active_components <- sum(gamma)
  log_prod_p_gamma <- n_active_components * log(prior_component_selection) +
    (r - n_active_components) * log(1-prior_component_selection)
  return(log_G + log_prod_p_gamma)
}

log_proposal_l_density <- function(gamma_l, eta_l, 
                                   gamma_l_prime, eta_l_prime, P_l) {
  if (gamma_l==1 & gamma_l_prime==0 & sum(eta_l_prime)==0) {
    return(0)
  } else if (gamma_l==0 & gamma_l_prime==1) {
    return(
      ifelse(eta_l_prime==1, log(P_l), 
             1 - log(P_l)) %>% sum()
    )
  } else {
    stop("Error in log_proposal_l_density evaluation possibly due to issue in gamma and eta proposal.")
  }
}

# Begin MCMC
log_dmvnorm_vector <- get_log_dmvnorm_vector(gamma, x, U, sigma2, p_m, n)

log_target <- log_target_density(gamma, eta, log_dmvnorm_vector, 
                                 prior_component_selection, 
                                 prior_variable_selection, r, p_m)

if (verbose==TRUE) {
  initial_conditions <- list(gamma=gamma, eta=eta, sigma2=sigma2, tau2=tau2, U=U,
                             log_dmvnorm_vector=log_dmvnorm_vector, log_target=log_target)
  # r x n_iteration proposals and accept/ rejects occur, so last two dimensions represent this
  log_proposal_forward_chain <- matrix(nrow = r, ncol = n_iterations)
  log_proposal_reverse_chain <- matrix(nrow = r, ncol = n_iterations)
  log_target_new_chain <- matrix(nrow = r, ncol = n_iterations)
  log_acceptance_ratio_chain <- matrix(nrow = r, ncol = n_iterations)
  acceptance_indicator_chain <- matrix(nrow = r, ncol = n_iterations)
  gamma_new_chain <- array(dim = c(r, r, n_iterations))
  eta_new_chain <- array(dim = c(r, p_m, r, n_iterations))
  P_chain <- array(dim = c(r, p_m, n_iterations))
}

for (iter in 1:n_iterations) {
  
  if (verbose==TRUE) {
    print("iter:")
    print(iter) 
    }

  for (l in 1:r) {
    
    if (verbose==TRUE) { 
      print("component:")
      print(l) 
      print("gamma previous:")
      print(gamma)
    }
    
    # Initialize from previous result
    gamma_new <- gamma
    eta_new <- eta
    
    if (gamma_new[l] == 1) {
      # Propose component deactivation
      gamma_new[l] <- 0
      log_dmvnorm_vector <- get_log_dmvnorm_vector(gamma_new, x, U, sigma2, p_m, n)
      # Propose feature deactivation
      eta_new[l,] <- rep(0, p_m)
    } else {
      # Propose component activation
      gamma_new[l] <- 1
      log_dmvnorm_vector <- get_log_dmvnorm_vector(gamma_new, x, U, sigma2, p_m, n)
      # Propose feature activation
      P_l <- get_P_l(l, log_dmvnorm_vector, eta_new, r, p_m)
      if (verbose == TRUE) {
        P_chain[l, , iter] <- P_l
      }
      eta_new[l, ] <- rbinom(n = p_m, size = 1, prob = P_l)
    }
    
    # Calculate log proposal densities
    log_proposal_forward <- log_proposal_l_density(gamma[l], eta[l, ], 
                                           gamma_new[l], eta_new[l, ], P_l)
    
    log_proposal_reverse <- log_proposal_l_density(gamma_new[l], eta_new[l, ], 
                           gamma[l], eta[l, ], P_l)
    
    # Calculate log target density
    log_target_new <- log_target_density(gamma_new, eta_new, log_dmvnorm_vector, 
                                     prior_component_selection, 
                                     prior_variable_selection, r, p_m)
    
    log_acceptance_ratio <- log_target_new + log_proposal_reverse - log_target - log_proposal_forward
    
    if (verbose==TRUE) {
      print("gamma new:")
      print(gamma_new)
      print("log_acceptance_ratio:")
      print(log_acceptance_ratio)
      log_proposal_reverse_chain[l, iter] <- log_proposal_reverse
      log_proposal_forward_chain[l, iter] <- log_proposal_forward
      log_target_new_chain[l,iter] <- log_target_new
      log_acceptance_ratio_chain[l,iter] <- log_acceptance_ratio
      gamma_new_chain[,l,iter] <- gamma_new
      eta_new_chain[,,l,iter] <- eta_new
    }
    
    if (log(runif(1)) < log_acceptance_ratio) {
      # Accept gamma, eta, log_target
      gamma <- gamma_chain[,iter] <- gamma_new
      eta <- eta_chain[,,iter] <- eta_new
      log_target <- log_target_chain[iter] <- log_target_new
      if (verbose == TRUE) {
        acceptance_indicator_chain[l,iter] <- 1
      }
    } else {
      # Reject by maintaining previous result
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
      save.image(file = paste0("data/", Sys.Date(), "_image_", seed, ".RData"))
    }
  }
  
}

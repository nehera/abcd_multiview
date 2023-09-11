library(tidyverse)

r <- 4
n_obs <- 200
p_m <- 10
prob_component_selection <- 0.5
prob_feature_selection <- 0.5

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

set.seed(1)
gamma <- initialize_gamma()
eta <- initialize_eta(gamma)

sigma2 <- rep(1, p_m)
tau2 <- matrix(1, nrow = r, ncol = p_m)
simulation_results <- readRDS("data/2023-09-08_simulation_results.rds")
U <- simulation_results$U # U is a fixed set of covariates

m <- 1 # Start with the 1st view
data_list <- readRDS("data/2023-09-08_simulation_data_list.rds")
x <- data_list[[m]]

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

calculate_P_lj <- function(l, gamma, eta_j, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection) {
  gamma_1 <- replace(gamma, l, 1)
  eta_1 <- replace(eta_j, l, 1)
  eta_0 <- replace(eta_j, l, 0)
  log_G_j_1 <- calculate_log_G_j(gamma_1, eta_1, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  log_G_j_0 <- calculate_log_G_j(gamma_1, eta_0, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  max_arg <- max(log_G_j_1, log_G_j_0) # Use logsumexp trick
  a <- log_G_j_1 - max_arg
  b <- log_G_j_0 - max_arg
  P_lj <- exp(a) / (exp(a)+exp(b))
  return(P_lj)
}

## -- Calculate log_target_density_l 
log_target_density_l <- function(l, gamma, eta, sigma2, tau2, U, n_obs, p_m, 
                                 x, prob_component_selection, prob_feature_selection) {
  sum_log_target_density_lj <- 0
  for (j in 1:p_m) {
    log_dmvnorm_j <- calculate_log_dmvnorm_j(eta[, j], sigma2[j], tau2[,j], U, n_obs, x[,j])
    # TODO understand if the prior probablity or P_lj should be used in place of prob_feature_selection
    sum_log_target_density_lj <- sum_log_target_density_lj + log_dmvnorm_j + 
      eta[l,j]*log(prob_feature_selection) + (1-eta[l,j])*log(1-prob_feature_selection)
  }
  log_target_density_l <- gamma[l]*log(prob_component_selection) + 
    (1-gamma[l])*log(1-prob_component_selection) + sum_log_target_density_lj
  return(log_target_density_l)
}


## Testing
# log_target_0 <- log_target_density_l(l, gamma, eta, sigma2, tau2, U, n_obs, p_m, x, prob_component_selection, prob_feature_selection)
# 
# gamma_1 <- replace(gamma, l, 1)
# eta_1 <- eta
# eta_1[l, ] <- rbinom(n = p_m, size = 1, prob = prob_feature_selection)
# 
# log_target_1 <- log_target_density_l(l, gamma_1, eta_1, sigma2, tau2, U, n_obs, p_m, x, prob_component_selection, prob_feature_selection)
# 
# log_target_0 - log_target_1 # Note, the absolute difference is MUCH smaller in this implementation :-) 

# TODO fix bug in log_proposal_backward
# log_proposal_backward <- log_proposal_l_density(l, gamma, eta, gamma_new, eta_new, sigma2, tau2, U, n_obs, p_m, x, prob_feature_selection)

## -- Calculate log_proposal_l
log_proposal_l_density <- function(l, gamma_prime, eta_prime, gamma, eta,
                                   sigma2, tau2, U, n_obs, p_m, x, prob_feature_selection) {
  if (gamma[l]==1 & gamma_prime[l]==0 & sum(eta_prime[l,])==0) {
    return(0)
  } else if (gamma[l]==0 & gamma_prime[l]==1) {
    log_proposal_l_density <- 0
    for (j in 1:p_m) {
      P_lj <- calculate_P_lj(l, gamma_prime, eta_prime[, j], sigma2[j], tau2[, j], U, n_obs, x[, j], prob_feature_selection)
      # TODO understand if any P_lj's should be 1 or 0 and if so, what a good value for c is
      c <- 10^-6
      if (P_lj==1) {
        P_lj <- P_lj - c
      } else if (P_lj==0) {
        P_lj <- P_lj + c
      }
      log_proposal_l_density <- log_proposal_l_density + eta_prime[l,j]*log(P_lj) + (1-eta_prime[l,j])*log(1-P_lj) 
    }
    return(log_proposal_l_density)
  } else {
    stop("Error in log_proposal_l_density evaluation.")
  }
}

## Testing
# log_proposal_forward <- log_proposal_l_density(l, gamma_1, eta_1, gamma, eta, sigma2, tau2, U, n_obs, p_m, x, prob_feature_selection)
# log_proposal_backward <- log_proposal_l_density(l, gamma, eta, gamma_1, eta_1, sigma2, tau2, U, n_obs, p_m, x, prob_feature_selection)
# log_proposal_backward - log_proposal_forward
# 
# log_A <- log_target_1 + log_proposal_backward - log_target_0 - log_proposal_forward
# log_A

## -- MCMC 2
n_sample <- 5000
n_burnin <- 1000
n_iterations <- n_sample + n_burnin

set.seed(1)
gamma_chain <- matrix(nrow = r, ncol = n_iterations)
eta_chain <- array(dim = c(r, p_m, n_iterations))
gamma <- initialize_gamma()
eta <- initialize_eta(gamma)

for (iter in 1:n_iterations) {
  print(paste("iter:", iter))
  for (l in 1:r) {
    # Propose values
    gamma_new <- replace(gamma, l, 1-gamma[l])
    eta_new <- eta
    if (gamma_new[l] == 0) {
      eta_new[l, ] <- rep(0, p_m)
    } else {
      eta_new[l, ] <- rbinom(n = p_m, size = 1, prob = prob_feature_selection) # Understand if P_lj's should be used
    }
    # Calculate log acceptance ratio
    log_target <- log_target_density_l(l, gamma, eta, sigma2, tau2, U, n_obs, p_m, x, prob_component_selection, prob_feature_selection)
    log_target_new <- log_target_density_l(l, gamma_new, eta_new, sigma2, tau2, U, n_obs, p_m, x, prob_component_selection, prob_feature_selection)
    log_proposal_forward <- log_proposal_l_density(l, gamma_new, eta_new, gamma, eta, sigma2, tau2, U, n_obs, p_m, x, prob_feature_selection)
    log_proposal_backward <- log_proposal_l_density(l, gamma, eta, gamma_new, eta_new, sigma2, tau2, U, n_obs, p_m, x, prob_feature_selection)
    log_acceptance_ratio <- log_target_new + log_proposal_backward - log_target - log_proposal_forward
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
      # TODO remove for loops
      P_l <- rep(NA, p_m) 
      for (j in 1:p_m) {
        P_l[j] <- calculate_P_lj(l, gamma, eta[, j], sigma2[j], tau2[, j], U, n_obs, x[, j], prob_feature_selection)
      }
      for (j in 1:p_m) {
        if (runif(1) < P_l[j]) {
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
}

# Save output
# TODO remove seed hardcoding
output_name <- paste0("data/", Sys.Date(), "_gamma_chain_seed_", 1, ".rds")
saveRDS(gamma_chain, file = output_name)
output_name <- paste0("data/", Sys.Date(), "_eta_chain_seed_", 1, ".rds")
saveRDS(eta_chain, file = output_name)
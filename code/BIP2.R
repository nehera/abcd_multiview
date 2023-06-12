#### ----- Development ---- ####
library(BIPnet)
## Set Inputs
set.seed(1)
simulation_results <- Simulate(setting=1)
data_list <- list(simulation_results$X1, 
                 simulation_results$X2, 
                 simulation_results$Y)
indic_var <- c(0, 0, 1) # metadata 
method <- "BIP" # method without grouping information to start
group_list <- NULL # default when grouping information not included
## Get Results for Testing
bip_0 <- BIP(dataList = data_list, IndicVar = indic_var, Method = method)
# MCMC Parameters
n_sample <- 5000
n_burnin <- 1000
n_max_models <- 50
# Model Hyper-parameters
r <- 4 # number of components
prior_component_selection <- c(1, 1) # Hyperparameters of a beta(a,b) distribution; prior distribution of the probability of selecting components
prior_response_selection <- c(1, 1) # Hyperparameters of a beta(a,b) distribution; prior distribution of the probability of selecting the response on a component.
prior_residual_variance <- c(0.01, 0.01) # Hyperparameters of residual variance sigma_j^2(m)
prior_b0 <- c(2, 2) # Hyperparameters of a gamma distribution; prior distribution of the paramater b_l0 that controls the shrinkage of loadings when there is no grouping information
prior_b <- c(1, 1) # Hyperparameters of a gamma distribution; prior distribution of group effect coefficients blj, j = 1, ..., K_l.
prior_group_selection <- c(1, 1) # Hyperparameters of a beta(a,b) distribution; prior distribution of the probability of selecting groups.
prior_variable_selection <- 0.05 # Prior probability for variable selection

#### ----- Begin Function ---- ####
# TODO: Initiate development from here on since .RData includes everything up to this line.
library(tidyverse)
library(Rcpp)
library(mvtnorm) # for dmvnorm

# BIP <- function(data_list = data_list, # TODO: This default might cause probs
#                 indic_var = indic_var, # TODO: This default might cause probs
#                 group_list = NULL,
#                 method = method, # TODO: This default might cause probs
#                 r = 4,
#                 n_sample = 5000,
#                 n_burnin = 1000,
#                 n_max_models = 50,
#                 prior_component_selection = c(1,1),
#                 prior_response_selection = c(1,1),
#                 prior_b0 = c(2,2),
#                 prior_b = c(1,1),
#                 prior_group_selection = c(1,1),
#                 prior_variable_selection = 0.05) {

##### ---- Initialize inputs/ outputs

# Define meta-parameters
omics_index <- which(indic_var == 0)
response_index <- which(indic_var == 1)
covariate_index <- which(indic_var == 2)
n_views <- length(data_list)
n_features <- sapply(data_list, ncol)

if (method=="BIP") {
  cat("Method applied is BIP, which does not incorporate grouping information.")
  group_list <- lapply(n_features, function(n) { matrix(1, n, 1) })
} else if (method=="BIPnet") {
  cat("Method applied is BIPnet, which incorporates grouping information by the argument group_list.")
  if (is.null(group_list)) {
    stop("For method BIPnet, you must pass a list of grouping information to argument group_list.")
  }
} else {
  stop("Argument method must be either BIP or BIPnet.")
}

n_groups <- sapply(group_list, ncol)
n_observations <- sapply(data_list, nrow) %>% unique
n_iterations <- n_sample + n_burnin

## Ensure parameters passed are appropriate
if (length(n_observations) > 1) {
  stop("Each view must contain the same n_observations i.e. the same number of rows.")
}
if (n_sample <= n_burnin) {
  stop("Argument n_burnin must be smaller than argument n_sample, the number of MCMC iterations.")
}
if (n_sample <= 100) {
  stop("Please specify a larger number of MCMC iterations.")
}

# TODO: Understand how r established when passed as NULL.
if (is.null(r)) {
  mysvd <- lapply(omics_index, function(i) svd(data_list[[i]]))
  mysumsvd <- lapply(1:length(omics_index), function(i) cumsum(mysvd[[i]]$d) / max(cumsum(mysvd[[i]]$d)) * 100)
  KMax <- max(unlist(lapply(1:length(omics_index), function(i) min(which(mysumsvd[[i]]>=80, arr.ind = TRUE)))))
  r <- min(KMax+1, 10)
}
## Scale data
# TODO: Understand which views should be scaled and if we might scale for n_features
# TODO: Add back mean_data and sd_data for reporting purposes
data_list[omics_index] <- lapply(data_list[omics_index], scale)

#### ---- MCMC functions

## Functions for initialization

init_matrix_list <- function(initial_value, n_views, n_row, n_col_vector) {
  lapply(1:n_views, function(m, n_row, n_col_vector) {
    matrix(initial_value, nrow = n_row, ncol = n_col_vector[m])
  }, n_row, n_col_vector)
}

init_vector_list <- function(initial_value, n_views, len_vector) {
  lapply(1:n_views, function(m, initial_value, len_vector) {
    rep(initial_value, len_vector[m])
  }, initial_value, len_vector)
}

## Functions for sampling posterior

sample_intercept <- function(Y, A_outcome, U, Sigma2_outcome, Sigma2_0 = 100) {
  n <- nrow(Y)
  Y_star <- matrix(nrow = n, ncol = 1)
  # TODO: Remove for loop
  for (i in 1:n_observations) {
    U_star_i <- U[i,] %*% A_outcome
    Y_star[i,1] <- Y[i,1] - U_star_i
  }
  Y_star_mean <- mean(Y_star)
  invSig2 <- n / Sigma2_outcome + 1 / Sigma2_0
  intercept <- (n * Y_star_mean) / (invSig2 * Sigma2_outcome) + sqrt(1 / invSig2) * rnorm(1)
  return(intercept)
}

sample_beta_posterior <- function(n_samples = 1, prior_alpha, prior_beta, n, summation) {
  posterior_alpha <- prior_alpha + summation
  posterior_beta <- prior_beta + n - summation
  return(rbeta(n_samples, posterior_alpha, posterior_beta))
}

#### ---- MCMC initialization

# Step 0. Initialize parameters to estimate
# TODO: Initialize from priors (rather than fixed)

## Factor analysis parameters

U <- matrix(rnorm(n_observations*r),
            nrow = n_observations, ncol = r)

A <- init_matrix_list(initial_value = 0, n_views, r, n_features)

Sigma2 <- init_vector_list(initial_value = 1, n_views, n_features)

## Variable selection parameters

# Component level 
q_component_level <- rbeta(n_views, prior_component_selection[1], prior_component_selection[2]) # Probability for mth data type, a component is active

Gamma <- init_vector_list(initial_value = 1, n_views, rep(r, n_views))

# Variable level 

q_variable_level <- init_vector_list(initial_value = prior_variable_selection, n_views, rep(r, n_views)) # Probability for the mth data type, for the lth active component, a variable is active

# Fix Covariable Component and Variable Selection Probabilities to 1
if (length(covariate_index) > 0) {
  q_component_level[covariate_index] <- 1 
  q_variable_level[[covariate_index]] <- rep(1, r)
}

Eta <- init_matrix_list(initial_value = 1, n_views, r, n_features) # Thierry's code seemingly refers to eta parameters as "rho"

Tau <- init_matrix_list(initial_value = 1, n_views, r, n_features)

## Grouping information parameters

Lambda <- init_matrix_list(initial_value = 1, n_views, r, n_features)

b_0 <- init_vector_list(initial_value = 0.1, n_views, rep(r, n_views))

b_l_dot <- init_matrix_list(initial_value = 0.1, n_views, r, n_groups) # Group effect on component

R <- init_matrix_list(initial_value = 1, n_views, r, n_groups) # Binary matrices that denotes if group k contributes to active component l

# TODO: Consider initiation from rbeta
q_r <- init_vector_list(initial_value = 0.5, n_views, rep(r, n_views)) # Probability for the mth data type, for the lth component, a group is active

t <- 1 # TODO: Remove post-development

#### ---- MCMC sample posteriors

for (t in 1:n_iterations) {
  
  # Step 0. Sample intercept for response variable
  
  intercept <- sample_intercept(Y = data_list[[response_index]],
                   A_outcome = A[[response_index]], U,
                   Sigma2_outcome = Sigma2[[response_index]])
  
  # Adjust Y for subsequent sampling
  Y_prime <- data_list[[response_index]] - intercept
  
  # Step 1. Sample component and variable selection probabilities
  Gamma_sums <- sapply(Gamma, sum)
  for (m in c(omics_index, response_index)) {
    if (m %in% omics_index) {
      q_component_level[m] <- sample_beta_posterior(prior_alpha = prior_component_selection[1],
                                                    prior_beta = prior_component_selection[2],
                                                    n = r, summation = Gamma_sums[m])
      
      # Probability of variable selection fixed for omics data by prior_variable_selection
    } else if (m %in% response_index) {
      q_component_level[m] <- sample_beta_posterior(prior_alpha = prior_response_selection[1],
                                                    prior_beta = prior_response_selection[2],
                                                    n = r, summation = Gamma_sums[m])
      q_variable_level[[m]] <- rep(q_component_level[m], r)
    } 
  }
  
  m <- 1 # TODO: Remove post-dev
  
  for (m in omics_index) {

    X_m <- data_list[[m]]
    l <- 1 # TODO: Remove post-dev
  
  # Step 2. Sample component, feature activation parameters (gamma, eta respectively) by Metropolis-Hastings
    
    for (l in 1:r) {
      gamma_old <- Gamma[[m]][l]
      eta_old <- Eta[[m]][l,]
      # Propose new values
      if (gamma_old == 1) {
        # Propose deactivation
        gamma_prime <- 0
        eta_prime <- rep(0, length(eta_old))
      } else if (gamma_old == 0) {
        # Propose activation
        gamma_prime <- 1
        # Sample eta_prime
        get_logProd_mlj <- function(l, q_variable_level_m, Eta_mj.) {
          active_probs <- q_variable_level_m
          inactive_probs <- 1 - q_variable_level_m
          probs <- c(active_probs[Eta_mj.], inactive_probs[1 - Eta_mj.])
          return(sum(log(probs)))
        }
        get_U_active_m <- function(m, Gamma, U) {
          U_active_index <- which(Gamma_m == 1)
          U[, U_active_index]
        }
        calculate_mvnorm_variance_mj <- function(m, j, U_active_m, Sigma2, Tau) {
          n_observations <- nrow(U_active_m)
          Sigma2[[m]][j] * U_active_m %*% diag(Tau[[m]][, j]) %*% t(U_active_m) + diag(n_observations)
        }
        get_log_dmvnorm_mj <- function(X_mj, mvnorm_variance_mj, n_observations) {
          mvtnorm::dmvnorm(x = X_mj, mean = rep(0, n_observations), sigma = mvnorm_variance_mj, log = TRUE)
        }
        get_P_lj <- function(logG_lj_0, logG_lj_1) {
          x_0 <- max(logG_lj_0, logG_lj_1)
          x <- logG_lj_1 - x_0
          y <- logG_lj_0 - x_0
          exp(x) / (exp(x) + exp(y))
        }
        eta_l_prime <- numeric(length(eta_old))
        for (j in 1:n_features[m]) {
          logProd_m <- q_variable_level[[m]] %>% calculate_logProd_m()
          U_active_m <- get_U_active_m(m, Gamma, U)
          mvnorm_variance_mj <- calculate_mvnorm_variance_mj(m, j, U_active_m, Sigma2, Tau)
          logdmvnorm_mj <- get_log_dmvnorm_mj(X_mj = data_list[[m]][, j],
                             mvnorm_variance_mj,
                             n_observations)
          Eta_mj <- Eta[[m]][, j]
          Eta_mj0 <- Eta_mj 
          Eta_mj0[l] <- 0
          Eta_mj1 <- Eta_mj
          Eta_mj1[l] <- 1
          logProd_lj_0 <- get_logProd_mlj(l, q_variable_level[[m]], Eta_mj0)
          logProd_lj_1 <- get_logProd_mlj(l, q_variable_level[[m]], Eta_mj1)
          logG_lj_0 <- logdmvnorm_mj + logProd_lj_0
          logG_lj_1 <- logdmvnorm_mj + logProd_lj_1
          P_lj <- get_P_lj(logG_lj_0, logG_lj_1)
          eta_l_prime[j] <- rbinom(1, 1, prob = P_lj)
        }
      }
      
      # Compute the acceptance probability
      
      # alpha <- min(1, target_density(x_star) / target_density(x[i - 1]))
      
      # Accept or reject the proposed values
      
      # if (runif(1) < alpha) {
      #   x[i] <- x_star
      # } else {
      #   x[i] <- x[i - 1]
      # }
      
    }
    
    # Step 3. Sample variance parameters sigma, tau, lambda
    
    active_components <- which(Gamma[[m]] == 1)
    active_features <- which(Eta[[m]] == 1)
    
    Sigma2_m_prime <- Sigma2[[m]]
    Tau_m_prime <- Tau[[m]]
    Lambda_m_prime <- Lambda[[m]]
    
    for (l in 1:r) {
    
      for (j in 1:n_features[m]) {
        
        if (j %in% active_features) {
          
          ## Sample from conditionals
          alpha <- prior_residual_variance[1] + n_observations/2
          sigma_j <- U[, active_components] %*% 
            diag(Tau[[m]][,j]) %*% 
            t(U[, active_components]) + 
            diag(n_observations)
          sigma_j_inverse <- sigma_j %>% solve # TODO: implement shortcut for getting sigma_j_inverse
          beta <- 1/2 * t(X_m[,j]) %*% sigma_j_inverse %*% X_m[,j] + prior_residual_variance[2]
          Sigma2_m_prime[j] <- 1/ rgamma(n = 1, alpha, beta)
          
          mu_prime <- sqrt(2 * Lambda[[m]][l,j] * Sigma2[[m]][j] / A[[m]][l,j]^2)
          lambda_prime <- 2 * Lambda[[m]][l,j]
          # TODO: Ensure inverse gaussian simulation appropriate
          set.seed(1)
          r_inverseGaussian <- function(mu, lambda) {
            v <- rnorm(1)
            y <- v^2
            x <- mu + (mu^2 * y)/(2*lambda) - (mu/(2*lambda)) * sqrt(4*mu*lambda*y + mu^2*y^2)
            test <- runif(1)
            if (test <= (mu)/(mu + x))
              return(x)
            else
              return((mu^2)/x)
          }
          Tau_m_prime[l,j] <- 1 / r_inverseGaussian(mu_prime, lambda_prime)
          
          # Lambda_m_prime[l,j] <- rgamma(n = 1, 
          #                               shape = ,
          #                               rate = )
          
        } else {
          
          ## Sample from pseudo-priors
          
        }
      
      }
      
    }
    
    # Step 4. Sample A
      
    active_components <- which(Gamma[[m]] == 1)
    active_features <- which(Eta[[m]] == 1)
    for (j in 1:n_features[m]) {
      sigma_a_j <- 1/Sigma2[[m]][j] * 
        (t(U[,active_components]) %*% U[,active_components] + diag(length(active_components))) %>%
        solve # Get the inverse
      mu_a_j <- sigma_a_j %*% 
        t(U[,active_components]) %*% 
        data_list[[m]][,j]
      a_prime <- MASS::mvrnorm(n = n_features[[m]],
                               mu = mu_a_j,
                               Sigma = sigma_a_j) %>% t
      
    }
    
    # Step 5. Sample U
    
    # Step 6. Sample group effect parameters r, b by Metropolis-Hastings
    
    }
    
  }
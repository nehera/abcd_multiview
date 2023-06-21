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
prior_residual_variance <- c(1, 1) # Hyperparameters of inverseGamma residual variance sigma_j^2(m)
prior_b0 <- c(2, 2) # Hyperparameters of a gamma distribution; prior distribution of the paramater b_l0 that controls the shrinkage of loadings when there is no grouping information
prior_b <- c(1, 1) # Hyperparameters of a gamma distribution; prior distribution of group effect coefficients blj, j = 1, ..., K_l.
prior_lambda_alpha <- 1 # TODO: Understand what this should be in practice
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

get_U_active_m <- function(Gamma_m, U) {
  U_active_index <- which(Gamma_m == 1)
  U[, U_active_index]
}

get_logProd_mlj <- function(q_variable_level_m, Eta_mj.) {
  active_probs <- q_variable_level_m
  inactive_probs <- 1 - q_variable_level_m
  probs <- c(active_probs[Eta_mj.], inactive_probs[1 - Eta_mj.])
  return(sum(log(probs)))
}

calculate_Sigma_j <- function(j, U_active_m, Sigma2_m, Tau2_m) {
  n <- nrow(U_active_m)
  if (ncol(U_active_m) == 0) {
    return(diag(n))
  } else {
    return(U_active_m %*% diag(Tau2_m[, j]) %*% t(U_active_m) + diag(n))
  }
}

calculate_mvnorm_variance_mj <- function(j, U_active_m, Sigma2_m, Tau2_m) {
  Sigma2_m[j] * calculate_Sigma_j(j, U_active_m, Sigma2_m, Tau2_m)
}

get_log_dmvnorm_mj <- function(X_mj, mvnorm_variance_mj) {
  n <- length(X_mj)
  mvtnorm::dmvnorm(x = X_mj, mean = rep(0, n), sigma = mvnorm_variance_mj, log = TRUE)
}

get_logG <- function(j, U_active_m, Sigma2_m, Tau2_m, 
                        X_mj, Eta_mj, q_variable_level_m) {
  mvnorm_variance_mj <- calculate_mvnorm_variance_mj(j, U_active_m, Sigma2_m, Tau2_m)
  logdmvnorm_mj <- get_log_dmvnorm_mj(X_mj, mvnorm_variance_mj)
  logProd_lj <- get_logProd_mlj(q_variable_level_m, Eta_mj)
  return(logdmvnorm_mj + logProd_lj)
}

sample_eta_prime_ml <- function(l, q_variable_level_m, Gamma_m, Eta_m,
                                X_m, U_active_m, Sigma2_m, Tau2_m) {
  
  n <- nrow(X_m)
  p_m <- ncol(X_m)
  
  get_P_lj <- function(logG_lj_0, logG_lj_1) {
    x_0 <- max(logG_lj_0, logG_lj_1)
    x <- logG_lj_1 - x_0
    y <- logG_lj_0 - x_0
    exp(x) / (exp(x) + exp(y))
  }
  
  eta_old <- Eta_m[l, ]
  eta_l_prime <- numeric(length(eta_old))
  
  for (j in 1:p_m) {
    
    Eta_mj <- Eta_m[, j]
    Eta_mj0 <- Eta_mj 
    Eta_mj0[l] <- 0
    Eta_mj1 <- Eta_mj
    Eta_mj1[l] <- 1
    
    logG_lj_1 <- get_logG(j, U_active_m, Sigma2[[m]], Tau2[[m]], 
                data_list[[m]][, j], Eta_mj1, 
                q_variable_level[[m]])
    logG_lj_0 <- get_logG(j, U_active_m, Sigma2[[m]], Tau2[[m]], 
                data_list[[m]][, j], Eta_mj0, 
                q_variable_level[[m]])
    
    P_lj <- get_P_lj(logG_lj_0, logG_lj_1)
    
    eta_l_prime[j] <- rbinom(1, 1, prob = P_lj)
    
  }
  
  return(eta_l_prime)
  
}

#### ---- MCMC initialization

# Step 0. Initialize parameters to estimate
# TODO: Initialize from priors (rather than fixed)

## Latent factor

U <- matrix(rnorm(n_observations*r), nrow = n_observations, ncol = r)

# Component activation

q_component_level <- rbeta(n_views, prior_component_selection[1], prior_component_selection[2]) # Probability for mth data type, a component is active

init_Gamma_m <- function(m, q_component_level, covariate_index) {
  if (length(covariate_index) > 0) {
    if (m %in% covariate_index) {
      rep(1, r) # Covariates activation forced to 1
    }
  } else {
    rbinom(n = r, size = 1, q_component_level[m])
  }
}

Gamma <- lapply(1:n_views, init_Gamma_m, 
                q_component_level, covariate_index)

init_Eta_m <- function(m, Gamma, prior_variable_selection, 
                       n_features, covariate_index, response_index) {
  if (length(covariate_index) > 0) {
    if (m %in% covariate_index) {
      matrix(1, nrow = r, ncol = n_features[m]) # Covariates activation forced to 1
    }
  } else if (m %in% response_index) {
    matrix(Gamma[[response_index]], ncol = 1) # Response variable activation forced to match component activation
  } else {
    
    init_omics_variable_activation_l <- function(gamma_l, prior_variable_selection, n_features_m) {
      if (gamma_l == 1) {
        rbinom(n = n_features[m], size = 1, prior_variable_selection)
      } else {
        rep(0, n_features[m])
      }
    }
    
    t(sapply(Gamma[[m]], init_omics_variable_activation_l, 
             prior_variable_selection, n_features[m]))
    
  }
}

Eta <- lapply(1:n_views, init_Eta_m,
              Gamma, prior_variable_selection,
              n_features, covariate_index, response_index)

init_Sigma2_m <- function(m, prior_residual_variance, n_features) {
  1/ rgamma(n = n_features[m], shape = prior_residual_variance[1], rate = prior_residual_variance[2])
}

Sigma2 <- lapply(1:n_views, init_Sigma2_m,
                 prior_residual_variance, n_features)

init_Tau2_m <- function(m, r, n_features, Tau2_prior_fixed = 1) {
  matrix(Tau2_prior_fixed, nrow = r, ncol = n_features[m])
} 

Tau2 <- lapply(1:n_views, init_Tau2_m, r, n_features)

# TODO: Remove loops in function
init_A_m <- function(m, r, n_features, Gamma, Eta, Sigma2, Tau2) {
  A_m <- matrix(nrow = r, ncol = n_features[m])
  for (l in 1:r) {
    a_l. <- rep(0, n_features[m])
    if (Gamma[[m]][l] == 1) {
      for (j in 1:n_features[m]) {
        if (Eta[[m]][l,j] == 1) {
          a_l.[j] <- rnorm(1, mean = 0, sd = sqrt(Sigma2[[m]][j] * Tau2[[m]][l,j]))
        }
      }
    }
    A_m[l, ] <- a_l.
  }
  return(A_m)
}

A <- lapply(1:n_views, init_A_m, r, n_features, 
            Gamma, Eta, Sigma2, Tau2)

## Grouping information parameters

init_vector_list <- function(initial_value, n_views, len_vector) {
  lapply(1:n_views, function(m, initial_value, len_vector) {
    rep(initial_value, len_vector[m])
  }, initial_value, len_vector)
}

b_0 <- init_vector_list(initial_value = 0.1, n_views, rep(r, n_views)) # b_0 fixed

q_r <- rbeta(1, prior_group_selection[1], prior_group_selection[2]) # TODO: Update in MCMC

K <- sapply(group_list, ncol) # number of groups in each view

init_R_m <- function(m, q_r, r, K) {
  matrix(rbinom(n = r * K[m], size = 1, prob = q_r),
         nrow = r, ncol = K[m])
}

R <- lapply(1:n_views, init_R_m, q_r, r, K)

init_B_m <- function(m, R, prior_b) {
  dims <- dim(R[[m]])
  r_vector <- as.vector(R[[m]])
  b_vector <- rgamma(sum(r_vector), prior_b[1], prior_b[2])
  r_vector[which(r_vector==1)] <- b_vector
  matrix(r_vector, nrow = dims[1], ncol = dims[2])
}

B <- lapply(1:n_views, init_B_m, R, prior_b)

init_Lambda2_m <- function(m, prior_lambda_alpha, b_0, group_list, B) {
  b_0_m <- b_0[[m]]; P_m <- group_list[[m]]; B_m <- B[[m]]
  n_l <- nrow(B_m); n_j <- nrow(P_m)
  get_beta_m_lj <- function(l, j, b_0_m, P_m, B_m) {
    b_0_m[l] + as.vector( t(P_m[j, ]) %*% B_m[l, ] )
  }
  beta_m <- matrix(nrow = n_l, ncol = n_j)
  for (l in 1:n_l) {
    for (j in 1:n_j) {
      beta_m[l,j] <- get_beta_m_lj(l, j, b_0_m, P_m, B_m)
    }
  }
  beta_m_vector <- as.vector(beta_m)
  Lambda2_m <- rgamma(length(beta_m_vector), prior_lambda_alpha, beta_m_vector) %>%
    matrix(nrow = nrow(beta_m), ncol = ncol(beta_m))
  return(Lambda2_m)
  
}

Lambda2 <- lapply(1:n_views, init_Lambda2_m,
                  prior_lambda_alpha, b_0, group_list, B)

t <- 1 # TODO: Remove post-development

#### ---- MCMC sample posteriors

for (t in 1:n_iterations) {
  
  # Step 0. Sample intercept for response variable
  
  intercept <- sample_intercept(Y = data_list[[response_index]],
                   A_outcome = A[[response_index]], U,
                   Sigma2_outcome = Sigma2[[response_index]])
  
  # Adjust Y for subsequent sampling
  Y_prime <- data_list[[response_index]] - intercept
  X <- data_list
  X[[response_index]] <- Y_prime
  
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

    l <- 1 # TODO: Remove post-dev
  
  # Step 2. Sample component, feature activation parameters (gamma, eta respectively) by Metropolis-Hastings
    
    U_active_m <- get_U_active_m(Gamma[[m]], U)
    Gamma_prime <- Gamma[[m]]
    Eta_prime <- Eta[[m]]
    for (l in 1:r) {

      # Propose new values
      if (Gamma[[m]][l] == 1) {
        # Propose deactivation
        Gamma_prime[l] <- 0
        Eta_prime[l, ] <- rep(0, ncol(Eta_prime))
      } else if (Gamma[[m]][l] == 0) {
        # Propose activation
        Gamma_prime[l] <- 1
        Eta_prime[l, ] <- sample_eta_prime_ml(l, q_variable_level[[m]], Gamma[[m]], Eta[[m]],
                                            X[[m]], U_active_m, Sigma2[[m]], Tau2[[m]])
        
      }
      
      # Compute the acceptance probability
      
      log_target_density <- function(Gamma_star, Eta_star) {
        r <- nrow(Eta_star)
        p_m <- ncol(Eta_star)
        logG <- numeric(p_m)
        for (j in 1:p_m) {
          logG[j] <- get_logG(j, U_active_m, Sigma2[[m]], Tau2[[m]], 
                              X[[m]][, j], Eta[[m]][, j], 
                              q_variable_level[[m]])
        }
        return(sum(logG))
      }
      
      alpha <- min(1, log_target_density(Gamma_prime, Eta_prime) - 
        log_target_density(Gamma[[m]], Eta[[m]]))
      
      # Accept or reject the proposed values
      
      if (runif(1) < alpha) {
        Gamma[[m]] <- Gamma_prime
        Eta[[m]] <- Eta_prime
      } 
      
    }
    
    # Step 3. Sample variance parameters sigma, Tau2, lambda
    
    sample_Sigma2_m <- function(m, prior_residual_variance, X, Tau2) {
      X_m <- X[[m]]
      n <- nrow(X_m); n_j <- ncol(X_m)
      alpha <- prior_residual_variance[1] + n/ 2
      betas <- numeric(n_j)
      for (j in 1:n_j) {
        Sigma_j_inverse <- calculate_Sigma_j(j, U_active_m, Sigma2_m, Tau2_m) %>% solve() # TODO: Consider Woodbury Identity for Speeding inversion
        betas[j] <- 1/ 2 * t(X_m[, j]) %*% Sigma_j_inverse %*% X_m[, j] + prior_residual_variance[2]
      }
      1/ rgamma(n_j, alpha, betas)
    }
    
    Sigma2_prime <- lapply(1:n_views, sample_Sigma2_m, prior_residual_variance, X, Tau2)
    

    sample_Tau2_Lambda2_m <- function(m, prior_lambda_alpha, Eta, group_list, 
                                      Tau2, Lambda2, Sigma2, A, b_0, B) {
      
      dims <- dim(Eta[[m]])
      active_index <- which(Eta[[m]] == 1)
      P_m <- group_list[[m]]

      Tau2_m_prime <- Tau2[[m]]
      Lambda2_m_prime <- Lambda2[[m]]
      
      for (l in 1:dims[1]) {
        for (j in 1:dims[2]) {
          if ( (l*j) %in% active_index ) {
            mu_prime <- sqrt( 2 * Lambda2[[m]][l,j] * Sigma2[[m]][j] / A[[m]][l,j]^2 )
            lambda2_prime <- 2 * Lambda2[[m]][l,j]
            Tau2_m_prime[l,j] <- 1/ statmod::rinvgauss(1, mean = mu_prime, shape = lambda2_prime)
            Lambda2_m_prime[l,j] <- rgamma(1, prior_lambda_alpha + 1, 
                                 b_0[[m]][l] + t(P_m[j, ]) %*% B[[m]][l, ] + Tau2[[m]][l,j])
            
          } else {
            # Pseudo-priors
            Tau2_m_prime[l,j] <- rexp(1, rate = Lambda2[[m]][l,j])
            Lambda2_m_prime[l,j] <- rgamma(prior_lambda_alpha, b_0[[m]][l])
          }
        }
      }
      
      return(list(Tau2_m_prime, Lambda2_m_prime))
    }

    Tau2_Lambda2 <- lapply(1:n_views, sample_Tau2_Lambda2_m, prior_lambda_alpha,
           Eta, group_list, Tau2, Lambda2, Sigma2, A, b_0, B)
    # Split Tau2_Lambda2 list
    Tau2 <- lapply(Tau2_Lambda2, `[[`, 1)
    Lambda2 <- lapply(Tau2_Lambda2, `[[`, 2)
    
    # Step 4. Sample A
    
    sample_A_m <- function(m, A, U, Sigma2, Gamma, Eta, X) {
      U_active_m <- get_U_active_m(Gamma[[m]], U)
      n_components_active <- sum(Gamma[[m]])
      
      if (n_components_active==0) {
        return(matrix(0, nrow = nrow(A[[m]]), ncol = ncol(A[[m]])))
      }
      
      Sigma_a_inverse_part <- t(U_active_m) %*% U_active_m + diag(n_components_active)
      
      Sigma_a_list <- lapply(Sigma2[[m]], function(s2_j, S_inverse = Sigma_a_inverse_part) { 
        solve(s2_j * S_inverse) # TODO: Consider Woodbury Identity for Matrix Inversion
        })
      p_m <- ncol(X[[m]])
      Mu_a_list <- lapply(1:p_m, function(j, Sigma_a = Sigma_a_list, U_active = U_active_m, X_m = X[[m]]) {
        Sigma_a[[j]] %*% t(U_active) %*% X_m[, j]
      })
      

      A_prime_active <- sapply(1:p_m, function(j, Mu_a = Mu_a_list, Sigma_a = Sigma_a_list) {
        rmvnorm(1, mean = as.vector(Mu_a[[j]]), sigma = Sigma_a[[j]])
      } )
      
      A_prime <- A[[m]]
      active_component_index <- which(Gamma[[m]]==1)
      A_prime[active_component_index, ] <- A_prime_active
      # replace with 0's if Eta[l,j] is 0
      A_prime[which(Eta[[m]] == 0)] <- 0
      
      return(A_prime)
    }
    
    A <- lapply(1:n_views, sample_A_m, A, U, Sigma2, Gamma, Eta, X) 
    
    # Step 5. Sample U
    
    sample_U <- function(A, X, Sigma2, r, n_observations) {
      a_matrix <- rlist::list.cbind(A)
      x_matrix <- rlist::list.cbind(X)
      sigma2_vector <- unlist(Sigma2)
      sigma_u_i <- a_matrix %*% diag(sigma2_vector) %*% t(a_matrix) + diag(r)
      
      mu_u_list <- lapply(1:n_observations, function(i, Sigma_u_i = sigma_u_i, 
                                                     A_matrix = a_matrix, 
                                                     Sigma2_vector = sigma2_vector, 
                                                     X_matrix = x_matrix) {
        Sigma_u_i %*% A_matrix %*% diag(Sigma2_vector) %*% X_matrix[i, ]
      })
      
      U_prime <- sapply(1:n_observations, function(i, Mu_u_list = mu_u_list, Sigma_u_i = sigma_u_i) {
        rmvnorm(1, mean = as.vector(Mu_u_list[[i]]), sigma = sigma_u_i)
      }) %>% t()
      return(U_prime) 
    }

    U_prime <- sample_U(A, X, Sigma2, r, n_observations)
    
    # Step 6. Sample group effect parameters r, b by Metropolis-Hastings
    
    # TODO: Implement Step B5
    k=1
    R_m <- R[[m]]
    B_m <- B[[m]]
    propose_br_m <- function(m, ) {
      
      for (l in 1:r) {
        
        R_m_l <- R_m[l, ]
        K <- length(R_m_l)
        
        for (k in K) {
          
          propose_RB_m_lk <- function(k, R_m_l, prior_b) {
            if (R_m_l[k] == 1) {
              R_m_lk_prime <- 0
              B_m_lk_prime <- 0
              return(list(r=R_m_lk_prime, b=B_m_lk_prime))
            } else {
              R_m_lk_prime <- 1
              B_m_lk_prime <- rgamma(1, prior_b[1], prior_b[2])
              return(list(r=R_m_lk_prime, b=B_m_lk_prime))
            }
          }
          
          # Propose new r_lk, b_lk
          
          rb_prime <- propose_RB_m_lk(k, R_m_l, priorb)
          rb <- list(r=R_m[l,k], b=B_m[l,k])
          
          # TODO: Compute the acceptance probability
        
          rb_log_target_density <- function(rb) {
            
            
            
          }

          alpha <- min(1, rb_log_target_density(rb_prime) - rb_log_target_density(rb))
          
          # Accept or reject the proposed values
          
          if (runif(1) < alpha) {
            R_m[l,k] <- rb_prime$r
            B_m[l,k] <- rb_prime$b
          } 
          
        }
        
      }
      

        
        
    }
    
    }
    
  }
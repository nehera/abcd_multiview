## Get CLI arguments
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}
seed <- as.numeric(args[1])
set.seed(seed)

## Load Packages
library(tidyverse)
library(Rcpp)
library(mvtnorm)
library(parallel) 
RNGkind("L'Ecuyer-CMRG") # for assigning separate, reproducible RNG streams to parallel computations

## Load Scaled Data
data_list <- readRDS("data/2023-06-30_simulated_data.rds")
indic_var <- c(0, 0, 1) 
method <- "BIP" # method without grouping information to start
group_list <- NULL # default when grouping information not included

## Set MCMC Parameters
n_sample <- 5000
n_burnin <- 1000
n_max_models <- 50

## Set Hyper-Parameters
r <- 4 # number of components
prior_component_selection <- c(1, 1) # Hyperparameters of a beta(a,b) distribution; prior distribution of the probability of selecting components
prior_response_selection <- c(1, 1) # Hyperparameters of a beta(a,b) distribution; prior distribution of the probability of selecting the response on a component.
prior_residual_variance <- c(1, 1) # Hyperparameters of inverseGamma residual variance sigma_j^2(m)
# TODO: Assess whether prior_b0 helpful
prior_b0 <- c(2, 2) # Hyperparameters of a gamma distribution; prior distribution of the paramater b_l0 that controls the shrinkage of loadings when there is no grouping information
prior_b <- c(1, 1) # Hyperparameters of a gamma distribution; prior distribution of group effect coefficients blj, j = 1, ..., K_l.
prior_lambda_alpha <- 1 # TODO: Understand what this should be in practice
prior_group_selection <- c(1, 1) # Hyperparameters of a beta(a,b) distribution; prior distribution of the probability of selecting groups.
prior_variable_selection <- 0.05 # Prior probability for variable selection

#### ----- Begin BIP2 ---- ####

# BIP2 <- function(data_list = data_list, # TODO: This default might cause probs
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

#### ---- Initialize MCMC meta-parameters

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
n_observations <- sapply(data_list, nrow) %>% unique()
n_iterations <- n_sample + n_burnin

# Ensure parameters passed are appropriate
if (length(n_observations) > 1) {
  stop("Each view must contain the same n_observations i.e. the same number of rows.")
}
# if (n_sample <= n_burnin) {
#   stop("Argument n_burnin must be smaller than argument n_sample, the number of MCMC iterations.")
# }
# if (n_sample <= 100) {
#   stop("Please specify a larger number of MCMC iterations.")
# }

# TODO: Understand how r established when passed as NULL.
if (is.null(r)) {
  mysvd <- lapply(omics_index, function(i) svd(data_list[[i]]))
  mysumsvd <- lapply(1:length(omics_index), function(i) cumsum(mysvd[[i]]$d) / max(cumsum(mysvd[[i]]$d)) * 100)
  KMax <- max(unlist(lapply(1:length(omics_index), function(i) min(which(mysumsvd[[i]]>=80, arr.ind = TRUE)))))
  r <- min(KMax+1, 10)
}

#### ---- Initialize MCMC data structures

intercept <- numeric(n_iterations)

U <- array(data = rnorm(n_observations*r), dim = c(n_observations, r, n_iterations))

# TODO: Unfix q_component_level
q_component_level <- matrix(0.5, nrow = n_views, ncol = n_iterations) # Probability for mth data type, a component is active

if (length(covariate_index) > 0) {
  q_component_level[, covariate_index] <- 1
}

init_Gamma_m <- function(m, q_component_level, covariate_index, n_iterations) {
  if (length(covariate_index) > 0) {
    if (m %in% covariate_index) {
      matrix(1, nrow = r, ncol = n_iterations) # Covariates activation forced to 1
    }
  } else {
    n_active <- numeric(1)
    while (n_active == 0) {
      Gamma_m <- matrix(rbinom(n = r, size = 1, q_component_level[m, 1]), 
                        nrow = r, ncol = n_iterations)
      n_active <- sum(Gamma_m[, 1])
    }
    return(Gamma_m)
  }
}

Gamma <- lapply(1:n_views, init_Gamma_m, 
                q_component_level, covariate_index, n_iterations)

init_Eta_m <- function(m, Gamma, prior_variable_selection, 
                       n_features, covariate_index, response_index, n_iterations) {
  if (length(covariate_index) > 0) {
    if (m %in% covariate_index) {
      # Covariate activation forced to 1
      array(1, dim = c(r, n_features[m], n_iterations)) 
    }
  } else if (m %in% response_index) {
    # Response variable activation forced to match component activation
    Gamma[[response_index]][, 1] %>%
      matrix(nrow = r, ncol = n_features[m]) %>% 
      array(dim = c(r, n_features[m], n_iterations)) 
  } else {
    # Omics variable activation initiated from prior distribution
    init_omics_variable_activation_l <- function(gamma_l, prior_variable_selection, n_features_m) {
      if (gamma_l == 1) {
        rbinom(n = n_features[m], size = 1, prior_variable_selection)
      } else {
        rep(0, n_features[m])
      }
    }
    sapply(Gamma[[m]][, 1], init_omics_variable_activation_l, prior_variable_selection, n_features[m]) %>%
      t() %>% array(dim = c(r, n_features[m], n_iterations))
  }
}

Eta <- lapply(1:n_views, init_Eta_m,
              Gamma, prior_variable_selection,
              n_features, covariate_index, response_index, n_iterations)

init_Sigma2_m <- function(m, prior_residual_variance, n_features, n_iterations) {
  1/ rgamma(n = n_features[m], shape = prior_residual_variance[1], rate = prior_residual_variance[2]) %>%
    matrix(nrow = n_features[m], ncol = n_iterations)
}

Sigma2 <- lapply(1:n_views, init_Sigma2_m,
                 prior_residual_variance, n_features, n_iterations)

init_Tau2_m <- function(m, r, n_features, Tau2_prior_fixed = 1) {
  array(data = Tau2_prior_fixed, dim = c(r, n_features[m], n_iterations))
} 

Tau2 <- lapply(1:n_views, init_Tau2_m, r, n_features)

init_A_m <- function(m, r, n_features, Gamma, Eta, Sigma2, Tau2, n_iterations) {
  A_m_1 <- matrix(0, nrow = r, n_features[m])
  active_index <- which(Eta[[m]][,,1] == 1)
  sigma2_temp <- Sigma2[[m]][,1] %>% 
    matrix(nrow = r, ncol = n_features[m], byrow = TRUE)
  active_sigmas <- sigma2_temp[active_index]
  active_taus <- Tau2[[m]][,,1][active_index]
  active_sds <- sqrt(active_sigmas * active_taus)
  # TODO: Remove loops in function
  for (i in 1:length(active_index)) {
    A_m_1[active_index[i]] <- rnorm(n = 1, mean = 0, sd = active_sds[i])
  }
  array(A_m_1, dim = c(r, n_features[m], n_iterations))
}

A <- lapply(1:n_views, init_A_m, r, n_features, 
            Gamma, Eta, Sigma2, Tau2, n_iterations)

b_0 <- rep(list(rep(0.1, r)), n_views) # TODO: Assess initialization from prior, but fix for now

q_r <- rbeta(1, prior_group_selection[1], prior_group_selection[2]) # TODO: Update in MCMC

K <- sapply(group_list, ncol) # number of groups in each view

init_R_m <- function(m, q_r, r, K, n_iterations) {
  array(rbinom(n = r * K[m], size = 1, prob = q_r), 
        dim = c(r, K[m], n_iterations))
}

R <- lapply(1:n_views, init_R_m, q_r, r, K, n_iterations)

init_B_m <- function(m, R, prior_b, r, K, n_iterations) {
  group_activation <- as.vector(R[[m]][,,1])
  coefficients <- rgamma(sum(group_activation), prior_b[1], prior_b[2])
  active_index <- which(group_activation==1)
  group_activation[active_index] <- coefficients
  array(group_activation, dim = c(r, K[m], n_iterations))
}

B <- lapply(1:n_views, init_B_m, R, prior_b, r, K, n_iterations)

init_Lambda2_m <- function(m, prior_lambda_alpha, group_list, b_0, B, 
                           r, n_features, n_iterations) {
  P_m <- group_list[[m]]
  b_0_m <- b_0[[m]]
  B_m_1 <- B[[m]][,,1]
  Lambda2_m1 <- matrix(0, nrow = r, ncol = n_features[m])
  for (l in 1:r) {
    for (j in 1:n_features[m]) {
      beta <- b_0_m[l] + as.vector( t(P_m[j, ]) %*% B_m_1[l] )
      Lambda2_m1[l,j] <- rgamma(1, prior_lambda_alpha, beta)
    }
  }
  Lambda2_m <- Lambda2_m1 %>%
    array(dim = c(r, n_features[m], n_iterations))
}

Lambda2 <- lapply(1:n_views, init_Lambda2_m, prior_lambda_alpha, group_list, b_0, B, 
                  r, n_features, n_iterations)

#### ---- Initialize MCMC sampling functions

reshape_views <- function(Gamma, Eta, X, Sigma2, Tau2) {
  n_views <- length(Gamma)
  view_list <- list()
  for (m in 1:n_views) {
    view_list[[m]] <- list(Gamma_m = Gamma[[m]], Eta_m = Eta[[m]], X_m = X[[m]],
                           Sigma2_m = Sigma2[[m]], Tau2_m = Tau2[[m]], m=m)
  }
  return(view_list)
}

extract_view_iter <- function(view_list, iter) {
  # view_m <- view_list[[1]]
  extract_view_m_iter <- function(view_m, iter) {
    list2env(view_m, envir = environment())
    Gamma_m_iter <- Gamma_m[, iter]
    Eta_m_iter <- Eta_m[,, iter]
    Sigma2_m_iter <- Sigma2_m[, iter]
    Tau2_m_iter <- Tau2_m[,, iter]
    return(list(Gamma_m = Gamma_m_iter, Eta_m = Eta_m_iter,
                X_m = X_m, Sigma2_m = Sigma2_m_iter, Tau2_m = Tau2_m_iter, m=m))
  }
  return(lapply(view_list, extract_view_m_iter, iter))
}

store_view_iter <- function(view_list, view_list_iter, iter) {
  for (m in 1:n_views) {
    view_list[[m]]$Gamma_m[, iter] <- view_list_iter[[m]]$Gamma_m # Store only gamma for now
  }    
  # TODO: Store other vars too
  return(view_list)
}

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


get_Sigma_j <- function(j, Gamma_m, U_iter, Sigma2_m, Tau2_m) {
  active_index <- which(Gamma_m == 1)
  U_active <- U_iter[, active_index] %>%
    matrix(ncol = length(active_index))
  n <- nrow(U_active)
  if (ncol(U_active) == 0) {
    return(diag(n))
  } else {
    return(U_active %*% diag(Tau2_m[active_index, j]) %*% t(U_active) + diag(n))
  }
}

get_logG <- function(j, Gamma_m, Eta_mj, U_iter, Sigma2_m, Tau2_m, 
                     X_m, prior_variable_selection) {
  n <- nrow(X_m)
  mvnorm_variance_mj <- Sigma2_m[j] * get_Sigma_j(j, Gamma_m, U_iter, Sigma2_m, Tau2_m)
  logdmvnorm_mj <- mvtnorm::dmvnorm(x = X_m[,j], mean = rep(0, n), sigma = mvnorm_variance_mj, log = TRUE)
  n_active <- sum(unlist(Eta_mj))
  logprod_lj <- n_active*log(prior_variable_selection) + 
    (r-n_active)*log(1-prior_variable_selection)
  return(logdmvnorm_mj + logprod_lj)
}

get_P_lj <- function(logG_lj_0, logG_lj_1) {
  x_0 <- max(logG_lj_0, logG_lj_1)
  x <- logG_lj_1 - x_0
  y <- logG_lj_0 - x_0
  exp(x) / (exp(x) + exp(y))
}

propose_eta_prime_m_l <- function(l, prior_variable_selection, 
                                  Gamma_m, Eta_m, X_m, U_iter,
                                  Sigma2_m, Tau2_m) {
  n <- nrow(X_m)
  p_m <- ncol(X_m)
  eta_m_l <- Eta_m[l, ]
  eta_m_l_prime <- numeric(length(eta_m_l))
  
  for (j in 1:p_m) {
    
    Eta_mj <- Eta_m[, j]
    Eta_mj0 <- Eta_mj 
    Eta_mj0[l] <- 0
    Eta_mj1 <- Eta_mj
    Eta_mj1[l] <- 1
    
    logG_lj_1 <- get_logG(j, Gamma_m, Eta_mj1, U_iter, Sigma2_m, Tau2_m, 
                          X_m, prior_variable_selection)
    
    logG_lj_0 <- get_logG(j, Gamma_m, Eta_mj0, U_iter, Sigma2_m, Tau2_m, 
                          X_m, prior_variable_selection)
    
    P_lj <- get_P_lj(logG_lj_0, logG_lj_1)
    
    eta_m_l_prime[j] <- rbinom(1, 1, prob = P_lj)
    
  }
  return(eta_m_l_prime)
}


propose_GammaEta_prime_m_l <- function(l, Gamma_m, Eta_m,
                                       X_m, U_iter, 
                                       Sigma2_m, Tau2_m,
                                       prior_variable_selection,
                                       m, covariate_index) {
  # Propose new values
  if (Gamma_m[l] == 1) {
    
    # If component active, propose deactivation
    return(list(gamma_prime_m_l = 0, eta_prime_m_l = rep(0, n_features[m])))
    
  } else if (m %in% covariate_index) {
    
    # If component inactive, propose activation
    return(list(gamma_prime_m_l = 1, eta_prime_m_l = rep(1, n_features[m])))
    
  } else {
    eta_prime_m_l <- propose_eta_prime_m_l(l, prior_variable_selection, 
                                           Gamma_m, Eta_m,
                                           X_m, U_iter, 
                                           Sigma2_m, Tau2_m)
    return(list(gamma_prime_m_l = 1, eta_prime_m_l = eta_prime_m_l))
  }
}

log_target_density <- function(Gamma_prime, Eta_prime, U_iter, Sigma2_m, Tau2_m, 
                               X_m, prior_variable_selection) {
  r <- nrow(Eta_prime)
  p_m <- ncol(Eta_prime)
  logG <- numeric(p_m)
  for (j in 1:p_m) { # TODO: Remove for loop
    logG[j] <- get_logG(j, Gamma_prime, Eta_prime[,j], U_iter, 
                        Sigma2_m, Tau2_m, X_m, prior_variable_selection)
  }
  return(sum(logG))
}

sample_GammaEta <- function(view_m_iter, U_iter, prior_variable_selection, covariate_index) {
  list2env(view_m_iter, envir = environment())
  # get meta-parameters
  r <- length(Gamma_m)
  #### YOU ARE HERE
  GammaEta_prime_list <- lapply(1:r, propose_GammaEta_prime_m_l,
                                Gamma_m, Eta_m, X_m, U_iter,
                                Sigma2_m, Tau2_m,
                                prior_variable_selection,
                                m, covariate_index)
  
  Gamma_prime <- sapply(GammaEta_prime_list, `[[`, 1)
  Eta_prime <- sapply(GammaEta_prime_list, `[[`, 2) %>% t()
  
  # Compute the acceptance probability
  alpha <- min(1, log_target_density(Gamma_prime, Eta_prime, U_iter, Sigma2_m, Tau2_m, 
                                     X_m, prior_variable_selection) - 
                 log_target_density(Gamma_m, Eta_m, U_iter, Sigma2_m, Tau2_m, 
                                    X_m, prior_variable_selection))
  
  # Accept or reject the proposed values
  if (runif(1) < alpha) {
    view_m_iter$Gamma_m <- Gamma_prime
    view_m_iter$Eta_m <- Eta_prime
    return(view_m_iter)
  } else {
    return(view_m_iter)
  }
  
}

#### ---- MCMC Posterior Sampling

iter <- 2 # TODO: Remove post-development and uncomment for iter in 2:n_iterations

# for (iter in 2:n_iterations) {
  
  # Step 0. Sample intercept for response variable
  
  intercept <- sample_intercept(Y = data_list[[response_index]],
                   A_outcome = A[[response_index]][,, iter-1], U[,, iter-1],
                   Sigma2_outcome = Sigma2[[response_index]][, iter-1]) %>%
    rep(n_iterations) # TODO: Sample intercept at each iteration
  
  # Adjust Y for subsequent sampling
  X <- data_list
  X[[response_index]] <- data_list[[response_index]] - intercept[iter]
  
  ## -- Transform data structures for sampling steps 1-5
  
  view_list <- reshape_views(Gamma, Eta, X, Sigma2, Tau2)
  
  # Step 0 continued. Sample component and variable selection probabilities
  # TODO: Assess whether or not to update component and variable selection probs
  
  # sample_beta_posterior <- function(n_samples = 1, prior_alpha, prior_beta, n, summation) {
  #   posterior_alpha <- prior_alpha + summation
  #   posterior_beta <- prior_beta + n - summation
  #   return(rbeta(n_samples, posterior_alpha, posterior_beta))
  # }
  
  # Gamma_sums <- sapply(Gamma, sum)
  # for (m in c(omics_index, response_index)) {
  #   if (m %in% omics_index) {
  #     q_component_level[m] <- sample_beta_posterior(prior_alpha = prior_component_selection[1],
  #                                                   prior_beta = prior_component_selection[2],
  #                                                   n = r, summation = Gamma_sums[m])
  #     
  #     # Probability of variable selection fixed for omics data by prior_variable_selection
  #   } else if (m %in% response_index) {
  #     q_component_level[m] <- sample_beta_posterior(prior_alpha = prior_response_selection[1],
  #                                                   prior_beta = prior_response_selection[2],
  #                                                   n = r, summation = Gamma_sums[m])
  #     q_variable_level[[m]] <- rep(q_component_level[m], r)
  #   } 
  # }
    
    for (iter in 2:n_iterations) {
      
      # Step 0. Extract last iteration's data
      
      view_list_iter <- extract_view_iter(view_list, iter-1)
      
      # Step 2. Sample component, feature activation parameters (gamma, eta respectively) by Metropolis-Hastings
      
      view_list_iter[-response_index] <- mclapply(view_list_iter[-response_index], sample_GammaEta, U_iter = U[,, iter],
                                                prior_variable_selection, covariate_index, mc.cores = 2)
      
      Gamma_response_iter <- view_list_iter[[response_index]]$Gamma_m
      if (sum(Gamma_response_iter)==r) {
        n_active <- numeric(1)
        while (n_active == 0) {
          Gamma_response_prime <- rbinom(n = r, size = 1, q_component_level[response_index, iter])
          n_active <- sum(Gamma_prime)
        }
      } else {
        Gamma_response_prime <- ifelse(Gamma_response_iter==1, 0, 1)
      }
      
      view_list_iter[[response_index]]$Gamma_m <- Gamma_response_prime
      view_list_iter[[response_index]]$Eta_m <- Gamma_response_prime
      
      # Step Final. Store iteration's data
      # TODO: Storing only new Gammas for now. Want to store more.
      view_list <- store_view_iter(view_list, view_list_iter, iter)
      
      if (iter %% 500 == 0) {
        # Save progress
        fileConn <- file(paste0("iter_", seed, ".txt"))
        writeLines(paste("Seed", seed, "is on iter", iter), fileConn)
        close(fileConn)
        # Save output
        output_name <- paste0("data/", Sys.Date(), "_view_list_seed_", seed, ".rds")
        saveRDS(view_list, file = output_name)
      }
      
    }
  
  #   # Step 3. Sample variance parameters sigma, Tau2, lambda
  #   
  #   sample_Sigma2_m <- function(m, prior_residual_variance, X, Tau2) {
  #     X_m <- X[[m]]
  #     n <- nrow(X_m); n_j <- ncol(X_m)
  #     alpha <- prior_residual_variance[1] + n/ 2
  #     betas <- numeric(n_j)
  #     for (j in 1:n_j) {
  #       Sigma_j_inverse <- get_Sigma_j(j, U_active_m, Sigma2_m, Tau2_m) %>% solve() # TODO: Consider Woodbury Identity for Speeding inversion
  #       betas[j] <- 1/ 2 * t(X_m[, j]) %*% Sigma_j_inverse %*% X_m[, j] + prior_residual_variance[2]
  #     }
  #     1/ rgamma(n_j, alpha, betas)
  #   }
  #   
  #   Sigma2_prime <- lapply(1:n_views, sample_Sigma2_m, prior_residual_variance, X, Tau2)
  #   
  # 
  #   sample_Tau2_Lambda2_m <- function(m, prior_lambda_alpha, Eta, group_list, 
  #                                     Tau2, Lambda2, Sigma2, A, b_0, B) {
  #     
  #     dims <- dim(Eta[[m]])
  #     active_index <- which(Eta[[m]] == 1)
  #     P_m <- group_list[[m]]
  # 
  #     Tau2_m_prime <- Tau2[[m]]
  #     Lambda2_m_prime <- Lambda2[[m]]
  #     
  #     for (l in 1:dims[1]) {
  #       for (j in 1:dims[2]) {
  #         if ( (l*j) %in% active_index ) {
  #           mu_prime <- sqrt( 2 * Lambda2[[m]][l,j] * Sigma2[[m]][j] / A[[m]][l,j]^2 )
  #           lambda2_prime <- 2 * Lambda2[[m]][l,j]
  #           Tau2_m_prime[l,j] <- 1/ statmod::rinvgauss(1, mean = mu_prime, shape = lambda2_prime)
  #           Lambda2_m_prime[l,j] <- rgamma(1, prior_lambda_alpha + 1, 
  #                                b_0[[m]][l] + t(P_m[j, ]) %*% B[[m]][l, ] + Tau2[[m]][l,j])
  #           
  #         } else {
  #           # Pseudo-priors
  #           Tau2_m_prime[l,j] <- rexp(1, rate = Lambda2[[m]][l,j])
  #           Lambda2_m_prime[l,j] <- rgamma(prior_lambda_alpha, b_0[[m]][l])
  #         }
  #       }
  #     }
  #     
  #     return(list(Tau2_m_prime, Lambda2_m_prime))
  #   }
  # 
  #   Tau2_Lambda2 <- lapply(1:n_views, sample_Tau2_Lambda2_m, prior_lambda_alpha,
  #          Eta, group_list, Tau2, Lambda2, Sigma2, A, b_0, B)
  #   # Split Tau2_Lambda2 list
  #   Tau2 <- lapply(Tau2_Lambda2, `[[`, 1)
  #   Lambda2 <- lapply(Tau2_Lambda2, `[[`, 2)
  #   
  #   # Step 4. Sample A
  #   
  #   sample_A_m <- function(m, A, U, Sigma2, Gamma, Eta, X) {
  #     U_active_m <- get_U_active_m(Gamma[[m]], U)
  #     n_components_active <- sum(Gamma[[m]])
  #     
  #     if (n_components_active==0) {
  #       return(matrix(0, nrow = nrow(A[[m]]), ncol = ncol(A[[m]])))
  #     }
  #     
  #     Sigma_a_inverse_part <- t(U_active_m) %*% U_active_m + diag(n_components_active)
  #     
  #     Sigma_a_list <- lapply(Sigma2[[m]], function(s2_j, S_inverse = Sigma_a_inverse_part) { 
  #       solve(s2_j * S_inverse) # TODO: Consider Woodbury Identity for Matrix Inversion
  #       })
  #     p_m <- ncol(X[[m]])
  #     Mu_a_list <- lapply(1:p_m, function(j, Sigma_a = Sigma_a_list, U_active = U_active_m, X_m = X[[m]]) {
  #       Sigma_a[[j]] %*% t(U_active) %*% X_m[, j]
  #     })
  #     
  # 
  #     A_prime_active <- sapply(1:p_m, function(j, Mu_a = Mu_a_list, Sigma_a = Sigma_a_list) {
  #       rmvnorm(1, mean = as.vector(Mu_a[[j]]), sigma = Sigma_a[[j]])
  #     } )
  #     
  #     A_prime <- A[[m]]
  #     active_component_index <- which(Gamma[[m]]==1)
  #     A_prime[active_component_index, ] <- A_prime_active
  #     # replace with 0's if Eta[l,j] is 0
  #     A_prime[which(Eta[[m]] == 0)] <- 0
  #     
  #     return(A_prime)
  #   }
  #   
  #   A <- lapply(1:n_views, sample_A_m, A, U, Sigma2, Gamma, Eta, X) 
  #   
  #   # Step 5. Sample U
  #   
  #   sample_U <- function(A, X, Sigma2, r, n_observations) {
  #     a_matrix <- rlist::list.cbind(A)
  #     x_matrix <- rlist::list.cbind(X)
  #     sigma2_vector <- unlist(Sigma2)
  #     sigma_u_i <- a_matrix %*% diag(sigma2_vector) %*% t(a_matrix) + diag(r)
  #     
  #     mu_u_list <- lapply(1:n_observations, function(i, Sigma_u_i = sigma_u_i, 
  #                                                    A_matrix = a_matrix, 
  #                                                    Sigma2_vector = sigma2_vector, 
  #                                                    X_matrix = x_matrix) {
  #       Sigma_u_i %*% A_matrix %*% diag(Sigma2_vector) %*% X_matrix[i, ]
  #     })
  #     
  #     U_prime <- sapply(1:n_observations, function(i, Mu_u_list = mu_u_list, Sigma_u_i = sigma_u_i) {
  #       rmvnorm(1, mean = as.vector(Mu_u_list[[i]]), sigma = sigma_u_i)
  #     }) %>% t()
  #     return(U_prime) 
  #   }
  # 
  #   U_prime <- sample_U(A, X, Sigma2, r, n_observations)
  #   
  #   # Step 6. Sample group effect parameters r, b by Metropolis-Hastings
  #   
  #   # TODO: Implement Step B5
  #   k=1
  #   R_m <- R[[m]]
  #   B_m <- B[[m]]
  #   propose_br_m <- function(m, ) {
  #     
  #     for (l in 1:r) {
  #       
  #       R_m_l <- R_m[l, ]
  #       K <- length(R_m_l)
  #       
  #       for (k in K) {
  #         
  #         propose_RB_m_lk <- function(k, R_m_l, prior_b) {
  #           if (R_m_l[k] == 1) {
  #             R_m_lk_prime <- 0
  #             B_m_lk_prime <- 0
  #             return(list(r=R_m_lk_prime, b=B_m_lk_prime))
  #           } else {
  #             R_m_lk_prime <- 1
  #             B_m_lk_prime <- rgamma(1, prior_b[1], prior_b[2])
  #             return(list(r=R_m_lk_prime, b=B_m_lk_prime))
  #           }
  #         }
  #         
  #         # Propose new r_lk, b_lk
  #         
  #         rb_prime <- propose_RB_m_lk(k, R_m_l, priorb)
  #         rb <- list(r=R_m[l,k], b=B_m[l,k])
  #         
  #         # TODO: Compute the acceptance probability
  #       
  #         rb_log_target_density <- function(rb) {
  #           eta_m_l <- Eta[[m]][l, ]
  #           eta_m_l[1] <- 1 # TODO remove post-dev
  #           active_feature_index <- which(eta_m_l == 1)
  #           if (length(active_feature_index) == 0) {
  #             "Fix this bug!"
  #           }
  #           
  #           b_0_l <- b_0[[m]][l]
  #           B_m_l <- B_m[l, ]
  #           R_m_l <- R_m[l, ]
  #           lambda_m_l <- Lambda2[[m]][l, ]
  # 
  #           logsum <- numeric(1)
  #           for (j in 1:length(eta_m_l)) {
  #             
  #             P_j <- group_list[[m]][j, ]
  #             
  #             summ <- numeric(1)
  #             for (k in 1:K) {
  #               summ <- summ + P_j[k] * R_m_l[k] * B_m_l[k]
  #             }
  #             beta <- b_0_l + summ
  # 
  #             p_j <- prior_lambda_alpha * log(beta) +
  #               -beta * lambda_m_l[j] +
  #               (prior_lambda_alpha-1) * log(lambda_m_l[j]) 
  #             # TODO: Complete summation while accounting for issues in notation.
  #             # + log((1-))
  #             
  #           }
  #           
  #         }
  # 
  #         alpha <- min(1, rb_log_target_density(rb_prime) - rb_log_target_density(rb))
  #         
  #         # Accept or reject the proposed values
  #         
  #         if (runif(1) < alpha) {
  #           R_m[l,k] <- rb_prime$r
  #           B_m[l,k] <- rb_prime$b
  #         } 
  #         
  #       }
  #       
  #     }
  #     
  # 
  #       
  #       
  #   }
  #   
  #   # }
  #   
  # # }
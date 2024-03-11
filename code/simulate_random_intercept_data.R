library(tidyverse)

## -- Define custom functions

# simulate_A assumes all features important across active components
simulate_A <- function(n_important_components, n_important_features) {
  A <- matrix(0, nrow = r, ncol = p_m) 
  if (n_important_components > 0) {
    index_important_components <- seq(to = n_important_components)
    index_important_features <- seq(to = n_important_features)
    n_nonzero_a <- n_important_components * n_important_features 
    nonzero_a <- matrix(rnorm(n_nonzero_a), 
                        nrow = n_important_components, 
                        ncol = n_important_features)
    A[index_important_components, index_important_features] <- nonzero_a
  }
  return(A)
}

simulate_omics_data <- function(n_views=2, n_obs=200, p_m=10, r=4,
                              prob_feature_importance=0.5, 
                              prob_component_importance=0.5) {
  
  n_important_features <- floor(prob_feature_importance*p_m)
  n_important_components <- floor(prob_component_importance*r)
  
  index_important_components <- seq(to = n_important_components)
  index_important_features <- seq(to = n_important_features) 
  
  X_list <- list()
  U <- matrix(data = rnorm(n_obs*r), nrow = n_obs, ncol = r)
  A_list <- list()
  E_list <- list()
  
  for (m in 1:n_views) {
    A_list[[m]] <- simulate_A(n_important_components, n_important_features)
    E_list[[m]] <- matrix(data = rnorm(n_obs*p_m), nrow = n_obs, ncol = p_m)
    X_list[[m]] <- U %*% A_list[[m]] + E_list[[m]]
  }
  
  omics_results <- list(X=X_list, U=U, A=A_list,
                        index_important_components= index_important_components,
                        index_important_features= index_important_features)
  
  return(omics_results)
}

simulate_re_data <- function(n_views=2, n_obs=200, p_m=10, r=4,
                          prob_feature_importance=0.5, 
                          prob_component_importance=0.5,
                          nu2=1, n_sites=5, n_covars=1) {
  
  omics_data <- simulate_omics_data(n_views, n_obs, p_m, r, prob_feature_importance, prob_component_importance)
  
  # Outcome model
  alpha_0 <- rnorm(1)
  alpha <- matrix(0, nrow = r, ncol = 1)
  alpha[index_important_components, ] <- rnorm(length(index_important_components))
  
  # Covariates
  W <- matrix(rnorm(n_obs), nrow = n_obs, ncol = n_covars)
  beta <- rnorm(n_covars) %>% matrix(ncol = 1)
  
  # Sample random intercept design matrix
  Z <- rmultinom(n_obs, size = 1, prob = rep(1/ n_sites, n_sites)) %>% t() # convert into n x n_sites design matrix
  n_per_site <- colSums(Z)
  xi_s <- rnorm(n_sites, sd = sqrt(nu2)) %>% matrix(ncol = 1) # make into column vector
  
  # Observation-level residusal
  epsilon <- matrix(rnorm(n_obs), nrow = n_obs)
  
  Y <- alpha_0 + Z %*% xi_s + omics_data$U %*% alpha + W %*% beta + epsilon
  
  return(list(Y=Y, alpha_0=alpha_0, alpha=alpha, W=W, beta=beta, 
              nu2=nu2, n_sites=n_sites, Z=Z, n_per_site=n_per_site,
              xi_s=xi_s, X=omics_data$X, U=omics_data$U, A=omics_data$A,
              index_important_components=index_important_components, 
              index_important_features=index_important_features))
}

## -- Simulate data
set.seed(473)
data <- simulate_data()

# # Calculate outcome
# y_tilde <- y - alpha_0 - U%*%A - W%*%beta # Essentially Z%*%xi_s + epsilon
# 
# est_conditional_params <- function(prior_mu, prior_var, post_var, x) {
#   n <- nrow(x)
#   conditional_var <- (1/ prior_var + n/ post_var)^(-1)
#   conditional_mu <- conditional_var * (prior_mu/ prior_var + sum(x)/ post_var)
#   return(c(conditional_mu, conditional_var))
# }
# 
# # Get conditional distribution paramters for each site
# conditional_params <- matrix(NA, nrow = n_sites, ncol = 2)
# for (s in 1:n_sites) {
#   site_index <- which(Z[, s]==1)
#   x_site <- y_tilde[site_index, , drop = FALSE]
#   conditional_params[s, ] <- est_conditional_params(prior_mu = 0, prior_var = 1, post_var = 1, x = x_site)
# }
# 
# apply(conditional_params, 2, mean)

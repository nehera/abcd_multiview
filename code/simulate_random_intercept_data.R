# Load packages
library(tidyverse)

# Define default arguments in a list
default_args <- list(
  n_views = 2,
  n_obs = 200,
  p_m = 10,
  r = 4,
  prob_feature_importance = 0.5,
  prob_component_importance = 0.5,
  nu2_site = 1,
  nu2_family = 0.5,
  n_sites = 5,
  n_families_per_site = 3,
  n_individs_per_family = 2,
  n_covars = 1,
  alpha_0=1, 
  sigma2=1,
  seed = 1
)

# Import arguments into the global environment
list2env(default_args, envir = .GlobalEnv)

## -- Define custom functions
# simulate_A assumes all features important across active components
simulate_A <- function(r, p_m, n_important_components, n_important_features) {
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

# Simulates omics data assuming features are active in balanced fashion i.e.
# activation pattern is same across views
# TODO consider scenarios where gamma^(m) and Eta^(m) != for all views
simulate_omics_data <- function(n_views=2, n_obs=200, p_m=10, r=4,
                              prob_feature_importance=0.5, 
                              prob_component_importance=0.5,
                              sigma2=1) {
  
  n_important_features <- floor(prob_feature_importance*p_m)
  n_important_components <- floor(prob_component_importance*r)
  
  index_important_components <- seq(to = n_important_components)
  index_important_features <- seq(to = n_important_features)
  
  gamma <- rep(0, r)
  gamma[index_important_components] <- 1
  Eta <- matrix(0, nrow = r, ncol = p_m)
  Eta[index_important_components, index_important_features] <- 1
  
  X_list <- list()
  U <- matrix(data = rnorm(n_obs*r), nrow = n_obs, ncol = r)
  A_list <- list()
  E_list <- list()
  
  for (m in 1:n_views) {
    A_list[[m]] <- simulate_A(r, p_m, n_important_components, n_important_features)
    E_list[[m]] <- matrix(data = rnorm(n_obs*p_m, sd = sqrt(sigma2)), nrow = n_obs, ncol = p_m)
    X_list[[m]] <- U %*% A_list[[m]] + E_list[[m]]
  }
  
  omics_results <- list(X=X_list, U=U, A=A_list,
                        index_important_components=index_important_components,
                        index_important_features=index_important_features,
                        gamma=gamma, Eta=Eta)
  
  return(omics_results)
}

# simulate_re_data <- function(n_views=2, n_obs=200, p_m=10, r=4,
#                           prob_feature_importance=0.5, 
#                           prob_component_importance=0.5,
#                           nu2=1, n_sites=5, n_covars=1, seed=1) {
#   
#   set.seed(seed)
#   
#   omics_data <- simulate_omics_data(n_views, n_obs, p_m, r, prob_feature_importance, prob_component_importance)
#   
#   # Outcome model
#   alpha_0 <- 1 # Grand intercept fixed
#   alpha <- matrix(0, nrow = r, ncol = 1)
#   
#   alpha[omics_data$index_important_components, ] <- rnorm(length(omics_data$index_important_components))
#   
#   # TODO Covariates
#   # W <- matrix(rnorm(n_obs), nrow = n_obs, ncol = n_covars)
#   # beta <- rnorm(n_covars) %>% matrix(ncol = 1)
#   
#   # Sample random intercept design matrix
#   Z <- rmultinom(n_obs, size = 1, prob = rep(1/ n_sites, n_sites)) %>% t() # convert into n x n_sites design matrix
#   n_per_site <- colSums(Z)
#   xi_s <- rnorm(n_sites, sd = sqrt(nu2)) %>% matrix(ncol = 1) # make into column vector
#   
#   # Observation-level residusal
#   epsilon <- matrix(rnorm(n_obs), nrow = n_obs)
#   
#   Y <- alpha_0 + Z %*% xi_s + omics_data$U %*% alpha + epsilon # TODO W %*% beta 
#   
#   return(list(Y=Y, alpha_0=alpha_0, alpha=alpha, # TODO W=W, beta=beta, 
#               nu2=nu2, n_sites=n_sites, Z=Z, n_per_site=n_per_site,
#               xi_s=xi_s, X=omics_data$X, U=omics_data$U, A=omics_data$A,
#               index_important_components=omics_data$index_important_components, 
#               index_important_features=omics_data$index_important_features))
# }

## -- Simulate data with random effects that are nested
# We start with the balanced case i.e. the same number of families per site, 
# roughly the same n of observations per family, and 
# roughly the same n of observations per site.
simulate_re_data_nested <- function(n_views=2, n_obs=200, p_m=10, r=4,
                                    prob_feature_importance=0.5, 
                                    prob_component_importance=0.5,
                                    nu2_site=1, nu2_family=0.5,
                                    n_sites=5, n_families_per_site=3, 
                                    n_individs_per_family = 2,
                                    n_covars=1, alpha_0=1, sigma2=1, seed=1) {
  
  
  set.seed(seed)
  omics_data <- simulate_omics_data(n_views, n_obs, p_m, r, prob_feature_importance, prob_component_importance)
  
  # Sample latent factor loadings
  alpha <- matrix(0, nrow = r, ncol = 1)
  alpha[omics_data$index_important_components, ] <- rnorm(length(omics_data$index_important_components))
  
  # Sample random intercept design matrix for sites
  Z_site <- rmultinom(n_obs, size = 1, prob = rep(1 / n_sites, n_sites)) %>% t()
  xi_sites <- rnorm(n_sites, sd = sqrt(nu2_site)) %>% matrix(ncol = 1)
  
  # Sample random intercept design matrix for families nested within sites
  Z_family <- matrix(0, nrow = n_obs, ncol = n_sites * n_families_per_site)
  site_assignments <- apply(Z_site, 1, which.max) # Identify the site each observation belongs to
  for (i in 1:n_obs) {
    site_index <- site_assignments[i]
    family_index <- sample(n_families_per_site, 1) # Randomly assign an observation to a family within the site
    Z_family[i, (site_index - 1) * n_families_per_site + family_index] <- 1
  }
  
  # Sample xi_familes centered at xi_sites
  xi_families <- matrix(0, nrow = n_sites, ncol = n_families_per_site)
  for (site in 1:n_sites) {
    xi_families[site, ] <- xi_sites[site] + rnorm(n_families_per_site, mean = 0, sd = sqrt(nu2_family))
  }
  
  temp <- xi_families
  
  print(temp)
  
  # Family effects as a single vector
  xi_families <- as.vector(t(xi_families)) %>% matrix(ncol = 1)
  
  # Sample residuals
  epsilon <- matrix(rnorm(n_obs, sd = sqrt(sigma2)), nrow = n_obs)
  
  # Combine effects
  Y <- alpha_0 + Z_family %*% xi_families + omics_data$U %*% alpha + epsilon
  
  return(list(Y=Y, Z_site=Z_site, Z_family=Z_family, xi_sites=xi_sites, xi_families=xi_families, 
              X=omics_data$X, U=omics_data$U, A=omics_data$A, alpha_0=alpha_0, alpha=alpha,
              gamma=omics_data$gamma, Eta=omics_data$Eta, nu2 = matrix(c(nu2_family, nu2_site), ncol = 1)))
}

################################################################################
# Major mods
################################################################################
simulate_re_data_nested <- function(n_views=2, n_obs=200, p_m=10, r=4,
                                    prob_feature_importance=0.5, 
                                    prob_component_importance=0.5,
                                    nu2_site=1, nu2_family=0.5,
                                    n_sites=5, n_families_per_site=3, 
                                    n_individs_per_family = 2,
                                    n_covars=1, alpha_0=1, sigma2=1, seed=1) {
  
  
  n_obs <- n_sites*n_families_per_site*n_individs_per_family # Remove?
  
  set.seed(seed)
  omics_data <- simulate_omics_data(n_views, n_obs, p_m, r, prob_feature_importance, prob_component_importance)
  
  # Sample latent factor loadings
  alpha <- matrix(0, nrow = r, ncol = 1)
  alpha[omics_data$index_important_components, ] <- rnorm(length(omics_data$index_important_components))
  
  # Specify design matrix for sites
  Z_site <- kronecker(diag(n_sites), rep(1, n_families_per_site*n_individs_per_family))
  xi_sites <- rnorm(n_sites, sd = sqrt(nu2_site)) %>% matrix(ncol = 1)
  
  # Specify design matrix for families nested within sites
  Z_family <- kronecker(diag(n_sites*n_families_per_site), rep(1, n_individs_per_family))
  
  
  # Sample xi_familes centered at xi_sites
  xi_families <- matrix(0, nrow = n_sites, ncol = n_families_per_site)
  for (site in 1:n_sites) {
    xi_families[site, ] <- xi_sites[site] + rnorm(n_families_per_site, mean = 0, sd = sqrt(nu2_family))
  }
  
  # Family effects as a single vector
  xi_families <- as.vector(t(xi_families)) %>% matrix(ncol = 1)
  
  # Sample residuals
  epsilon <- matrix(rnorm(n_obs, sd = sqrt(sigma2)), nrow = n_obs)
  
  # Combine effects
  Y <- alpha_0 + Z_family %*% xi_families + omics_data$U %*% alpha + epsilon
  
  return(list(Y=Y, Z_site=Z_site, Z_family=Z_family, xi_sites=xi_sites, xi_families=xi_families, 
              X=omics_data$X, U=omics_data$U, A=omics_data$A, alpha_0=alpha_0, alpha=alpha,
              gamma=omics_data$gamma, Eta=omics_data$Eta, nu2 = matrix(c(nu2_family, nu2_site), ncol = 1)))
}

## -- Example Usage
data_nested <- simulate_re_data_nested(n_obs=4, n_sites=4, n_families_per_site=4, n_individs_per_family = 4)

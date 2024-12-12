generate_unif_samples <- function(n) {
  # Randomly decide which interval each sample should come from
  # Note, this assumes that each unif interval is of equivalent length
  choose_interval <- sample(c(1, 2), n, replace = TRUE)
  # Pre-allocate the result vector
  sampled_values <- numeric(n)
  # Generate the samples based on the chosen interval
  for (i in 1:n) {
    if (choose_interval[i] == 1) {
      sampled_values[i] <- runif(1, min = -0.5, max = -0.3)
    } else {
      sampled_values[i] <- runif(1, min = 0.3, max = 0.5)
    }
  }
  return(sampled_values)
}

# Simulate loadings
simulate_A <- function(r, features_per_view, n_shared_components, n_important_features) {
  A <- matrix(0, nrow = r, ncol = features_per_view)
  if (n_shared_components > 0) {
    nonzero_a = matrix(generate_unif_samples(n_shared_components * n_important_features), 
                       nrow = n_shared_components, ncol = n_important_features)
    index_shared_components <- seq(to = n_shared_components)
    index_important_features <- seq(to = n_important_features)
    A[index_shared_components, index_important_features] <- nonzero_a
    A[abs(A)<=10^-8] <- 0
  }
  return(A)
}

# Simulates omics data assuming features are active in balanced fashion 
# i.e. activation pattern is same across views
simulate_omics_data <- function(n_views, N_obs, features_per_view, r,
                                prob_feature_importance,
                                prob_component_shared, sigma2) {
  
  n_important_features <- floor(prob_feature_importance*features_per_view)
  n_shared_components <- floor(prob_component_shared*r)
  
  index_shared_components <- seq(to = n_shared_components)
  index_important_features <- seq(to = n_important_features)
  
  gamma <- rep(0, r)
  gamma[index_shared_components] <- 1
  Eta <- matrix(0, nrow = r, ncol = features_per_view)
  Eta[index_shared_components, index_important_features] <- 1
  
  X_list <- list()
  U <- matrix(data = rnorm(N_obs*r), nrow = N_obs, ncol = r)
  A_list <- list()
  E_list <- list()
  
  for (m in 1:n_views) {
    A <- simulate_A(r, features_per_view, n_shared_components, n_important_features)
    epsilon <- rnorm(N_obs*features_per_view, 0, sd = sqrt(sigma2)) %>% matrix(nrow = N_obs, ncol = features_per_view)
    X_list[[m]] <- U %*% A + epsilon
    A_list[[m]] <- A
  }
  
  omics_results <- list(X=X_list, U=U, A=A_list,
                        index_shared_components=index_shared_components,
                        index_important_features=index_important_features,
                        gamma=gamma, Eta=Eta)
  
  return(omics_results)
}

# # # For dev
# seed = 1
# sigma2_ksi_true = 1.5
# sigma2_theta_true = 0.75
# n_views = 4
# features_per_view = 150
# r = 6

# Simulates omics data and then outcome data with random effects
simulation_study_data_generator <- function(seed = 1,
                                            sigma2_ksi_true, 
                                            sigma2_theta_true,
                                            covars = F,
                                            n_views = 4,
                                            features_per_view = 150,
                                            r = 6,
                                            train_set = T,
                                            dev_set = F) {
  
  set.seed(seed)
  
  # For development, we overwrite certain args
  if (dev_set) {
    features_per_view <- 10
  }
  
  # Fix parameters
  prob_feature_importance <- 0.2
  prob_component_important_to_outcome <- 0.5
  prob_component_important_to_omics_views <- 1

  N_sites <- 30
  if (train_set) {
    n_families_per_site <- 100
    if (dev_set) {
      n_families_per_site <- 10
    } 
  } else {
    n_families_per_site <- 25
    if (dev_set) {
      n_families_per_site <- 2 
    }
  }
  n_individs_per_family <- 2
  
  N_families <- N_sites*n_families_per_site
  N_obs <- N_sites*n_families_per_site*n_individs_per_family
  # Specify design matrix for sites
  Z_site <- kronecker(diag(N_sites), rep(1, n_families_per_site*n_individs_per_family))
  # Specify design matrix for families nested within sites
  Z_family <- kronecker(diag(N_families), rep(1, n_individs_per_family))
  
  mu <- 1
  sigma2_ksi <- sigma2_ksi_true # Site variance.
  sigma2_theta <- sigma2_theta_true # Family:site variance (vec of len N_sites)
  sigma2 <- 1 # Residual variance fixed
  
  # Simulate random effects
  ksi <- rnorm(N_sites, mu, sd = sqrt(sigma2_ksi)) %>% matrix(ncol = 1)
  
  if (covars == T) {
    # Simulate covariates
    n_covars <- 2
    W <- matrix(rnorm(N_obs*n_covars), nrow = N_obs, ncol = n_covars)
    beta <- matrix(rep(1, n_covars), ncol = 1)
  } else {
    n_covars <- 1
    W <- matrix(1, nrow = N_obs, ncol = n_covars)
    # Covariates don't have an effect
    beta <- matrix(rep(0, n_covars), ncol = n_covars)
  }
  
  if(length(sigma2_theta) != N_sites) {
    stop("Length of sigma2_theta must equal N_sites.")
  }
  
  # Simulate omics views
  omics_data <- simulate_omics_data(n_views, N_obs, features_per_view, r, 
                                    prob_feature_importance, 
                                    prob_component_shared = prob_component_important_to_omics_views, 
                                    sigma2)
  U <- omics_data$U
  
  # Fix latent factor loadings to 1 in important components
  alpha <- matrix(0, nrow = r, ncol = 1)
  n_important_components <- floor(prob_component_important_to_outcome*r)
  alpha[omics_data$index_important_components, ] <- rep(1, n_important_components) 
  
  # Sample theta_sf|ksi_s ~ N(ksi_s, sigma2_theta_s)
  theta <- matrix(0, nrow = N_sites, ncol = n_families_per_site)
  for (s in 1:N_sites) {
    theta[s, ] <- rnorm(n_families_per_site, mean = ksi[s], sd = sqrt(sigma2_theta[s]))
  }
  
  # Family effects as a single vector
  theta <- as.vector(t(theta)) %>% matrix(ncol = 1)
  
  epsilon <- rnorm(N_obs, 0, sd = sqrt(sigma2)) %>% matrix(ncol = 1)
  
  # Combine effects
  Y <- U %*% alpha + Z_family %*% theta + epsilon
  
  if (covars == T) {
    # Add in the covariate effect
    Y <- Y + W %*% beta
  }
  
  return(list(Y=Y, Z_site=Z_site, Z_family=Z_family, ksi=ksi, theta=theta, 
              X=omics_data$X, U=omics_data$U, A=omics_data$A, 
              mu=mu, alpha=alpha, W=W, beta=beta, gamma=omics_data$gamma, Eta=omics_data$Eta, 
              sigma2 = sigma2, nu2 = list(sigma2_ksi=sigma2_ksi, sigma2_theta=sigma2_theta)))
}
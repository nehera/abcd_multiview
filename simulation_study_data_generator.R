# Simulate loadings
simulate_A <- function(r, p_m, n_important_components, n_important_features) {
  A <- matrix(0, nrow = r, ncol = p_m)
  if (n_important_components > 0) {
    a1=-.5;
    b1=-.3;
    V1a=matrix(0, n_important_components, n_important_features) 
    for (l in 1:n_important_components){
      # V1a[l,]=runif(n_important_features, a1, b1) 
      V1a[l,]=1
    }
    nonzero_a = V1a
    index_important_components <- seq(to = n_important_components)
    index_important_features <- seq(to = n_important_features)
    A[index_important_components, index_important_features] <- nonzero_a
    A[abs(A)<=10^-8] <- 0
  }
  return(A)
}

# Simulates omics data assuming features are active in balanced fashion 
# i.e. activation pattern is same across views
simulate_omics_data <- function(n_views, N_obs, p_m, r,
                                prob_feature_importance,
                                prob_component_importance, sigma2) {
  
  n_important_features <- floor(prob_feature_importance*p_m)
  n_important_components <- floor(prob_component_importance*r)
  
  index_important_components <- seq(to = n_important_components)
  index_important_features <- seq(to = n_important_features)
  
  gamma <- rep(0, r)
  gamma[index_important_components] <- 1
  Eta <- matrix(0, nrow = r, ncol = p_m)
  Eta[index_important_components, index_important_features] <- 1
  
  X_list <- list()
  U <- matrix(data = rnorm(N_obs*r), nrow = N_obs, ncol = r)
  A_list <- list()
  E_list <- list()
  
  for (m in 1:n_views) {
    A <- simulate_A(r, p_m, n_important_components, n_important_features)
    epsilon <- rnorm(N_obs*p_m, 0, sd = sqrt(sigma2)) %>% matrix(nrow = N_obs, ncol = p_m)
    X_list[[m]] <- U %*% A + epsilon
    A_list[[m]] <- A
  }
  
  omics_results <- list(X=X_list, U=U, A=A_list,
                        index_important_components=index_important_components,
                        index_important_features=index_important_features,
                        gamma=gamma, Eta=Eta)
  
  return(omics_results)
}

# Simulates omics data and then outcome data with random effects
simulation_study_data_generator <- function(seed = 1, 
                                            high_signal_to_noise = T,
                                            high_kurtosis = T,
                                            covars = F,
                                            train_set = T) {
  
  set.seed(seed)
  
  # Fix parameters
  N_sites <- 30
  if (train_set) {
    n_families_per_site <- 100
  } else {
    n_families_per_site <- 25
  }
  n_individs_per_family <- 2
  
  N_families <- N_sites*n_families_per_site
  N_obs <- N_sites*n_families_per_site*n_individs_per_family
  # Specify design matrix for sites
  Z_site <- kronecker(diag(N_sites), rep(1, n_families_per_site*n_individs_per_family))
  # Specify design matrix for families nested within sites
  Z_family <- kronecker(diag(N_families), rep(1, n_individs_per_family))
  
  # Mapping of families to sites
  Z_family_to_site <- t(Z_family) %*% Z_site
  
  n_views <- 4
  p_m <- 100
  r <- 4
  prob_feature_importance <- 0.05
  prob_component_importance <- 0.5
  
  mu <- 1
  sigma2_ksi <- 0.8 # Site variance. Note, this is not applicable when kurtosis is set = "High"
  sigma2_theta <- rep(0.2, N_sites) # Family:site variance
  
  if (high_signal_to_noise) {
    sigma2 <- 0.1 # Noise variance is low relative to random effect variances
  } else {
    sigma2 <- 1 # Noise variance is same as sum of random effect variances
  }
  
  if (high_kurtosis == T) {
    # Simulate site effects from t distribution
    # ksi_s ~ mu + t(3)
    ksi <- mu + rt(N_sites, df = 3) %>% matrix(ncol = 1)
  } else {
    # Simulate site effects from normal distrution
    # ksi_s ~ N(mu, sigma2_ksi)
    ksi <- rnorm(N_sites, mu, sd = sqrt(sigma2_ksi)) %>% matrix(ncol = 1)
  }
  
  if (covars == T) {
    # Simulate covariates
    n_covars <- 2
    W <- matrix(rnorm(N_obs*n_covars), nrow = N_obs, ncol = n_covars)
    beta <- matrix(rep(1, n_covars), ncol = 1)
  } 
  
  if(length(sigma2_theta) != N_sites) {
    stop("Length of sigma2_theta must equal N_sites.")
  }
  
  # Simulate omics views
  omics_data <- simulate_omics_data(n_views, N_obs, p_m, r, 
                                    prob_feature_importance, 
                                    prob_component_importance, sigma2)
  U <- omics_data$U
  
  # Fix latent factor loadings to 1 in important components
  alpha <- matrix(0, nrow = r, ncol = 1)
  n_important_components <- floor(prob_component_importance*r)
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
  
  return(list(Y=Y, Z_site=Z_site, Z_family=Z_family, Z_family_to_site=Z_family_to_site, 
              ksi=ksi, theta=theta, X=omics_data$X, U=omics_data$U, A=omics_data$A, 
              mu=mu, alpha=alpha, W=W, beta=beta, gamma=omics_data$gamma, Eta=omics_data$Eta, 
              sigma2 = sigma2, nu2 = list(sigma2_ksi=sigma2_ksi, sigma2_theta=sigma2_theta)))
}
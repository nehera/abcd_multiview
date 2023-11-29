library(tidyverse)

simulate_iid_data <- function(n_views=2, n_obs=200, p_m=10, r=4,
                              prob_feature_importance=0.5, 
                              prob_component_importance=0.5) {

  n_important_features <- floor(prob_feature_importance*p_m)
  n_important_components <- floor(prob_component_importance*r)
  
  index_important_components <- seq(to = n_important_components)
  index_important_features <- seq(to = n_important_features) 

  # simulate_A assumes all features important across active components
  simulate_A <- function(n_important_components, n_important_features) {
    index_important_components <- seq(to = n_important_components)
    index_important_features <- seq(to = n_important_features)
    n_nonzero_a <- n_important_components * n_important_features 
    nonzero_a <- matrix(rnorm(n_nonzero_a), 
                        nrow = n_important_components, 
                        ncol = n_important_features)
    A <- matrix(0, nrow = r, ncol = p_m) 
    A[index_important_components, index_important_features] <- nonzero_a
    return(A)
  }
  
  A_list <- list()
  for (m in 1:n_views) {
    A_list[[m]] <- simulate_A(n_important_components, n_important_features)
  }
  
  U <- matrix(data = rnorm(n_obs*r), nrow = n_obs, ncol = r)
  
  E_list <- list()
  for (m in 1:n_views) {
    E_list[[m]] <- matrix(data = rnorm(n_obs*p_m), nrow = n_obs, ncol = p_m)
  }
  
  X_list <- list()
  for (m in 1:n_views) {
    X_list[[m]] <- U %*% A_list[[m]] + E_list[[m]]
  }

  a <- ifelse(1:r %in% index_important_components, 1, 0) %>%
    matrix(nrow = r)
  e <- matrix(rnorm(n_obs), nrow = n_obs)
    
  Y <- U %*% a + e
  
  simulation_results <- list(X_list=X_list, Y=Y,
                             index_important_components=seq(n_important_components),
                             index_important_features=seq(n_important_features),
                             U=U, A_list=A_list)
  
  return(simulation_results)
}


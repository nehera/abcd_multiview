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
# omics_data <- simulate_omics_data(n_views, N_obs, features_per_view, r, 
#                                   n_important_features, 
#                                   prob_component_shared = prob_component_important_to_omics_views, 
#                                   sigma2)
simulate_omics_data <- function(n_views, N_obs, features_per_view, r,
                                n_important_features,
                                prob_component_shared, sigma2) {
  
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
                                            dev_set = F) {
  
  set.seed(seed)
  
  # Fix parameters
  n_important_features <- 30
  # For development, we overwrite certain args
  if (dev_set) {
    features_per_view <- 10
    n_important_features <- 2
  }
  prob_component_important_to_outcome <- 0.5
  prob_component_important_to_omics_views <- 1

  N_sites <- 30
  if (dev_set) {
      n_families_per_site <- 10
      } else {
        n_families_per_site <- 120
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
  # sigma2 <- .1 # Residual variance fixed
  
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
                                    n_important_features, 
                                    prob_component_shared = prob_component_important_to_omics_views, 
                                    sigma2)
  U <- omics_data$U
  
  # Fix latent factor loadings to 1 in important components
  alpha <- matrix(0, nrow = r, ncol = 1)
  n_important_components <- floor(prob_component_important_to_outcome*r)
  alpha[1:n_important_components, ] <- 1
  
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

split_data_by_seed <- function(data_set, train_frac = 0.8, seed = 123) {
  set.seed(seed)  # Set the seed for reproducibility
  
  src_subject_id <- 1:nrow(data_set$Z_site)
  site_id_l <- apply(data_set$Z_site, 1, function(r) which(r == 1))
  rel_family_id <- apply(data_set$Z_family, 1, function(r) which(r == 1))
  
  cluster_data <- data.frame(
    src_subject_id = src_subject_id,
    site_id_l = site_id_l,
    rel_family_id = rel_family_id
  )
  
  # Function to create train/test split for each site
  split_train_test <- function(data, train_frac) {
    unique_families <- unique(data$rel_family_id)
    train_families <- sample(unique_families, size = floor(train_frac * length(unique_families)))
    train_data <- data %>% filter(rel_family_id %in% train_families)
    test_data <- data %>% filter(!rel_family_id %in% train_families)
    list(train = train_data, test = test_data)
  }
  
  # Split data by site
  split_by_site <- split(cluster_data, cluster_data$site_id_l)
  
  # Apply split function to each site
  split_data <- lapply(split_by_site, split_train_test, train_frac = train_frac)
  
  # Combine train and test sets
  train_data <- bind_rows(lapply(split_data, `[[`, "train"))
  test_data <- bind_rows(lapply(split_data, `[[`, "test"))
  
  # Create train and test sets by subsetting relevant elements
  train_indices_combined <- train_data$src_subject_id
  test_indices_combined <- test_data$src_subject_id
  
  train_set <- list(
    Y = matrix(data_set$Y[train_indices_combined], ncol = 1),
    Z_site = data_set$Z_site[train_indices_combined, ],
    Z_family = data_set$Z_family[train_indices_combined, ],
    X = lapply(data_set$X, function(mat) mat[train_indices_combined, ]),
    U = data_set$U[train_indices_combined, ],
    W = data_set$W[train_indices_combined, ]
  )
  
  test_set <- list(
    Y = matrix(data_set$Y[test_indices_combined], ncol = 1),
    Z_site = data_set$Z_site[test_indices_combined, ],
    Z_family = data_set$Z_family[test_indices_combined, ],
    X = lapply(data_set$X, function(mat) mat[test_indices_combined, ]),
    U = data_set$U[test_indices_combined, ],
    W = data_set$W[test_indices_combined, ]
  )
  
  # Copy unchanged elements to both train and test sets
  elements_to_copy <- c("mu", "alpha", "beta", "gamma", "Eta", "sigma2", "nu2", "ksi", "theta", "A")
  for (element in elements_to_copy) {
    train_set[[element]] <- data_set[[element]]
    test_set[[element]] <- data_set[[element]]
  }
  
  return(list(train_set = train_set, test_set = test_set))
}

# # Example usage
# result <- split_data_by_seed(data_set, train_frac = 0.8, seed = 42)
# train_set <- result$train_set
# test_set <- result$test_set

reshape_for_BIP <- function(train_set, n_views, covars_flag = TRUE) {
  # Initialize the list with the response variable
  trainList <- list(train_set$Y)
  
  # Add views to the trainList
  for (m in 1:n_views) {
    trainList[[m+1]] <- train_set$X[[m]]
  }
  
  # Initialize IndicVar
  IndicVar <- c(1, rep(0, n_views))
  
  # Add covariates if covars_flag is TRUE
  if (covars_flag == TRUE) {
    covar_index <- length(trainList) + 1
    trainList[[covar_index]] <- train_set$W
    IndicVar[covar_index] <- 2
  }
  
  return(list(trainList = trainList, IndicVar = IndicVar))
}

check_BIPmixed_coverage <- function(BIPmixed_result, train_set, n_iter, n_burnin, n_views) {
  # Combine samples into a 3D array
  combine_samples_nested <- function(samples_list, n_iter, n_chains) {
    n_params <- ncol(samples_list[[1]]$mu_samples) +
      ncol(samples_list[[1]]$beta_samples) +  
      ncol(samples_list[[1]]$ksi_samples) +
      ncol(samples_list[[1]]$theta_samples) +
      ncol(samples_list[[1]]$sigma2_ksi_samples) +
      ncol(samples_list[[1]]$sigma2_samples) +
      ncol(samples_list[[1]]$sigma2_theta_samples)
    
    combined_array <- array(NA, dim = c(n_iter, n_chains, n_params))
    
    for (chain in 1:n_chains) {
      samples <- samples_list[[chain]]
      combined_array[, chain, 1] <- samples$mu
      combined_array[, chain, 2:(1 + ncol(samples$beta_samples))] <- samples$beta_samples
      combined_array[, chain, (1 + ncol(samples$beta_samples) + 1):(1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples))] <- samples$ksi_samples
      combined_array[, chain, (1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + 1):(1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples))] <- samples$theta_samples
      combined_array[, chain, 1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + 1] <- samples$sigma2_ksi_samples
      combined_array[, chain, (1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + 2):(1 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + 1 + ncol(samples$sigma2_theta_samples))] <- samples$sigma2_theta_samples
      combined_array[, chain, (2 + ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + ncol(samples$sigma2_theta_samples) + 1)] <- samples$sigma2_samples
    }
    
    return(combined_array)
  }
  
  # Sample list for BIPmixed
  samples_list <- list(BIPmixed_result)
  n_chains <- length(samples_list)
  combined_samples <- combine_samples_nested(samples_list, n_iter, n_chains)
  
  # Parameter names
  param_names <- c("mu",
                   paste0("beta_", 1:ncol(samples_list[[1]]$beta_samples)),
                   paste0("ksi_", 1:ncol(samples_list[[1]]$ksi_samples)),
                   paste0("theta_", 1:ncol(samples_list[[1]]$theta_samples)),
                   "sigma2_ksi",
                   paste0("sigma2_theta_", 1:ncol(samples_list[[1]]$sigma2_theta_samples)),
                   "sigma2")
  
  # Assign parameter names
  dimnames(combined_samples) <- list(
    iterations = NULL,
    chains = NULL,
    parameters = param_names
  )
  
  # Use rstan::monitor to focus on variance parameters
  mcmc_summary <- monitor(combined_samples)
  
  # Extract summary statistics
  mean_values <- round(mcmc_summary$mean, 4)
  median_values <- round(mcmc_summary$`50%`, 4)
  lower_bounds <- round(mcmc_summary$`2.5%`, 4)
  upper_bounds <- round(mcmc_summary$`97.5%`, 4)
  
  # Combine the true values into a single vector, ordered according to the MCMC parameters
  mu_true <- train_set$mu
  ksi_true <- train_set$ksi
  theta_true <- train_set$theta
  beta_true <- train_set$beta
  sigma2_ksi_true <- train_set$nu2$sigma2_ksi
  sigma2_theta_true <- train_set$nu2$sigma2_theta
  sigma2_true <- train_set$sigma2
  true_values <- c(mu_true, beta_true, ksi_true, theta_true, sigma2_ksi_true, sigma2_theta_true, sigma2_true)
  
  # Initial values
  init_values <- samples_list[[1]]$initial_values
  mu_init <- init_values$mu_init
  beta_init <- init_values$beta_init
  ksi_init <- init_values$ksi_init
  theta_init <- init_values$theta_init
  sigma2_ksi_init <- init_values$sigma2_ksi_init
  sigma2_theta_init <- init_values$sigma2_theta_init
  sigma2_init <- init_values$sigma2_init
  initial_values <- c(mu_init, beta_init, ksi_init, theta_init,
                      sigma2_ksi_init, sigma2_theta_init, sigma2_init)
  
  # Create comparison dataframe
  comparison <- data.frame(
    param_name = param_names,
    lower = lower_bounds,
    mean = mean_values,
    median = median_values,
    upper = upper_bounds,
    true_value = round(true_values, 4),
    initial = round(initial_values, 4)
  )
  
  # Check credible interval coverage
  comparison$within_credible_interval <- with(comparison, true_value >= lower & true_value <= upper)
  
  # Filter and analyze variance parameters
  variance_comparison <- comparison %>% filter(grepl("^sigma2$|^sigma2_ksi$|^sigma2_theta_", param_name))
  variance_param_coverage <- data.frame(
    param_name = variance_comparison$param_name,
    within_credible_interval = variance_comparison$within_credible_interval
  ) %>%
    mutate(param_type = case_when(
      param_name == "sigma2_ksi" ~ "sigma2_ksi",
      grepl("^sigma2_theta_", param_name) ~ "sigma2_theta",
      param_name == "sigma2" ~ "sigma2"
    )) %>%
    group_by(param_type) %>%
    summarise(
      total = n(),
      within_ci = sum(within_credible_interval),
      proportion_within_ci = within_ci / total
    )
  
  print(variance_param_coverage)
  
  # Create variance trace plot
  variance_df <- data.frame(
    sigma2 = as.vector(BIPmixed_result$sigma2_samples),
    sigma2_ksi = as.vector(BIPmixed_result$sigma2_ksi_samples)
  ) %>%
    mutate(iter = 1:n()) %>%
    gather(key = "parameter", value = "value", -iter)
  
  variance_trace_plot <- ggplot(variance_df, aes(x = iter, y = value, color = parameter)) +
    geom_line(alpha = 0.7) +
    geom_vline(xintercept = n_burnin, linetype = "dashed", color = "black") + 
    labs(x = "Iteration", y = "Value", color = "Parameter") +
    theme_minimal()
  
  print(variance_trace_plot)
  
  return(list(
    comparison = comparison,
    variance_param_coverage = variance_param_coverage,
    variance_trace_plot = variance_trace_plot
  ))
}

check_global_variable_selection <- function(result, IndicVar, train_set, method_name) {
  # Estimate variable selection performance globally
  eta_true_global <- train_set$TrueVar1 # as.numeric(colSums(train_set$Eta) > 0)
  omics_index <- which(IndicVar == 0)
  
  # Compute variable selection criteria for each omics index
  VarSelGlobalPerformance_list <- lapply(result$VarSelMeanGlobal[omics_index], function(pred.prob) {
    BIPnet::ComputeVarCriteria(pred.prob, eta_true_global, thres = 0.5) %>%
      as.data.frame() %>%
      mutate(Method = method_name)
  })
  
  # Combine list of data frames and add index column
  combine_dataframes_with_index <- function(df_list, index_col_name = "View") {
    for (i in seq_along(df_list)) {
      df_list[[i]][[index_col_name]] <- i
    }
    combined_df <- do.call(rbind, df_list)
    return(combined_df)
  }
  
  # Combine the performance results into a single data frame
  combined_df <- combine_dataframes_with_index(VarSelGlobalPerformance_list)
  
  return(combined_df)
}

reshape_for_BIPpredict <- function(test_set, n_views, covars_flag = TRUE) {
  # Initialize the test list with the views (X)
  testList <- test_set$X
  
  # Initialize IndicVar
  IndicVar <- rep(0, n_views)
  
  # Add covariates if covars_flag is TRUE
  if (covars_flag == TRUE) {
    covar_index <- length(testList) + 1
    testList[[covar_index]] <- test_set$W
  }
  
  return(list(testList = testList, IndicVar = IndicVar))
}

# Define a function to calculate MSE, bias squared, prediction variance
calculate_prediction_metrics <- function(Y, y_preds) {
  # # Scale Y and y_preds
  # Y <- (Y - mean(Y)) / sd(Y)
  # y_preds <- (y_preds - mean(y_preds)) / sd(y_preds)
  # Calculate performance metrics
  mse <- mean((Y - y_preds)^2)
  bias2 <- (mean(y_preds) - mean(Y))^2
  var_pred <- var(y_preds)
  mean_pred <- mean(y_preds)
  corr_pred <- cor(Y, y_preds)
  data.frame(MSE = mse, Bias2 = bias2, 
             Variance = var_pred, 
             Mean = mean_pred,
             Correlation = corr_pred)
}
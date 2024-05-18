## -- Start with Unnested

library(parallel)
library(MASS)
library(rstan)
library(lme4)
library(dplyr)

# Modified Gibbs sampler to accept a seed
gibbs_sampler <- function(y, X, study_site, n_iter, priors, seed) {
  set.seed(seed)
  
  n <- length(y)
  J <- length(unique(study_site))
  
  # Priors
  beta_prior_mean <- priors$beta_prior_mean
  beta_prior_var <- priors$beta_prior_var
  sigma_u_prior_a <- priors$sigma_u_prior_a
  sigma_u_prior_b <- priors$sigma_u_prior_b
  sigma_prior_a <- priors$sigma_prior_a
  sigma_prior_b <- priors$sigma_prior_b
  
  # Initialize parameters
  beta <- rep(0, 2)
  u <- rep(0, J)
  sigma_u2 <- 1
  sigma2 <- 1
  
  # Storage for samples
  samples <- list(beta = matrix(NA, nrow = n_iter, ncol = 2),
                  u = matrix(NA, nrow = n_iter, ncol = J),
                  sigma_u2 = rep(NA, n_iter),
                  sigma2 = rep(NA, n_iter))
  
  for (iter in 1:n_iter) {
    # Sample beta
    V_beta <- solve(t(X) %*% X / sigma2 + diag(1 / beta_prior_var, 2))
    m_beta <- V_beta %*% (t(X) %*% (y - u[study_site]) / sigma2 + beta_prior_mean / beta_prior_var)
    beta <- mvrnorm(1, mu = m_beta, Sigma = V_beta)
    
    # Sample u
    for (j in 1:J) {
      idx <- which(study_site == j)
      y_j <- y[idx]
      X_j <- X[idx, ]
      V_u <- 1 / (length(y_j) / sigma2 + 1 / sigma_u2)
      m_u <- V_u * sum(y_j - X_j %*% beta) / sigma2
      u[j] <- rnorm(1, mean = m_u, sd = sqrt(V_u))
    }
    
    # Sample sigma_u2
    a_sigma_u <- sigma_u_prior_a + J / 2
    b_sigma_u <- sigma_u_prior_b + sum(u^2) / 2
    sigma_u2 <- 1 / rgamma(1, shape = a_sigma_u, rate = b_sigma_u)
    
    # Sample sigma2
    a_sigma <- sigma_prior_a + n / 2
    b_sigma <- sigma_prior_b + sum((y - X %*% beta - u[study_site])^2) / 2
    sigma2 <- 1 / rgamma(1, shape = a_sigma, rate = b_sigma)
    
    # Store samples
    samples$beta[iter, ] <- beta
    samples$u[iter, ] <- u
    samples$sigma_u2[iter] <- sigma_u2
    samples$sigma2[iter] <- sigma2
  }
  
  return(samples)
}

# Function to combine samples into a 3D array for rstan::monitor
combine_samples <- function(samples_list) {
  n_iter <- nrow(samples_list[[1]]$beta)
  n_chains <- length(samples_list)
  n_params <- ncol(samples_list[[1]]$beta) + ncol(samples_list[[1]]$u) + 2  # beta, u, sigma_u2, sigma2
  
  # Create an empty array
  combined_array <- array(NA, dim = c(n_iter, n_chains, n_params))
  
  # Fill in the array
  for (chain in 1:n_chains) {
    samples <- samples_list[[chain]]
    combined_array[, chain, 1:ncol(samples$beta)] <- samples$beta
    combined_array[, chain, (ncol(samples$beta) + 1):(ncol(samples$beta) + ncol(samples$u))] <- samples$u
    combined_array[, chain, ncol(samples$beta) + ncol(samples$u) + 1] <- samples$sigma_u2
    combined_array[, chain, ncol(samples$beta) + ncol(samples$u) + 2] <- samples$sigma2
  }
  
  return(combined_array)
}

# Example data
set.seed(123)
n_iter <- 10000
n_chains <- 4
n <- 300
J <- 30
study_site <- rep(1:J, each = n / J)
X <- cbind(1, rnorm(n))
beta_true <- c(1, 2)
sigma_u_true <- 2
u_true <- rnorm(J, sd = sigma_u_true)
sigma_true <- 1
y <- X %*% beta_true + u_true[study_site] + rnorm(n, sd = sigma_true)

# Priors
priors <- list(beta_prior_mean = c(0, 0),
               beta_prior_var = c(100, 100),
               sigma_u_prior_a = 2,
               sigma_u_prior_b = 1,
               sigma_prior_a = 2,
               sigma_prior_b = 1)

# Run Gibbs sampler in parallel
seeds <- 1:n_chains
samples_list <- mclapply(seeds, function(seed) gibbs_sampler(y, X, study_site, n_iter, priors, seed), mc.cores = n_chains)

# Combine samples
combined_samples <- combine_samples(samples_list)

# Assign parameter names
dimnames(combined_samples) <- list(
  iterations = NULL,
  chains = NULL,
  parameters = c(paste0("beta_", 1:ncol(samples_list[[1]]$beta)),
                 paste0("u_", 1:ncol(samples_list[[1]]$u)),
                 "sigma_u2", "sigma2")
)

# Use rstan::monitor
mcmc_summary <- monitor(combined_samples)

# Extract summary statistics
mean_values <- round(mcmc_summary$mean, 4)
median_values <- round(mcmc_summary$`50%`, 4)
lower_bounds <- round(mcmc_summary$`2.5%`, 4)
upper_bounds <- round(mcmc_summary$`97.5%`, 4)

# Combine the true values into a single vector, ordered according to the MCMC parameters
true_values <- c(beta_true, u_true, sigma_u_true^2, sigma_true^2)

# Parameter names
param_names <- c(paste0("beta_", 1:ncol(samples_list[[1]]$beta)),
                 paste0("u_", 1:ncol(samples_list[[1]]$u)),
                 "sigma_u2", "sigma2")

# Initialize a dataframe to hold the comparisons
comparison <- data.frame(
  param_name = param_names,
  lower = lower_bounds,
  mean = mean_values,
  median = median_values,
  upper = upper_bounds,
  true_value = round(true_values, 4)
)

# Add a logical vector to see if true values are within the credible intervals
comparison$within_credible_interval <- with(comparison, true_value >= lower & true_value <= upper)

# Print the result to check each parameter
print(comparison)

# % of all parameters within credible interval
cat("% of all parameters within credible interval: ", mean(comparison$within_credible_interval), "\n")

# Filter for random intercept parameters (those starting with "u_")
u_comparison <- comparison %>% filter(grepl("^u_", param_name))

# % of random intercept parameters within credible interval
cat("% of random intercept parameters within credible interval: ", mean(u_comparison$within_credible_interval), "\n")

# % of random intercept parameters with correct sign
cat("% of sampled RE parameters with correct sign: ", mean(sign(u_comparison$mean) == sign(u_comparison$true_value)), "\n")

# Frequentist model
df_temp <- data.frame(y = y, x = X[, 2], site = factor(study_site))
freq_model <- lmer(y ~ x + (1 | site), data = df_temp)
summary(freq_model)

## -- Nested with N Study Sites Worth of Family:Site Variance Parameters

library(parallel)
library(MASS)
library(rstan)
library(lme4)
library(dplyr)

# Modified Gibbs sampler to accept a seed and handle nested random effects
gibbs_sampler_nested <- function(y, X, study_site, family, n_iter, priors, seed) {
  set.seed(seed)
  
  n <- length(y)
  J <- length(unique(study_site))
  F_count <- length(unique(family))
  
  # Priors
  beta_prior_mean <- priors$beta_prior_mean
  beta_prior_var <- priors$beta_prior_var
  sigma_u_prior_a <- priors$sigma_u_prior_a
  sigma_u_prior_b <- priors$sigma_u_prior_b
  sigma_v_prior_a <- priors$sigma_v_prior_a
  sigma_v_prior_b <- priors$sigma_v_prior_b
  sigma_prior_a <- priors$sigma_prior_a
  sigma_prior_b <- priors$sigma_prior_b
  
  # Initialize parameters
  beta <- rep(0, 2)
  u <- rep(0, J)
  v <- rep(0, F_count)
  sigma_u2 <- 1
  sigma_v2 <- rep(1, J)  # One variance per study site
  sigma2 <- 1
  
  # Storage for samples
  samples <- list(beta = matrix(NA, nrow = n_iter, ncol = 2),
                  u = matrix(NA, nrow = n_iter, ncol = J),
                  v = matrix(NA, nrow = n_iter, ncol = F_count),
                  sigma_u2 = rep(NA, n_iter),
                  sigma_v2 = matrix(NA, nrow = n_iter, ncol = J),
                  sigma2 = rep(NA, n_iter))
  
  for (iter in 1:n_iter) {
    # Sample beta
    V_beta <- solve(t(X) %*% X / sigma2 + diag(1 / beta_prior_var, 2))
    m_beta <- V_beta %*% (t(X) %*% (y - u[study_site] - v[family]) / sigma2 + beta_prior_mean / beta_prior_var)
    beta <- mvrnorm(1, mu = m_beta, Sigma = V_beta)
    
    # Sample u
    for (j in 1:J) {
      idx <- which(study_site == j)
      y_j <- y[idx]
      X_j <- X[idx, ]
      v_j <- v[family[idx]]
      V_u <- 1 / (length(y_j) / sigma2 + 1 / sigma_u2)
      m_u <- V_u * sum(y_j - X_j %*% beta - v_j) / sigma2
      u[j] <- rnorm(1, mean = m_u, sd = sqrt(V_u))
    }
    
    # Sample v
    for (f in 1:F_count) {
      idx <- which(family == f)
      y_f <- y[idx]
      X_f <- X[idx, ]
      u_f <- u[study_site[idx]]
      site <- study_site[idx[1]]  # All indices in idx have the same site
      V_v <- 1 / (length(y_f) / sigma2 + 1 / sigma_v2[site])
      m_v <- V_v * sum(y_f - X_f %*% beta - u_f) / sigma2
      v[f] <- rnorm(1, mean = m_v, sd = sqrt(V_v))
    }
    
    # Sample sigma_u2
    a_sigma_u <- sigma_u_prior_a + J / 2
    b_sigma_u <- sigma_u_prior_b + sum(u^2) / 2
    sigma_u2 <- 1 / rgamma(1, shape = a_sigma_u, rate = b_sigma_u)
    
    # Sample sigma_v2 for each study site
    for (j in 1:J) {
      families_in_site <- unique(family[study_site == j])
      a_sigma_v <- sigma_v_prior_a + length(families_in_site) / 2
      b_sigma_v <- sigma_v_prior_b + sum(v[families_in_site]^2) / 2
      sigma_v2[j] <- 1 / rgamma(1, shape = a_sigma_v, rate = b_sigma_v)
    }
    
    # Sample sigma2
    a_sigma <- sigma_prior_a + n / 2
    b_sigma <- sigma_prior_b + sum((y - X %*% beta - u[study_site] - v[family])^2) / 2
    sigma2 <- 1 / rgamma(1, shape = a_sigma, rate = b_sigma)
    
    # Store samples
    samples$beta[iter, ] <- beta
    samples$u[iter, ] <- u
    samples$v[iter, ] <- v
    samples$sigma_u2[iter] <- sigma_u2
    samples$sigma_v2[iter, ] <- sigma_v2
    samples$sigma2[iter] <- sigma2
  }
  
  return(samples)
}

# Function to combine samples into a 3D array for rstan::monitor
combine_samples_nested <- function(samples_list) {
  n_iter <- nrow(samples_list[[1]]$beta)
  n_chains <- length(samples_list)
  n_params <- ncol(samples_list[[1]]$beta) + ncol(samples_list[[1]]$u) + ncol(samples_list[[1]]$v) + 2 + ncol(samples_list[[1]]$sigma_v2)  # beta, u, v, sigma_u2, sigma2, sigma_v2
  
  # Create an empty array
  combined_array <- array(NA, dim = c(n_iter, n_chains, n_params))
  
  # Fill in the array
  for (chain in 1:n_chains) {
    samples <- samples_list[[chain]]
    combined_array[, chain, 1:ncol(samples$beta)] <- samples$beta
    combined_array[, chain, (ncol(samples$beta) + 1):(ncol(samples$beta) + ncol(samples$u))] <- samples$u
    combined_array[, chain, (ncol(samples$beta) + ncol(samples$u) + 1):(ncol(samples$beta) + ncol(samples$u) + ncol(samples$v))] <- samples$v
    combined_array[, chain, ncol(samples$beta) + ncol(samples$u) + ncol(samples$v) + 1] <- samples$sigma_u2
    combined_array[, chain, ncol(samples$beta) + ncol(samples$u) + ncol(samples$v) + 2] <- samples$sigma2
    combined_array[, chain, (ncol(samples$beta) + ncol(samples$u) + ncol(samples$v) + 3):(ncol(samples$beta) + ncol(samples$u) + ncol(samples$v) + 2 + ncol(samples$sigma_v2))] <- samples$sigma_v2
  }
  
  return(combined_array)
}

# Example data
set.seed(123)
n_iter <- 10000
n_chains <- 4
n <- 900
J <- 30
F_count <- 30
study_site <- rep(1:J, each = n / J)
family <- rep(1:F_count, each = n / F_count)
X <- cbind(1, rnorm(n))
beta_true <- c(1, 2)
sigma_u_true <- 2
sigma_v_true <- 1.5 # Assume each f:s variance parameter is same across sites
u_true <- rnorm(J, sd = sigma_u_true)
v_true <- rnorm(F_count, sd = sigma_v_true)
sigma_true <- 1
y <- X %*% beta_true + u_true[study_site] + v_true[family] + rnorm(n, sd = sigma_true)

# Priors
priors <- list(beta_prior_mean = c(0, 0),
               beta_prior_var = c(100, 100),
               sigma_u_prior_a = 2,
               sigma_u_prior_b = 1,
               sigma_v_prior_a = 2,
               sigma_v_prior_b = 1,
               sigma_prior_a = 2,
               sigma_prior_b = 1)

# Run Gibbs sampler in parallel
seeds <- 1:n_chains
samples_list <- mclapply(seeds, function(seed) gibbs_sampler_nested(y, X, study_site, family, n_iter, priors, seed), mc.cores = n_chains)

# Combine samples
combined_samples <- combine_samples_nested(samples_list)

# Assign parameter names
dimnames(combined_samples) <- list(
  iterations = NULL,
  chains = NULL,
  parameters = c(paste0("beta_", 1:ncol(samples_list[[1]]$beta)),
                 paste0("u_", 1:ncol(samples_list[[1]]$u)),
                 paste0("v_", 1:ncol(samples_list[[1]]$v)),
                 "sigma_u2", "sigma2",
                 paste0("sigma_v2_", 1:ncol(samples_list[[1]]$sigma_v2)))
)

# Use rstan::monitor
mcmc_summary <- monitor(combined_samples)

# Extract summary statistics
mean_values <- round(mcmc_summary$mean, 4)
median_values <- round(mcmc_summary$`50%`, 4)
lower_bounds <- round(mcmc_summary$`2.5%`, 4)
upper_bounds <- round(mcmc_summary$`97.5%`, 4)

# Combine the true values into a single vector, ordered according to the MCMC parameters
true_values <- c(beta_true, u_true, v_true, sigma_u_true^2, sigma_true^2, rep(sigma_v_true^2, J))

# Parameter names
param_names <- c(paste0("beta_", 1:ncol(samples_list[[1]]$beta)),
                 paste0("u_", 1:ncol(samples_list[[1]]$u)),
                 paste0("v_", 1:ncol(samples_list[[1]]$v)),
                 "sigma_u2", "sigma2",
                 paste0("sigma_v2_", 1:ncol(samples_list[[1]]$sigma_v2)))

# Initialize a dataframe to hold the comparisons
comparison <- data.frame(
  param_name = param_names,
  lower = lower_bounds,
  mean = mean_values,
  median = median_values,
  upper = upper_bounds,
  true_value = round(true_values, 4)
)

# Add a logical vector to see if true values are within the credible intervals
comparison$within_credible_interval <- with(comparison, true_value >= lower & true_value <= upper)

# % of all parameters within credible interval
cat("% of all parameters within credible interval: ", mean(comparison$within_credible_interval), "\n")

print("Fixed Effect Estimation vs. Truth:")

# Filter for fixed effect parameters
fixed_effect_comparison <- comparison %>% filter(grepl("^beta_", param_name))

# Print the result to check fixed effect parameters
print(fixed_effect_comparison)

print("Random Intercept Estimation vs. Truth:")

# Filter for random intercept parameters (those starting with "u_")
random_intercept_comparison <- comparison %>% filter(grepl("^u_|^v_", param_name))

# % of random intercept parameters within credible interval
cat("% of random intercept parameters within credible interval: ", mean(random_intercept_comparison$within_credible_interval), "\n")

# % of random intercept parameters with correct sign
random_intercept_comparison$correct_sign <- sign(random_intercept_comparison$mean) == sign(random_intercept_comparison$true_value)
cat("% of random intercept parameters with correct sign: ", mean(random_intercept_comparison$correct_sign), "\n")

# Print the result to check random intercept parameters
print(random_intercept_comparison)

print("Variance Parameter Estimation vs. Truth:")

# Filter the comparison dataframe for variance parameters
variance_comparison <- comparison %>% filter(grepl("^sigma2$|^sigma_u2$|^sigma_v2_", param_name))

# Print the filtered result to check each variance parameter
print(variance_comparison)

# % of variance parameters within credible interval
cat("% of variance parameters within credible interval: ", mean(variance_comparison$within_credible_interval), "\n")

# Make traceplots

library(ggplot2)
library(gridExtra)

# Set the number of chains to plot
n_chains_to_plot <- 2

# Define the parameters to plot
parameters_to_plot <- c(
  paste0("beta_", 1:ncol(samples_list[[1]]$beta)),
  paste0("u_", 1),
  paste0("v_", 1),
  "sigma_u2",
  "sigma2",
  paste0("sigma_v2_", 1:2)
)

# True values for the parameters
true_values <- c(
  beta_true, 
  u_true[1], 
  v_true[1], 
  sigma_u_true^2, 
  sigma_true^2, 
  rep(sigma_v_true^2, 2)
)

# Prepare a list to store ggplot objects
plot_list <- list()

# Create trace plots for each parameter in each chain
for (i in seq_along(parameters_to_plot)) {
  param <- parameters_to_plot[i]
  true_value <- true_values[i]
  
  for (chain in 1:n_chains_to_plot) {
    # Extract the values for the parameter and chain
    param_values <- combined_samples[, chain, which(dimnames(combined_samples)[[3]] == param)]
    
    # Create a dataframe for ggplot
    df <- data.frame(
      iter = 1:n_iter,
      value = param_values,
      chain = factor(chain)
    )
    
    # Create the trace plot
    p <- ggplot(df, aes(x = iter, y = value, color = chain)) +
      geom_line() +
      geom_hline(yintercept = true_value, linetype = "dashed", color = "black") +
      geom_vline(xintercept = n_iter / 2, linetype = "dashed", color = "red") +
      labs(x = ifelse(param == tail(parameters_to_plot, 1), "Iteration", ""), y = ifelse(chain == 1, param, "")) +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Remove y-axis label for chains > 1
    if (chain > 1) {
      p <- p + theme(axis.title.y = element_blank())
    }
    
    # Add the plot to the list
    plot_list[[paste(param, chain, sep = "_")]] <- p
  }
}

# Arrange the plots into a grid
grid_plots <- do.call(grid.arrange, c(plot_list, ncol = n_chains_to_plot))

# Make Density Plots

library(ggplot2)
library(gridExtra)
library(scales)

# Set the number of chains to plot
n_chains_to_plot <- 2

# Define the parameters to plot
parameters_to_plot <- c(
  paste0("beta_", 1:ncol(samples_list[[1]]$beta)),
  paste0("u_", 1),
  paste0("v_", 1),
  "sigma_u2",
  "sigma2",
  paste0("sigma_v2_", 1:2)
)

# True values for the parameters
true_values <- c(
  beta_true, 
  u_true[1], 
  v_true[1], 
  sigma_u_true^2, 
  sigma_true^2, 
  rep(sigma_v_true^2, 2)
)

# Determine the x-axis limits for each class of parameters
beta_limits <- range(combined_samples[(n_iter / 2 + 1):n_iter, , grep("^beta_", dimnames(combined_samples)[[3]])])
u_limits <- range(combined_samples[(n_iter / 2 + 1):n_iter, , grep("^u_", dimnames(combined_samples)[[3]])])
v_limits <- range(combined_samples[(n_iter / 2 + 1):n_iter, , grep("^v_", dimnames(combined_samples)[[3]])])
sigma_u2_limits <- range(combined_samples[(n_iter / 2 + 1):n_iter, , grep("^sigma_u2", dimnames(combined_samples)[[3]])])
sigma2_limits <- range(combined_samples[(n_iter / 2 + 1):n_iter, , grep("^sigma2", dimnames(combined_samples)[[3]])])
sigma_v2_limits <- c(0, 25)  # Truncate to a maximum of 25

# Determine the y-axis limits to be (0, 1) for all density plots
y_limits <- c(0, 1)

# Prepare a list to store ggplot objects
plot_list <- list()

# Create density plots for each parameter in each chain after burn-in
for (i in seq_along(parameters_to_plot)) {
  param <- parameters_to_plot[i]
  true_value <- true_values[i]
  
  # Determine x-axis limits based on parameter class
  if (grepl("^beta_", param)) {
    x_limits <- beta_limits
  } else if (grepl("^u_", param)) {
    x_limits <- u_limits
  } else if (grepl("^v_", param)) {
    x_limits <- v_limits
  } else if (grepl("^sigma_u2", param)) {
    x_limits <- sigma_u2_limits
  } else if (grepl("^sigma2", param)) {
    x_limits <- sigma2_limits
  } else if (grepl("^sigma_v2_", param)) {
    x_limits <- sigma_v2_limits
  }
  
  for (chain in 1:n_chains_to_plot) {
    # Extract the values for the parameter and chain after burn-in
    param_values <- combined_samples[(n_iter / 2 + 1):n_iter, chain, which(dimnames(combined_samples)[[3]] == param)]
    
    # Create a dataframe for ggplot
    df <- data.frame(
      value = param_values,
      chain = factor(chain)
    )
    
    # Compute the density
    density_values <- density(param_values)
    df_density <- data.frame(
      x = density_values$x,
      y = density_values$y / max(density_values$y),  # Normalize the density to unit scale
      chain = factor(chain)
    )
    
    # Create the density plot
    p <- ggplot(df_density, aes(x = x, y = y, fill = chain)) +
      geom_line(aes(color = chain), size = 1) +
      geom_vline(xintercept = true_value, linetype = "dashed", color = "black") +
      scale_x_continuous(limits = x_limits, labels = label_number(accuracy = 0.1)) +
      scale_y_continuous(limits = y_limits, labels = label_number(accuracy = 0.1)) +
      labs(x = ifelse(param == tail(parameters_to_plot, 1), "Value", ""), y = ifelse(chain == 1, param, "")) +
      theme_minimal() +
      theme(legend.position = "none", axis.title.y = element_text(size = 8))
    
    # Remove y-axis label for chains > 1
    if (chain > 1) {
      p <- p + theme(axis.title.y = element_blank())
    }
    
    # Add the plot to the list
    plot_list[[paste(param, chain, sep = "_")]] <- p
  }
}

# Arrange the plots into a grid
grid_plots <- do.call(grid.arrange, c(plot_list, ncol = n_chains_to_plot))

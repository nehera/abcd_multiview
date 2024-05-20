# Install Rcpp and RcppArmadillo if you haven't already
# install.packages("Rcpp")
# install.packages("RcppArmadillo")

library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(rstan)

# Source the C++ code
sourceCpp("gibbs_sampler_nested.cpp")

# Combine samples into a 3D array for rstan::monitor
combine_samples_nested <- function(samples_list) {
  n_iter <- nrow(samples_list[[1]]$beta_samples)
  n_chains <- length(samples_list)
  n_params <- ncol(samples_list[[1]]$beta_samples) + ncol(samples_list[[1]]$u_samples) + ncol(samples_list[[1]]$v_samples) + 2 + ncol(samples_list[[1]]$sigma_v2_samples)  # beta, u, v, sigma_u2, sigma2, sigma_v2
  
  # Create an empty array
  combined_array <- array(NA, dim = c(n_iter, n_chains, n_params))
  
  # Fill in the array
  for (chain in 1:n_chains) {
    samples <- samples_list[[chain]]
    combined_array[, chain, 1:ncol(samples$beta_samples)] <- samples$beta_samples
    combined_array[, chain, (ncol(samples$beta_samples) + 1):(ncol(samples$beta_samples) + ncol(samples$u_samples))] <- samples$u_samples
    combined_array[, chain, (ncol(samples$beta_samples) + ncol(samples$u_samples) + 1):(ncol(samples$beta_samples) + ncol(samples$u_samples) + ncol(samples$v_samples))] <- samples$v_samples
    combined_array[, chain, ncol(samples$beta_samples) + ncol(samples$u_samples) + ncol(samples$v_samples) + 1] <- samples$sigma_u2_samples
    combined_array[, chain, ncol(samples$beta_samples) + ncol(samples$u_samples) + ncol(samples$v_samples) + 2] <- samples$sigma2_samples
    combined_array[, chain, (ncol(samples$beta_samples) + ncol(samples$u_samples) + ncol(samples$v_samples) + 3):(ncol(samples$beta_samples) + ncol(samples$u_samples) + ncol(samples$v_samples) + 2 + ncol(samples$sigma_v2_samples))] <- samples$sigma_v2_samples
  }
  
  return(combined_array)
}

# Example data
set.seed(123)
n_iter <- 6000
n_chains <- 4
n <- 30 * 30 * 3
J <- 30 # Number of Sites
F_count <- 30 * 30 # Number of Families in Total
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
# Note, seed cannot be set in C++
samples_list <- mclapply(seeds, function(seed) { 
  gibbs_sampler_nested(y, X, study_site, family, n_iter, priors$beta_prior_mean, priors$beta_prior_var, priors$sigma_u_prior_a, priors$sigma_u_prior_b, priors$sigma_v_prior_a, priors$sigma_v_prior_b, priors$sigma_prior_a, priors$sigma_prior_b)
}, mc.cores = n_chains)

# Combine samples
combined_samples <- combine_samples_nested(samples_list)

# Assign parameter names
dimnames(combined_samples) <- list(
  iterations = NULL,
  chains = NULL,
  parameters = c(paste0("beta_", 1:ncol(samples_list[[1]]$beta_samples)),
                 paste0("u_", 1:ncol(samples_list[[1]]$u_samples)),
                 paste0("v_", 1:ncol(samples_list[[1]]$v_samples)),
                 "sigma_u2", "sigma2",
                 paste0("sigma_v2_", 1:ncol(samples_list[[1]]$sigma_v2_samples)))
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
param_names <- c(paste0("beta_", 1:ncol(samples_list[[1]]$beta_samples)),
                 paste0("u_", 1:ncol(samples_list[[1]]$u_samples)),
                 paste0("v_", 1:ncol(samples_list[[1]]$v_samples)),
                 "sigma_u2", "sigma2",
                 paste0("sigma_v2_", 1:ncol(samples_list[[1]]$sigma_v2_samples)))

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

# Set the number of chains to plot
n_chains_to_plot <- 2

# Define the parameters to plot
parameters_to_plot <- c(
  paste0("beta_", 1:ncol(samples_list[[1]]$beta_samples)),
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

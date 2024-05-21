#include <RcppArmadillo.h>
#include <random>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// Function to sample from a normal distribution
vec rnorm_cpp(int n, double mean, double sd) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> d(mean, sd);
  
  vec samples(n);
  for (int i = 0; i < n; ++i) {
    samples[i] = d(gen);
  }
  return samples;
}

// Function to sample from an inverse gamma distribution
double rinvgamma_cpp(double shape, double scale) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::gamma_distribution<> d(shape, 1.0 / scale);
  return 1.0 / d(gen);
}

// [[Rcpp::export]]
List gibbs_sampler_nested(const vec& y, const mat& W, uvec study_site, uvec family, int n_iter, const vec& beta_prior_mean, const vec& beta_prior_var, double sigma_ksi_prior_a, double sigma_ksi_prior_b, double sigma_theta_prior_a, double sigma_theta_prior_b, double sigma_prior_a, double sigma_prior_b) {
  
  // Adjust for 0-based indexing in C++
  study_site -= 1;
  family -= 1;
  
  int n = y.n_elem;
  int S_count = max(study_site) + 1;
  int F_count = max(family) + 1;
  int n_beta = W.n_cols;
  
  // Initialize parameters
  vec beta = zeros<vec>(n_beta);
  vec ksi = zeros<vec>(S_count);
  vec theta = zeros<vec>(F_count);
  double sigma_ksi2 = 1.0;
  vec sigma_theta2 = ones<vec>(S_count);
  double sigma2 = 1.0;
  
  // Storage for samples
  mat beta_samples(n_iter, n_beta, fill::zeros);
  mat ksi_samples(n_iter, S_count, fill::zeros);
  mat theta_samples(n_iter, F_count, fill::zeros);
  vec sigma_ksi2_samples(n_iter, fill::zeros);
  mat sigma_theta2_samples(n_iter, S_count, fill::zeros);
  vec sigma2_samples(n_iter, fill::zeros);
  
  for (int iter = 0; iter < n_iter; ++iter) {
    // Sample beta
    mat WtW = W.t() * W;
    mat V_beta = inv(WtW / sigma2 + diagmat(1.0 / beta_prior_var));
    vec m_beta = V_beta * (W.t() * (y - ksi(study_site) - theta(family)) / sigma2 + beta_prior_mean / beta_prior_var);
    beta = mvnrnd(m_beta, V_beta);
    
    // Sample ksi
    for (int s = 0; s < S_count; ++s) {
      uvec idx = find(study_site == s);
      vec y_s = y(idx);
      mat W_s = W.rows(idx);
      vec v_s = theta(family(idx));
      double V_ksi = 1.0 / (y_s.n_elem / sigma2 + 1.0 / sigma_ksi2);
      double m_ksi = V_ksi * sum(y_s - W_s * beta - v_s) / sigma2;
      ksi[s] = rnorm_cpp(1, m_ksi, sqrt(V_ksi))[0];
    }
    
    // Sample theta
    for (int f = 0; f < F_count; ++f) {
      uvec idx = find(family == f);
      vec y_f = y(idx);
      mat W_f = W.rows(idx);
      vec ksi_f = ksi(study_site(idx));
      int site = study_site(idx[0]);
      double V_theta = 1.0 / (y_f.n_elem / sigma2 + 1.0 / sigma_theta2[site]);
      double m_theta = V_theta * sum(y_f - W_f * beta - ksi_f) / sigma2;
      theta[f] = rnorm_cpp(1, m_theta, sqrt(V_theta))[0];
    }
    
    // Sample sigma_ksi2
    double a_sigma_ksi = sigma_ksi_prior_a + S_count / 2.0;
    double b_sigma_ksi = sigma_ksi_prior_b + sum(square(ksi)) / 2.0;
    sigma_ksi2 = rinvgamma_cpp(a_sigma_ksi, b_sigma_ksi);
    
    // Sample sigma_theta2 for each study site
    for (int s = 0; s < S_count; ++s) {
      uvec families_in_site = find(study_site == s);
      families_in_site = unique(family(families_in_site));
      if (!families_in_site.is_empty()) {
        for (size_t idx = 0; idx < families_in_site.n_elem; ++idx) {
          if (families_in_site[idx] < 0 || families_in_site[idx] >= F_count) {
            Rcpp::Rcerr << "Invalid family index: " << families_in_site[idx] << std::endl;
            Rcpp::stop("Index out of bounds");
          }
        }
        double a_sigma_theta = sigma_theta_prior_a + families_in_site.n_elem / 2.0;
        double b_sigma_theta = sigma_theta_prior_b + sum(square(theta(families_in_site))) / 2.0;
        sigma_theta2[s] = rinvgamma_cpp(a_sigma_theta, b_sigma_theta);
      }
    }
    
    // Sample sigma2
    double a_sigma = sigma_prior_a + n / 2.0;
    double b_sigma = sigma_prior_b + sum(square(y - W * beta - ksi(study_site) - theta(family))) / 2.0;
    sigma2 = rinvgamma_cpp(a_sigma, b_sigma);
    
    // Store samples
    beta_samples.row(iter) = beta.t();
    ksi_samples.row(iter) = ksi.t();
    theta_samples.row(iter) = theta.t();
    sigma_ksi2_samples[iter] = sigma_ksi2;
    sigma_theta2_samples.row(iter) = sigma_theta2.t();
    sigma2_samples[iter] = sigma2;
  }
  
  // Return samples as a list
  return List::create(
    Named("beta_samples") = beta_samples,
    Named("ksi_samples") = ksi_samples,
    Named("theta_samples") = theta_samples,
    Named("sigma_ksi2_samples") = sigma_ksi2_samples,
    Named("sigma_theta2_samples") = sigma_theta2_samples,
    Named("sigma2_samples") = sigma2_samples
  );
}


/*** R
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
library(scales)

# Source the C++ code
# sourceCpp("gibbs_sampler_nested.cpp")

# Combine samples into a 3D array for rstan::monitor
combine_samples_nested <- function(samples_list) {
  n_iter <- nrow(samples_list[[1]]$beta_samples)
  n_chains <- length(samples_list)
  n_params <- ncol(samples_list[[1]]$beta_samples) + 
    ncol(samples_list[[1]]$ksi_samples) + 
    ncol(samples_list[[1]]$theta_samples) + 
    2 + ncol(samples_list[[1]]$sigma_theta2_samples)  # beta, ksi, theta, sigma_ksi2, sigma2, sigma_theta2
  
  # Create an empty array
  combined_array <- array(NA, dim = c(n_iter, n_chains, n_params))
  
  # Fill in the array
  for (chain in 1:n_chains) {
    samples <- samples_list[[chain]]
    combined_array[, chain, 1:ncol(samples$beta_samples)] <- samples$beta_samples
    combined_array[, chain, (ncol(samples$beta_samples) + 1):(ncol(samples$beta_samples) + ncol(samples$ksi_samples))] <- samples$ksi_samples
    combined_array[, chain, (ncol(samples$beta_samples) + ncol(samples$ksi_samples) + 1):(ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples))] <- samples$theta_samples
    combined_array[, chain, ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + 1] <- samples$sigma_ksi2_samples
    combined_array[, chain, ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + 2] <- samples$sigma2_samples
    combined_array[, chain, (ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + 3):(ncol(samples$beta_samples) + ncol(samples$ksi_samples) + ncol(samples$theta_samples) + 2 + ncol(samples$sigma_theta2_samples))] <- samples$sigma_theta2_samples
  }
  
  return(combined_array)
}

# Example data
set.seed(123)
n_iter <- 5000
n_chains <- 4
n <- 30 * 30 * 3
S_count <- 30 # Number of Sites
F_count <- 30 * 30 # Number of Families in Total
study_site <- rep(1:S_count, each = n / S_count)
family <- rep(1:F_count, each = n / F_count)
W <- cbind(1, rnorm(n))
beta_true <- c(1, 2)
sigma_ksi_true <- 2
sigma_theta_true <- 1.5 # Assume each f:s variance parameter is same across sites
ksi_true <- rnorm(S_count, sd = sigma_ksi_true)
theta_true <- rnorm(F_count, sd = sigma_theta_true)
sigma_true <- 1
y <- W %*% beta_true + ksi_true[study_site] + theta_true[family] + rnorm(n, sd = sigma_true)

# Priors
priors <- list(beta_prior_mean = c(0, 0),
               beta_prior_var = c(100, 100),
               sigma_ksi_prior_a = 2,
               sigma_ksi_prior_b = 1,
               sigma_theta_prior_a = 2,
               sigma_theta_prior_b = 1,
               sigma_prior_a = 2,
               sigma_prior_b = 1)

# Run Gibbs sampler in parallel
seeds <- 1:n_chains
# Note, seed cannot be set in C++
samples_list <- mclapply(seeds, function(seed) { 
  gibbs_sampler_nested(y, W, study_site, family, n_iter, priors$beta_prior_mean, priors$beta_prior_var, 
                       priors$sigma_ksi_prior_a, priors$sigma_ksi_prior_b, 
                       priors$sigma_theta_prior_a, priors$sigma_theta_prior_b, 
                       priors$sigma_prior_a, priors$sigma_prior_b)
}, mc.cores = n_chains)

# Combine samples
combined_samples <- combine_samples_nested(samples_list)

# Assign parameter names
dimnames(combined_samples) <- list(
  iterations = NULL,
  chains = NULL,
  parameters = c(paste0("beta_", 1:ncol(samples_list[[1]]$beta_samples)),
                 paste0("ksi_", 1:ncol(samples_list[[1]]$ksi_samples)),
                 paste0("theta_", 1:ncol(samples_list[[1]]$theta_samples)),
                 "sigma_ksi2", "sigma2",
                 paste0("sigma_theta2_", 1:ncol(samples_list[[1]]$sigma_theta2_samples)))
)

# Use rstan::monitor
mcmc_summary <- monitor(combined_samples)

# Extract summary statistics
mean_values <- round(mcmc_summary$mean, 4)
median_values <- round(mcmc_summary$`50%`, 4)
lower_bounds <- round(mcmc_summary$`2.5%`, 4)
upper_bounds <- round(mcmc_summary$`97.5%`, 4)

# Combine the true values into a single vector, ordered according to the MCMC parameters
true_values <- c(beta_true, ksi_true, theta_true, sigma_ksi_true^2, sigma_true^2, rep(sigma_theta_true^2, S_count))

# Parameter names
param_names <- c(paste0("beta_", 1:ncol(samples_list[[1]]$beta_samples)),
                 paste0("ksi_", 1:ncol(samples_list[[1]]$ksi_samples)),
                 paste0("theta_", 1:ncol(samples_list[[1]]$theta_samples)),
                 "sigma_ksi2", "sigma2",
                 paste0("sigma_theta2_", 1:ncol(samples_list[[1]]$sigma_theta2_samples)))

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

# Filter for random intercept parameters (those starting with "ksi_")
random_intercept_comparison <- comparison %>% filter(grepl("^ksi_|^theta_", param_name))

# % of random intercept parameters within credible interval
cat("% of random intercept parameters within credible interval: ", mean(random_intercept_comparison$within_credible_interval), "\n")

# % of random intercept parameters with correct sign
random_intercept_comparison$correct_sign <- sign(random_intercept_comparison$mean) == sign(random_intercept_comparison$true_value)
cat("% of random intercept parameters with correct sign: ", mean(random_intercept_comparison$correct_sign), "\n")

# Print the result to check random intercept parameters
print(random_intercept_comparison)

print("Variance Parameter Estimation vs. Truth:")

# Filter the comparison dataframe for variance parameters
variance_comparison <- comparison %>% filter(grepl("^sigma2$|^sigma_ksi2$|^sigma_theta2_", param_name))

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
  paste0("ksi_", 1),
  paste0("theta_", 1),
  "sigma_ksi2",
  "sigma2",
  paste0("sigma_theta2_", 1:2)
)

# True values for the parameters
true_values <- c(
  beta_true, 
  ksi_true[1], 
  theta_true[1], 
  sigma_ksi_true^2, 
  sigma_true^2, 
  rep(sigma_theta_true^2, 2)
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
u_limits <- range(combined_samples[(n_iter / 2 + 1):n_iter, , grep("^ksi_", dimnames(combined_samples)[[3]])])
v_limits <- range(combined_samples[(n_iter / 2 + 1):n_iter, , grep("^theta_", dimnames(combined_samples)[[3]])])
sigma_ksi2_limits <- range(combined_samples[(n_iter / 2 + 1):n_iter, , grep("^sigma_ksi2", dimnames(combined_samples)[[3]])])
sigma2_limits <- range(combined_samples[(n_iter / 2 + 1):n_iter, , grep("^sigma2", dimnames(combined_samples)[[3]])])
sigma_theta2_limits <- c(0, 25)  # Truncate to a maximum of 25

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
  } else if (grepl("^ksi_", param)) {
    x_limits <- u_limits
  } else if (grepl("^theta_", param)) {
    x_limits <- v_limits
  } else if (grepl("^sigma_ksi2", param)) {
    x_limits <- sigma_ksi2_limits
  } else if (grepl("^sigma2", param)) {
    x_limits <- sigma2_limits
  } else if (grepl("^sigma_theta2_", param)) {
    x_limits <- sigma_theta2_limits
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

# Observation: Fixed effect intercept & Site-level random intercepts 
# seem to take longer to converge. This is likely attributable to 
# anti-correlation between the Fixed effect intercept & 
# Site-level random intercepts, which is apparent in the traceplots. 
*/

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
List gibbs_sampler_nested(vec y, mat W, arma::mat Z_family, arma::mat Z_site, 
                          arma::mat Z_family_to_site, int n_iter, double mu_prior_var, 
                          const vec& mu_beta, const vec& beta_prior_var, 
                          double sigma_ksi_prior_a, double sigma_ksi_prior_b, 
                          double sigma_theta_prior_a, double sigma_theta_prior_b, 
                          double sigma_prior_a, double sigma_prior_b, 
                          int r, mat U, vec alpha) {
  
  // Set GSL random number generator seed
  // long seed=1;
  // gsl_rng * rr = gsl_rng_alloc (gsl_rng_rand48);
  // gsl_rng_set (rr, seed);
  
  // Calculates the number of clusters
  int N_sites = Z_site.n_cols;
  int N_families = Z_family.n_cols;
  int N_obs = y.n_elem;
  int n_beta = W.n_cols;
  vec y_tilde(N_obs);
  
  //
  //arma::mat Sigma_0_inv = inv(U*eye(r, r)*trans(U) + eye(N_obs, N_obs));
  arma::mat D = eye(r, r);
  arma::mat Sigma_0_inv = eye(N_obs, N_obs) - U*inv(inv(D) + U.t()*U)*U.t();
  
  // Initialize parameters
  double mu = as_scalar(rnorm_cpp(1, 0, std::sqrt(mu_prior_var)));
  double sigma2_ksi = rinvgamma_cpp(sigma_ksi_prior_a, sigma_ksi_prior_b);
  vec sigma2_theta = zeros(N_sites);
  
  //double sigma2 = rinvgamma_cpp(sigma_prior_a, sigma_prior_b);
  
  // vec beta = mvnrnd(mu_beta, diagmat(beta_prior_var));
  
  vec ksi = rnorm_cpp(N_sites, 0, std::sqrt(sigma2_ksi));
  
  vec theta = zeros(N_families);  // Changed to initialize theta to zero
  for(int s = 0; s < N_sites; s++) {
    sigma2_theta(s) = rinvgamma_cpp(sigma_theta_prior_a, sigma_theta_prior_b);
    arma::uvec indices_family_in_site = arma::find(Z_family_to_site.col(s) != 0); // Get indices of families that belong to site s
    for (int f : indices_family_in_site) {
      theta(f) = as_scalar(rnorm_cpp(1, ksi(s), std::sqrt(sigma2_theta(s))));
    }
  }
  
  y_tilde = y - Z_family*theta - U*alpha;
  vec beta = inv(trans(W)*W)*trans(W)*y_tilde;
  double sigma2 = arma::dot(trans(y-W*beta), (y-W*beta))/(N_obs - W.n_cols);
    
  // Store initial values
  double mu_init = mu;
  vec beta_init = beta;
  vec ksi_init = ksi;
  vec theta_init = theta;
  double sigma2_ksi_init = sigma2_ksi;
  vec sigma2_theta_init = sigma2_theta;
  double sigma2_init = sigma2;  
  
  // Storage for samples
  vec mu_samples(n_iter, fill::zeros);
  mat beta_samples(n_iter, n_beta, fill::zeros);
  mat ksi_samples(n_iter, N_sites, fill::zeros);
  mat theta_samples(n_iter, N_families, fill::zeros);
  vec sigma2_ksi_samples(n_iter, fill::zeros);
  mat sigma2_theta_samples(n_iter, N_sites, fill::zeros);
  vec sigma2_samples(n_iter, fill::zeros);
  
  for (int iter = 0; iter < n_iter; iter++) {
    // Sample beta
    y_tilde = y - Z_family*theta - U*alpha; // R_beta
    mat V_beta = inv(trans(W)*W/sigma2 + inv(diagmat(beta_prior_var)));
    vec m_beta = V_beta * (W.t()*y_tilde/sigma2 + inv(diagmat(beta_prior_var))*mu_beta);
    
    beta = mvnrnd(m_beta, V_beta);
    
    // Sample random effects
    y_tilde = y - W*beta - U*alpha; // R_theta
    for (int s = 0; s < N_sites; s++) {
      arma::uvec families_in_s = arma::find(Z_family_to_site.col(s) != 0); // Get indices of families that belong
      
      // Sample family-level intercepts theta
      for (int f : families_in_s) {
        arma::uvec individuals_in_f = arma::find(Z_family.col(f) == 1); // Get indices of observations that belong
        double sum_y_tilde = sum(y_tilde.elem(individuals_in_f));
        int n_sf = individuals_in_f.n_elem;
        double V_theta = 1.0 / (n_sf/sigma2 + 1.0/sigma2_theta(s));
        double m_theta = V_theta*(sum_y_tilde/sigma2 + ksi(s)/sigma2_theta(s));  // Corrected the mean calculation
        theta(f) = rnorm_cpp(1, m_theta, sqrt(V_theta))[0];
      }
      
      // Sample site-level intercepts ksi
      double sum_theta = sum(theta.elem(families_in_s));
      int n_s = families_in_s.n_elem;
      double V_ksi = 1.0 / (n_s/sigma2_theta(s) + 1.0/sigma2_ksi);
      double m_ksi = V_ksi*(mu + sum_theta/sigma2_theta(s));
      ksi(s) = as_scalar(rnorm_cpp(1, m_ksi, sqrt(V_ksi)));
      
      // Sample sigma2_theta for s-th site
      double a_sigma_theta = sigma_theta_prior_a + n_s/2.0;
      double b_sigma_theta = sigma_theta_prior_b + sum(square(theta.elem(families_in_s) - ksi(s)*ones(n_s))) / 2.0;
      sigma2_theta[s] = rinvgamma_cpp(a_sigma_theta, b_sigma_theta);    
      
    }
    
    // Sample overall mean mu
    double V_mu = 1.0/(1.0/mu_prior_var + N_sites/sigma2_ksi);
    double m_mu = V_mu*(sum(ksi)/sigma2_ksi);
    mu = as_scalar(rnorm_cpp(1, m_mu, std::sqrt(V_mu)));
    
    // Sample sigma2_ksi
    double a_sigma_ksi = sigma_ksi_prior_a + N_sites/2.0;
    double b_sigma_ksi = sigma_ksi_prior_b + sum(square(ksi - mu))/2.0;
    sigma2_ksi = rinvgamma_cpp(a_sigma_ksi, b_sigma_ksi);    
    
    // Sample sigma
    y_tilde = y - W*beta - Z_family*theta; // R_sigma
    double a_sigma = sigma_prior_a + N_obs/2.0;
    double b_sigma = sigma_prior_b + arma::dot(y_tilde, Sigma_0_inv*y_tilde)/2.0;
    sigma2 = rinvgamma_cpp(a_sigma, b_sigma);
    
    // Store samples
    mu_samples(iter) = mu;
    beta_samples.row(iter) = trans(beta);
    ksi_samples.row(iter) = ksi.t();
    theta_samples.row(iter) = theta.t();
    sigma2_ksi_samples(iter) = sigma2_ksi;
    sigma2_theta_samples.row(iter) = sigma2_theta.t();
    sigma2_samples(iter) = sigma2;
  }
  
  // Return samples as a list
  return List::create(
    Named("mu_samples") = mu_samples,
    Named("beta_samples") = beta_samples,
    Named("ksi_samples") = ksi_samples,
    Named("theta_samples") = theta_samples,
    Named("sigma2_ksi_samples") = sigma2_ksi_samples,
    Named("sigma2_theta_samples") = sigma2_theta_samples,
    Named("sigma2_samples") = sigma2_samples,
    Named("initial_values") = List::create(
      Named("mu_init") = mu_init,
      Named("beta_init") = beta_init,
      Named("ksi_init") = ksi_init,
      Named("theta_init") = theta_init,
      Named("sigma2_ksi_init") = sigma2_ksi_init,
      Named("sigma2_theta_init") = sigma2_theta_init,
      Named("sigma2_init") = sigma2_init
    )
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
library(lme4)
library(coda)

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
simulate_omics_data <- function(n_views=2, N_obs=200, p_m=10, r=4,
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
  U <- matrix(data = rnorm(N_obs*r), nrow = N_obs, ncol = r)
  A_list <- list()
  E_list <- list()
  
  for (m in 1:n_views) {
    A_list[[m]] <- simulate_A(r, p_m, n_important_components, n_important_features)
    E_list[[m]] <- matrix(data = rnorm(N_obs*p_m, sd = sqrt(sigma2)), nrow = N_obs, ncol = p_m)
    X_list[[m]] <- U %*% A_list[[m]] + E_list[[m]]
  }
  
  omics_results <- list(X=X_list, U=U, A=A_list,
                        index_important_components=index_important_components,
                        index_important_features=index_important_features,
                        gamma=gamma, Eta=Eta)
  
  return(omics_results)
}

simulate_re_data_nested <- function(n_views=2, p_m=10, r=4,
                                    prob_feature_importance=0.5, 
                                    prob_component_importance=0.5,
                                    sigma2_ksi=1, sigma2_theta=rep(1, 5),
                                    N_sites=5, n_families_per_site=3, 
                                    n_individs_per_family = 2,
                                    n_covars=1, mu=1, sigma2=1, seed=1,
                                    balanced=T) {
  
  # Define default arguments in a list
  default_args <- list(
    n_views=2, 
    p_m=10, 
    r=4,
    prob_feature_importance=0.5, 
    prob_component_importance=0.5,
    sigma2_ksi=1, 
    sigma2_theta=rep(1, 5),
    N_sites=5, 
    n_families_per_site=3, 
    n_individs_per_family = 2,
    n_covars=1, 
    mu=1, 
    sigma2=1, 
    seed=1
  )
  
  # Import arguments into the global environment
  #list2env(default_args, envir = .GlobalEnv)

  # Outcome model
  set.seed(seed)
  
  if(length(sigma2_theta) != N_sites) {
    sigma2_theta <- rep(1, N_sites)
  }
  
  if(balanced == T) {
    N_families <- N_sites*n_families_per_site
    N_obs <- N_sites*n_families_per_site*n_individs_per_family
    
    # Specify design matrix for sites
    Z_site <- kronecker(diag(N_sites), rep(1, n_families_per_site*n_individs_per_family))
    
    # Specify design matrix for families nested within sites
    Z_family <- kronecker(diag(N_families), rep(1, n_individs_per_family))
  }
  
  # Mapping of families to sites
  Z_family_to_site <- t(Z_family) %*% Z_site
  
  omics_data <- simulate_omics_data(n_views, N_obs, p_m, r, prob_feature_importance, prob_component_importance, sigma2)
  U <- omics_data$U
  
  # Sample latent factor loadings
  alpha <- matrix(0, nrow = r, ncol = 1)
  #alpha[omics_data$index_important_components, ] <- rnorm(length(omics_data$index_important_components))
  alpha <- rnorm(r, 0, sqrt(sigma2))
  
  # Simulate ksi_s ~ N(0, sigma2_ksi)
  ksi <- rnorm(N_sites, mu, sd = sqrt(sigma2_ksi)) %>% matrix(ncol = 1)
  
  # Sample theta_sf|ksi_s ~ N(ksi_s, sigma2_theta_s)
  theta <- matrix(0, nrow = N_sites, ncol = n_families_per_site)
  for (s in 1:N_sites) {
    theta[s, ] <- rnorm(n_families_per_site, mean = ksi[s], sd = sqrt(sigma2_theta[s]))
  }
  
  # Family effects as a single vector
  theta <- as.vector(t(theta)) %>% matrix(ncol = 1)
  
  W <- NULL
  
  if(n_covars > 0) {
    for(k in 1:n_covars) {
      W <- cbind(W, rnorm(N_obs))
    }
  }
  
  beta <- matrix(rnorm(n_covars), ncol = 1)
  
  # Sample residuals
  epsilon <- matrix(rnorm(N_obs, sd = sqrt(sigma2)), nrow = N_obs)
  
  # Combine effects
  Y <- W%*%beta + Z_family%*%theta + U%*%alpha + epsilon # Add omics data
  
  return(list(Y=Y, Z_site=Z_site, Z_family=Z_family, Z_family_to_site=Z_family_to_site, ksi=ksi, theta=theta, 
              X=omics_data$X, U=omics_data$U, A=omics_data$A, mu=mu, alpha=alpha, W=W, beta=beta,
              gamma=omics_data$gamma, Eta=omics_data$Eta, sigma2 = sigma2,
              nu2 = list(sigma2_ksi=sigma2_ksi_true, sigma2_theta=sigma2_theta_true)))
}

# Combine samples into a 3D array for rstan::monitor
combine_samples_nested <- function(samples_list, n_iter, n_chains) {
  n_params <- ncol(samples_list[[1]]$mu_samples) + 
    ncol(samples_list[[1]]$beta_samples) +  # FIGURE OUT WHAT TO DO ABOUT THIS IF THERE'S NO FIXED EFFECTS
    ncol(samples_list[[1]]$ksi_samples) + 
    ncol(samples_list[[1]]$theta_samples) + 
    ncol(samples_list[[1]]$sigma2_ksi_samples) + 
    ncol(samples_list[[1]]$sigma2_samples) + 
    ncol(samples_list[[1]]$sigma2_theta_samples)  # beta, ksi, theta, sigma2_ksi, sigma2, sigma2_theta
  
  # Create an empty array
  combined_array <- array(NA, dim = c(n_iter, n_chains, n_params))
  
  # Fill in the array
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

n_covars <- 2
N_sites <- 30
n_families_per_site <- 30
n_individs_per_family <- 10

# Priors
priors <- list(mu_prior_var = 100,
               mu_beta = rep(0, n_covars),
               beta_prior_var = rep(100, n_covars),
               sigma_ksi_prior_a = 3,
               sigma_ksi_prior_b = 2,
               sigma_theta_prior_a = 3,
               sigma_theta_prior_b = 2,
               sigma_prior_a = 3,
               sigma_prior_b = 2)

#sigma2_ksi_true <- 1
#sigma2_theta_true <- rep(1, N_sites)
#sigma2_true <- 1

sigma2_ksi_true <- 1/rgamma(1, priors$sigma_ksi_prior_a, priors$sigma_ksi_prior_b)
sigma2_theta_true <- 1/rgamma(N_sites, priors$sigma_theta_prior_a, priors$sigma_theta_prior_b)
sigma2_true <- 1/rgamma(1, priors$sigma_prior_a, priors$sigma_prior_b)
#sigma2_true <- 1

r <- 4

simulation_results <- simulate_re_data_nested(N_sites = N_sites, 
                                              n_families_per_site = n_families_per_site,
                                              n_individs_per_family = n_individs_per_family,
                                              sigma2_ksi = sigma2_ksi_true, 
                                              sigma2_theta = sigma2_theta_true,
                                              sigma2 = sigma2_true,
                                              n_covars = n_covars, n_views = 1)
# dataList <- list(simulation_results$X[[1]],
#                  simulation_results$X[[2]],
#                  simulation_results$Y)

mu_true <- simulation_results$mu
ksi_true <- simulation_results$ksi
theta_true <- simulation_results$theta
beta_true <- simulation_results$beta
U_true <- simulation_results$U
alpha_true <- simulation_results$alpha
W <- simulation_results$W

W[,-1] <- scale(W[,-1])

Z_site <- simulation_results$Z_site
Z_family <- simulation_results$Z_family
Z_family_to_site <- simulation_results$Z_family_to_site
y <- simulation_results$Y 

N_sites <- ncol(Z_site)
N_families <- ncol(Z_family)
N_obs <- nrow(Z_family)

n_chains <- 1
n_iter <- 5000
n_burnin <- floor(n_iter*0.5)

RE_df <- data.frame(y = y,
                    site = which(Z_site == 1, arr.ind = T)[,2],
                    family = which(Z_family == 1, arr.ind = T)[,2])
# Source the C++ code
# sourceCpp("gibbs_sampler_nested.cpp")

# Run Gibbs sampler in parallel
seeds <- 1:n_chains

start_time <- Sys.time()

# Note, seed cannot be set in C++
samples_list <- mclapply(seeds, function(seed) { 
  gibbs_sampler_nested(y, W, Z_family, Z_site, Z_family_to_site, n_iter, priors$mu_prior_var, priors$mu_beta, priors$beta_prior_var, 
                       priors$sigma_ksi_prior_a, priors$sigma_ksi_prior_b, 
                       priors$sigma_theta_prior_a, priors$sigma_theta_prior_b, 
                       priors$sigma_prior_a, priors$sigma_prior_b,
                       r, U_true, alpha_true)
}, mc.cores = n_chains)

end_time <- Sys.time()

# Extract initial values
init_values <- samples_list[[1]]$initial_values
mu_init <- init_values$mu_init
beta_init <- init_values$beta_init
ksi_init <- init_values$ksi_init
theta_init <- init_values$theta_init
sigma2_ksi_init <- init_values$sigma2_ksi_init
sigma2_theta_init <- init_values$sigma2_theta_init
sigma2_init <- init_values$sigma2_init

# Combine samples
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

# Use rstan::monitor
mcmc_summary <- monitor(combined_samples)


# Extract summary statistics
mean_values <- round(mcmc_summary$mean, 4)
median_values <- round(mcmc_summary$`50%`, 4)
lower_bounds <- round(mcmc_summary$`2.5%`, 4)
upper_bounds <- round(mcmc_summary$`97.5%`, 4)

# Combine the true values into a single vector, ordered according to the MCMC parameters
true_values <- c(mu_true, beta_true, ksi_true, theta_true, sigma2_ksi_true, sigma2_theta_true, sigma2_true)


# Initial values
initial_values <- c(mu_init, beta_init, ksi_init, theta_init, 
                    sigma2_ksi_init, sigma2_theta_init, sigma2_init)

# Initialize a dataframe to hold the comparisons
comparison <- data.frame(
  param_name = param_names,
  lower = lower_bounds,
  mean = mean_values,
  median = median_values,
  upper = upper_bounds,
  true_value = round(true_values, 4),
  initial = round(initial_values, 4)
)

# Add a logical vector to see if true values are within the credible intervals
comparison$within_credible_interval <- with(comparison, true_value >= lower & true_value <= upper)

print("Fixed Effect Estimation vs. Truth:")

# Filter for fixed effect parameters
fixed_effect_comparison <- comparison %>% filter(grepl("^mu|^beta_", param_name))

# Print the result to check fixed effect parameters
print(fixed_effect_comparison)

print("Random Intercept Estimation vs. Truth:")

# Filter for random intercept parameters (those starting with "ksi_")
random_intercept_comparison <- comparison %>% filter(grepl("^ksi_|^theta_", param_name))

# % of random intercept parameters with correct sign
random_intercept_comparison$correct_sign <- sign(random_intercept_comparison$mean) == sign(random_intercept_comparison$true_value)

# Print the result to check random intercept parameters
print(random_intercept_comparison)

print("Variance Parameter Estimation vs. Truth:")

# Filter the comparison dataframe for variance parameters
variance_comparison <- comparison %>% filter(grepl("^sigma2$|^sigma2_ksi$|^sigma2_theta_", param_name))

# Print the filtered result to check each variance parameter
print(variance_comparison)

# Set the number of chains to plot
n_chains_to_plot <- min(n_chains, 1)

################################################################################
parameter_groups <- list(
  mu = "mu",
  FE = c(paste0("beta_", 1:ncol(samples_list[[1]]$beta_samples))),
  site_intercepts = c(paste0("ksi_", 1:N_sites)),
  family_intercepts = c(paste0("theta_", 1:N_families)),
  site_variance = "sigma2_ksi",
  family_variance = c(paste0("sigma2_theta_", 1:N_sites)),
  error_variance = "sigma2"
)

parameter_group_true_values <- list(
  mu = mu_true,
  FE = beta_true, 
  site_intercepts = ksi_true, 
  family_intercepts = theta_true, 
  site_variance = sigma2_ksi_true, 
  family_variance = sigma2_theta_true,
  error_variance = sigma2_true
)

CI_list <- list(
  mu = comparison %>% filter(grepl("^mu", param_name)) %>% select(lower, upper), 
  FE = comparison %>% filter(grepl("^beta_", param_name)) %>% select(lower, upper), 
  site_intercepts = comparison %>% filter(grepl("^ksi_", param_name)) %>% select(lower, upper), 
  family_intercepts = comparison %>% filter(grepl("^theta_", param_name)) %>% select(lower, upper), 
  site_variance = comparison %>% filter(grepl("^sigma2_ksi", param_name)) %>% select(lower, upper), 
  family_variance = comparison %>% filter(grepl("^sigma2_theta_", param_name)) %>% select(lower, upper),
  error_variance = comparison %>% filter(grepl("^sigma2$", param_name)) %>% select(lower, upper)
)

################################################################################

# Determine the y-axis limits to be (0, 1) for all density plots
y_limits <- c(0, 1)

for(g in seq_along(parameter_groups)) {
  number_displayed <- 8

  if(length(parameter_groups[[g]]) > number_displayed) {
    display_index <- sort(sample(seq_along(parameter_groups[[g]]), number_displayed))
  } else {
    display_index <- seq_along(parameter_groups[[g]])
  }

  # Determine x-axis limits based on parameter class
  x_limits <- combined_samples[, chain, which(dimnames(combined_samples)[[3]] %in% parameter_groups[[g]][display_index])]
  
  # Prepare a list to store ggplot objects
  plot_list <- list()

  # Create trace plots for each parameter in each chain
  for (i in display_index) {
    param <- parameter_groups[[g]][i]
    true_value <- parameter_group_true_values[[g]][i]
    initial <- init_values[[g]][i]
    CI <- as.data.frame(CI_list[[g]])[i,]

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
        geom_hline(yintercept = true_value, linetype = "dashed", color = "black",
                   linewidth = (1 + 1/(8*number_displayed))) +
        geom_vline(xintercept = n_burnin, linetype = "dashed", color = "darkred",
                   linewidth = (1 + 1/(8*number_displayed))) +
        annotate("segment", x = 0, xend = 0,
                 y = CI$lower, yend = CI$upper,
                 color = "black", linewidth = (1 + 1/(8*number_displayed)),
                 arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm"))) +
        labs(x = ifelse(param == tail(parameter_groups[[g]], 1), "Iteration", ""), y = ifelse(chain == 1, bquote(.(param)), "")) +
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

  # Arrange the trace plots into a grid
  grid_plots <- do.call(grid.arrange, c(plot_list, ncol = n_chains_to_plot))

  # Create density plots for each parameter in each chain after burn-in
  for (i in display_index) {
    param <- parameter_groups[[g]][i]
    true_value <- parameter_group_true_values[[g]][i]
    CI <- as.data.frame(CI_list[[g]])[i,]

    # Determine x-axis limits based on parameter class
    if (grepl("^mu", param)) {
      x_limits <- mu_limits
    } else if (grepl("^beta_", param)) {
      x_limits <- beta_limits
    } else if (grepl("^ksi_", param)) {
      x_limits <- ksi_limits
    } else if (grepl("^theta_", param)) {
      x_limits <- theta_limits
    } else if (grepl("^sigma2_ksi", param)) {
      x_limits <- sigma2_ksi_limits
    } else if (grepl("^sigma2$", param)) {
      x_limits <- sigma2_limits
    } else if (grepl("^sigma2_theta_", param)) {
      x_limits <- sigma2_theta_limits
    }

    for (chain in 1:n_chains_to_plot) {
      # Extract the values for the parameter and chain after burn-in
      param_values <- combined_samples[(n_burnin + 1):n_iter, chain, which(dimnames(combined_samples)[[3]] == param)]

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
        geom_line(aes(color = chain), linewidth = 1) +
        geom_vline(xintercept = true_value, linetype = "dashed", color = "black", linewidth = 1) +
        annotate("segment", x = CI$lower, xend = CI$upper,
                 y = 0.5, yend = 0.5,
                 color = "black", linewidth = (1 + 1/(8*number_displayed)),
                 arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm"))) +
        scale_x_continuous(limits = x_limits, labels = label_number(accuracy = 0.01)) +
        scale_y_continuous(limits = y_limits, labels = label_number(accuracy = 0.01)) +
        labs(x = ifelse(param == tail(parameter_groups[[g]], 1), "Value", ""), y = ifelse(chain == 1, param, "")) +
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

}

# Prepare a list to store ggplot objects
plot_list <- list()

# Traceplot of joint ksi prior
library(mvtnorm)
#ldmvnorm()

RE_df$w1 <- W[,1]
RE_df$w2 <- W[,2]

summary(lmer(y~(1|site) + (1|family:site) + w1 + w2, data = RE_df))

# Observation: Fixed effect intercept & Site-level random intercepts 
# seem to take longer to converge. This is likely attributable to 
# anti-correlation between the Fixed effect intercept & 
# Site-level random intercepts, which is apparent in the traceplots. 

sprintf("BIP duration: %f", end_time-start_time)
cat("% of all parameters within credible interval: ", mean(comparison$within_credible_interval), "\n")
cat("% of random intercept parameters within credible interval: ", mean(random_intercept_comparison$within_credible_interval), "\n")
cat("% of random intercept parameters with correct sign: ", mean(random_intercept_comparison$correct_sign), "\n")
cat("% of variance parameters within credible interval: ", mean(variance_comparison$within_credible_interval), "\n")

simulation_settings <- data.frame(
  variable = c("# of Sites", "# of Families", "# of Observations", "# of Iterations"),
  value = c(N_sites, N_families, N_obs, n_iter)
)

print(simulation_settings)
*/

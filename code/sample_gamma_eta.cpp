#include <cstdlib>
#include <iostream>
#include <string>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

#include <utils.h>

// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::export]]
vec initialize_gamma(int r=4, double prior_component_selection=0.5) {
  //vec gamma = rbinom(r, 1, prior_component_selection);
  vec gamma = randbin(r, 1, prior_component_selection);
  return gamma;
}

// [[Rcpp::export]]
mat initialize_Eta(vec gamma, int r=4, int p_m=10, double prior_feature_selection=0.5) {
  mat Eta(r, p_m, fill::value(datum::nan));
  
  for(int l = 0; l < r; l++) {
    if(gamma(l) == 1) {
      vec Eta_l = randbin(p_m, 1, prior_feature_selection);
      Eta.row(l) = Eta_l.t();
    } else {
      Eta.row(l).zeros();
    }
  }
  
  return Eta;
}

// [[Rcpp::export]]
mat calculate_mvnorm_var_j(vec eta_j, double sigma2_j, vec tau2_j, mat U, int n_obs) {
  // Where's tau2_j used?
  int n_components_active_in = sum(eta_j);
  mat Sigma2_j(n_obs, n_obs, fill::value(datum::nan));
  
  if (n_components_active_in > 0) {
    uvec components_active_in = get_indices_of_ones(eta_j);
    Sigma2_j = U.cols(components_active_in) *
      eye(n_components_active_in, n_components_active_in) *
      U.cols(components_active_in).t() +
      eye(n_obs, n_obs);
  } else {
    Sigma2_j = eye(n_obs, n_obs);
  }
  mat mvnorm_var = sigma2_j * Sigma2_j;
  return mvnorm_var;
}

// [[Rcpp::export]]
double calculate_log_dmvnorm_j(vec x_j, vec mu, mat U, vec eta_j, vec tau2_j, double sigma2_j, int n_obs) {
  int n_active_components = sum(eta_j);
  double Sigma_det = 0;
  mat Sigma_inv(n_obs, n_obs, fill::value(datum::nan));
  
  if (n_active_components > 0) {
    uvec active_components = get_indices_of_ones(eta_j);
    
    mat D(n_active_components, n_active_components, fill::value(0));
    D.diag() = tau2_j.elem(active_components);
    U = U.cols(active_components);
    
    Sigma_det = pow(sigma2_j, n_obs)*det(inv(D) + U.t()*U)*det(D);
    Sigma_inv = (1.0/sigma2_j)*get_woodbury_inv(U, D);
  } else {
    // Sigma2_j = I, an n-by-n identity matrix
    Sigma_det = pow(sigma2_j, n_obs);
    Sigma_inv = (1.0/sigma2_j)*eye(n_obs, n_obs);
  }
  
  double log_dens = -0.5*n_obs*log(2*datum::pi) - 0.5*log(Sigma_det) - 0.5*as_scalar(trans(x_j-mu)*Sigma_inv*(x_j-mu));
  
  return log_dens;
}

// [[Rcpp::export]]
double calculate_log_G_j(vec mu, vec gamma, vec eta_j, double sigma2_j, vec tau2_j, mat U, int n_obs, vec x_j, double prob_feature_selection) {
  double log_dmvnorm_j = calculate_log_dmvnorm_j(x_j, mu, U, eta_j, tau2_j, sigma2_j, n_obs);
  uvec active_gamma = get_indices_of_ones(gamma);
  int n_active_gamma = active_gamma.n_elem; // TODO account for all inactive
  int n_active_eta_given_gamma_1 = sum(eta_j.elem(active_gamma));
  double log_G_j = log_dmvnorm_j + 
    n_active_eta_given_gamma_1 * log(prob_feature_selection) + 
    (n_active_gamma - n_active_eta_given_gamma_1) * log(1 - prob_feature_selection);
  return log_G_j;
}

// [[Rcpp::export]]
vec calculate_log_PQ_lj(int l, vec mu, vec gamma, vec eta_j, double sigma2_j, vec tau2_j, mat U, int n_obs, vec x_j, double prob_feature_selection) {
  vec log_PQ(2);
  vec gamma_1 = gamma;
  vec eta_1 = eta_j;
  vec eta_0 = eta_j;
  gamma_1(l) = 1;
  eta_1(l) = 1;
  eta_0(l) = 0;
  double log_G_j_1 = calculate_log_G_j(mu, gamma_1, eta_1, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection);
  double log_G_j_0 = calculate_log_G_j(mu, gamma_1, eta_0, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection);
  
  // Use logsumexp trick
  double max_arg = std::max(log_G_j_1, log_G_j_0); 
  double a = log_G_j_1 - max_arg;
  double b = log_G_j_0 - max_arg;

  double log_P_lj = a - log(exp(a)+exp(b)); 
  double log_Q_lj = b - log(exp(a)+exp(b));
  log_PQ(0) = log_P_lj;
  log_PQ(1) = log_Q_lj;
  
  return log_PQ;
}

// [[Rcpp::export]]
double calculate_eta_lj_threshold(int l, vec mu, vec gamma, vec eta_j, double sigma2_j, vec tau2_j, mat U, int n_obs, vec x_j, double prob_feature_selection) {
  vec gamma_1 = gamma;
  vec eta_1 = eta_j;
  vec eta_0 = eta_j;
  gamma_1(l) = 1;
  eta_1(l) = 1;
  eta_0(l) = 0;
  
  double log_G_j_1 = calculate_log_G_j(mu, gamma_1, eta_1, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection);
  double log_G_j_0 = calculate_log_G_j(mu, gamma_1, eta_0, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection);
  return log_G_j_1 - log_G_j_0;
}

// [[Rcpp::export]]
double log_target_density_l(int l, vec mu, vec gamma, mat Eta, vec sigma2, mat tau2, mat U, int n_obs, int p_m, 
                            mat X, double prob_component_selection, double prob_feature_selection) {
  double sum_log_target_density_lj = 0;
  for (int j=0; j<p_m; j++) {
    // TODO understand if log_dmvnorm_j should only use data within. 
    // Note, this should cancel in the log_acceptance_ratio subtraction.
    double log_dmvnorm_j = calculate_log_dmvnorm_j(X.col(j), mu, U, Eta.col(j), tau2.col(j), sigma2(j), n_obs);
    sum_log_target_density_lj = sum_log_target_density_lj + log_dmvnorm_j + 
      Eta(l,j)*log(prob_feature_selection) + (1-Eta(l,j))*log(1-prob_feature_selection);
  }
  double log_target_density_l = gamma(l)*log(prob_component_selection) + 
    (1-gamma(l))*log(1-prob_component_selection) + sum_log_target_density_lj;
  return log_target_density_l;
}

// [[Rcpp::export]]
double log_proposal_l_density(int l, vec mu, vec gamma_prime, mat Eta_prime, 
                              vec gamma, mat Eta, vec sigma2, mat tau2, mat U, 
                              int n_obs, int p_m, mat X, double prob_feature_selection) {
  if ((gamma(l)==1) & (gamma_prime(l)==0) & (sum(Eta_prime.row(l))==0)) {
    return 0;
  } else if ((gamma(l)==0) & (gamma_prime(l)==1)) {
    double log_proposal_l_density = 0;
    for (int j=0; j<p_m; j++) {
      vec PQ_lj = calculate_log_PQ_lj(l, mu, gamma_prime, Eta_prime.col(j), sigma2(j), 
                                      tau2.col(j), U, n_obs, X.col(j), prob_feature_selection);
      log_proposal_l_density = log_proposal_l_density + Eta_prime(l,j) * PQ_lj(0) + (1-Eta_prime(l,j)) * PQ_lj(1);
    } 
    return log_proposal_l_density;
  } else { 
    throw std::invalid_argument("Error in log_proposal_l_density evaluation.");
    return datum::nan;
  }
} 

// Main function

// [[Rcpp::export]]
Rcpp::List main_sample_gamma_Eta(int n_iterations, int n_burnin,
                          int r, int n_obs, int p_m, 
                          double prob_component_selection, 
                          double prob_feature_selection, 
                          arma::mat X, arma::vec sigma2, 
                          arma::mat tau2, arma::mat U) {
  
  // Fix mu to a vec of zeros. Note, the functions that use mu likely don't require it as an argument. 
  arma::vec mu = arma::zeros<arma::vec>(n_obs);
  
  // Set frequency of result printing
  int k = 500;
  
  // Display MCMC parameters
  std::cout << "n_iterations:" << std::endl;
  std::cout << n_iterations << std::endl;
  std::cout << "n_burnin:" << std::endl;
  std::cout << n_burnin << std::endl;
  
  // Initialize intermediate data structures
  std::cout << "Initializing intermediate data structures..." << std::endl;
  vec gamma = initialize_gamma(r, prob_component_selection);
  mat Eta = initialize_Eta(gamma, r, p_m, prob_feature_selection);
  
  std::cout << "Initial gamma: " << gamma << std::endl;
  
  // Initialize posterior chain structures 
  std::cout << "Initializing posterior summary structures..." << std::endl;
  mat gamma_chain = arma::zeros(r, n_iterations);
  arma::cube Eta_chain(r, p_m, n_iterations, arma::fill::zeros);
  
  // Sample for n_iterations
  std::cout << "Starting MCMC sampling..." << std::endl;
  for (int iter = 0; iter < n_iterations; iter++) {
    
    if (iter % k == 0) {
    std::cout << "iter:" << std::endl;
    std::cout << iter << std::endl;
    std::cout << "n Selected Components: " << sum(gamma) << std::endl;
    std::cout << "n Selected Components x Features:" << sum(Eta) << std::endl;
    }
    
    // Sample gamma_Eta
  
    for (int l = 0; l < r; l++) {
      
      if (iter % k == 0) {
      std::cout << "l:" << std::endl;
      std::cout << l << std::endl;
      std::cout << "gamma:" << std::endl;
      std::cout << gamma << std::endl;
      }
      
      vec gamma_new = gamma; 
      mat Eta_new = Eta;
      gamma_new[l] = 1 - gamma_new[l];
      
      if (iter % k == 0) {
      std::cout << "gamma_new:" << std::endl;
      std::cout << gamma_new << std::endl;
      }
      
      if (gamma_new[l] == 0) {
        for (int j = 0; j < p_m; j++) {
          Eta_new(l, j) = 0;
        }
      } else {
        for (int j = 0; j < p_m; j++) {
          // Calculate eta_lj_threshold
          double eta_lj_threshold = calculate_eta_lj_threshold(l, mu, gamma_new, Eta.col(j), sigma2[j], tau2.col(j), U, n_obs, X.col(j), prob_feature_selection);
          // Generate a random number between 0 and 1
          double random_u = arma::randu();
          // Turn on/ turn off Eta_new_lj
          double logit_random_u = log(random_u/(1.0-random_u));
          if (logit_random_u < eta_lj_threshold) {
            Eta_new(l, j) = 1;
          } else {
            Eta_new(l, j) = 0;
          }
        }
      }
      
      // Calculate log acceptance ratio
      double log_target_new = log_target_density_l(l, mu, gamma_new, Eta_new, sigma2, tau2, U, 
                                                   n_obs, p_m, X, prob_component_selection, prob_feature_selection);
      double log_target = log_target_density_l(l, mu, gamma_new, Eta_new, sigma2, tau2, U, 
                                                   n_obs, p_m, X, prob_component_selection, prob_feature_selection);
      double log_proposal_backward = log_proposal_l_density(l, mu, gamma, Eta, gamma_new, Eta_new, sigma2, tau2, U, 
                                                            n_obs, p_m, X, prob_feature_selection);
      double log_proposal_forward = log_proposal_l_density(l, mu, gamma_new, Eta_new, gamma, Eta, sigma2, tau2, U, 
                                                            n_obs, p_m, X, prob_feature_selection);
      
      double log_acceptance_ratio = log_target_new - log_target + log_proposal_backward - log_proposal_forward;
      
      // Accept/ reject proposed gamma and eta
      double random_u = arma::randu();
      double log_random_u = log(random_u);
      if (log_random_u < log_acceptance_ratio) {
        gamma[l] = gamma_new[l];
        Eta.row(l) = Eta_new.row(l);
      }
      
      // Gibb's sample to mix feature activation parameters
      if (gamma[l]==1) {

        for (int j = 0; j < p_m; j++) {
          // Calculate eta_lj_threshold
          double eta_lj_threshold = calculate_eta_lj_threshold(l, mu, gamma, Eta.col(j), sigma2[j], tau2.col(j), U, n_obs, X.col(j), prob_feature_selection);
          
          // Turn on/ turn off eta_lj
          double random_u = arma::randu();
          double logit_random_u = log(random_u/(1.0-random_u));
          if (logit_random_u < eta_lj_threshold) {
            Eta(l, j) = 1;
          } else {
            Eta(l, j) = 0;
          }
        }
      }
      
    }
    
    gamma_chain.col(iter) = gamma;
    Eta_chain.slice(iter) = Eta;
    
  }

  return Rcpp::List::create(
    Rcpp::Named("gamma_chain") = gamma_chain,
    Rcpp::Named("Eta_chain") = Eta_chain
  );

}


/*** R
library(tidyverse)

# Define parameters
seed <- 1
r <- 4
n_obs <- 200
p_m <- 10
prob_component_selection <- 0.25
prob_feature_selection <- 0.5
n_iterations <- 5000
n_burnin <- 1000

# Generate data
source("simulate_simple_data.R")

set.seed(1)
simulation_results <- simulate_iid_data(prob_component_importance = prob_component_selection)
data_list <- list(simulation_results$X_list[[1]], 
                  simulation_results$X_list[[2]], 
                  simulation_results$Y)
data_list <- lapply(data_list, scale)

# Start with the 1st view
m <- 1 
X <- data_list[[m]]

# Fix covariates
sigma2 <- rep(1, p_m)
tau2 <- matrix(1, nrow = r, ncol = p_m)
U <- simulation_results$U

# Define custom functions
get_gamma_df <- function(gamma_chain) {
  n_iterations <- ncol(gamma_chain)
  gamma_df <- gamma_chain[, 1:n_iterations] %>% 
    apply(MARGIN = 1, FUN = cummean) %>% 
    as.data.frame() %>% mutate(iteration = 1:n_iterations) %>% 
    gather(key = "gamma", value = "MPP", -iteration) %>%
    mutate(gamma = gamma %>% as.factor()) 
  return(gamma_df)
}

# Get MCMC results
start_time <- Sys.time()
gamma_Eta <- main_sample_gamma_Eta(n_iterations, n_burnin, r, n_obs, p_m,
                      prob_component_selection, prob_feature_selection,
                      X, sigma2, tau2, U)
end_time <- Sys.time()
print("MCMC Duration:")
print(end_time - start_time)

gamma_chain <- gamma_Eta$gamma_chain
Eta_chain <- gamma_Eta$Eta_chain

print("Component Selection Mean:")
gamma_chain[, (n_burnin+1):n_iterations] %>% apply(MARGIN = 1, FUN = mean)

print("Variable Selection Mean:")
Eta_chain[,,(n_burnin+1):n_iterations] %>% apply(MARGIN = c(1,2), FUN = mean)

gamma_df <- get_gamma_df(gamma_chain)
gamma_plot <- ggplot(gamma_df, aes(x = iteration, y = MPP, color = gamma)) +
  geom_line() + geom_vline(xintercept = n_burnin, linetype = "dashed", color = "red") +
  labs(x = "iteration", y = "MPP", title = "Trace plot for gamma")
print(gamma_plot)

# Analyze feature selection in component 1
features_of_interest <- c(1, 2, 6, 7)
feature_names <- paste("j", features_of_interest, sep = "_")
eta_df <- Eta_chain[1, features_of_interest, 1:n_iterations] %>%
  apply(MARGIN = 1, FUN = cummean) %>% 
  as.data.frame() %>% rename_at(vars(names(.)), ~ feature_names) %>% mutate(iteration = 1:n_iterations) %>% 
  gather(key = "eta", value = "MPP", -iteration) %>%
  mutate(eta = eta %>% as.factor()) 

eta_1_plot <- ggplot(eta_df, aes(x = iteration, y = MPP, color = eta)) + 
  geom_line() + geom_vline(xintercept = n_burnin, linetype = "dashed", color = "red") +
  labs(x = "iteration", y = "MPP", title = "Trace plot for eta_1.")

print(eta_1_plot)
*/
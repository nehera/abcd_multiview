#include <cstdlib>
#include <iostream>
#include <string>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

#include <utils.h>

// // Set parameters
// int r = 4;
// int n_obs = 200;
// int p_m = 10;
// double prob_component_selection = 0.5;
// double prob_feature_selection = 0.5;
// int n_sample = 5000;
// int n_burnin = 1000;
// int n_iterations = n_sample + n_burnin;
// 
// // Start with the 1st view
// int m = 1;
// 
// // Fix covariates
// vec sigma2 = ones(p_m);
// mat tau2 = ones(r, p_m);
// //mat U = simulation_results$U;

// [[Rcpp::export]]
vec initialize_gamma(int r=4, double prior_component_selection=0.5) {
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

// Custom struct: A custom struct can hold multiple objects and return that container. 
// This approach allows you to group related objects together and return them as a single unit. 

// // Define custom structs to hold gamma and Eta from the mth view
// struct struct_gamma_Eta_combined {
//   arma::vec gamma;
//   arma::mat Eta;
// };
// 
// // Define custom struct for sample_gamma_eta_l to return
// struct gamma_eta_l_combined {
//   int gamma_l; 
//   arma::vec eta_l;
// };
// 
// // A function that returns a dummy gamma_Eta_combined struct
// struct_gamma_Eta_combined combine_gamma_Eta(vec gamma, mat Eta) {
//   return {gamma, Eta};
// }

// // Function for exporting custom struct to R
// SEXP gamma_Eta_combined_to_R(struct_gamma_Eta_combined gamma_Eta) {
//   // Create a List to hold gamma and Eta
//   Rcpp::List list;
//   list["gamma"] = gamma_Eta.gamma;
//   list["Eta"] = gamma_Eta.Eta;
//   return Rcpp::wrap(list);
// }

// The main function should eventually return Rcpp::List

// [[Rcpp::export]]
int main_sample_gamma_Eta(int n_iterations, int n_burnin,
                          int r, int n_obs, int p_m, 
                          double prob_component_selection, 
                          double prob_feature_selection, 
                          arma::mat X, arma::vec sigma2, 
                          arma::mat tau2, arma::mat U) {
  
  // Fix mu to a vec of zeros. Note, the functions that use mu likely don't require it as an argument. 
  arma::vec mu = arma::zeros<arma::vec>(n_obs);
  
  // Display MCMC parameters
  std::cout << "n_iterations:" << std::endl;
  std::cout << n_iterations << std::endl;
  std::cout << "n_burnin:" << std::endl;
  std::cout << n_burnin << std::endl;
  
  // Initialize intermediate data structures
  std::cout << "Initializing intermediate data structures..." << std::endl;
  vec gamma = initialize_gamma(r, prob_component_selection);
  mat Eta = initialize_Eta(gamma, r, p_m, prob_feature_selection);
  // Initialize posterior chain structures 
  std::cout << "Initializing posterior summary structures..." << std::endl;
  mat gamma_chain = arma::zeros(r, n_iterations);
  arma::cube Eta_chain(r, p_m, n_iterations, arma::fill::zeros);
  
  // Sample for n_iterations
  std::cout << "Starting MCMC sampling..." << std::endl;
  for (int iter = 0; iter < n_iterations; iter++) {
    
    std::cout << "iter:" << std::endl;
    std::cout << iter << std::endl;
    std::cout << "n Selected Components:" << std::endl;
    std::cout << sum(gamma) << std::endl;
    std::cout << "n Selected Components x Features:" << std::endl;
    std::cout << sum(Eta) << std::endl;
    
    // TODO Summarize the MPPs and print
    
    // Sample gamma_Eta
  
    for (int l = 0; l < r; l++) {
      
      std::cout << "l:" << std::endl;
      std::cout << l << std::endl;
      
      // Sample gamma_Eta_l
      
      // Propose new values for the lth component

      std::cout << "gamma:" << std::endl;
      std::cout << gamma << std::endl;
      vec gamma_new = gamma; // Does changing gamma_new impact gamma?
      mat Eta_new = Eta;
      gamma_new[l] = 1 - gamma_new[l];  
      std::cout << "gamma_new:" << std::endl;
      std::cout << gamma_new << std::endl;
      
      if (gamma_new[l] == 0) {
        for (int j = 0; j < p_m; j++) {
          Eta_new(l, j) = 0;
        }
      } else {
        for (int j = 0; j < p_m; j++) {
          // Generate a random number between 0 and 1
          double random_u = arma::randu();
          // Perform Bernoulli draw
          int eta_new_lj = (random_u < prob_feature_selection) ? 1 : 0;
          Eta_new(l, j) = eta_new_lj;
        }
      }
      
      // Calculate log acceptance ratio
      // Dummy vars result in log_acceptance_ratio ~ 0.5
      // double log_target_new = 1.0;
      // double log_target = 1.6931472;
      // double log_proposal_backward = 1.0;
      // double log_proposal_forward = 1.0;
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
        
        std::cout << "Proposal accepted." << std::endl;
        gamma[l] = gamma_new[l];
        Eta.row(l) = Eta_new.row(l);
      }
      
      // Gibb's sample to mix feature activation parameters
      if (gamma[l]==1) {
        std::cout << "Gibb's sampling feature activation parameters..." << std::endl;
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
    
    // Store posterior sample
    std::cout << "Storing posterior sample..." << std::endl;
    
    gamma_chain.col(iter) = gamma;
    Eta_chain.slice(iter) = Eta;
    
    }
  
  // Write Eta_chain and gamma_chain to files
  gamma_chain.save("gamma_chain.txt", arma::raw_ascii);
  Eta_chain.save("Eta_chain.txt", arma::raw_ascii);
  
  std::cout << "Files written successfully." << std::endl;
  
  return 0;
}


/*** R
# R functions
initialize_gamma_R <- function(r=4, prior_component_selection=0.5) {
  gamma <- rbinom(n = r, size = 1, prob = prior_component_selection)
  return(gamma)
}

initialize_eta_R <- function(gamma, r=4, p_m=10, prior_feature_selection=0.5) {
  eta <- matrix(nrow = r, ncol = p_m)
  for (l in 1:r) {
    if (gamma[l] == 1) {
      eta[l, ] <- rbinom(n = p_m, size = 1, prior_feature_selection)
    } else {
      eta[l, ] <- rep(0, p_m)
    }
  }
  return(eta)
}

calculate_mvnorm_var_j_R <- function(eta_j, sigma2_j, tau2_j, U, n_obs) {
  n_components_active_in <- sum(eta_j)
  if (n_components_active_in > 0) {
    components_active_in <- which(eta_j==1)
    Sigma2_j <- U[, components_active_in, drop = FALSE] %*% 
      diag(n_components_active_in) %*% 
      t(U[, components_active_in, drop = FALSE]) +
      diag(n_obs) 
  } else {
    Sigma2_j <- diag(n_obs)
  }
  mvnorm_var <- sigma2_j * Sigma2_j
  return(mvnorm_var)
}

calculate_log_dmvnorm_j_R <- function(eta_j, sigma2_j, tau2_j, U, n_obs, x_j) {
  mvnorm_var_j <- calculate_mvnorm_var_j(eta_j, sigma2_j, tau2_j, U, n_obs)
  # TODO use woodbury matrix identity for log_dmvnorm evaluation for sigma inversion
  # TODO compare the results of the woodbury identity function with those of dmvnorm
  log_dmvnorm_j <- mvtnorm::dmvnorm(x = x_j, mean = rep(0, n_obs), sigma = mvnorm_var_j, log = TRUE) 
  return(log_dmvnorm_j)
}

calculate_log_G_j_R <- function(gamma, eta_j, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection) {
  log_dmvnorm_j <- calculate_log_dmvnorm_j_R(eta_j, sigma2_j, tau2_j, U, n_obs, x_j)
  active_gamma <- which(gamma == 1)
  n_active_gamma <- length(active_gamma) # TODO account for all inactive
  n_active_eta_given_gamma_1 <- eta_j[active_gamma] %>% sum()
  log_G_j <- log_dmvnorm_j + 
    n_active_eta_given_gamma_1 * log(prob_feature_selection) + 
    (n_active_gamma - n_active_eta_given_gamma_1) * log(1 - prob_feature_selection)
  return(log_G_j) 
}

calculate_log_PQ_lj_R <- function(l, gamma, eta_j, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection) {
  gamma_1 <- replace(gamma, l, 1)
  eta_1 <- replace(eta_j, l, 1)
  eta_0 <- replace(eta_j, l, 0)
  log_G_j_1 <- calculate_log_G_j_R(gamma_1, eta_1, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  log_G_j_0 <- calculate_log_G_j_R(gamma_1, eta_0, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  # Use logsumexp trick
  max_arg <- max(log_G_j_1, log_G_j_0) 
  a <- log_G_j_1 - max_arg
  b <- log_G_j_0 - max_arg
  # TODO understand if the log(exp(...)) is ok.
  log_P_lj <- a - log(exp(a)+exp(b)) 
  log_Q_lj <- b - log(exp(a)+exp(b))
  return(c(log_P_lj, log_Q_lj))
}

calculate_eta_lj_threshold_R <- function(l, gamma, eta_j, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection) {
  gamma_1 <- replace(gamma, l, 1)
  eta_1 <- replace(eta_j, l, 1)
  eta_0 <- replace(eta_j, l, 0)
  log_G_j_1 <- calculate_log_G_j_R(gamma_1, eta_1, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  log_G_j_0 <- calculate_log_G_j_R(gamma_1, eta_0, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  return(log_G_j_1 - log_G_j_0)
}

log_target_density_l_R <- function(l, gamma, eta, sigma2, tau2, U, n_obs, p_m, 
                                   x, prob_component_selection, prob_feature_selection) {
  sum_log_target_density_lj <- 0
  for (j in 1:p_m) {
    # TODO understand if log_dmvnorm_j should only use data within. Note, this should cancel in the log_acceptance_ratio subtraction.
    log_dmvnorm_j <- calculate_log_dmvnorm_j_R(eta[, j], sigma2[j], tau2[,j], U, n_obs, x[,j]) 
    sum_log_target_density_lj <- sum_log_target_density_lj + log_dmvnorm_j + 
      eta[l,j]*log(prob_feature_selection) + (1-eta[l,j])*log(1-prob_feature_selection)
  }
  log_target_density_l <- gamma[l]*log(prob_component_selection) + 
    (1-gamma[l])*log(1-prob_component_selection) + sum_log_target_density_lj
  return(log_target_density_l)
}

log_proposal_l_density_R <- function(l, gamma_prime, eta_prime, gamma, eta,
                                     sigma2, tau2, U, n_obs, p_m, x, prob_feature_selection) {
  if (gamma[l]==1 & gamma_prime[l]==0 & sum(eta_prime[l,])==0) {
    return(0)
  } else if (gamma[l]==0 & gamma_prime[l]==1) {
    log_proposal_l_density <- 0
    for (j in 1:p_m) {
      PQ_lj <- calculate_log_PQ_lj_R(l, gamma_prime, eta_prime[, j], sigma2[j], tau2[, j], U, n_obs, x[, j], prob_feature_selection)
      log_proposal_l_density <- log_proposal_l_density + eta_prime[l,j] * PQ_lj[1] + (1-eta_prime[l,j]) * PQ_lj[2]
    }
    return(log_proposal_l_density)
  } else {
    stop("Error in log_proposal_l_density evaluation.")
  }
}

############################################

# Generating data
source("simulate_simple_data.R")
job <- 1
simulation_results <- simulate_iid_data(seed=job)
data_list <- list(simulation_results$X_list[[1]], 
                  simulation_results$X_list[[2]], 
                  simulation_results$Y)
data_list <- lapply(data_list, scale)
m <- 1 
X <- data_list[[m]]

reindex_r_cpp <- function(x) return(x-1)
reindex_cpp_r <- function(x) return(x+1)

set.seed(123)
r <- 4
n_obs <- 200
p_m <- 10
prob_component_selection <- 0.5
prob_feature_selection <- 0.5

# Fix covariates
sigma2 <- rep(1, p_m)
tau2 <- matrix(1, nrow = r, ncol = p_m)
U <- simulation_results$U

I <- diag(n_obs)
mu <- rep(0, n_obs)


gamma <- initialize_gamma(r = 4, prior_component_selection = 0.5)
gamma
initialize_gamma_R(r = 4, prior_component_selection = 0.5)

Eta <- initialize_eta(gamma, r, p_m, 0.5)
Eta

j <- 1
l <- 2

gamma_new <- replace(gamma, l, 1-gamma[l])
Eta_new <- Eta
if (gamma_new[l] == 0) {
  Eta_new[l, ] <- rep(0, p_m)
} else {
  Eta_new[l, ] <- rbinom(n = p_m, size = 1, prob = prob_feature_selection)
}

Var_j <- calculate_mvnorm_var_j(Eta[,j], sigma2[j], tau2[,j], U, n_obs)
Var_j_R <- calculate_mvnorm_var_j_R(Eta[,j], sigma2[j], tau2[j], U, n_obs)

which(round(Var_j, 6) != round(Var_j_R, 6))

calculate_log_dmvnorm_j(X[,j], mu, U, Eta[,j], tau2[,j], sigma2[j], n_obs)
calculate_log_dmvnorm_j_R(Eta[,j], sigma2[j], tau2[,j], U, n_obs, X[,j])

calculate_log_G_j(mu, gamma, Eta[,j], sigma2[j], tau2[,j], U, n_obs, X[,j], 
                  prob_feature_selection)
calculate_log_G_j_R(gamma, Eta[,j], sigma2[j], tau2[,j], U, n_obs, X[,j], prob_feature_selection)

calculate_log_PQ_lj(reindex_r_cpp(l), mu,  gamma, Eta[,j], sigma2[j], tau2[,j], U, n_obs, 
                    X[,j], prob_feature_selection)
calculate_log_PQ_lj_R(l, gamma, Eta[,j], sigma2[j], tau2[,j], U, n_obs, X[,j], prob_feature_selection)

calculate_eta_lj_threshold(reindex_r_cpp(l), mu,  gamma, Eta[,j], sigma2[j], tau2[,j], U, 
                           n_obs, X[,j], prob_feature_selection)
calculate_eta_lj_threshold_R(l, gamma, Eta[,j], sigma2[j], tau2[,j], U, n_obs, X[,j], prob_feature_selection)

log_target_density_l(reindex_r_cpp(l), mu, gamma, Eta, sigma2, tau2, U, n_obs, p_m, 
                     X, prob_component_selection, prob_feature_selection)
log_target_density_l_R(l, gamma, Eta, sigma2, tau2, U, n_obs, p_m, 
                       X, prob_component_selection, prob_feature_selection)

log_proposal_l_density(reindex_r_cpp(l), mu, gamma_new, Eta_new, gamma, Eta, sigma2, tau2, 
                       U, n_obs, p_m, X, prob_feature_selection)
log_proposal_l_density_R(l, gamma_new, Eta_new, gamma, eta, sigma2, tau2, U, n_obs, p_m, X, prob_feature_selection)

##################################################################

set.seed(1)
start_time <- Sys.time()
main_sample_gamma_Eta(n_iterations=5000, n_burnin = 1000, 
                      r=4, n_obs, p_m, 
                      prob_component_selection, 
                      prob_feature_selection, 
                      X, sigma2, tau2, U) 
end_time <- Sys.time()
print("MCMC Duration:")
print(end_time - start_time)
*/
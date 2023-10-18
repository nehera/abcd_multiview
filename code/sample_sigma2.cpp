#include <cstdlib>
#include <iostream>
#include <string>
#include <random>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

#include "utils.cpp" 

// [[Rcpp::export]]
float sample_gamma(float alpha, float beta) {
  
  if (alpha <= 0 || beta <= 0) {
    throw std::runtime_error("Error: Alpha and Beta must be positive.");
  }
  
  // Generate a random sample from a gamma distribution
  float sample = arma::randg<arma::vec>(1, arma::distr_param(alpha, beta))(0);
  
  return sample;
}

// [[Rcpp::export]]
float sample_sigma2_j(float alpha_0, float beta_0, vec x_j, mat Sigma_j_inv, int n_obs) {

  float n_float = static_cast<float>(n_obs); // Convert int to float for division
  float alpha = alpha_0 + n_obs/2;
  vec beta_vec = 1/2 * x_j.t() * Sigma_j_inv * x_j + beta_0; // This returns a vec of len 1
  float beta = static_cast<float>(beta_vec(0)); // convert to float for sample_gamma
  
  float sigma2_j = 1.0/ sample_gamma(alpha, beta);
  return sigma2_j;
}

// [[Rcpp::export]]
vec sample_sigma2(float alpha_0, float beta_0, mat X, 
                  mat U, mat Tau, mat Eta, int n_obs, int p_m) {
  
  vec sigma2;
  
  for (int j=0; j < n_obs; j++) {
    
    // get_Sigma_j_inv
    vec tau_j = Tau.col(j);
    vec eta_j_vec = Eta.col(j);
    uvec eta_j = conv_to<uvec>::from(eta_j_vec);
    mat Sigma_j_inv = get_Sigma_j_inv(U, tau_j, eta_j, n_obs);
    
    // sample sigma2_j
    vec x_j = X.col(j);
    float sigma2_j = sample_sigma2_j(alpha_0, beta_0, x_j, Sigma_j_inv, n_obs);
    
    // store sigma2_j
    sigma2[j] = sigma2_j;
    
  }
  
  return sigma2;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
library(tidyverse)

set.seed(1)
n_samples <- 10000
alpha <- 1
beta <- 1
gamma_R <- rgamma(n_samples, alpha, beta)
gamma_cpp <- sample_gamma(n_samples, alpha, beta)

# Sample data (replace with your actual data)
# set.seed(123)
# gamma_R <- rnorm(100, mean = 0, sd = 1)
# gamma_cpp <- rnorm(100, mean = 2, sd = 1)

# Combine the vectors and create a data frame
data <- data.frame(
  Value = c(gamma_R, gamma_cpp),
  Group = rep(c("gamma_R", "gamma_cpp"), each = 100)
)

# Create a density plot
p <- ggplot(data, aes(x = Value, fill = Group, alpha = Group)) +
  geom_density(adjust = 1/4) +
  scale_fill_manual(values = c("gamma_R" = "blue", "gamma_cpp" = "red")) +
  scale_alpha_manual(values = c("gamma_R" = 0.5, "gamma_cpp" = 0.5)) +
  labs(title = "Density Plot of gamma_R and gamma_cpp", x = "Value") +
  theme_minimal()

# Display the plot
print(p)
*/

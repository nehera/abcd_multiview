#include <cstdlib>
#include <iostream>
#include <string>
#include <random>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec sample_gamma(int n, float alpha, float beta) {
  if (alpha <= 0 || beta <= 0) {
    throw std::runtime_error("Error: Alpha and Beta must be positive.");
  }
  
  // Create a random number generator engine
  std::random_device rd;
  std::mt19937 gen(rd());
  
  // Create an inverse gamma distribution
  std::gamma_distribution<float> rgamma_distribution(alpha, beta);
  
  arma::vec sample(n);
  for (int i = 0; i < n; i++) {
    // Generate a random sample from the distribution
    sample[i] = rgamma_distribution(gen);
  }
  
  return sample;
}

// [[Rcpp::export]]
float sample_sigma2_j(float alpha_0, float beta_0, int n_obs, vec x_j, mat Sigma_j) {
  // Check if the length of x_j is equal to n_obs
  if (x_j.n_elem != n_obs) {
    // Throw an exception and stop the function
    throw std::runtime_error("Length of x_j is not equal to n_obs.");
  }
  
  
  float n_float = static_cast<float>(n_obs); // Convert int to float for division
  float alpha = alpha_0 + n_obs/2;
  double beta = 1/2 * x_j.t() * Sigma_j.i() * x_j + beta_0;
  
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

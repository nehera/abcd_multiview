#include <cstdlib>
#include <iostream>
#include <string>
#include <random>
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// Function to get indices of elements equal to 1
uvec get_indices_of_ones(vec x) {
  int n = x.size();
  uvec indices;
  for (int i = 0; i < n; i++) {
    if (x[i] == 1) {
      // Create a new uvector one element longer and copy the elements
      arma::uvec newIndices(indices.n_elem + 1);
      newIndices.head(indices.n_elem) = indices;
      // Append the new element to the end
      newIndices(indices.n_elem) = i;
      // Replace the original indices with the newIndices if needed
      indices = newIndices;
    }
  }
  return indices;
}

// Function to sample from a multivariate normal distribution
// In this code: mean is the mean vector of the multivariate normal distribution.
// covariance is the covariance matrix of the multivariate normal distribution.
// The function performs the following steps:
// Generates random samples z from a standard normal distribution (mean = 0, variance = 1).
// Performs the Cholesky decomposition of the covariance matrix covariance to obtain the lower triangular matrix L.
// Samples from the multivariate normal distribution by multiplying the Cholesky-decomposed L with the random samples z and adding it to the mean vector.
// You can call the sample_mvnorm_example function with your mean vector and covariance matrix to sample from the multivariate normal distribution. 
arma::vec sample_mvnorm(const arma::vec mean, const arma::mat covariance) {
  // Generate random samples from a standard normal distribution
  arma::vec z = arma::randn<arma::vec>(mean.n_elem);
  // Perform Cholesky decomposition of the covariance matrix
  arma::mat L = arma::chol(covariance, "lower");
  // Sample from the multivariate normal distribution
  arma::vec sample = mean + L * z;
  return sample;
}

mat get_woodbury_inv(mat U, mat D) {
  mat I = eye(U.n_rows, U.n_rows);
  mat Sigma_inv = I - U*inv(inv(D) + U.t()*U)*U.t();
  return Sigma_inv;
}

arma::vec randbin(int n, int size=1, double prob=0.5) {
  vec randbin_result = arma::zeros<arma::vec>(n);
  int bern = 0;
  // Draw n observations from binomial distribution
  for (int i = 0; i < n; i++) {
    // Draw one observation from binomial distribution
    int obs_tally = 0;
    for (int j = 0; j < size; j++) {
      // Generate a random number between 0 and 1
      double random_u = arma::randu();
      // Perform Bernoulli draw
      bern = (random_u < prob) ? 1 : 0;
      obs_tally = obs_tally + bern;
    }
    randbin_result[i] = obs_tally;
  }
  return randbin_result;
}


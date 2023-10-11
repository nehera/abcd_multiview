#include <cstdlib>
#include <iostream>
#include <string>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Function to get indices of elements equal to 1
uvec getIndicesOfOnes(uvec x) {
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

// Function to extract columns from a matrix based on 'x'
arma::mat extractColumns(arma::mat m, uvec x) {
  // Get indices of elements equal to 1
  uvec indices = getIndicesOfOnes(x);
  // Extract the selected columns
  arma::mat selectedColumns = m.cols(indices);
  return selectedColumns;
}

// Function for taking the inverse of a matrix using the Woodbury Identity
mat woodburyInv(mat U, mat D) {
  mat I = eye(U.n_rows, U.n_rows);
  mat sigma_inv = I - U*inv(inv(D) + U.t()*U)*U.t();
  return sigma_inv;
}

// Function to sample from a multivariate normal distribution
// In this code: mean is the mean vector of the multivariate normal distribution.
// covariance is the covariance matrix of the multivariate normal distribution.
// The function performs the following steps:
// Generates random samples z from a standard normal distribution (mean = 0, variance = 1).
// Performs the Cholesky decomposition of the covariance matrix covariance to obtain the lower triangular matrix L.
// Samples from the multivariate normal distribution by multiplying the Cholesky-decomposed L with the random samples z and adding it to the mean vector.
// You can call the sampleMultivariateNormalExample function with your mean vector and covariance matrix to sample from the multivariate normal distribution. 
arma::vec sampleMultivariateNormal(const arma::vec mean, const arma::mat covariance) {
  // Generate random samples from a standard normal distribution
  arma::vec z = arma::randn<arma::vec>(mean.n_elem);
  // Perform Cholesky decomposition of the covariance matrix
  arma::mat L = arma::chol(covariance, "lower");
  // Sample from the multivariate normal distribution
  arma::vec sample = mean + L * z;
  return sample;
}

// Example usage
// [[Rcpp::export]]
arma::vec sampleMultivariateNormalExample(const arma::vec mean, const arma::mat covariance) {
  return sampleMultivariateNormal(mean, covariance);
}

// [[Rcpp::export]]
arma::vec sample_a_j(vec x_j, float sigma2_j, mat U, 
                     uvec gamma, uvec eta_j, int r) {
  
  int n_components_active = sum(gamma);
  
  mat Sigma_a(r, r, arma::fill::zeros);
  vec mu_a(r, arma::fill::zeros);
  vec a_j(r, arma::fill::zeros);
  
  if (n_components_active > 0) {
    
    uvec indices = getIndicesOfOnes(gamma);

    mat activeU = U.cols(indices);
    
    mat Sigma_a_inv = (1/ sigma2_j) * (activeU.t() * activeU + eye(n_components_active, n_components_active));
    Sigma_a = Sigma_a_inv.i();
    
    mu_a = Sigma_a * activeU.t() * x_j;
    
    vec active_a_j = sampleMultivariateNormal(mu_a, Sigma_a);
    
    // If eta_lj = 0, a_lj = 0
    for (int i = 0; i < n_components_active; i++) {
      int index = indices[i];
      if (eta_j[index] == 1) {
        a_j[index] = active_a_j[i];
      }
    }
  } 
  
  // Note, if n_components_active == 0, we return a vec of zeros.
  return a_j;
}

// [[Rcpp::export]]
arma::vec sample_A(mat X, vec sigma2, mat U, 
                   uvec gamma, mat eta, int r, int p_m) {
  
  arma::mat A(r, p_m, arma::fill::zeros);
  
  for (int j = 0; j < p_m; j++) {
    std::cout << j << std::endl;
    vec x_j = X.col(j);
    float sigma2_j = sigma2[j];
    // eta_j needs to be a uvec
    arma::vec doubleVec = eta.col(j);
    arma::uvec eta_j = arma::conv_to<arma::uvec>::from(doubleVec.head(r)); // Use .head(r) to ensure a consistent size
    vec a_j = sample_a_j(x_j, sigma2_j, U, gamma, eta_j, r);
    std::cout << a_j << std::endl;
    A.col(j) = a_j;
    std::cout << A.col(j) << std::endl;
  }
  
  return A;
  
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
source("simulate_simple_data.R")
simulation_results = simulate_iid_data(n_obs = 100)
m = 1 # start with first view
j = 1
sigma2_j = 1
U = simulation_results$U
gamma = c(1,1,0,0)
eta_j = c(1,1,0,0)
r = 4
X = simulation_results$X_list[[m]]
x_j = X[, j]
print("Sampled a_j:")
a_j = sample_a_j(x_j, sigma2_j, U, gamma, eta_j, r)
a_j
print("True a_j:")
simulation_results$A_list[[m]][, j]

# Now we sample the whole A matrix
p_m = 10
# TODO remove hardcoding
Eta = matrix(c(rep(eta_j, 5), rep(0, r*5)), nrow = r, ncol = p_m) 
sigma2 = rep(sigma2_j, p_m)

print("Sampled A:")
sample_A(X, sigma2, U, gamma, Eta, r, p_m)
print("True A:")
simulation_results$A_list[[m]]
*/

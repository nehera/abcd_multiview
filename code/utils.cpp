#include <cstdlib>
#include <iostream>
#include <string>
#include <random>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// Function to get indices of elements equal to 1
uvec get_indices_of_ones(uvec x) {
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

// [[Rcpp::export]]
mat get_woodbury_inv(mat U, mat D) {
  mat I = eye(U.n_rows, U.n_rows);
  mat sigma_inv = I - U*inv(inv(D) + U.t() * U) * U.t();
  return sigma_inv;
}

// [[Rcpp::export]]
arma::mat get_Sigma_j(mat U, vec tau_j, uvec eta_j, int n_obs) {
  
  mat Sigma_j;
  
  int n_components_active_in = arma::sum(eta_j);
  
  arma::uvec components_active_in;
  arma::mat Sigma2_j;
  
  if (n_components_active_in > 0) {
    
    components_active_in = find(eta_j == 1);
    mat U_active = U.cols(components_active_in);
    vec tau_active = tau_j(components_active_in);
    mat tau_diag = diagmat(tau_active);
    Sigma2_j = U_active * tau_diag * U_active.t() + eye(n_obs, n_obs);
    
  } else {
    
    Sigma2_j = eye(n_obs, n_obs);
    
  }
  
  return Sigma_j;
  
}

// [[Rcpp::export]]
arma::mat get_Sigma_j_inv(mat U, vec tau_j, uvec eta_j, int n_obs) {
  
  mat Sigma_j_inv;
  
  int n_components_active_in = sum(eta_j);
  
  if (n_components_active_in > 0) {
    
    uvec components_active_in = find(eta_j == 1);
    mat U_active = U.cols(components_active_in);
    vec tau_active = tau_j(components_active_in);
    mat tau_diag = diagmat(tau_active);
    Sigma_j_inv = get_woodbury_inv(U_active, tau_diag);
    
  } else {
    
    Sigma_j_inv = eye(n_obs, n_obs);
    
  }
  
  return Sigma_j_inv;
  
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/

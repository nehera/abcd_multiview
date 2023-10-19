#include <cstdlib>
#include <iostream>
#include <string>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// Set parameters
int r = 4;
int n_obs = 200;
int p_m = 10;
double prob_component_selection = 0.5;
double prob_feature_selection = 0.5;
int n_sample = 5000;
int n_burnin = 1000;
int n_iterations = n_sample + n_burnin;

// Start with the 1st view
int m = 1;
//x <- data_list[[m]

// Fix covariates
vec sigma2 = ones(p_m);
mat tau2 = ones(r, p_m);
//mat U = simulation_results$U;

// [[Rcpp::export]]
int main() {
  arma_rng::set_seed(123);
  // TO-DO: add code for submitting MSI job
  cout << n_iterations;
}

// Function to get indices of elements equal to 1
// [[Rcpp::export]]
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
vec initialize_gamma(int r=4, double prior_component_selection=0.5) {
  vec gamma = Rcpp::rbinom(r, 1, prior_component_selection);
  
  return gamma;
}

// [[Rcpp::export]]
mat initialize_eta(vec gamma, int r=4, int p_m=10, double prior_feature_selection=0.5) {
  mat Eta(r, p_m, fill::value(datum::nan));
  
  for(int l = 0; l < r; l++) {
    if(gamma(l) == 1) {
      rowvec Eta_l = Rcpp::rbinom(p_m, 1, prior_feature_selection);
      Eta.row(l) = Eta_l;
    } else {
      Eta.row(l).zeros();
    }
  }
  
  return Eta;
}

// [[Rcpp::export]]
mat calculate_mvnorm_var_j(uvec eta_j, double sigma2_j, vec tau2_j, mat U, int n_obs) {
  // Where's tau2_j used?
  int n_components_active_in = sum(eta_j);
  mat Sigma2_j(n_obs, n_obs, fill::value(datum::nan));
  
  if (n_components_active_in > 0) {
    // Utils file?
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
mat woodburyInv(mat U, mat D) {
  mat I = eye(U.n_rows, U.n_rows);
  mat Sigma_inv = I - U*inv(inv(D) + U.t()*U)*U.t();
  
  return Sigma_inv;
}

// [[Rcpp::export]]
double calculate_log_dmvnorm_j(vec x, vec mu, mat U, mat D) {
  // This function takes two matrices U and D as inputs instead of a covariance matrix,
  // such that UDU'+I = Sigma, the covariance matrix
  double k = x.n_elem;
  int r = eta_j.n_elem;
  mat D(r, r, fill::value(0));
  D.diag() = tau2_j;
  
  double Sigma_det = det(inv(D) + U.t()*U)*det(D);
  mat Sigma_j = calculate_mvnorm_var_j(eta_j, sigma2_j, tau2_j, U, n_obs);
  mat Sigma_inv = woodburyInv(U, D);
  
  double log_dens = -0.5*k*log(2*datum::pi) - 0.5*log(Sigma_det) - 0.5*as_scalar(trans(x-mu)*Sigma_inv*(x-mu));
  
  return log_dens;
}

calculate_log_dmvnorm_j <- function(eta_j, sigma2_j, tau2_j, U, n_obs, x_j) {
  mvnorm_var_j <- calculate_mvnorm_var_j(eta_j, sigma2_j, tau2_j, U, n_obs)
// TODO use woodbury matrix identity for log_dmvnorm evaluation for sigma inversion
// TODO compare the results of the woodbury identity function with those of dmvnorm
  log_dmvnorm_j <- mvtnorm::dmvnorm(x = x_j, mean = rep(0, n_obs), sigma = mvnorm_var_j, log = TRUE) 
  return(log_dmvnorm_j)
}



/*** R
set.seed(123)
r <- 4
n_obs <- 200
p_m <- 10
prob_component_selection <- 0.5
prob_feature_selection <- 0.5

# Fix covariates
sigma2 <- rep(1, p_m)
tau2 <- matrix(1, nrow = r, ncol = p_m)
U <- matrix(rnorm(n_obs*r), ncol = r)
D <- diag(r)
I <- diag(n_obs)
sigma <- U %*% D %*% t(U) + I

gamma <- initialize_gamma(r = 4, prior_component_selection = 0.5)
gamma

Eta <- initialize_eta(gamma, r, p_m, 0.5)
Eta

j <- 1
Sigma_j <- calculate_mvnorm_var_j(Eta[,j], sigma2[j], tau2[,j], U, n_obs)
Sigma_j[1:3, 1:3]
sigma[1:3, 1:3]
*/

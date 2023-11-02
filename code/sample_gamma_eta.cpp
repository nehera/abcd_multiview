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
double calculate_log_dmvnorm_j(vec x_j, vec mu_j, mat U, vec eta_j, vec tau2_j, double sigma2_j, uvec active_components, int n_obs) {
  // This function takes two matrices U and D as inputs instead of a covariance matrix,
  // such that UDU'+I = Sigma, the covariance matrix
  int r = U.n_cols;
  mat D(r, r, fill::value(0));
  D.diag() = tau2_j.elem(active_components);
  
  double Sigma_det = det(inv(D) + U.t()*U)*det(D);
  mat Sigma_inv = woodburyInv(U, D);
   
  mat Sigma = U*D*U.t() + eye(n_obs, n_obs);
  
  double log_dens = -0.5*n_obs*log(2*datum::pi) - 0.5*log(Sigma_det) - 0.5*as_scalar(trans(x_j-mu_j)*Sigma_inv*(x_j-mu_j));
  
  return log_dens;
}

// [[Rcpp::export]]
double calculate_log_G_j(vec mu_j, vec gamma, vec eta_j, double sigma2_j, vec tau2_j, mat U, uvec active_components, int n_obs, vec x_j, double prob_feature_selection) {
  double log_dmvnorm_j = calculate_log_dmvnorm_j(x_j, mu_j, U, eta_j, tau2_j, sigma2_j, active_components, n_obs);
  int n_active_gamma = active_components.n_elem; // TODO account for all inactive
  int n_active_eta_given_gamma_1 = sum(eta_j.elem(active_components));
  double log_G_j = log_dmvnorm_j + 
      n_active_eta_given_gamma_1 * log(prob_feature_selection) + 
      (n_active_gamma - n_active_eta_given_gamma_1) * log(1 - prob_feature_selection);
  return log_G_j;
}

// TO-DO - test starting here
// [[Rcpp::export]]
vec calculate_log_PQ_lj(int l, vec mu_j, vec gamma, vec eta_j, double sigma2_j, vec tau2_j, mat U, uvec active_components, int n_obs, vec x_j, double prob_feature_selection) {
  vec log_PQ(2, fill::value(datum::nan));
  vec gamma_1 = gamma;
  vec eta_1 = eta_j;
  vec eta_0 = eta_j;
  gamma_1(l) = 1;
  eta_1(l) = 1;
  eta_0(l) = 0;
  double log_G_j_1 = calculate_log_G_j(mu_j, gamma_1, eta_1, sigma2_j, tau2_j, U, active_components, n_obs, x_j, prob_feature_selection);
  double log_G_j_0 = calculate_log_G_j(mu_j, gamma_1, eta_0, sigma2_j, tau2_j, U, active_components, n_obs, x_j, prob_feature_selection);
  // Use logsumexp trick
  double max_arg = std::max(log_G_j_1, log_G_j_0); 
  double a = log_G_j_1 - max_arg;
  double b = log_G_j_0 - max_arg;
  // TODO understand if the log(exp(...)) is ok.
  double log_P_lj = a - log(exp(a)+exp(b)); 
  double log_Q_lj = b - log(exp(a)+exp(b));
  log_PQ(0) = log_P_lj;
  log_PQ(1) = log_Q_lj;
  
  cout << eta_1 << ": "<< eta_0 << "\n";
  cout << log_G_j_1<< ": " << log_G_j_0 << "\n";
  return log_PQ;
}

// [[Rcpp::export]]
double calculate_eta_lj_threshold(int l, vec mu_j, vec gamma, vec eta_j, double sigma2_j, vec tau2_j, mat U, uvec active_components, int n_obs, vec x_j, double prob_feature_selection) {
  vec gamma_1 = gamma;
  vec eta_1 = eta_j;
  vec eta_0 = eta_j;
  gamma_1(l) = 1;
  eta_1(l) = 1;
  eta_0(l) = 0;
  
  double log_G_j_1 = calculate_log_G_j(mu_j, gamma_1, eta_1, sigma2_j, tau2_j, U, active_components, n_obs, x_j, prob_feature_selection);
  double log_G_j_0 = calculate_log_G_j(mu_j, gamma_1, eta_0, sigma2_j, tau2_j, U, active_components, n_obs, x_j, prob_feature_selection);
  return log_G_j_1 - log_G_j_0;
}

// [[Rcpp::export]]
double log_target_density_l(int l, vec mu, vec gamma, vec eta, double sigma2, vec tau2, mat U, uvec active_components, int n_obs, int p_m, 
                                 vec x_j, double prob_component_selection, double prob_feature_selection) {
  double sum_log_target_density_lj = 0;
  for (int j=0; j<p_m; j++) {
// TODO understand if log_dmvnorm_j should only use data within. Note, this should cancel in the log_acceptance_ratio subtraction.
    double log_dmvnorm_j = calculate_log_dmvnorm_j(x_j, mu, U, eta, tau2, sigma2, active_components, n_obs);
    double sum_log_target_density_lj = sum_log_target_density_lj + log_dmvnorm_j + 
      eta(l,j)*log(prob_feature_selection) + (1-eta(l,j))*log(1-prob_feature_selection);
  }
  double log_target_density_l = gamma(l)*log(prob_component_selection) + 
    (1-gamma(l))*log(1-prob_component_selection) + sum_log_target_density_lj;
  return(log_target_density_l);
    
}

// Main function



/*** R
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
U <- matrix(rnorm(n_obs*r), ncol = r)

I <- diag(n_obs)
x_j <- rep(0, n_obs)
mu_j <- rep(0, n_obs)


gamma <- initialize_gamma(r = 4, prior_component_selection = 0.5)
gamma

Eta <- initialize_eta(gamma, r, p_m, 0.5)
Eta

j <- 1
l <- 2
n_components_active <- sum(Eta[,j])
active_components <- which(Eta[,j] == 1)
U_gamma <- U[,active_components]

sigma <- U_gamma %*% diag(tau2[active_components,j]) %*% t(U_gamma) + I
Sigma_j <- calculate_mvnorm_var_j(Eta[,j], sigma2[j], tau2[,j], U, n_obs)

active_components <- reindex_r_cpp(active_components)
l <- reindex_r_cpp(l)

calculate_log_dmvnorm_j(x_j, mu_j, U_gamma, Eta[,j], tau2[,j], sigma2[j], active_components, n_obs)
mvtnorm::dmvnorm(x_j, mu_j, Sigma_j, log=T)

#(mu_j, gamma, eta_j, sigma2_j,  tau2_j, U, active_components, n_obs, x_j, prob_feature_selection)
calculate_log_G_j(mu_j, gamma, Eta[,j], sigma2[j], tau2[,j], U_gamma, active_components, n_obs, x_j, prob_feature_selection)
calculate_log_PQ_lj(l, mu_j,  gamma, Eta[,j], sigma2[j], tau2[,j], U_gamma, active_components, n_obs, x_j, prob_feature_selection)
calculate_eta_lj_threshold(l, mu_j,  gamma, Eta[,j], sigma2[j], tau2[,j], U_gamma, active_components, n_obs, x_j, prob_feature_selection)
log_target_density_l(l, mu_j, gamma, Eta[,j], sigma2[j], tau2[,j], U_gamma, active_components, n_obs, p_m, 
                     x, prob_component_selection, prob_feature_selection)
*/

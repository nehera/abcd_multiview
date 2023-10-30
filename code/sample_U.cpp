#include <cstdlib>
#include <iostream>
#include <string>
#include <random>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

#include <utils.h>

// [[Rcpp::export]]
arma::mat sample_U(int n_obs, int p, int r, mat X_combined,
                   mat A_combined, vec sigma2_combined) {
  
  // Initialize U
  arma::mat U(n_obs, r, arma::fill::zeros);

  // Calculate Sigma_u
  arma::mat Sigma_u_inv = A_combined * diagmat(sigma2_combined) * A_combined.t() + eye(r,r);
  arma::mat Sigma_u = Sigma_u_inv.i();

  for (int i = 0; i < n_obs; i++) {

    // Extract the ith observation's data
    // .row returns a row vector, so we take t() to ensure we have a column vector
    arma::mat X_i(p, 1, arma::fill::zeros);
    X_i = X_combined.row(i).t(); 
    
    // Calculate mu_u_i
    vec mu_u_i = Sigma_u * A_combined * diagmat(sigma2_combined) * X_i;
    
    // Sample u_i from rmvnorm
    vec u_i = sample_mvnorm(mu_u_i, Sigma_u);
    
    // Repace ith row of U with u_i row vector
    U.row(i) = u_i.t(); 
    
  }
  
  return U;
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
source("simulate_simple_data.R")

n_obs <- 100
r <- 4

simulation_results <- simulate_iid_data(n_obs = 100)
attach(simulation_results)

X_combined <- cbind(Y, do.call(cbind, X_list))

p <- ncol(X_combined)

# TODO Understand what A_0 should be?
A_0 <- matrix(rnorm(r), nrow = r, ncol = 1)
A_combined <- cbind(A_0, do.call(cbind, A_list))

# TODO Understand if sigma2_combined should be fixed to 1 for now. 
sigma2_combined <- rep(1, p)

U_sample <- sample_U(n_obs, p, r, X_combined, A_combined, sigma2_combined)

print("i:")
i = 1
print(i)

print("Sampled u_i:")
U_sample[i, ] %>% head()

print("True u_i:")
U[i, ] %>% head()

# Check U_mean with all else set to truth
n_sample <- 5000
n_burnin <- 1000
n_iterations <- n_sample + n_burnin
U_chain <- array(NA, dim = c(n_obs, r, n_iterations))
for (iter in 1:n_iterations) {
  U_chain[,,iter] <- sample_U(n_obs, p, r, X_combined, A_combined, sigma2_combined)
}
# Note, we don't take after burn-in since no mixing required
U_mean <- apply(U_chain, MARGIN = c(1,2), mean)
print("Mean U across iterations:")
U_mean
print("Difference between U_mean and U_truth:")
U_diff <- U_mean - U
U_diff
print("Frobenius Norm of Difference between U_mean and U_truth:")
norm(U_diff, type = "F")
*/

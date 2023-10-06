#include <cstdlib>
#include <iostream>
#include <string>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
mat woodburyInv(mat U, mat D) {
  mat I = eye(U.n_rows, U.n_rows);
  mat sigma_inv = I - U*inv(inv(D) + U.t()*U)*U.t();
  
  return sigma_inv;
}

// [[Rcpp::export]]
double log_dmvnorm(vec x, vec mu, mat U, mat D) {
  // This function takes two matrices U and D as inputs instead of a covariance matrix,
  // such that UDU'+I = Sigma, the covariance matrix
  double k = x.n_elem;
  double sigma_det = det(inv(D) + U.t()*U)*det(D);
  mat sigma_inv = woodburyInv(U, D);
  
  double log_dens = -0.5*k*log(2*datum::pi) - 0.5*log(sigma_det) - 0.5*as_scalar(trans(x-mu)*sigma_inv*(x-mu));
  
  return log_dens;
}

/*** R
library(mvtnorm)

set.seed(123)
n <- 10
r <- 4
I <- diag(n)
U <- matrix(rnorm(n*r), ncol = r)
D <- diag(r)
sigma <- U %*% D %*% t(U) + I

which(round(woodburyInv(U, D), 6) != round(solve(sigma), 6))
which(round(det(solve(D) + t(U)%*%U)*det(D), 6) != round(det(sigma), 6))
which(round(exp(log_dmvnorm(rep(0, n), rep(0, n), U, D)), 6) != round(dmvnorm(rep(0, n), rep(0, n), sigma), 6))

log_dmvnorm(rep(0, n), rep(0, n), U, D)
dmvnorm(rep(0, n), rep(0, n), sigma, log = T)
*/

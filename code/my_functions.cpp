#include <cstdlib>
#include <iostream>
#include <string>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
// TODO revise!!!
double sample_intercept_cpp(arma::mat Y, mat A_outcome, mat U, double Sigma2_outcome, double Sigma2_0 = 100) {
  int n = Y.n_rows;
  mat Y_star(n, 1); Y_star.fill(0);
  for (int i=0; i<n; i++) {
    mat U_star_i = U.row(i) * A_outcome;
    Y_star(i, 0) = Y(i, 0) - U_star_i(0,0); // consider matrix/ vector workaround?
  }
  vec Y_star_mean = arma::mean(Y_star); // one element
  double invSig2 = n/ Sigma2_outcome + 1/ Sigma2_0;
  double u = Rcpp::rnorm(1)(0);
  double intercept = (n * Y_star_mean(0)) / (invSig2 * Sigma2_outcome) + sqrt(1/ invSig2) * u; //Rcpp::rnorm(1);
  return intercept;
}

// [[Rcpp::export]]
double get_P_lj_cpp(double logG_lj_0, double logG_lj_1) {
  double x_0 = std::max(logG_lj_0, logG_lj_1);
  double x = logG_lj_1 - x_0;
  double y = logG_lj_0 - x_0;
  return (exp(x) / (exp(x) + exp(y)));
}

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

/*** R
timesTwo(42)
# get_P_lj_cpp(-3, -4)
# get_P_lj(-3, -4)
*/

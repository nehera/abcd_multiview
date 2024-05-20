#include <RcppArmadillo.h>
#include <random>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// Function to sample from a normal distribution
vec rnorm_cpp(int n, double mean, double sd) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> d(mean, sd);
  
  vec samples(n);
  for (int i = 0; i < n; ++i) {
    samples[i] = d(gen);
  }
  return samples;
}

// Function to sample from an inverse gamma distribution
double rinvgamma_cpp(double shape, double scale) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::gamma_distribution<> d(shape, 1.0 / scale);
  return 1.0 / d(gen);
}

// [[Rcpp::export]]
List gibbs_sampler_nested(const vec& y, const mat& X, uvec study_site, uvec family, int n_iter, const vec& beta_prior_mean, const vec& beta_prior_var, double sigma_u_prior_a, double sigma_u_prior_b, double sigma_v_prior_a, double sigma_v_prior_b, double sigma_prior_a, double sigma_prior_b) {
  
  // Adjust for 0-based indexing in C++
  study_site -= 1;
  family -= 1;
  
  int n = y.n_elem;
  int J = max(study_site) + 1;
  int F_count = max(family) + 1;
  
  // Initialize parameters
  vec beta = zeros<vec>(2);
  vec u = zeros<vec>(J);
  vec v = zeros<vec>(F_count);
  double sigma_u2 = 1.0;
  vec sigma_v2 = ones<vec>(J);
  double sigma2 = 1.0;
  
  // Storage for samples
  mat beta_samples(n_iter, 2, fill::zeros);
  mat u_samples(n_iter, J, fill::zeros);
  mat v_samples(n_iter, F_count, fill::zeros);
  vec sigma_u2_samples(n_iter, fill::zeros);
  mat sigma_v2_samples(n_iter, J, fill::zeros);
  vec sigma2_samples(n_iter, fill::zeros);
  
  for (int iter = 0; iter < n_iter; ++iter) {
    // Sample beta
    mat XtX = X.t() * X;
    mat V_beta = inv(XtX / sigma2 + diagmat(1.0 / beta_prior_var));
    vec m_beta = V_beta * (X.t() * (y - u(study_site) - v(family)) / sigma2 + beta_prior_mean / beta_prior_var);
    beta = mvnrnd(m_beta, V_beta);
    
    // Sample u
    for (int j = 0; j < J; ++j) {
      uvec idx = find(study_site == j);
      vec y_j = y(idx);
      mat X_j = X.rows(idx);
      vec v_j = v(family(idx));
      double V_u = 1.0 / (y_j.n_elem / sigma2 + 1.0 / sigma_u2);
      double m_u = V_u * sum(y_j - X_j * beta - v_j) / sigma2;
      u[j] = rnorm_cpp(1, m_u, sqrt(V_u))[0];
    }
    
    // Sample v
    for (int f = 0; f < F_count; ++f) {
      uvec idx = find(family == f);
      vec y_f = y(idx);
      mat X_f = X.rows(idx);
      vec u_f = u(study_site(idx));
      int site = study_site(idx[0]);
      double V_v = 1.0 / (y_f.n_elem / sigma2 + 1.0 / sigma_v2[site]);
      double m_v = V_v * sum(y_f - X_f * beta - u_f) / sigma2;
      v[f] = rnorm_cpp(1, m_v, sqrt(V_v))[0];
    }
    
    // Sample sigma_u2
    double a_sigma_u = sigma_u_prior_a + J / 2.0;
    double b_sigma_u = sigma_u_prior_b + sum(square(u)) / 2.0;
    sigma_u2 = rinvgamma_cpp(a_sigma_u, b_sigma_u);
    
    // Sample sigma_v2 for each study site
    for (int j = 0; j < J; ++j) {
      uvec families_in_site = find(study_site == j);
      families_in_site = unique(family(families_in_site));
      if (!families_in_site.is_empty()) {
        for (size_t idx = 0; idx < families_in_site.n_elem; ++idx) {
          if (families_in_site[idx] < 0 || families_in_site[idx] >= F_count) {
            Rcpp::Rcerr << "Invalid family index: " << families_in_site[idx] << std::endl;
            Rcpp::stop("Index out of bounds");
          }
        }
        double a_sigma_v = sigma_v_prior_a + families_in_site.n_elem / 2.0;
        // std::cout << "v(families_in_site): " << v(families_in_site).t() << std::endl;
        double b_sigma_v = sigma_v_prior_b + sum(square(v(families_in_site))) / 2.0;
        // std::cout << "b_sigma_v: " << b_sigma_v << std::endl;
        sigma_v2[j] = rinvgamma_cpp(a_sigma_v, b_sigma_v);
      }
    }
    
    // Sample sigma2
    double a_sigma = sigma_prior_a + n / 2.0;
    double b_sigma = sigma_prior_b + sum(square(y - X * beta - u(study_site) - v(family))) / 2.0;
    sigma2 = rinvgamma_cpp(a_sigma, b_sigma);
    
    // Store samples
    beta_samples.row(iter) = beta.t();
    u_samples.row(iter) = u.t();
    v_samples.row(iter) = v.t();
    sigma_u2_samples[iter] = sigma_u2;
    sigma_v2_samples.row(iter) = sigma_v2.t();
    sigma2_samples[iter] = sigma2;
  }
  
  // Return samples as a list
  return List::create(
    Named("beta_samples") = beta_samples,
    Named("u_samples") = u_samples,
    Named("v_samples") = v_samples,
    Named("sigma_u2_samples") = sigma_u2_samples,
    Named("sigma_v2_samples") = sigma_v2_samples,
    Named("sigma2_samples") = sigma2_samples
  );
}


/*** R
# # Example data
# set.seed(123)
# n_iter <- 1000
# n <- 900
# J <- 30
# F_count <- 30
# study_site <- rep(1:J, each = n / J)
# family <- rep(1:F_count, each = n / F_count)
# X <- cbind(1, rnorm(n))
# beta_true <- c(1, 2)
# sigma_u_true <- 2
# sigma_v_true <- 1.5
# u_true <- rnorm(J, sd = sigma_u_true)
# v_true <- rnorm(F_count, sd = sigma_v_true)
# sigma_true <- 1
# y <- X %*% beta_true + u_true[study_site] + v_true[family] + rnorm(n, sd = sigma_true)
# 
# # Priors
# priors <- list(beta_prior_mean = c(0, 0),
#                beta_prior_var = c(100, 100),
#                sigma_u_prior_a = 2,
#                sigma_u_prior_b = 1,
#                sigma_v_prior_a = 2,
#                sigma_v_prior_b = 1,
#                sigma_prior_a = 2,
#                sigma_prior_b = 1)
# 
# # Run Gibbs sampler
# samples <- gibbs_sampler_nested(y, X, study_site, family, n_iter,
#                                 priors$beta_prior_mean, priors$beta_prior_var,
#                                 priors$sigma_u_prior_a, priors$sigma_u_prior_b,
#                                 priors$sigma_v_prior_a, priors$sigma_v_prior_b,
#                                 priors$sigma_prior_a, priors$sigma_prior_b)
# 
# # Summary of results
# summary(samples$beta_samples)
# summary(samples$u_samples)
# summary(samples$v_samples)
# summary(samples$sigma_u2_samples)
# summary(samples$sigma_v2_samples)
# summary(samples$sigma2_samples)
*/

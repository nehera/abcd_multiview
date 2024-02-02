#include <cstdlib>
#include <iostream>
#include <string>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

//#include <sample_sigma2.h>
#include <utils.h>


// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::export]]
vec initialize_gamma(int r=4, double prior_component_selection=0.5) {
  //vec gamma = rbinom(r, 1, prior_component_selection);
  vec gamma = randbin(r, 1, prior_component_selection);
  return gamma;
}

// [[Rcpp::export]]
mat initialize_Eta(vec gamma, int r=4, int p_m=10, double prior_feature_selection=0.5) {
  mat Eta(r, p_m, fill::value(datum::nan));
  
  for(int l = 0; l < r; l++) {
    if(gamma(l) == 1) {
      vec Eta_l = randbin(p_m, 1, prior_feature_selection);
      Eta.row(l) = Eta_l.t();
    } else {
      Eta.row(l).zeros();
    }
  }
  
  return Eta;
}

// [[Rcpp::export]]
double calculate_log_dmvnorm_j(vec x_j, vec mu, mat U, vec eta_j, vec tau2_j, double sigma2_j, int n_obs) {
  int n_active_components = sum(eta_j);
  double Sigma_det = 0;
  mat Sigma_inv(n_obs, n_obs, fill::value(datum::nan));
  
  if (n_active_components > 0) {
    uvec active_components = get_indices_of_ones(eta_j);
    
    //mat D(n_active_components, n_active_components, fill::value(0));
    U = U.cols(active_components);
    // Get the dimensions of U_active
    int numRows = U.n_rows;
    int numCols = U.n_cols;
    // Create a new matrix D of the same size
    arma::mat D(numCols, numCols, arma::fill::zeros);
    //D.diag() = tau2_j.elem(active_components);

    // Assign selected elements from tau2_j to the diagonal of D
    for (int i = 0; i < numCols; ++i) {
      D(i, i) = tau2_j(active_components(i));
    }
    
    //Sigma_det = pow(sigma2_j, n_obs)*det(inv(D) + U.t()*U)*det(D);
    
    //mat Sigma_j = U*D*U.t() + eye(n_obs, n_obs);
    //Sigma_det = pow(sigma2_j, n_obs)*det(Sigma_j);
    
    Sigma_inv = (1.0/sigma2_j)*get_woodbury_inv(U, D);
    Sigma_det = pow(sigma2_j, n_obs)/ det(Sigma_inv);
    
  } else {
    // Sigma2_j = I, an n-by-n identity matrix
    Sigma_det = pow(sigma2_j, n_obs);
    Sigma_inv = (1.0/sigma2_j)*eye(n_obs, n_obs);
  }
  
  std::cout << "Sigma_det: " << Sigma_det << endl;
  
  double log_dens = -0.5*n_obs*log(2*datum::pi) - 0.5*log(Sigma_det) - 0.5*as_scalar(trans(x_j-mu)*Sigma_inv*(x_j-mu));
  
  return log_dens;
}

// calculate_quad_form_and_log_dmvnorm_j returns a vec 
// 0th index: quadratic form t(x_j)*Sigma_j_inv*x_j
// 1st index: log_dmvnorm_j

// [[Rcpp::export]]
vec calculate_quad_form_and_log_dmvnorm_j(vec x_j, vec mu, mat U, vec eta_j, vec tau2_j, double sigma2_j, int n_obs) {
  int n_active_components = sum(eta_j);
  double Sigma_det = 0;
  mat Sigma_inv(n_obs, n_obs, fill::value(datum::nan));
  if (n_active_components > 0) {
    uvec active_components = get_indices_of_ones(eta_j);
    mat D(n_active_components, n_active_components, fill::value(0));
    D.diag() = tau2_j.elem(active_components);
    U = U.cols(active_components);
    Sigma_det = pow(sigma2_j, n_obs)*det(inv(D) + U.t()*U)*det(D);
    Sigma_inv = (1.0/sigma2_j)*get_woodbury_inv(U, D);
  } else {
    // Sigma2_j = I, an n-by-n identity matrix
    Sigma_det = pow(sigma2_j, n_obs);
    Sigma_inv = (1.0/sigma2_j)*eye(n_obs, n_obs);
  }
  double quad_form = as_scalar(trans(x_j-mu)*Sigma_inv*(x_j-mu)); 
  double log_dens = -0.5*n_obs*log(2*datum::pi) - 0.5*log(Sigma_det) - 0.5*quad_form;
  vec quad_form_and_log_dmvnorm_j = join_cols(vec{quad_form}, vec{log_dens});
  return quad_form_and_log_dmvnorm_j;
}

// [[Rcpp::export]]
double calculate_log_G_j(vec mu, vec gamma, vec eta_j, double sigma2_j, vec tau2_j, mat U, int n_obs, vec x_j, double prob_feature_selection) {
  double log_dmvnorm_j = calculate_log_dmvnorm_j(x_j, mu, U, eta_j, tau2_j, sigma2_j, n_obs);
  uvec active_gamma = get_indices_of_ones(gamma);
  int n_active_gamma = active_gamma.n_elem; // TODO account for all inactive
  int n_active_eta_given_gamma_1 = sum(eta_j.elem(active_gamma));
  double log_G_j = log_dmvnorm_j + 
    n_active_eta_given_gamma_1 * log(prob_feature_selection) + 
    (n_active_gamma - n_active_eta_given_gamma_1) * log(1 - prob_feature_selection);
  return log_G_j;
}

// [[Rcpp::export]]
vec calculate_log_PQ_lj(int l, vec mu, vec gamma, vec eta_j, double sigma2_j, vec tau2_j, mat U, int n_obs, vec x_j, double prob_feature_selection) {
  vec log_PQ(2);
  vec gamma_1 = gamma;
  vec eta_1 = eta_j;
  vec eta_0 = eta_j;
  gamma_1(l) = 1;
  eta_1(l) = 1;
  eta_0(l) = 0;
  double log_G_j_1 = calculate_log_G_j(mu, gamma_1, eta_1, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection);
  double log_G_j_0 = calculate_log_G_j(mu, gamma_1, eta_0, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection);
  
  // Use logsumexp trick
  double max_arg = std::max(log_G_j_1, log_G_j_0); 
  double a = log_G_j_1 - max_arg;
  double b = log_G_j_0 - max_arg;
  
  double log_P_lj = a - log(exp(a)+exp(b)); 
  double log_Q_lj = b - log(exp(a)+exp(b));
  log_PQ(0) = log_P_lj;
  log_PQ(1) = log_Q_lj;
  
  return log_PQ;
}

// [[Rcpp::export]]
double calculate_eta_lj_threshold(int l, vec mu, vec gamma, vec eta_j, double sigma2_j, vec tau2_j, mat U, int n_obs, vec x_j, double prob_feature_selection) {
  vec gamma_1 = gamma;
  vec eta_1 = eta_j;
  vec eta_0 = eta_j;
  gamma_1(l) = 1;
  eta_1(l) = 1;
  eta_0(l) = 0;
  
  double log_G_j_1 = calculate_log_G_j(mu, gamma_1, eta_1, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection);
  double log_G_j_0 = calculate_log_G_j(mu, gamma_1, eta_0, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection);
  return log_G_j_1 - log_G_j_0;
}

// // [[Rcpp::export]]
// double log_target_density_l(int l, vec mu, vec gamma, mat Eta, vec sigma2, mat Tau, mat U, int n_obs, int p_m, 
//                             mat X, double prob_component_selection, double prob_feature_selection) {
//   double sum_log_target_density_lj = 0;
//   for (int j=0; j<p_m; j++) {
//     // TODO understand if log_dmvnorm_j should only use data within. 
//     // Note, this should cancel in the log_acceptance_ratio subtraction.
//     double log_dmvnorm_j = calculate_log_dmvnorm_j(X.col(j), mu, U, Eta.col(j), Tau.col(j), sigma2(j), n_obs);
//     sum_log_target_density_lj = sum_log_target_density_lj + log_dmvnorm_j + 
//       Eta(l,j)*log(prob_feature_selection) + (1-Eta(l,j))*log(1-prob_feature_selection);
//   }
//   double log_target_density_l = gamma(l)*log(prob_component_selection) + 
//     (1-gamma(l))*log(1-prob_component_selection) + sum_log_target_density_lj;
//   return log_target_density_l;
// }

// [[Rcpp::export]]
double log_target_density_l(int l, vec mu, vec gamma, mat Eta, vec sigma2, mat Tau, mat U, int n_obs, int p_m, 
                            mat X, double prob_component_selection, double prob_feature_selection) {
  double sum_log_target_density_lj = 0;
  for (int j=0; j<p_m; j++) {
    // TODO understand if log_dmvnorm_j should only use data within. 
    // Note, this should cancel in the log_acceptance_ratio subtraction.
    double log_dmvnorm_j = calculate_log_dmvnorm_j(X.col(j), mu, U, Eta.col(j), Tau.col(j), sigma2(j), n_obs);
    
    if(gamma(l) == 0) {
      sum_log_target_density_lj = sum_log_target_density_lj + log_dmvnorm_j;
    } else {
      sum_log_target_density_lj = sum_log_target_density_lj + log_dmvnorm_j + 
        Eta(l,j)*log(prob_feature_selection) + (1-Eta(l,j))*log(1-prob_feature_selection);  
    }
  }
  double log_target_density_l = gamma(l)*log(prob_component_selection) + 
    (1-gamma(l))*log(1-prob_component_selection) + sum_log_target_density_lj;
  return log_target_density_l;
}

// [[Rcpp::export]]
double diff_target_density_l(int l, vec mu, 
                             vec gamma_prime, mat Eta_prime, 
                             vec gamma, mat Eta, 
                             vec sigma2, mat Tau, mat U, int n_obs, int p_m, 
                             mat X, double prob_component_selection, 
                             double prob_feature_selection) {
  
  double diff_log_prob_component = log(1-prob_component_selection)-log(prob_component_selection);
  if (gamma_prime(l)==1) {
    double diff_log_prob_component = -diff_log_prob_component;
  }
  
  // std::cout << "diff_log_prob_component: " << diff_log_prob_component << std::endl;
  
  double diff_log_dmvnorm = 0;
  
  int n_active = 0;
  
  for (int j=0; j<p_m; j++) {
    double log_dmvnorm_prime_j = calculate_log_dmvnorm_j(X.col(j), mu, U, Eta_prime.col(j), Tau.col(j), sigma2(j), n_obs);
    double log_dmvnorm_j = calculate_log_dmvnorm_j(X.col(j), mu, U, Eta.col(j), Tau.col(j), sigma2(j), n_obs);
    diff_log_dmvnorm = diff_log_dmvnorm + log_dmvnorm_prime_j - log_dmvnorm_j;
    if (gamma_prime(l)==0) {
      n_active = n_active + Eta.col(j)[l];
    } else {
      n_active = n_active + Eta_prime.col(j)[l];
    }
  }
  
  // std::cout << "diff_log_dmvnorm: " << diff_log_dmvnorm << std::endl;
  // std::cout << "n_active: " << n_active << std::endl;
  
  double diff_log_prob_feature = -log(prob_feature_selection)*n_active -log(1-prob_feature_selection)*(p_m-n_active);
  
  if (gamma_prime(l)==1) {
    diff_log_prob_feature = -diff_log_prob_feature;
    }
  
  // std::cout << "diff_log_prob_feature: " << diff_log_prob_feature << std::endl;
  
  double diff_log_target_density_l = diff_log_prob_component + diff_log_dmvnorm + diff_log_prob_feature;
  return diff_log_target_density_l;
}

// [[Rcpp::export]]
double log_proposal_l_density(int l, vec mu, vec gamma_prime, mat Eta_prime, 
                              vec gamma, mat Eta, vec sigma2, mat Tau, mat U, 
                              int n_obs, int p_m, mat X, double prob_feature_selection) {
  if ((gamma(l)==1) & (gamma_prime(l)==0) & (sum(Eta_prime.row(l))==0)) {
    return 0;
  } else if ((gamma(l)==0) & (gamma_prime(l)==1)) {
    double log_proposal_l_density = 0;
    for (int j=0; j<p_m; j++) {
      vec PQ_lj = calculate_log_PQ_lj(l, mu, gamma_prime, Eta_prime.col(j), sigma2(j), 
                                      Tau.col(j), U, n_obs, X.col(j), prob_feature_selection);
      log_proposal_l_density = log_proposal_l_density + Eta_prime(l,j) * PQ_lj(0) + (1-Eta_prime(l,j)) * PQ_lj(1);
    } 
    return log_proposal_l_density;
  } else { 
    throw std::invalid_argument("Error in log_proposal_l_density evaluation.");
    return datum::nan;
  }
}

// [[Rcpp::export]]
double diff_log_proposal_l(int l, vec mu, vec gamma_prime, mat Eta_prime, 
                           vec gamma, mat Eta, vec sigma2, mat Tau, mat U, 
                           int n_obs, int p_m, mat X, double prob_feature_selection) {
  double diff_log_proposal_l = 0;
  if(gamma_prime(l) == 0) {
    for (int j=0; j<p_m; j++) {
      vec PQ_lj = calculate_log_PQ_lj(l, mu, gamma_prime, Eta_prime.col(j), sigma2(j), 
                                      Tau.col(j), U, n_obs, X.col(j), prob_feature_selection);
      diff_log_proposal_l = diff_log_proposal_l + Eta(l,j) * PQ_lj(0) + (1-Eta(l,j)) * PQ_lj(1);
    }
  } else {
    for (int j=0; j<p_m; j++) {
      vec PQ_lj = calculate_log_PQ_lj(l, mu, gamma, Eta.col(j), sigma2(j), 
                                      Tau.col(j), U, n_obs, X.col(j), prob_feature_selection);
      diff_log_proposal_l = diff_log_proposal_l - Eta_prime(l,j) * PQ_lj(0) - (1-Eta_prime(l,j)) * PQ_lj(1);
    }
  }
  
  return diff_log_proposal_l;
} 

// [[Rcpp::export]]
float sample_gamma(float alpha, float beta) {
  
  if (alpha <= 0 || beta <= 0) {
    throw std::runtime_error("Error: Alpha and Beta must be positive.");
  }
  
  // Generate a random sample from a gamma distribution
  float sample = arma::randg<arma::vec>(1, arma::distr_param(alpha, beta))(0);
  
  return sample;
}

// [[Rcpp::export]]
arma::mat get_Sigma_j(mat U, vec tau_j, vec eta_j, int n_obs) {
  
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
arma::mat get_Sigma_j_inv(mat U, vec tau_j, vec eta_j, int n_obs) {
  
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

// [[Rcpp::export]]
float sample_sigma2_j(float alpha_0, float beta_0, vec x_j, mat Sigma_j_inv, int n_obs) {
  
  float n_float = static_cast<float>(n_obs); // Convert int to float for division
  float alpha = alpha_0 + n_float/2;
  vec beta_vec = 1/2 * x_j.t() * Sigma_j_inv * x_j + beta_0; // This returns a vec of len 1
  float beta = static_cast<float>(beta_vec(0)); // convert to float for sample_gamma
  
  float sigma2_j = 1.0/sample_gamma(alpha, 1.0/beta);
  
  return sigma2_j;
}

// [[Rcpp::export]]
vec sample_sigma2(float alpha_0, float beta_0, mat X, 
                  mat U, mat Tau, mat Eta, int n_obs, int p_m) {
  vec sigma2(p_m);
  vec tau_j;
  vec eta_j;
  mat Sigma_j_inv;
  vec x_j;
  float sigma2_j = 0;
  
  for (int j=0; j < p_m; j++) {
    // get_Sigma_j_inv
    tau_j = Tau.col(j);
    eta_j = Eta.col(j);
    //uvec eta_j = conv_to<uvec>::from(eta_j_vec);
    
    //Sigma_j_inv = eye(n_obs, n_obs);
    Sigma_j_inv = get_Sigma_j_inv(U, tau_j, eta_j, n_obs);
    
    // sample sigma2_j
    x_j = X.col(j);
    sigma2_j = sample_sigma2_j(alpha_0, beta_0, x_j, Sigma_j_inv, n_obs);
    
    // store sigma2_j
    sigma2(j) = sigma2_j;
    
  }
  
  return sigma2;
}

// [[Rcpp::export]]
arma::mat sample_U(int n_obs, int p, int r, mat X_combined,
                   mat A_combined, vec sigma2_combined, arma::mat &U) {
  
  // Initialize U
  //arma::mat U(n_obs, r, arma::fill::zeros);
  
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
    
    // Sample u_i from rmvnorm. Note: Returned u_i is column vector. 
    vec u_i = sample_mvnorm(mu_u_i, Sigma_u);
    
    // Repace ith row of U with u_i row vector
    U.row(i) = u_i.t(); 
    
  }
  
  return U;
  
}

// Functions for sampling A matrix
// [[Rcpp::export]]
arma::vec sample_a_j(vec x_j, float sigma2_j, mat U, 
                     vec gamma, vec eta_j, int r) {
  
  int n_components_active = sum(gamma);
  
  mat Sigma_a(r, r, arma::fill::zeros);
  vec mu_a(r, arma::fill::zeros);
  vec a_j(r, arma::fill::zeros);
  
  if (n_components_active > 0) {
    
    uvec indices = get_indices_of_ones(gamma);
    
    mat activeU = U.cols(indices);
    
    mat Sigma_a_inv = (1.0/sigma2_j) * (activeU.t() * activeU + eye(n_components_active, n_components_active));
    Sigma_a = Sigma_a_inv.i();
    
    mu_a = Sigma_a * activeU.t() * x_j;
    
    vec active_a_j = sample_mvnorm(mu_a, Sigma_a);
    
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
arma::mat sample_A(mat X, vec sigma2, mat U, 
                   vec gamma, mat eta, int r, int p_m) {
  
  arma::mat A(r, p_m, arma::fill::zeros);
  
  for (int j = 0; j < p_m; j++) {
    
    vec x_j = X.col(j);
    float sigma2_j = sigma2[j];
    // eta_j needs to be a uvec
    arma::vec eta_j = eta.col(j);
    //arma::uvec eta_j = arma::conv_to<arma::uvec>::from(doubleVec.head(r)); // Use .head(r) to ensure a consistent size
    vec a_j = sample_a_j(x_j, sigma2_j, U, gamma, eta_j, r);
    A.col(j) = a_j; 
    
  }
  
  return A;
  
}

// [[Rcpp::export]]
void main_sample_gamma_Eta(int iter, int n_burnin,
                           int r, int n_obs, int p_m, 
                           double prob_component_selection, 
                           double prob_feature_selection, 
                           arma::mat X, arma::vec sigma2, 
                           arma::mat Tau, arma::mat U,
                           arma::vec &gamma, arma::mat &Eta, arma::vec mu,
                           arma::mat &gamma_chain, arma::cube &Eta_chain) {
  
  int k = 3000;
  vec gamma_new = gamma; // Does changing gamma_new impact gamma?
  mat Eta_new = Eta;
  
  // Sample for n_iterations
  if (iter % k == 0) {
    std::cout << "iter: " << iter << std::endl;
    std::cout << "n Selected Components: " << sum(gamma) << std::endl;
    std::cout << "n Selected Components x Features: " << std::endl;
    std::cout << sum(Eta) << std::endl;
  }
  // TODO Summarize the MPPs and print
  
  // Sample gamma_Eta
  for (int l = 0; l < r; l++) {
    if(iter % k == 0) {
      std::cout << "l:" << std::endl;
      std::cout << l << std::endl;
      
      // Sample gamma_Eta_l
      
      // Propose new values for the lth component
      
      //std::cout << "gamma:" << std::endl;
      //std::cout << gamma << std::endl;
    }
    
    gamma_new[l] = 1 - gamma_new[l];
    
    if (iter % k == 0) {
      std::cout << "gamma_new:" << std::endl;
      std::cout << gamma_new << std::endl;
    }
    
    if (gamma_new[l] == 0) {
      for (int j = 0; j < p_m; j++) {
        Eta_new(l, j) = 0;
      }
    } else {
      for (int j = 0; j < p_m; j++) {
        double eta_lj_threshold = calculate_eta_lj_threshold(l, mu, gamma, Eta.col(j), sigma2[j], Tau.col(j), U, n_obs, X.col(j), prob_feature_selection);
        
        // Generate a random number between 0 and 1
        double random_u = arma::randu();
        double logit_random_u = log(random_u/(1.0-random_u));
        
        int eta_new_lj = (logit_random_u < eta_lj_threshold) ? 1 : 0;
        Eta_new(l, j) = eta_new_lj;
      }
    }
    
    // Calculate log acceptance ratio
    // Dummy vars result in log_acceptance_ratio ~ 0.5
    // double log_target_new = 1.0;
    // double log_target = 1.6931472;
    // double log_proposal_backward = 1.0;
    // double log_proposal_forward = 1.0;
    // double log_target_new = log_target_density_l(l, mu, gamma_new, Eta_new, sigma2, Tau, U, 
    //                                              n_obs, p_m, X, prob_component_selection, prob_feature_selection);
    // 
    // std::cout << "log_target_new: " << log_target_new << std::endl;
    // 
    // double log_target = log_target_density_l(l, mu, gamma, Eta, sigma2, Tau, U, 
    //                                          n_obs, p_m, X, prob_component_selection, prob_feature_selection);
    // 
    // std::cout << "log_target: " << log_target << std::endl;
    // 
    // std::cout << "diff_log_target_0: " << log_target_new - log_target << std::endl;
    
    double diff_log_target = diff_target_density_l(l, mu, gamma_new, Eta_new, 
                                                   gamma, Eta, 
                                                   sigma2, Tau, U, n_obs, p_m, 
                                                   X, prob_component_selection, 
                                                   prob_feature_selection);
    
    std::cout << "diff_log_target: " << diff_log_target << std::endl;

    double diff_log_proposal = diff_log_proposal_l(l, mu, gamma_new, Eta_new, gamma, Eta, sigma2, Tau, U, 
                                                   n_obs, p_m, X, prob_feature_selection);
    
    std::cout << "diff_log_proposal: " << diff_log_proposal << std::endl;
    
    //double log_acceptance_ratio = log_target_new - log_target + diff_log_proposal;
    double log_acceptance_ratio = diff_log_target + diff_log_proposal;
    
    std::cout << "log_acceptance_ratio: " << log_acceptance_ratio << std::endl;
    
    // Accept/ reject proposed gamma and eta
    double random_u = arma::randu();
    double log_random_u = log(random_u);
    if (log_random_u < log_acceptance_ratio) {
      
      //std::cout << "Proposal accepted." << std::endl;
      gamma[l] = gamma_new[l];
      Eta.row(l) = Eta_new.row(l);
      
      //n_accepted(l, iter) = 1;
      
    } else {
      //std::cout << "Proposal rejected." << std::endl;
    }
    
    // Gibb's sample to mix feature activation parameters
    if (gamma[l]==1) {
      //std::cout << "Gibb's sampling feature activation parameters..." << std::endl;
      for (int j = 0; j < p_m; j++) {
        // Calculate eta_lj_threshold
        double eta_lj_threshold = calculate_eta_lj_threshold(l, mu, gamma, Eta.col(j), sigma2[j], Tau.col(j), U, n_obs, X.col(j), prob_feature_selection);
        
        // Turn on/ turn off eta_lj
        double random_u = arma::randu();
        double logit_random_u = log(random_u/(1.0-random_u));
        if (logit_random_u < eta_lj_threshold) {
          Eta(l, j) = 1;
        } else {
          Eta(l, j) = 0;
        }
      }
    }
    
  }
  
  // Store posterior sample
  //std::cout << "Storing posterior sample..." << std::endl;
  
  gamma_chain.col(iter) = gamma;
  Eta_chain.slice(iter) = Eta;
  
  
  
  // Write Eta_chain and gamma_chain to files
  // gamma_chain.save("gamma_chain.txt", arma::raw_ascii);
  // Eta_chain.save("Eta_chain.txt", arma::raw_ascii);
  // 
  // std::cout << "Files written successfully." << std::endl;
  
  // return Rcpp::List::create(
  //   Rcpp::Named("gamma_chain") = gamma_chain,
  //   Rcpp::Named("Eta_chain") = Eta_chain,
  //   Rcpp::Named("n_accepted") = n_accepted
  // );
}

// [[Rcpp::export]]
Rcpp::List BIP(int n_iterations, int n_burnin, int r, int n_obs, int p_m,
               double prob_component_selection, double prob_feature_selection,
               arma::mat X, arma::vec sigma2, arma::mat Tau, arma::mat U,
               int p, arma::mat X_combined, arma::mat A_combined, arma::vec sigma2_combined) {
  // Fix mu to a vec of zeros. Note, the functions that use mu likely don't require it as an argument.
  arma::vec mu = arma::zeros<arma::vec>(n_obs);
  
  //
  mat n_accepted(r, n_iterations, fill::value(datum::nan));
  
  // Display MCMC parameters
  std::cout << "n_iterations:" << std::endl;
  std::cout << n_iterations << std::endl;
  std::cout << "n_burnin:" << std::endl;
  std::cout << n_burnin << std::endl;
  
  // Initialize intermediate data structures
  std::cout << "Initializing intermediate data structures..." << std::endl;
  vec gamma = initialize_gamma(r, prob_component_selection);
  mat Eta = initialize_Eta(gamma, r, p_m, prob_feature_selection);
  
  std::cout << "Initial gamma: " << std::endl; 
  std::cout << gamma << std::endl;
  
  std::cout << "Initial Eta: " << std::endl; 
  std::cout << Eta << std::endl;
  
  // Initialize posterior chain structures
  std::cout << "Initializing posterior summary structures..." << std::endl;
  mat gamma_chain = arma::zeros(r, n_iterations);
  arma::cube Eta_chain(r, p_m, n_iterations, arma::fill::zeros);
  mat sigma2_chain = arma::zeros(p_m, n_iterations);
  arma::cube U_chain(n_obs, r, n_iterations, arma::fill::zeros);
  
  //int a_0 = 1;
  //int b_0 = 1;
  
  //List sigma2_views = List::create();
  //List Tau_views = List::create();
  //vec sigma2_combined();
  
  std::cout << "Starting MCMC sampling..." << std::endl;
  for (int iter = 0; iter < n_iterations; iter++) {
    // for (int m = 0; m < M+1; m++) {
    //     p_m = P[m];
    //     sigma2_views[[m]] = rep(1, p_m);
    //     Tau_views[[m]] = matrix(rep(1, r*p_m), nrow = r, ncol = p_m);
    // }
    
    //Tau = Tau_views[1];
    main_sample_gamma_Eta(iter, n_burnin,
                          r=4, n_obs, p_m,
                          prob_component_selection,
                          prob_feature_selection,
                          X, sigma2, Tau, U,
                          gamma, Eta, mu, gamma_chain, Eta_chain);
    
    if (iter%500 == 0) {
      std::cout << "Test" << std::endl;
      std::cout << gamma << std::endl;
    }
    
    // What about A?
    
    //sigma2 = sample_sigma2(a_0, b_0, X, U, Tau, Eta, n_obs, p_m);
    sigma2_chain.col(iter) = sigma2;
    
    //U = sample_U(n_obs, p, r, X_combined, A_combined, sigma2_combined, U);
    U_chain.slice(iter) = U;
  }
  
  
  
  return Rcpp::List::create(
    Rcpp::Named("gamma_chain") = gamma_chain,
    Rcpp::Named("Eta_chain") = Eta_chain,
    Rcpp::Named("sigma2_chain") = sigma2_chain,
    Rcpp::Named("U_chain") = U_chain
  );
}

/*** R
# R functions
initialize_gamma_R <- function(r=4, prior_component_selection=0.5) {
  gamma <- rbinom(n = r, size = 1, prob = prior_component_selection)
  return(gamma)
}

initialize_eta_R <- function(gamma, r=4, p_m=10, prior_feature_selection=0.5) {
  eta <- matrix(nrow = r, ncol = p_m)
  for (l in 1:r) {
    if (gamma[l] == 1) {
      eta[l, ] <- rbinom(n = p_m, size = 1, prior_feature_selection)
    } else {
      eta[l, ] <- rep(0, p_m)
    }
  }
  return(eta)
}

calculate_mvnorm_var_j_R <- function(eta_j, sigma2_j, tau2_j, U, n_obs) {
  n_components_active_in <- sum(eta_j)
  if (n_components_active_in > 0) {
    components_active_in <- which(eta_j==1)
    Sigma2_j <- U[, components_active_in, drop = FALSE] %*% 
      diag(n_components_active_in) %*% 
      t(U[, components_active_in, drop = FALSE]) +
      diag(n_obs) 
  } else {
    Sigma2_j <- diag(n_obs)
  }
  mvnorm_var <- sigma2_j * Sigma2_j
  return(mvnorm_var)
}

calculate_log_dmvnorm_j_R <- function(eta_j, sigma2_j, tau2_j, U, n_obs, x_j) {
  mvnorm_var_j <- calculate_mvnorm_var_j_R(eta_j, sigma2_j, tau2_j, U, n_obs)
  # TODO use woodbury matrix identity for log_dmvnorm evaluation for sigma inversion
  # TODO compare the results of the woodbury identity function with those of dmvnorm
  log_dmvnorm_j <- mvtnorm::dmvnorm(x = x_j, mean = rep(0, n_obs), sigma = mvnorm_var_j, log = TRUE) 
  return(log_dmvnorm_j)
}

calculate_log_G_j_R <- function(gamma, eta_j, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection) {
  log_dmvnorm_j <- calculate_log_dmvnorm_j_R(eta_j, sigma2_j, tau2_j, U, n_obs, x_j)
  active_gamma <- which(gamma == 1)
  n_active_gamma <- length(active_gamma) # TODO account for all inactive
  n_active_eta_given_gamma_1 <- eta_j[active_gamma] %>% sum()
  log_G_j <- log_dmvnorm_j + 
    n_active_eta_given_gamma_1 * log(prob_feature_selection) + 
    (n_active_gamma - n_active_eta_given_gamma_1) * log(1 - prob_feature_selection)
  return(log_G_j) 
}

calculate_log_PQ_lj_R <- function(l, gamma, eta_j, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection) {
  gamma_1 <- replace(gamma, l, 1)
  eta_1 <- replace(eta_j, l, 1)
  eta_0 <- replace(eta_j, l, 0)
  log_G_j_1 <- calculate_log_G_j_R(gamma_1, eta_1, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  log_G_j_0 <- calculate_log_G_j_R(gamma_1, eta_0, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  # Use logsumexp trick
  max_arg <- max(log_G_j_1, log_G_j_0) 
  a <- log_G_j_1 - max_arg
  b <- log_G_j_0 - max_arg
  # TODO understand if the log(exp(...)) is ok.
  log_P_lj <- a - log(exp(a)+exp(b)) 
  log_Q_lj <- b - log(exp(a)+exp(b))
  return(c(log_P_lj, log_Q_lj))
}

calculate_eta_lj_threshold_R <- function(l, gamma, eta_j, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection) {
  gamma_1 <- replace(gamma, l, 1)
  eta_1 <- replace(eta_j, l, 1)
  eta_0 <- replace(eta_j, l, 0)
  log_G_j_1 <- calculate_log_G_j_R(gamma_1, eta_1, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  log_G_j_0 <- calculate_log_G_j_R(gamma_1, eta_0, sigma2_j, tau2_j, U, n_obs, x_j, prob_feature_selection)
  return(log_G_j_1 - log_G_j_0)
}

log_target_density_l_R <- function(l, gamma, eta, sigma2, Tau, U, n_obs, p_m, 
                                   x, prob_component_selection, prob_feature_selection) {
  sum_log_target_density_lj <- 0
  for (j in 1:p_m) {
    # TODO understand if log_dmvnorm_j should only use data within. Note, this should cancel in the log_acceptance_ratio subtraction.
    log_dmvnorm_j <- calculate_log_dmvnorm_j_R(eta[, j], sigma2[j], Tau[,j], U, n_obs, x[,j]) 
    sum_log_target_density_lj <- sum_log_target_density_lj + log_dmvnorm_j + 
      eta[l,j]*log(prob_feature_selection) + (1-eta[l,j])*log(1-prob_feature_selection)
  }
  log_target_density_l <- gamma[l]*log(prob_component_selection) + 
    (1-gamma[l])*log(1-prob_component_selection) + sum_log_target_density_lj
  return(log_target_density_l)
}

log_proposal_l_density_R <- function(l, gamma_prime, eta_prime, gamma, eta,
                                     sigma2, Tau, U, n_obs, p_m, x, prob_feature_selection) {
  if (gamma[l]==1 & gamma_prime[l]==0 & sum(eta_prime[l,])==0) {
    return(0)
  } else if (gamma[l]==0 & gamma_prime[l]==1) {
    log_proposal_l_density <- 0
    for (j in 1:p_m) {
      PQ_lj <- calculate_log_PQ_lj_R(l, gamma_prime, eta_prime[, j], sigma2[j], Tau[, j], U, n_obs, x[, j], prob_feature_selection)
      log_proposal_l_density <- log_proposal_l_density + eta_prime[l,j] * PQ_lj[1] + (1-eta_prime[l,j]) * PQ_lj[2]
    }
    return(log_proposal_l_density)
  } else {
    stop("Error in log_proposal_l_density evaluation.")
  }
}

get_gamma_df <- function(gamma_chain) {
  n_iterations <- ncol(gamma_chain)
  gamma_df <- gamma_chain[, 1:n_iterations] %>% 
    apply(MARGIN = 1, FUN = cummean) %>% 
    as.data.frame() %>% mutate(iteration = 1:n_iterations) %>% 
    gather(key = "gamma", value = "MPP", -iteration) %>%
    mutate(gamma = gamma %>% as.factor()) 
  return(gamma_df)
}

# Trace plot function
trace_plot <- function(df, parameter, prob_component_selection, prob_feature_selection, n_burnin) {
  ggplot(df,
         aes(x = iteration, y = MPP)) +
    geom_line(aes_string(color = parameter)) + geom_vline(xintercept = n_burnin,
                                                          linetype = "dashed", color = "red") +
    labs(x = "iteration", y = "MPP",
         title = sprintf("Trace plot for %s", parameter),
         subtitle = sprintf("prob_comp_selection = %.2f, prob_feat_selection = %.2f", 
                            prob_component_selection, 
                            prob_feature_selection))
}

############################################

# Generating data
source("00_simulate_simple_data.R")
job <- 1
simulation_results <- simulate_iid_data()

attach(simulation_results)
M <- length(X_list)
X_combined <- cbind(Y, do.call(cbind, X_list))
p <- ncol(X_combined)

r <- 4
n_obs <- 200
n_iterations <- 2000
n_burnin <- 1000
p_m <- 10
P <- c(1, rep(p_m, M))
prob_component_selection <- 0.5
prob_feature_selection <- 0.5
a_0 <- 1
b_0 <- 1
sigma2_views <- list()
Tau_views <- list()

# TODO Understand what A_0 should be?
A_0 <- matrix(rnorm(r), nrow = r, ncol = 1)
A_combined <- cbind(A_0, do.call(cbind, A_list))

data_list <- list(Y, X_list[[1]], X_list[[2]])
#data_list <- lapply(data_list, scale)
m <- 2 
X <- data_list[[m]]

# Scale data
#data_list <- lapply(data_list, scale)
indic_var <- c(0, 0, 1) 
method <- "BIP" # method without grouping information to start
group_list <- NULL # default when grouping information not included
#bip_0 <- BIPnet::BIP(dataList = data_list, IndicVar = indic_var, Method = method, probvarsel = prob_feature_selection)

reindex_r_cpp <- function(x) return(x-1)
reindex_cpp_r <- function(x) return(x+1)

# Set R and Rcpp seeds; NOTE: they are not generating the same random numbers
#RNGkind(sample.kind = "Rounding")
set.seed(job)
#set_seed(job)

# Some other things I tried
#rngCpp(5)
#cbind( runif(5), rnorm(5), rt(5, 5), rbeta(5, 1, 1))
#setArmaRNGseed(1)

# Fix covariates
# sigma2_j=20
for(m in 1:(M+1)) {
  p_m <- P[m]
  # sigma2_j < 0.0245 causes issues
  sigma2_views[[m]] <- rep(1, p_m)
  Tau_views[[m]] <- matrix(rep(1, r*p_m), nrow = r, ncol = p_m)
}

sigma2_combined <- unlist(sigma2_views)
p <- ncol(X_combined)

#sigma2 <- rep(1, p_m)
#Tau <- matrix(1, nrow = r, ncol = p_m)
Tau <- Tau_views[[m]]
sigma2 <- sigma2_views[[m]]

U <- simulation_results$U

I <- diag(n_obs)
mu <- rep(0, n_obs)

gamma <- initialize_gamma(r = 4, prior_component_selection = 0.5)
gamma
initialize_gamma_R(r = 4, prior_component_selection = 0.5)

Eta <- initialize_Eta(gamma, r, p_m, 0.5)
Eta

j <- 1
l <- 1

gamma_new <- replace(gamma, l, 1-gamma[l])
Eta_new <- Eta
if (gamma_new[l] == 0) {
  Eta_new[l, ] <- rep(0, p_m)
} else {
  Eta_new[l, ] <- rbinom(n = p_m, size = 1, prob = prob_feature_selection)
}

#Var_j <- calculate_mvnorm_var_j(Eta[,j], sigma2[j], Tau[,j], U, n_obs)
#Var_j_R <- calculate_mvnorm_var_j_R(Eta[,j], sigma2[j], Tau[j], U, n_obs)

#which(round(Var_j, 6) != round(Var_j_R, 6))

# calculate_log_dmvnorm_j(X[,j], mu, U, Eta[,j], Tau[,j], sigma2[j], n_obs)
# calculate_log_dmvnorm_j_R(Eta[,j], sigma2[j], Tau[,j], U, n_obs, X[,j])
# 
# calculate_log_G_j(mu, gamma, Eta[,j], sigma2[j], Tau[,j], U, n_obs, X[,j], 
#                   prob_feature_selection)
# calculate_log_G_j_R(gamma, Eta[,j], sigma2[j], Tau[,j], U, n_obs, X[,j], prob_feature_selection)
# 
# calculate_log_PQ_lj(reindex_r_cpp(l), mu,  gamma, Eta[,j], sigma2[j], Tau[,j], U, n_obs, 
#                     X[,j], prob_feature_selection)
# calculate_log_PQ_lj_R(l, gamma, Eta[,j], sigma2[j], Tau[,j], U, n_obs, X[,j], prob_feature_selection)
# 
# calculate_eta_lj_threshold(reindex_r_cpp(l), mu,  gamma, Eta[,j], sigma2[j], Tau[,j], U, 
#                            n_obs, X[,j], prob_feature_selection)
# calculate_eta_lj_threshold_R(l, gamma, Eta[,j], sigma2[j], Tau[,j], U, n_obs, X[,j], prob_feature_selection)
# 
# log_target_density_l(reindex_r_cpp(l), mu, gamma, Eta, sigma2, Tau, U, n_obs, p_m, 
#                      X, prob_component_selection, prob_feature_selection)
# log_target_density_l_R(l, gamma, Eta, sigma2, Tau, U, n_obs, p_m, 
#                        X, prob_component_selection, prob_feature_selection)
# 
# log_proposal_l_density(reindex_r_cpp(l), mu, gamma_new, Eta_new, gamma, Eta, sigma2, Tau, 
#                        U, n_obs, p_m, X, prob_feature_selection)
# log_proposal_l_density_R(l, gamma_new, Eta_new, gamma, eta, sigma2, Tau, U, n_obs, p_m, X, prob_feature_selection)

##################################################################

start_time <- Sys.time()
chains <- BIP(n_iterations=n_iterations, n_burnin = n_burnin,
              r=4, n_obs, p_m,
              prob_component_selection,
              prob_feature_selection,
              X, sigma2, Tau, U,
              p, X_combined, A_combined, sigma2_combined)
end_time <- Sys.time()
print("MCMC Duration:")
print(end_time - start_time)

attach(chains)
#n_accepted[is.nan(n_accepted)] <- 0
#apply(n_accepted, 1, mean)

gamma_df <- get_gamma_df(gamma_chain)
trace_plot(gamma_df, "gamma", prob_component_selection, prob_feature_selection, 
           n_burnin)

# component_selection_probabilities <- c(0.1, 0.5, 0.9)
# feature_selection_probabilities <- c(0.1, 0.5, 0.9)
# plots <- list()
# pct_accepted <- matrix(nrow = length(component_selection_probabilities)*length(feature_selection_probabilities),
#                        ncol = 4)
# setting <- 1
# 
# for(prob_component_selection in component_selection_probabilities) {
#   for(prob_feature_selection in feature_selection_probabilities) {
#     gamma_Eta <- main_sample_gamma_Eta(n_iterations=n_iterations, n_burnin = n_burnin,
#                                        r=4, n_obs, p_m,
#                                        prob_component_selection,
#                                        prob_feature_selection,
#                                        X, sigma2, Tau, U)
#     
#     attach(gamma_Eta)
#     n_accepted[is.nan(n_accepted)] <- 0
#     pct_accepted[setting, ] <- apply(n_accepted, 1, mean)
#     
#     gamma_df <- get_gamma_df(gamma_chain)
#     
#     plots[[setting]] <- trace_plot(gamma_df, "gamma", prob_component_selection, 
#                                    prob_feature_selection, n_burnin)
#     setting <- setting + 1
#   }
# }

features_of_interest <- c(1, 2, 6, 7)
feature_names <- paste("j", features_of_interest, sep = "_")
# TODO make trace plot function to reduce copy paste going forward
Eta_df <- Eta_chain[1, features_of_interest, 1:n_iterations] %>%
  apply(MARGIN = 1, FUN = cummean) %>% 
  as.data.frame() %>% rename_at(vars(names(.)), ~ feature_names) %>% mutate(iteration = 1:n_iterations) %>% 
  gather(key = "Eta", value = "MPP", -iteration) %>%
  mutate(eta = Eta %>% as.factor()) 

Eta_1_plot <- ggplot(Eta_df, 
                     aes(x = iteration, y = MPP, color = eta)) + 
  geom_line() + geom_vline(xintercept = n_burnin, 
                           linetype = "dashed", color = "red") +
  labs(x = "iteration", y = "MPP", 
       title = "Trace plot for eta_1.")

#print(Eta_1_plot)

temp <- t(sigma2_chain)
temp <- cbind(temp, 1:2000)
ggplot(as.data.frame(temp[seq(1, 2000, 10),]), aes(x=V11)) + geom_line(aes(y=V7))
*/
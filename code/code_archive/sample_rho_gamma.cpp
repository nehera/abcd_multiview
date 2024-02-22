#include <cstdlib>
#include <iostream>
#include <string>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

#include "utils.h"

// [[Rcpp::export]]
arma::vec get_loggauss_quadform_j(arma::vec Xj, double s2j, arma::mat U, arma::vec Tauj, arma::vec Etaj) {
  int n = Xj.n_elem;
  int n_active_components = sum(Etaj);
  double Sigma_det = 0;
  mat Sigma_inv(n, n, fill::value(datum::nan));
  if (n_active_components > 0) {
    uvec active_components = get_indices_of_ones(Etaj);
    mat D(n_active_components, n_active_components, fill::value(0));
    D.diag() = Tauj.elem(active_components);
    U = U.cols(active_components);
    Sigma_det = pow(s2j, n)*det(inv(D) + U.t()*U)*det(D);
    Sigma_inv = (1.0/s2j)*get_woodbury_inv(U, D);
  } else {
    // Sigma2_j = I, an n-by-n identity matrix
    Sigma_det = pow(s2j, n);
    Sigma_inv = (1.0/s2j)*eye(n, n);
  }
  arma::rowvec transXj = arma::conv_to<arma::rowvec>::from(Xj);
  double quadformj = as_scalar(transXj*Sigma_inv*Xj); 
  double loggaussj = -0.5*n*log(2*datum::pi) - 0.5*log(Sigma_det) - 0.5*quadformj;
  arma::vec loggauss_quadform_j = join_cols(vec{loggaussj}, vec{quadformj});
  return loggauss_quadform_j;
}

// [[Rcpp::export]]
Rcpp::List sample_gamma_Eta(int r, int n, int IndVar, int p,
                      arma::uvec gamma, arma::mat Tau, arma::mat U, arma::mat X,
                      arma::vec q1, // probabilities of component selection
                      double q2, // prob_feature_selection
                      arma::vec s2, // vec of sigma_j's
                      arma::mat Eta) {

  // The core logic of the function remains largely the same as Samplerhogamma,
  // but with variable names updated to match your new naming conventions (rho -> gamma, Gam -> Eta).
  // Implement the logic as per your original function, ensuring that all references to rho and Gam
  // are updated to gamma and Eta, respectively.

  // Initialize data structures
  int l,j;
  arma::uvec gammanew = arma::zeros<arma::uvec>(r);
  arma::mat Etanew = arma::zeros<arma::mat>(r, p);
  arma::vec quadForm = arma::zeros(p);
  arma::vec loggauss = arma::zeros(p);
  for (l=0; l<r; l++){
    for (j=0; j<p; j++){
      Etanew(l,j) = Eta(l,j);
    }
    gammanew[l] = gamma[l];
  }

  // IndVar=1 denotes outcome view
  if (IndVar==1) {
    int j = 0; // only one outcome feature, so j fixed to 0
    for (l=0; l<r; l++){
      
      arma::vec logGausQuadFormj;
      
      Eta(l,0)=0;
      logGausQuadFormj = get_loggauss_quadform_j(X.col(j), s2[j], U, Tau.col(j), Eta.col(j));
      double quadj = logGausQuadFormj[0];
      double loggaussj = logGausQuadFormj[1];
      
      Etanew(l,0)=1;
      logGausQuadFormj = get_loggauss_quadform_j(X.col(j), s2[j], U, Tau.col(j), Etanew.col(j));
      double quadnewj = logGausQuadFormj[0];
      double loggaussnewj = logGausQuadFormj[1];
      
      double rat=loggaussnewj+log(q2)-loggaussj-log(1-q2);
      double uni=arma::randu();
      if (log(uni/(1-uni))<rat){
        Etanew(l,0) = Eta(l,0) = gamma[l] = 1;
        loggauss[0]=loggaussnewj;
        quadForm[0]=quadnewj;
      } else {
        Etanew(l,0) = Eta(l,0) = gamma[l] = 0;
        loggauss[0]=loggaussj;
        quadForm[0]=quadj;
      }
    }
  } else { 
    // Sample component and feature activation indicators for omics/ covariate views
    for (l=0; l<r; l++){
      
      double logq=0; // What is logq again? 
      
      if (gamma[l]==0){
        
        // Propose to turn the component on
        gammanew[l]=1;
        double logpostnew=0;
        double logpostold=0;
        
        arma::vec quadForm1 = arma::zeros(p);
        arma::vec loggausnew1 = arma::zeros(p);
        
        double quadj=0;
        double quadnewj=0;
        double loggaussj=0;
        double loggaussnewj=0;
        
        // Propose to turn on/ off features
        for (j=0;j<p;j++){

          double logqj=0;
          double logmqj=0;
          arma::vec logGausQuadFormj;
          
          logpostold+=loggauss[j];
          
          Eta(l,j)=0;
          logGausQuadFormj = get_loggauss_quadform_j(X.col(j), s2[j], U, Tau.col(j), Eta.col(j));
          double quadj = logGausQuadFormj[0];
          double loggaussj = logGausQuadFormj[1];
          
          Etanew(l,j)=1;
          logGausQuadFormj = get_loggauss_quadform_j(X.col(j), s2[j], U, Tau.col(j), Etanew.col(j));
          double quadnewj = logGausQuadFormj[0];
          double loggaussnewj = logGausQuadFormj[1];
    
          double rat=-loggaussj+loggaussnewj+log(q1[l])-log(1-q1[l]); // log(P_lj / Q_lj) i.e. log(G_1)-log(G_0)
          double x1=loggaussnewj+log(q1[l]); // log G_1
          double x2=loggaussj+log(1-q1[l]); // log G_0
          double maxq = std::max(x1,x2);
          double den=maxq+log(exp(x1-maxq)+exp(x2-maxq)); // log(G_0 + G_1)
          logqj=x1-den; // log P_lj
          logmqj=x2-den; // log Q_lj = log (1-P_lj)
          double uni=arma::randu();
          if ((log(uni/(1-uni))<rat)){
            // Turn the feature on
            Eta(l,j)=1;
            Etanew(l,j)=1;
            logpostnew+=loggaussnewj+log(q1[l]);
            loggausnew1[j]=loggaussnewj;
            // log proposal difference- add the probability that a feature is on given a component is on
            logq+=logqj;
            quadForm1[j]=quadnewj;
          } else {
            Eta(l,j)=0;
            Etanew(l,j)=0;
            logq+=logmqj;
            // log proposal difference- add the probability that a feature is off given a component is on
            logpostnew+=loggaussj+log(1-q1[l]); // loggauss doesn't change from the last iteration when we are proposing the component be on and are keeping the feature off.
            quadForm1[j]=quadj;
            loggausnew1[j]=loggaussj;
          }
        }
        // Add the log probability that the component is on
        logpostnew+=log(q2);
        logpostold+=log(1-q2);
        double un=arma::randu();
        double rat1=logpostnew-logpostold-logq; // log acceptance ratio
        if (log(un)<rat1){
          // accept having the component on
          gamma[l]=1;
          // store the new log gauss and quadForms
          for (j=0; j<p; j++){
            quadForm[j] = quadForm1[j];
            loggauss[j] = loggausnew1[j];
          }
        } else {
          // stay off
          gamma[l]=0;
          for (j=0; j<p; j++) {
            Eta(l,j) = Etanew(l,j)=0;
          }
        }
        // initialize gamma new for the next iteration & remove quadform and loggaussnewj since unneccessary for the next step?
        gammanew[l] = gamma[l];
      } else {
        // gamma is on and we are proposing to turn it off
        gammanew[l]=0;
        // initialize proposal data structures
        double logpostnew=0;
        double logpostold=0;
        arma::vec quadForm1 = arma::zeros(p);
        arma::vec loggausnew1 = arma::zeros(p);

        double quadj=0;

        double loggaussnewj=0;
        double loggaussn=0;
        double logpq=0;

        for (j=0;j<p;j++){
          
          arma::vec logGausQuadFormj;
          
          logpostold += loggauss[j];
          
          Etanew(l,j) = 1;
          logGausQuadFormj = get_loggauss_quadform_j(X.col(j), s2[j], U, Tau.col(j), Etanew.col(j));
          quadj = logGausQuadFormj[0];
          loggaussn = logGausQuadFormj[1];

          Etanew(l,j) = 0;
          logGausQuadFormj = get_loggauss_quadform_j(X.col(j), s2[j], U, Tau.col(j), Etanew.col(j));
          quadj = logGausQuadFormj[0]; // quadj redundantly calculated. it could be re-used from above. 
          loggaussnewj = logGausQuadFormj[1];
          
          if (p!=1){
            double x1=loggaussn+log(q1[l]);
            double x2=loggaussnewj+log(1-q1[l]);
            double maxq=std::max(x1, x2);
            double den=maxq+log(exp(x1-maxq)+exp(x2-maxq));
            double logqj=x1-den;
            double logmqj=x2-den;
            logq += Eta(l,j)*log(q1[l])+(1-Eta(l,j))*log(1-q1[l]);
            logpq += Eta(l,j)*logqj+(1-Eta(l,j))*logmqj;
          }
          logpostnew+=loggaussnewj;
          loggausnew1[j]=loggaussnewj;
          quadForm1[j]=quadj;
          Etanew(l,j)=Eta(l,j);
        }
        logpostnew+=log(1-q2)+logpq;
        logpostold+=log(q2)+logq;
        double un=arma::randu();
        double rat1=logpostnew-logpostold;
        if (log(un)<rat1){
          gamma[l]=0;
          for (j=0;j<p;j++){
            Etanew(l,j)=Eta(l,j)=0;
            quadForm[j]=quadForm1[j];
            loggauss[j]=loggausnew1[j];
          }
        } else {
          for (j=0;j<p;j++){
            Etanew(l,j)=Eta(l,j);
          }
        }
        gammanew[l]=gamma[l];
      }
      /* Within move*/
      if (gamma[l] == 1) {

        double logpostold=0;
        arma::vec quadForm1 = arma::zeros(p);
        arma::vec loggausnew1 = arma::zeros(p);
        double quadj=0;
        double quadnewj=0;
        double loggaussj=0;
        double loggaussnewj=0;
        
        arma::vec logGausQuadFormj;

        for (j=0; j<p; j++){
          logpostold+=loggauss[j];
          
          logGausQuadFormj = get_loggauss_quadform_j(X.col(j), s2[j], U, Tau.col(j), Eta.col(j));
          quadj = logGausQuadFormj[0];
          loggaussj = logGausQuadFormj[1];
          
          Etanew(l,j)=1-Eta(l,j);
          logGausQuadFormj = get_loggauss_quadform_j(X.col(j), s2[j], U, Tau.col(j), Etanew.col(j));
          quadj = logGausQuadFormj[0];
          loggaussnewj = logGausQuadFormj[1];
          
          double rat=-loggaussj+loggaussnewj;
          double uni=arma::randu();
          if (Eta(l,j)==0){
            quadj=quadnewj;
            loggaussj=loggaussnewj;
            rat+=log(q1[l])-log(1-q1[l]);
          } else {
            quadnewj=quadj;
            loggaussnewj=loggaussj;
            rat=-rat+log(q1[l])-log(1-q1[l]);
          }
          if (log(uni/(1-uni))<rat){
            Eta(l,j)=1;
            quadForm[j]=quadj;
            loggauss[j]=loggaussj;
          } else {
            Eta(l,j)=0;
            quadForm[j]=quadnewj;
            loggauss[j]=loggaussnewj;
          }
          Etanew(l,j)=Eta(l,j);
        }
      }
    }
  }

  // After implementing the logic, return the result as a List
  return List::create(Named("gamma") = gamma,
                      Named("Eta") = Eta,
                      Named("quadForm") = quadForm,
                      Named("loggauss") = loggauss);
}

/*** R
# R code to test the function
# Generating data
set.seed(1)
source("00_simulate_simple_data.R")
simulation_results <- simulate_iid_data()

# Defining parameters
m = 1 # start with the first view
r = 4
n = 200
p = 10
Tau = matrix(1, nrow = r, ncol = p) # Tau fixed
U = simulation_results$U
# X = simulation_results$X_list[[m]]
X <- simulation_results$Y
s2 = rep(1, p)
prob_component_selection = 0.5
prob_feature_selection = 0.5

# Initializing gamma and Eta
gamma <- rbinom(n = r, size = 1, prob = prob_component_selection)
Eta <- matrix(NA, r, p)
for (l in 1:r) {
  if (gamma[l] == 1) {
    Eta[l,] <- rbinom(n = p, size = 1, prob_feature_selection)
  } else {
    Eta[l,] <- rep(0, p_m)
  }
}

j = 1

get_loggauss_quadform_j(X[,j], s2[j], U, Tau[,j], Eta[,j])

n_iter <- 2000
gamma_chain <- matrix(NA, r, n_iter)
Eta_chain <- array(NA, dim = c(p, r, n_iter))
for (i in 1:n_iter) {
  samp <- sample_gamma_Eta(r, n, IndVar=1, p, gamma, Tau, U, X,
                   q1= rep(prob_component_selection, r), # probabilities of component selection
                   q2= prob_feature_selection, # prob_feature_selection
                   s2, Eta)
  gamma <- samp$gamma
  Eta <- samp$Eta
  gamma_chain[, i] <- gamma
  Eta_chain[,, i] <- Eta
}
apply(gamma_chain, 1, mean)
apply(Eta_chain, c(1,2), mean)
*/

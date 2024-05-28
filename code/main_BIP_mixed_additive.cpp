#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <time.h>

#include "header.h"
#include "utils.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

#include <RcppGSL.h>
// [[Rcpp::depends(RcppGSL)]]

// Subroutines for extension

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

///// myfunction.c subroutines

void SampleIntercept(gsl_rng * rr, int n, int r, double * intercept, double * sigma2, double sigma20, double ** U, double ** A, double **y){
  int i,l;
  double meany=0;
  for (i = 0; i < n; i++){
    double ua=0;
    for (l=0;l<r;l++){
      ua+=U[i][l]*A[l][0];
    }
    meany+=y[i][0]-ua;
  }
  meany=meany/n;
  double invsig2=n/sigma2[0]+1/sigma20;
  *intercept=(n/(invsig2*sigma2[0]))*meany+sqrt(1/invsig2)*gsl_ran_ugaussian(rr);
}

void logPostGam(double *logpo, arma::vec IndVar, int Np, int r, int n, arma::vec P, 
                double *** Tau, double ** U,double *** X, double **s2, bool ** rho, 
                bool *** Gam, double** qv, double* q) {
  double logpost=0;
  int m,j,l;
  for (m=0;m<Np;m++){
    for (j=0;j<P[m];j++){
      double logp;double quad;
      logGausQuadForm(j,r, n,P[m], Tau[m], U,X[m], s2[m][j],&quad,Gam[m][j],&logp);
      logpost+=logp;
      for (l=0;l<r;l++){
        if (IndVar[m]!=2)
          logpost+=Gam[m][j][l]*log(qv[m][l])+(1-Gam[m][j][l])*log(1-qv[m][l]);
      }
    }

    if (IndVar[m]==0){
      for (l=0;l<r;l++){
        logpost+=rho[m][l]*log(q[m])+(1-rho[m][l])*log(1-q[m]);
      }
    }
  }
  *logpo=logpost;
}

void logGausQuadForm(int j, int r, int n,int p,double ** Tau, double ** U, double ** X, double s2, double* quadForm, bool * Gam, double * loggauss){
  int i,s,s1;
  int NZ1[r];
  int nx1=0;
  findc(r,Gam,0,NZ1, &nx1);
  
  double result=0;double quadF=0;
  if (nx1>0){
    gsl_vector *work1 =gsl_vector_alloc (nx1);
    double * Sigma1= static_cast<double*>(malloc(nx1*nx1*sizeof(double)));
    for (s=0;s<nx1;s++){
      for (s1=0;s1<=s;s1++){
        double a=0;
        for (i=0;i<n;i++){
          //if (Gam[NZ1[s]]*Gam[NZ1[s1]]==1)
          a+=U[i][NZ1[s]]*U[i][NZ1[s1]];
        }   
        Sigma1[s*nx1+s1]=Sigma1[s1*nx1+s]=a;
        
      }
      Sigma1[s*nx1+s]+=(1.0/Tau[NZ1[s]][j]);
      //printf("q1=%lf\n",Tau[NZ1[s]][j]);
    }
    gsl_matrix_view m1  = gsl_matrix_view_array (Sigma1, nx1,nx1);
    //printf("Error on Gam\n");
    
    gsl_linalg_cholesky_decomp (&m1.matrix);
    gsl_vector *xi =gsl_vector_alloc (n);
    for (i = 0; i < n; i++) {
      gsl_vector_set (xi, i, X[i][j]/sqrt(s2));
    }
    
    gsl_vector *y =gsl_vector_alloc (nx1);
    double sumT=0;
    for (s=0;s<nx1;s++){
      double a=0;                            
      //if (Gam[NZ1[s]]==1){
      sumT+=log(Tau[NZ1[s]][j]);
      for (i=0;i<n;i++){
        a+=U[i][NZ1[s]]*X[i][j]/sqrt(s2);
      }
      //}
      gsl_vector_set (y, s, a);
    }
    logmultigaussianT(xi, y,  &m1.matrix,&result,&quadF, work1);
    result=result-0.5*sumT-(n/2.0)*log(s2);
    gsl_vector_free (y);gsl_vector_free (work1);gsl_vector_free (xi);
    free(Sigma1); 
  } else {
    for (i = 0; i < n; i++){
      quadF+=pow(X[i][j],2);
    }
    quadF=quadF/s2;
    result=-(n/2.0)*log(s2)- 0.5*n*log(2.0*M_PI)-0.5*quadF;
  }
  
  *quadForm=quadF;
  *loggauss=result;
  
}

void sampleGam(gsl_rng * rr, int r, int n,int p, bool * rho,double ** Tau, double ** U,double ** X,double *s2,double* quadForm,bool ** Gam,double * loggauss,double q){
  int j,s,l;
  for (j=0;j<p;j++){
    for (l=0;l<r;l++) Gam[j][l]=rho[l]*Gam[j][l];
  }
  
  double phi=0.5;
  int NZ1[r];
  int nx1=0;
  findc(r,rho,0,NZ1, &nx1);
  if (nx1>0){
    for (j=0;j<p;j++){
      bool Gamnew2[r];
      for (l=0;l<r;l++) Gamnew2[l]=0;
      double logpostold=0;double logpostnew=0;
      bool Gamold1[nx1];
      for (s=0;s<nx1;s++){
        Gamold1[s]=Gam[j][NZ1[s]];
      }
      bool Gamnew1[nx1];
      proposal(nx1,Gamold1,Gamnew1,phi, rr);
      for (s=0;s<nx1;s++){
        Gamnew2[NZ1[s]]=Gamnew1[s];
        logpostold+=Gamold1[s]*log(q)+(1-Gamold1[s])*log(1-q);
        logpostnew+=Gamnew1[s]*log(q)+(1-Gamnew1[s])*log(1-q);
      }
      
      /*Logpost new*/
      double loggaussnew=0;double quadForm1=0;
      logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quadForm1,Gamnew2,&loggaussnew);
      logpostnew+=loggaussnew;
      logpostold+=loggauss[j];
      double u=gsl_ran_flat (rr, 0, 1);
      if ((log(u)<logpostnew-logpostold)){
        for (l=0;l<r;l++){
          Gam[j][l]=Gamnew2[l];
        }
        loggauss[j]=loggaussnew;
        quadForm[j]=quadForm1;
      }
    }
  }
}

double  logpost(int r, int n,int p, bool * rho,double ** Tau, double ** U,double ** X, double q,double* s2,double* quadForm,bool ** Gam,double * loggauss){
  double logpostrho=0;
  int l,j;
  for (j=0;j<p;j++){
    logGausQuadForm(j,r, n,p,Tau, U,X,s2[j],&quadForm[j],Gam[j],&loggauss[j]);
    logpostrho+= loggauss[j];
  }
  for (l=0;l<r;l++){
    logpostrho+=rho[l]*log(q)+(1-rho[l])*log(1-q);
  }
  return logpostrho;
}

void rho1(gsl_rng * rr,int r, int n,int p, bool * rho,double ** Tau, double ** U,double ** X, double q,double* s2,double* quadForm,bool** Gam,double *loggauss){
  int l,j;
  double logpostold=0;
  double logprior=0;
  for (l=0;l<r;l++){
    logprior+=rho[l]*log(q)+(1-rho[l])*log(1-q);
  }
  for (j=0;j<p;j++){
    logpostold+=loggauss[j];
  }
  logpostold+=logprior;
  double * quadForm1= static_cast<double*>(malloc(p*sizeof(double)));
  double * loggauss1= static_cast<double*>(malloc(p*sizeof(double)));

  double phi=0.5;
  bool rhonew[r];
  proposal(r,rho,rhonew,phi, rr);
  double logpostnew=logpost(r, n,p, rhonew,Tau,U,X, q,s2,quadForm1,Gam,loggauss1);
  double u=gsl_ran_flat (rr, 0, 1);
  if (log(u)<logpostnew-logpostold) {
    for (l=0;l<r;l++){
      rho[l]=rhonew[l];
    }
    for (j=0;j<p;j++){
      quadForm[j]=quadForm1[j];
      loggauss[j]=loggauss1[j];
    }
  } 
}

void SamplerhoGamma(gsl_rng * rr,int r, int n,int IndVar,int p, bool * rho,double ** Tau, double ** U,double ** X, double* q1,double q2,double* s2,double* quadForm,bool** Gam,double *loggauss){
  int l,j;
  bool *rhonew= static_cast<bool*>(malloc(r*sizeof(bool)));
  bool Gamnew[p][r];
  for (l=0;l<r;l++){
    for (j=0;j<p;j++){
      Gamnew[j][l]=Gam[j][l];
    } 
    rhonew[l]=rho[l];
  }
  
  if (IndVar==1){
    for (l=0;l<r;l++){
      double quad1=0;double quad2=0;
      double loggaussold=0; double loggaussnew=0;
      Gam[0][l]=0;
      logGausQuadForm(0,r, n,p, Tau,  U,X,s2[0],&quad1,Gam[0],&loggaussold);
      Gamnew[0][l]=1;
      logGausQuadForm(0,r, n,p,Tau,  U,X,s2[0],&quad2,Gamnew[0],&loggaussnew);
      double rat=loggaussnew+log(q2)-loggaussold-log(1-q2);
      double uni=gsl_ran_flat (rr, 0, 1);
      if (log(uni/(1-uni))<rat){
        Gamnew[0][l]=Gam[0][l]=rho[l]=1;
        loggauss[0]=loggaussnew;
        quadForm[0]=quad2;
      } else {
        Gamnew[0][l]=Gam[0][l]=rho[l]=0;
        loggauss[0]=loggaussold;
        quadForm[0]=quad1;
      }
    }
    
  } else {
    for (l=0;l<r;l++){
      double logq=0;
      if (rho[l]==0){
        // Propose to turn the component on
        rhonew[l]=1;
        double logpostnew=0;
        double logpostold=0;
        double * quadForm1= static_cast<double*>(malloc(p*sizeof(double)));
        double * loggausnew1= static_cast<double*>(malloc(p*sizeof(double)));
        double quad1=0; double quad2=0;
        double loggaussold=0; double loggaussnew=0;
        // Propose to turn on/ off features
        for (j=0;j<p;j++){
          double logqj=0;double logmqj=0;
          logpostold+=loggauss[j];
          Gam[j][l]=0;
          logGausQuadForm(j,r, n,p, Tau,  U,X,s2[j],&quad1,Gam[j],&loggaussold);
          Gamnew[j][l]=1;
          logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quad2,Gamnew[j],&loggaussnew);
          double rat=-loggaussold+loggaussnew+log(q1[l])-log(1-q1[l]); // log(P_lj / Q_lj)

          double x1=loggaussnew+log(q1[l]); // log G_1
          double x2=loggaussold+log(1-q1[l]); // log G_0
          double maxq=std::max(x1,x2);
          double den=maxq+log(exp(x1-maxq)+exp(x2-maxq)); // log(G_0 + G_1)
          logqj=x1-den; // log P_lj
          logmqj=x2-den; // log Q_lj = 1 - log P_lj
          
          double uni=gsl_ran_flat (rr, 0, 1);
          if ((log(uni/(1-uni))<rat)){
            // Turn the feature on
            Gam[j][l]=1;Gamnew[j][l]=1;
            logpostnew+=loggaussnew+log(q1[l]);
            loggausnew1[j]=loggaussnew;
            // log proposal difference- add the probability that a feature is on given a component is on
            logq+=logqj; 
            quadForm1[j]=quad2;
          } else {
            Gam[j][l]=0;Gamnew[j][l]=0;
            logq+=logmqj; 
            // log proposal difference- add the probability that a feature is off given a component is on
            logpostnew+=loggaussold+log(1-q1[l]);
            quadForm1[j]=quad1;
            loggausnew1[j]=loggaussold;
          }
        }
        // Add the log probability that the component is on
        logpostnew+=log(q2);
        logpostold+=log(1-q2);
        double un=gsl_ran_flat (rr, 0, 1);
        double rat1=logpostnew-logpostold-logq; // log acceptance ratio
        if (log(un)<rat1){
          // accept having the component on
          rho[l]=1;
          // store the new log gauss and quadForms
          for (j=0;j<p;j++){
            quadForm[j]=quadForm1[j];
            loggauss[j]=loggausnew1[j];
          }
        } else {
          // stay off
          rho[l]=0;
          for (j=0;j<p;j++) Gam[j][l]=Gamnew[j][l]=0;
        }
        // initialize gamma new for the next iteration & remove quadform and loggaussnew since unneccessary for the next step? 
        rhonew[l]=rho[l];
        free(quadForm1);free(loggausnew1);
      } else {
        // gamma is on and we are proposing to turn it off
        rhonew[l]=0;
        // initialize proposal data structures
        double logpostnew=0;
        double logpostold=0;
        double * quadForm1= static_cast<double*>(malloc(p*sizeof(double)));
        double * loggausnew1= static_cast<double*>(malloc(p*sizeof(double)));
        double quad1;
        double loggaussnew=0; double loggaussn=0;
        double logpq=0;
        for (j=0;j<p;j++){
          logpostold+=loggauss[j];
          Gamnew[j][l]=1;
          logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quad1,Gamnew[j],&loggaussn);
          Gamnew[j][l]=0;
          logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quad1,Gamnew[j],&loggaussnew);
          if (p!=1){
            double x1=loggaussn+log(q1[l]);
            double x2=loggaussnew+log(1-q1[l]);
            double maxq=std::max(x1,x2);
            double den=maxq+log(exp(x1-maxq)+exp(x2-maxq));
            double logqj=x1-den;
            double logmqj=x2-den;
            logq+=Gam[j][l]*log(q1[l])+(1-Gam[j][l])*log(1-q1[l]);
            logpq+=Gam[j][l]*logqj+(1-Gam[j][l])*logmqj;
          }
          logpostnew+=loggaussnew;
          loggausnew1[j]=loggaussnew;
          quadForm1[j]=quad1;
          Gamnew[j][l]=Gam[j][l];
        }
        logpostnew+=log(1-q2)+logpq;
        logpostold+=log(q2)+logq;
        double un=gsl_ran_flat (rr, 0, 1);
        double rat1=logpostnew-logpostold;
        if (log(un)<rat1){
          rho[l]=0;
          for (j=0;j<p;j++){
            Gamnew[j][l]=Gam[j][l]=0;
            quadForm[j]=quadForm1[j];
            loggauss[j]=loggausnew1[j];
          }
        } else {
          for (j=0;j<p;j++){
            Gamnew[j][l]=Gam[j][l];
          }
        }
        free(quadForm1);free(loggausnew1);
        rhonew[l]=rho[l];
      }
      /* Within move*/
      if (rho[l]==1){
        double logpostold=0;
        double * quadForm1= static_cast<double*>(malloc(p*sizeof(double)));
        double * loggausnew1= static_cast<double*>(malloc(p*sizeof(double)));
        double quad1,quad2;
        double loggaussold,loggaussnew;
        for (j=0;j<p;j++){
          logpostold+=loggauss[j];
          logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quad1,Gam[j],&loggaussold);
          Gamnew[j][l]=1-Gam[j][l];
          logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quad2,Gamnew[j],&loggaussnew);
          double rat=-loggaussold+loggaussnew;
          double uni=gsl_ran_flat (rr, 0, 1);
          if (Gam[j][l]==0){
            quad1=quad2;
            loggaussold=loggaussnew;
            rat+=log(q1[l])-log(1-q1[l]);
          } else {
            quad2=quad1;
            loggaussnew=loggaussold;
            rat=-rat+log(q1[l])-log(1-q1[l]);
          }
          if (log(uni/(1-uni))<rat){
            Gam[j][l]=1;
            quadForm[j]=quad1;
            loggauss[j]=loggaussold;
          } else {
            Gam[j][l]=0;
            quadForm[j]=quad2;
            loggauss[j]=loggaussnew;
          }
          Gamnew[j][l]=Gam[j][l];
        }
        free(quadForm1);free(loggausnew1);
      }
    }
  }
  free(rhonew); 
}

void sigma2(gsl_rng * rr, int p,int n,double a0, double b0, double* quadForm,double * s2){
  int j;
  for (j=0;j<p;j++){
    double inv=1/(b0+0.5*s2[j]*quadForm[j]);
    s2[j]=1/gsl_ran_gamma (rr, a0+n/2.0, inv);
  }
}

void EffectZero(int l,gsl_rng * rr,int K,int p,bool rho, bool * R,double* B,double *B0,bool ** Path,double alphab0, double betab0,double alpha, double * lambda2,double *AcceptB0,bool ** Gam){
  double alphasb0=1;double betasb0=1;
  int j,k1;
  double B0new=gsl_ran_gamma (rr, alphasb0, 1/betasb0);;
  double lograt=0;
  double kl=0;
  for (j=0;j<p;j++){
    double kb=0;double kb1=0;
    if (Gam[j][l]==1){
      for (k1=0;k1<K;k1++){
        kb+=Path[j][k1]*B[k1];
      }
      kl+=lambda2[j];
      kb=kb+*B0;
      kb1=kb+B0new;
      lograt+=alpha*(log(kb1)-log(kb));
    }
  }
  if (rho==1){
    lograt+=(betasb0-betab0-kl)*(B0new-*B0)+(-alphasb0+alphab0)*(log(B0new)-log(*B0));
    double uni=gsl_ran_flat (rr, 0, 1);
    if (log(uni)<lograt){
      *B0=B0new;
      *AcceptB0+=1;
    } 
  } else {
    *B0=gsl_ran_gamma (rr, alphab0, 1/betab0);
  } 
  
}

void TauLambda(int l,gsl_rng * rr,int K,int p,double *A,double* B,double B0,double * Tau,bool ** Path,double alpha, double * lambda2,double *s2,bool ** Gam){
  /* K is the number of pathways
   *  * R is the binary vector for pathway effect
   *   * Path is the grouping information; it is a binary matrix
   *    * B is the pathway effect
   *     * alphab and betab are priors for B effects
   *      * l is the component index
   *       */
  int j,k1;
  for (j=0;j<p;j++){

    if ((Gam[j][l]==1)&& (fabs(A[j])>0)){
      double mu=sqrt(2*lambda2[j]*s2[j]/pow(A[j],2.0));
      Tau[j]=1.0/inverseGaussian(rr,  mu, 2*lambda2[j]);
      if (Tau[j]<0)

        if (((Tau[j]-Tau[j])!=0)|(Tau[j]<0)){
          Tau[j]=1/mu+1/(2*lambda2[j]);
        }
        double kb=0;
        for (k1=0;k1<K;k1++){
          kb+=Path[j][k1]*B[k1];
        }
        lambda2[j]=gsl_ran_gamma (rr, alpha+1, 1/(B0+kb+Tau[j]));
    } else {
      Tau[j]=gsl_ran_exponential (rr, 1.0/lambda2[j]);
      lambda2[j]=gsl_ran_gamma (rr, alpha, 1/B0);
    }
  }
}

void GroupEffect(int l,gsl_rng * rr,bool rho,int K,int p, bool * R,double *A,double* B,double B0,double * Tau,bool ** Path,double alphab, double betab,double w,double alpha, double * lambda2,double *s2,double * AcceptR,bool **Gam){
  /* K is the number of pathways
   * R is the binary vector for pathway effect
   * Path is the grouping information; it is a binary matrix
   * B is the pathway effect
   * alphab and betab are priors for B effects
   * l is the component index
   */
  int j, k,k1;
  
  if (rho==0){
    for (k=0;k<K;k++){
      B[k]=0;R[k]=0;
    }
  } else {
    
    double alphas=2;double betas=2; // proposal parameters 
    for (k=0;k<K;k++){
      bool Rknew=0;double Bknew=0;
      /* Between model move*/
      if (R[k]==1){
        Rknew=0; Bknew=0;
        double lograt=0;
        double kl=0;
        for (j=0;j<p;j++){
          if (Gam[j][l]==1){
            double kb=0;double kb1=0;
            for (k1=0;k1<K;k1++){
              kb+=Path[j][k1]*B[k1];
            }
            kb1=kb-Path[j][k]*B[k];
            lograt+=alpha*(log(B0+kb1)-log(B0+kb));
            kl+=Path[j][k]*lambda2[j];
          }
        }
        lograt+=(betab+kl-betas)*B[k]+(alphas-alphab)*log(B[k])+gsl_sf_lngamma(alphab)-gsl_sf_lngamma(alphas)+alphas*log(betas)-alphab*log(betab)+log(1-w)-log(w);
        double uni=gsl_ran_flat (rr, 0, 1);
        if (log(uni)<lograt){
          B[k]=Bknew; R[k]=Rknew;
          AcceptR[k]+=1;
        }
      }
      else if (R[k]==0){
        Rknew=1; Bknew=gsl_ran_gamma (rr, alphas, 1/betas);;
        double lograt=0;
        double kl=0;
        for (j=0;j<p;j++){
          if (Gam[j][l]==1){
            double kb=0;double kb1=0;
            for (k1=0;k1<K;k1++){
              kb+=Path[j][k1]*B[k1];
            }
            kb1=kb-Path[j][k]*B[k];
            kb=kb1+Path[j][k]*Bknew;
            lograt+=alpha*(log(B0+kb)-log(B0+kb1));
            kl+=Path[j][k]*lambda2[j];
          }
        }
        lograt+=(betas-betab-kl)*Bknew+(-alphas+alphab)*log(Bknew)-gsl_sf_lngamma(alphab)+gsl_sf_lngamma(alphas)-alphas*log(betas)+alphab*log(betab)-log(1-w)+log(w);
        double uni=gsl_ran_flat (rr, 0, 1);
        if (log(uni)<lograt){
          B[k]=Bknew; R[k]=Rknew;
          AcceptR[k]+=1;
        }
      }
      /* Within model move*/
      
      if (R[k]==1){
        Bknew=gsl_ran_gamma (rr, alphas, 1/betas);;
        double lograt=0;
        double kl=0;
        for (j=0;j<p;j++){ 
          if (Gam[j][l]==1){
            double kb=0;double kb1=0;
            for (k1=0;k1<K;k1++){
              kb+=Path[j][k1]*B[k1];
            }
            kb1=kb-Path[j][k]*B[k]+Path[j][k]*Bknew;
            lograt+=alpha*(log(B0+kb1)-log(B0+kb));
            kl+=Path[j][k]*lambda2[j];
          }
        }
        lograt+=(betas-betab-kl)*(Bknew-B[k])+(-alphas+alphab)*(log(Bknew)-log(B[k]));
        double uni=gsl_ran_flat (rr, 0, 1);
        if (log(uni)<lograt){
          B[k]=Bknew;
        }
      }
    } // close for k
  } // End for else if (rho==0)
}

void LoadAOther(gsl_rng * rr,int r, int n,int p, bool * rho,double ** Tau,double ** A, double ** U,double ** X,double* s2,bool ** Gam){
  int s,s1,j,i;
  for (j=0;j<p;j++){
    for (s=0;s<r;s++){
      A[s][j]=0;
    }
    int NZ1[r];
    int nx1=0;
    findc(r,Gam[j],0,NZ1, &nx1);
    if (nx1>0){
      double * SigmaInv= static_cast<double*>(malloc(nx1*nx1*sizeof(double)));
      for (s=0;s<nx1;s++){
        for (s1=0;s1<=s;s1++){
          double a=0;
          for (i=0;i<n;i++){
            a+=U[i][NZ1[s]]*U[i][NZ1[s1]];
          }
          SigmaInv[s*nx1+s1]=SigmaInv[s1*nx1+s]=a/s2[j];
        }
        SigmaInv[s*nx1+s]+=1/(Tau[s][j]*s2[j]);
      }
      gsl_matrix_view m  = gsl_matrix_view_array (SigmaInv, nx1,nx1);
      gsl_linalg_cholesky_decomp (&m.matrix);
      gsl_vector *mu =gsl_vector_alloc (nx1);
      double xy[nx1];
      for (s=0; s<nx1;s++){
        xy[s]=0;
        for (i=0; i<n;i++){
          xy[s]+=U[i][NZ1[s]]*X[i][j]/s2[j];
        }
      }
      gsl_vector_view b= gsl_vector_view_array (xy, nx1);
      gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
      gsl_vector *result=gsl_vector_alloc (nx1);
      multivariate_gaussian (rr, mu, &m.matrix,result);
      for (s=0; s<nx1;s++){
        A[NZ1[s]][j]=gsl_vector_get(result,s);
      }
      free(SigmaInv); gsl_vector_free (mu);gsl_vector_free (result);
    }
  }
}

void EstimateLoad(gsl_rng * rr,int r, int n,int p, bool * rho,double ** Tau,double ** A, double ** U,double ** X,double* s2,bool ** Gam){
  int s,s1,j,i;
  for (j=0;j<p;j++){
    for (s=0;s<r;s++){
      A[s][j]=0;
    }
    int NZ1[r];
    int nx1=0;
    findc(r,Gam[j],0,NZ1, &nx1);
    if (nx1>0){
      double * SigmaInv= static_cast<double*>(malloc(nx1*nx1*sizeof(double)));
      for (s=0;s<nx1;s++){
        for (s1=0;s1<=s;s1++){
          double a=0;
          for (i=0;i<n;i++){
            a+=U[i][NZ1[s]]*U[i][NZ1[s1]];
          }
          SigmaInv[s*nx1+s1]=SigmaInv[s1*nx1+s]=a/s2[j];
        }
        SigmaInv[s*nx1+s]+=1/(Tau[s][j]*s2[j]);
      }
      gsl_matrix_view m  = gsl_matrix_view_array (SigmaInv, nx1,nx1);
      TRY
      {
        gsl_linalg_cholesky_decomp (&m.matrix);
      }
      CATCH
      {
        printf("\nError on the sampling of A.");
      }
      ETRY;
      gsl_vector *mu =gsl_vector_alloc (nx1);
      double xy[nx1];
      for (s=0; s<nx1;s++){
        xy[s]=0;
        for (i=0; i<n;i++){
          xy[s]+=U[i][NZ1[s]]*X[i][j]/s2[j];
        }
      }
      gsl_vector_view b= gsl_vector_view_array (xy, nx1);
      gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);

      for (s=0; s<nx1;s++){
        A[NZ1[s]][j]=gsl_vector_get(mu,s);
      }
      free(SigmaInv); gsl_vector_free (mu);
\
    }
  }
}

void proposal(int n,bool *R,bool *R1,float phi, gsl_rng * r) {
  int i;
  for (i=0; i<n;i++) R1[i]=R[i];
  
  int n1=0;//number of 1
  int n0=0;//number of zeros
  int v1[n];
  findc(n,R,1,v1,&n0);//find indices different from 1 i.e ==0;
  int v2[n-n0];
  findc(n,R,0,v2,&n1);//find indices different of zeros i.e ==1
  double u=gsl_ran_flat (r, 0, 1);

  if ((u < phi) || (n0 == 0) || (n1 == 0)) {
    int l= gsl_rng_uniform_int(r,n);
    
    R1[l] = 1 - R[l];
  } else {
    int l1=gsl_rng_uniform_int(r,n0);
    int l2=gsl_rng_uniform_int(r,n1);
    
    R1[v1[l1]] = R[v2[l2]];
    R1[v2[l2]] = R[v1[l1]];
  }
}

void SampleUU(gsl_rng * rr,int r, int n, int Np, arma::vec P, double *** A, 
              double ** U, double *** X, double** s2){
  
  int s,s1,j,i;
  int sumMark=0;
  for (i=0;i<Np;i++) sumMark+=P[i]; 
  
  double * SigmaInv= static_cast<double*>(malloc(r*r*sizeof(double)));
  for (s=0;s<r;s++){
    for (s1=0;s1<=s;s1++){
      double a=0;
      for (j=0;j<sumMark;j++){
        int k=0;
        for (i=0;i<Np;i++){
          if ((k<=j) && (j<k+P[i])) a+=A[i][s][j-k]*A[i][s1][j-k]/s2[i][j-k]; 
          k+=P[i];
        }
      }
      SigmaInv[s*r+s1]=SigmaInv[s1*r+s]=a;
    }
    SigmaInv[s*r+s]+=1;
  }
  
  gsl_matrix_view m  = gsl_matrix_view_array (SigmaInv, r,r);

  TRY
  {

    gsl_linalg_cholesky_decomp (&m.matrix);

  }
  CATCH
  {
    printf("\nError on the sampling of U.");
  }
  ETRY;
  int l;
  for (i=0;i<n;i++){
    gsl_vector *mu =gsl_vector_alloc (r);
    double Ax[r];
    for (s=0; s<r;s++){
      Ax[s]=0;
      for (j=0;j<sumMark;j++){
        int k=0;
        for (l=0;l<Np;l++){
          if ((k<=j) && (j<k+P[l])) Ax[s]+=X[l][i][j-k]*A[l][s][j-k]/s2[l][j-k]; 
          k+=P[l];
        }
      }
    }
    gsl_vector_view b= gsl_vector_view_array (Ax, r);
    gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
    gsl_vector *result=gsl_vector_alloc (r);
    multivariate_gaussian (rr, mu, &m.matrix,result);
    for (s=0; s<r;s++){
      U[i][s]=gsl_vector_get(result,s);
    }
    gsl_vector_free (mu);gsl_vector_free (result);
  }
  free(SigmaInv);
}

void SampleU(gsl_rng * rr,int r, int n,int p0,int p1,int p2,double *** A, double ** U,double *** X,double** s2){
  int s,s1,j,i;
  double * SigmaInv= static_cast<double*>(malloc(r*r*sizeof(double)));
  for (s=0;s<r;s++){
    for (s1=0;s1<=s;s1++){
      double a=0;
      for (j=0;j<p0+p1+p2;j++){
        if (j<p0) a+=A[0][s][j]*A[0][s1][j]/s2[0][j];
        else if (j<p0+p1) a+=A[1][s][j-p0]*A[1][s1][j-p0]/s2[1][j-p0];
        else a+=A[2][s][j-p1-p0]*A[2][s1][j-p1-p0]/s2[2][j-p1-p0];
      }
      SigmaInv[s*r+s1]=SigmaInv[s1*r+s]=a;
    }
    SigmaInv[s*r+s]+=1;
  }
  gsl_matrix_view m  = gsl_matrix_view_array (SigmaInv, r,r);
  TRY
  {
    gsl_linalg_cholesky_decomp (&m.matrix);
  }
  CATCH
  {
    printf("\nError on the sampling of U.");
  }
  ETRY;
  
  for (i=0;i<n;i++){
    gsl_vector *mu =gsl_vector_alloc (r);
    double Ax[r];
    for (s=0; s<r;s++){
      Ax[s]=0;
      for (j=0; j<p0+p1+p2;j++){
        if (j<p0)
          Ax[s]+=X[0][i][j]*A[0][s][j]/s2[0][j];
        else if (j<p0+p1)
          Ax[s]+=X[1][i][j-p0]*A[1][s][j-p0]/s2[1][j-p0];
        else Ax[s]+=X[2][i][j-p1-p0]*A[2][s][j-p1-p0]/s2[2][j-p1-p0];
      }
    }
    gsl_vector_view b= gsl_vector_view_array (Ax, r);
    gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
    gsl_vector *result=gsl_vector_alloc (r);
    multivariate_gaussian (rr, mu, &m.matrix,result);
    for (s=0; s<r;s++){
      U[i][s]=gsl_vector_get(result,s);
    }
    gsl_vector_free (mu);gsl_vector_free (result);
  }
  free(SigmaInv);
}

///// BIP.c main function

// [[Rcpp::export]]
Rcpp::List mainfunction(int Method, int n, arma::vec P, int r, int Np, arma::vec datasets,
                  arma::vec IndVar, arma::vec K, arma::vec Paths, int maxmodel, int nbrsample, // nbrsample=n_iter-n_burnin
                  int burninsample, arma::vec CompoSelMean, arma::vec VarSelMean,
                  arma::vec VarSelMeanGlobal, arma::vec GrpSelMean, arma::vec GrpEffectMean,
                  arma::vec IntGrpMean, arma::vec EstU, arma::vec EstSig2, double InterceptMean, arma::vec EstLoadMod,
                  arma::vec EstLoad, int nbrmodel1, arma::vec postgam, arma::vec priorcompsel, arma::vec priorcompselo,
                  arma::vec priorb0, arma::vec priorb, arma::vec priorgrpsel, double probvarsel,
                  // Arguments below are from gibbs_sampler_nested.cpp
                  const vec& y, const mat& W, uvec study_site, uvec family, int n_iter, 
                  const vec& beta_prior_mean, const vec& beta_prior_var, double sigma_ksi_prior_a, 
                  double sigma_ksi_prior_b, double sigma_theta_prior_a, double sigma_theta_prior_b, 
                  double sigma_prior_a, double sigma_prior_b) {
  
  // Set GSL random number generator seed
  long seed=1;
  gsl_rng * rr = gsl_rng_alloc (gsl_rng_rand48);
  gsl_rng_set (rr, seed);
  
  // Adjust for 0-based indexing in C++
  study_site -= 1;
  family -= 1;
  
  // int n = y.n_elem;
  int S_count = max(study_site) + 1;
  int F_count = max(family) + 1;
  int n_beta = W.n_cols;
  
  // Initialize parameters
  vec beta = zeros<vec>(n_beta);
  vec ksi = zeros<vec>(S_count);
  vec theta = zeros<vec>(F_count);
  double sigma_ksi2 = 1.0;
  vec sigma_theta2 = ones<vec>(S_count);
  double sigma2_0 = 1.0;
  
  // Storage for samples
  mat beta_samples(n_iter, n_beta, fill::zeros);
  mat ksi_samples(n_iter, S_count, fill::zeros);
  mat theta_samples(n_iter, F_count, fill::zeros);
  vec sigma_ksi2_samples(n_iter, fill::zeros);
  mat sigma_theta2_samples(n_iter, S_count, fill::zeros);
  vec sigma2_samples(n_iter, fill::zeros);
  
  // Initialize residualization intermediate
  mat y_tilde(n, 1, fill::zeros);
  
  setvbuf(stdout, NULL, _IONBF, 0);

  if (Method==1)
    printf("\nThe Method is GroupInfo\n");
  else if (Method==0) {
    printf("The Method is NoGroupInfo\n");
  }
  
  int i,l,j,k;
  int m;
  clock_t t1 = clock();

  printf("Number of MCMC samples after burn-in is %d\n",nbrsample);

  printf("Number of burn-in is %d\n",burninsample);

  printf("Number of samples is %d\n",n);

  for (m=0;m<Np;m++){
    int n_markers = P[m];
    printf("Number of markers in platform %d is %d\n", m, n_markers);
  }
  
  printf("Number of components is %d\n", r);

  double *** X= static_cast<double***>(malloc(Np*sizeof(double **)));
  double *** X1= static_cast<double***>(malloc(Np*sizeof(double **)));

  k=0;
  
  for (m=0;m<Np;m++){
    X[m]=dmatrix(0,n-1,0,P[m]-1);
    X1[m]=dmatrix(0,n-1,0,P[m]-1);
    for (i=0;i<n;i++){
      for (j=0;j<P[m];j++){
        X[m][i][j]=X1[m][i][j]=datasets[k];
        k++;
      }
    }
  }
  
  bool *** Path= static_cast<bool***>(malloc(Np*sizeof(bool **)));
  int m1=0; //bool pp=0;
  int kk=0;
  for (m=0;m<Np;m++){
    Path[m]=bmatrix(0,P[m]-1,0,K[m]-1);
    if (IndVar[m]==1){
      Path[m][0][0]=1;
    } else if (IndVar[m]==0) {
      for (j=0;j<P[m];j++){
        for (k=0;k<K[m];k++){
          Path[m][j][k]=Paths[kk];
          kk++;
        }
      }
      m1+=1;
    } else if (IndVar[m]==2){//covariates
      for (j=0;j<P[m];j++){
        for (k=0;k<K[m];k++){
          Path[m][j][k]=1;
        }
      }
    } 
  }
  

  double ** U=dmatrix(0,n-1,0,r-1);
  double ** meanU=dmatrix(0,n-1,0,r-1);
  for (i=0;i<n;i++){
    for (l=0;l<r;l++){
      U[i][l]=gsl_ran_ugaussian(rr);
      meanU[i][l]=0;
    }
  }
  
  double ***A= static_cast<double***>(malloc(Np*sizeof(double **)));
  for (m=0;m<Np;m++){
    A[m]=dmatrix(0,r-1,0,P[m]-1);
  }
  
  /* Hyperparameter*/
  double * q= static_cast<double*>(malloc(Np*sizeof(double)));
  double ** qv= static_cast<double**>(malloc(Np*sizeof(double*)));
  double ** wg= static_cast<double**>(malloc(Np*sizeof(double*)));
  for (m=0;m<Np;m++){
    q[m]=0.5;
    qv[m]= static_cast<double*>(malloc(r*sizeof(double)));
    wg[m]= static_cast<double*>(malloc(r*sizeof(double))); //proba. for group selection
  }

  // Define hyperparameters
  double a0=1; double b0=1;
  double alphab=priorb[0]; double betab=priorb[1];
  double al=priorcompsel[0]; double bl=priorcompsel[1]; // Hyper for q
  double al0=priorcompselo[0]; double bl0=priorcompselo[1]; // Hyper for q in the outcome
  double alphab0=priorb0[0]; double betab0=priorb0[1]; //Hyper for b0
  double alg=priorgrpsel[0]; double blg=priorgrpsel[1];
  double alpha=1;
  
  // Initialize parameters to estimate
  bool ** rhoest= static_cast<bool**>(malloc(Np*sizeof(bool*)));
  double ** rhomean= static_cast<double**>(malloc(Np*sizeof(double*)));
  double ** Gamvs= static_cast<double**>(malloc(Np*sizeof(double*)));
  bool *** R= static_cast<bool***>(malloc(Np*sizeof(bool**)));
  bool *** Gam= static_cast<bool***>(malloc(Np*sizeof(bool**)));
  double *** Gammean= static_cast<double***>(malloc(Np*sizeof(double**)));
  double *** AcceptR= static_cast<double***>(malloc(Np*sizeof(double**)));
  double *** Bmean= static_cast<double***>(malloc(Np*sizeof(double**)));
  double *** Rmean= static_cast<double***>(malloc(Np*sizeof(double**)));
  double *** B= static_cast<double***>(malloc(Np*sizeof(double**)));
  double *** lambda2= static_cast<double***>(malloc(Np*sizeof(double**)));
  double *** Tau= static_cast<double***>(malloc(Np*sizeof(double**)));
  double *** Taumean= static_cast<double***>(malloc(Np*sizeof(double**)));
  double ** B0= static_cast<double**>(malloc(Np*sizeof(double*)));
  double ** B0mean= static_cast<double**>(malloc(Np*sizeof(double*)));
  double ** AcceptB0= static_cast<double**>(malloc(Np*sizeof(double*)));
  double ** quadForm= static_cast<double**>(malloc(Np*sizeof(double*)));
  double ** loggauss= static_cast<double**>(malloc(Np*sizeof(double*)));
  double ** s2= static_cast<double**>(malloc(Np*sizeof(double*)));
  double ** s2Mean= static_cast<double**>(malloc(Np*sizeof(double*)));
  
  int mj=0;
  for (m=0;m<Np;m++){
    
    Gamvs[m]= static_cast<double*>(malloc(P[m]*sizeof(double)));
    Gam[m]=bmatrix(0,P[m]-1,0,r-1);
    Gammean[m]=dmatrix(0,P[m]-1,0,r-1);
    rhoest[m]= static_cast<bool*>(malloc(r*sizeof(bool)));
    rhomean[m]= static_cast<double*>(malloc(r*sizeof(double)));
    R[m]=bmatrix(0,r-1,0,K[m]-1);
    AcceptR[m]=dmatrix(0,r-1,0,K[m]-1);
    Bmean[m]=dmatrix(0,r-1,0,K[m]-1);
    Rmean[m]=dmatrix(0,r-1,0,K[m]-1);
    B[m]=dmatrix(0,r-1,0,K[m]-1);
    lambda2[m]=dmatrix(0,r-1,0,P[m]-1);
    B0[m]= static_cast<double*>(malloc(r*sizeof(double)));
    B0mean[m]= static_cast<double*>(malloc(r*sizeof(double)));
    AcceptB0[m]= static_cast<double*>(malloc(r*sizeof(double)));
    Tau[m]=dmatrix(0,r-1,0,P[m]-1);
    Taumean[m]=dmatrix(0,r-1,0,P[m]-1);
    
    for (l=0;l<r;l++){
      rhomean[m][l]=0; B0mean[m][l]=0; B0[m][l]=0.1; AcceptB0[m][l]=0;
      for (j=0;j<P[m];j++) Gammean[m][j][l]=0;
      
      double  uni=gsl_ran_flat (rr, 0, 1);
      qv[m][l]=probvarsel;
      wg[m][l]=0.5; //Prior prob group selection
      if (IndVar[m]==1) //Response
        qv[m][l]=0.5;
      else if (IndVar[m]==2) qv[m][l]=1;
      
      if (uni<0.5) {
        rhoest[m][l]=1;
        for (j=0;j<P[m];j++){
          if (gsl_ran_flat (rr, 0, 1)<qv[m][l]) {
            Gam[m][j][l]=1;A[m][l][j]=0.01;
          } else {
            Gam[m][j][l]=0;A[m][l][j]=0;
          }
        }
      } else {
        rhoest[m][l]=0;
        for (j=0;j<P[m];j++) {
          Gam[m][j][l]=0;A[m][l][j]=0;
        }
      }
      
      if (IndVar[m]==1){ //Response
        Gam[m][0][l]=rhoest[m][l];
      } else if (IndVar[m]==2){ //Clinical factors
        rhoest[m][l]=1;
        for (j=0;j<P[m];j++) {
          Gam[m][j][l]=1; A[m][l][j]=0.01;
        }
      }
      
      for (k=0;k<K[m];k++){
        Rmean[m][l][k]=0;Bmean[m][l][k]=0;AcceptR[m][l][k]=0;
        R[m][l][k]=0;B[m][l][k]=0;

        if ((Method==1)&& (IndVar[m]==0)){
          if (rhoest[m][l]==1){
            double  ui=gsl_ran_flat (rr, 0, 1);
            if (ui<wg[m][l]) {
              R[m][l][k]=1;B[m][l][k]=0.1;
            }
          }
        }
      }
      
      for (j=0;j<P[m];j++){
        Tau[m][l][j]=1; Taumean[m][l][j]=0;
        lambda2[m][l][j]=1;
        if (IndVar[m]==2){
          Tau[m][l][j]=100;
        }
      }
    }
    
    s2[m]= static_cast<double*>(malloc(P[m]*sizeof(double)));
    s2Mean[m]= static_cast<double*>(malloc(P[m]*sizeof(double))); 
    quadForm[m]= static_cast<double*>(malloc(P[m]*sizeof(double)));
    loggauss[m]= static_cast<double*>(malloc(P[m]*sizeof(double)));
    
    for (j=0;j<P[m];j++) {
      loggauss[m][j]=0;
      s2[m][j]=0.1; s2Mean[m][j]=0;
      quadForm[m][j]=Gamvs[m][j]=0; 
      mj+=1;
    }
  }
  
  // Dim of a model
  int dim=0;
  for (m=0;m<Np;m++){
    for (j=0;j<P[m];j++){
      dim+=P[m]*r;
    }
  }
  dim+=(Np-1)*r;

  int t;
  int N=nbrsample+burninsample;
  double intercp;
  bool ** rhomodel= static_cast<bool**>(malloc(nbrsample*sizeof(bool*)));
  for (t=0;t<nbrsample;t++){
    rhomodel[t]= static_cast<bool*>(malloc(dim*sizeof(bool)));
  }
  
  for (t=0;t<N;t++){
    for (m=0;m<Np;m++){
      if (IndVar[m]==1){
        
        // Sample beta
        mat WtW = W.t() * W;
        mat V_beta = inv(WtW / sigma2_0 + diagmat(1.0 / beta_prior_var));
        vec m_beta = V_beta * (W.t() * (y - ksi(study_site) - theta(family)) / sigma2_0 + beta_prior_mean / beta_prior_var);
        beta = mvnrnd(m_beta, V_beta);
        
        // Sample ksi
        for (int s = 0; s < S_count; ++s) {
          uvec idx = find(study_site == s);
          vec y_s = y(idx);
          mat W_s = W.rows(idx);
          vec v_s = theta(family(idx));
          double V_ksi = 1.0 / (y_s.n_elem / sigma2_0 + 1.0 / sigma_ksi2);
          double m_ksi = V_ksi * sum(y_s - W_s * beta - v_s) / sigma2_0;
          ksi[s] = rnorm_cpp(1, m_ksi, sqrt(V_ksi))[0];
        }
        
        // Sample theta
        for (int f = 0; f < F_count; ++f) {
          uvec idx = find(family == f);
          vec y_f = y(idx);
          mat W_f = W.rows(idx);
          vec ksi_f = ksi(study_site(idx));
          int site = study_site(idx[0]);
          double V_theta = 1.0 / (y_f.n_elem / sigma2_0 + 1.0 / sigma_theta2[site]);
          double m_theta = V_theta * sum(y_f - W_f * beta - ksi_f) / sigma2_0;
          theta[f] = rnorm_cpp(1, m_theta, sqrt(V_theta))[0];
        }
        
        // Sample sigma_ksi2
        double a_sigma_ksi = sigma_ksi_prior_a + S_count / 2.0;
        double b_sigma_ksi = sigma_ksi_prior_b + sum(square(ksi)) / 2.0;
        sigma_ksi2 = rinvgamma_cpp(a_sigma_ksi, b_sigma_ksi);
        
        // Sample sigma_theta2 for each study site
        for (int s = 0; s < S_count; ++s) {
          uvec families_in_site = find(study_site == s);
          families_in_site = unique(family(families_in_site));
          if (!families_in_site.is_empty()) {
            for (size_t idx = 0; idx < families_in_site.n_elem; ++idx) {
              if (families_in_site[idx] < 0 || families_in_site[idx] >= F_count) {
                Rcpp::Rcerr << "Invalid family index: " << families_in_site[idx] << std::endl;
                Rcpp::stop("Index out of bounds");
              }
            }
            double a_sigma_theta = sigma_theta_prior_a + families_in_site.n_elem / 2.0;
            double b_sigma_theta = sigma_theta_prior_b + sum(square(theta(families_in_site))) / 2.0;
            sigma_theta2[s] = rinvgamma_cpp(a_sigma_theta, b_sigma_theta);
          }
        }
        
        // Sample sigma2_0
        double a_sigma = sigma_prior_a + n / 2.0;
        double b_sigma = sigma_prior_b + sum(square(y - W * beta - ksi(study_site) - theta(family))) / 2.0;
        sigma2_0 = rinvgamma_cpp(a_sigma, b_sigma);
        
        // Store samples
        beta_samples.row(t) = beta.t();
        ksi_samples.row(t) = ksi.t();
        theta_samples.row(t) = theta.t();
        sigma_ksi2_samples[t] = sigma_ksi2;
        sigma_theta2_samples.row(t) = sigma_theta2.t();
        sigma2_samples[t] = sigma2_0;
        
        // Store values in Thierry's objects
        intercp = beta[0];
        s2[0] = &sigma2_0;
        
        // Residualize for downstream sampling
        y_tilde = y - W * beta - ksi(study_site) - theta(family);
        for (i=0;i<n;i++){
          for (j=0;j<P[m];j++){
            X1[m][i][j] = y_tilde(i,j);
          }
        }

      }
        
      double sumrho=0;
      
      for (l=0;l<r;l++){
        double sumeta=0;
        sumrho+=rhoest[m][l];
        for (j=0;j<P[m];j++){
          sumeta+=Gam[m][j][l];
        }
        if (IndVar[m]==2) qv[m][l]=1;
      }
      
      if (IndVar[m]==1){
        q[m]=gsl_ran_beta (rr, al0+sumrho, bl0+r-sumrho);
        for (l=0;l<r;l++){
          qv[m][l]=q[m];
        }
      } else if (IndVar[m]==0) {
        q[m]=gsl_ran_beta (rr, al+sumrho, bl+r-sumrho);
      }
      
      for (j=0;j<P[m];j++){
        logGausQuadForm(j,r, n,P[m], Tau[m], U,X1[m], s2[m][j],&quadForm[m][j],Gam[m][j],&loggauss[m][j]);
        if (t==0){
          loggauss[m][j]=-DBL_MAX;
        }
      }
      
      if (IndVar[m]!=2){
        SamplerhoGamma(rr,r, n,IndVar[m],P[m],rhoest[m],Tau[m], U,X1[m],qv[m],q[m],s2[m],quadForm[m],Gam[m],loggauss[m]);
      }
      
      if (IndVar[m]!=1){
        sigma2(rr,  P[m],n,a0, b0,quadForm[m],s2[m]);
      }

      LoadAOther(rr,r, n,P[m],rhoest[m],Tau[m],A[m],U,X1[m],s2[m],Gam[m]);
      
      for (l=0;l<r;l++){
        if (IndVar[m]!=2){
          TauLambda(l,rr,K[m],P[m],A[m][l],B[m][l],B0[m][l],Tau[m][l],Path[m],alpha,lambda2[m][l],s2[m],Gam[m]);
        }
        if ((Method==1)&&(IndVar[m]==0)){
          sumrho=0;
          for (k=0;k<K[m];k++) sumrho+=R[m][l][k];
          wg[m][l]=gsl_ran_beta(rr,alg+rhoest[m][l]*sumrho,blg+rhoest[m][l]*(K[m]-sumrho));
          GroupEffect(l,rr,rhoest[m][l],K[m],P[m], R[m][l],A[m][l],B[m][l],B0[m][l],Tau[m][l],Path[m],alphab, betab,wg[m][l],alpha, lambda2[m][l],s2[m],AcceptR[m][l],Gam[m]);
        }
        EffectZero(l,rr,K[m],P[m],rhoest[m][l], R[m][l],B[m][l],&B0[m][l],Path[m],alphab0, betab0,alpha, lambda2[m][l],&AcceptB0[m][l],Gam[m]);
      }
      
    }

    SampleUU(rr,r,n,Np,P,A,U,X1,s2);
    
    if (t>=burninsample){
      InterceptMean+=intercp/nbrsample;
      int rm=0;
      for (m=0;m<Np;m++){
        for (j=0;j<P[m];j++){
          for (l=0;l<r;l++){
            rhomodel[t-burninsample][rm]=Gam[m][j][l];
            rm++;
          }
        }
      }
      
      for (m=0;m<Np;m++){
        if (IndVar[m]!=1){
          for (l=0;l<r;l++){
            rhomodel[t-burninsample][rm]=rhoest[m][l];
            rm++;
          }
        }
      }
      
      for (l=0;l<r;l++){
        for (i=0;i<n;i++){
          meanU[i][l]+=U[i][l]/nbrsample;
        }
      }
      
      for (m=0;m<Np;m++){
        for (j=0;j<P[m];j++){
          s2Mean[m][j]+=s2[m][j]/nbrsample;
          int xx=1;
          for (l=0;l<r;l++){
            Gammean[m][j][l]+=(double) Gam[m][j][l]/nbrsample;
            Taumean[m][l][j]+=Tau[m][l][j]/nbrsample;
            xx*=1-Gam[m][j][l];
          }
          
          if (xx==0) Gamvs[m][j]+=1.0/nbrsample;
        }
        
        for (l=0;l<r;l++){
          B0mean[m][l]+=B0[m][l]/nbrsample; 
          rhomean[m][l]+=(double) rhoest[m][l]/nbrsample;
          for (k=0;k<K[m];k++){
            Rmean[m][l][k]+=(double) R[m][l][k]/nbrsample;
            Bmean[m][l][k]+=B[m][l][k]/nbrsample;
          }
        }
      }

      if (t%((N+5)/5)==1){

        printf("\n");
  
        printf("The number of iterations is  %d\n",t);
      }
    }
  }
  printf("\n");

  int sumP=0;int sumK=0;
  for (m=0;m<Np;m++){
    for (l=0;l<r;l++){
      CompoSelMean[m * r + l]=rhomean[m][l];
    }
    
    for (j=0;j<P[m];j++){
      for (l=0;l<r;l++){
        VarSelMean[sumP*r+j*r+l]=Gammean[m][j][l];
      }
    }
    
    for (l=0;l<r;l++){
      for (k=0;k<K[m];k++){
        GrpSelMean[sumK*r+l*K[m]+k]=Rmean[m][l][k];
        GrpEffectMean[sumK*r+l*K[m]+k]=Bmean[m][l][k];
      }
    }
    
    for (j=0;j<P[m];j++){
      VarSelMeanGlobal[sumP+j]=Gamvs[m][j];
    }
    
    for (l=0;l<r;l++){
      IntGrpMean[m*r+l]=B0mean[m][l];
    }
    sumP+=P[m]; sumK+=K[m];
  }
  
  /* Loading estimate for prediction using multiple models*/
  int u=0;
  for (i=0;i<n;i++){
    for (l=0;l<r;l++){
      EstU[u]=meanU[i][l];
      u+=1;
    }
  }
  
  int countmodel=0;
  int *modelidx= static_cast<int*>(malloc(nbrsample*sizeof(double)));
  bool **UniqModel=  UniqueModel(nbrsample, dim, rhomodel,modelidx,&countmodel);
  free_bmatrix(rhomodel,0,nbrsample-1,0,dim-1);
  double * logpo= static_cast<double*>(malloc(countmodel*sizeof(double)));
  for (t=0;t<countmodel;t++){
    int rm=0;
    int sumg=0;
    for (m=0;m<Np;m++){
      if (IndVar[m]==1){
        for (i=0;i<n;i++){
          for (j=0;j<P[m];j++){
            X1[m][i][j]=X[m][i][j]-InterceptMean;
          }
        }
      }
      
      for (j=0;j<P[m];j++){
        for (l=0;l<r;l++){
          Gam[m][j][l]=UniqModel[t][rm];rm++;
          sumg+= Gam[m][j][l];
        }
      }
    }
    for (m=0;m<Np;m++){
      if (IndVar[m]!=1){
        for (l=0;l<r;l++){
          rhoest[m][l]=UniqModel[t][rm];
          rm++;
        }
      } else {
        for (l=0;l<r;l++){
          rhoest[m][l]=Gam[m][0][l];
        }
      } 
    }
    
    logPostGam(&logpo[t],IndVar,Np,r, n,P,Tau, meanU,X1,s2Mean,rhoest,Gam,qv,q);
  }
  int * highmodelidx= static_cast<int*>(malloc(countmodel*sizeof(int)));
  sort(countmodel,logpo,highmodelidx);
  double maxlogpost=logpo[0];

  int nbrmax=std::min(maxmodel,countmodel);
  nbrmodel1 = nbrmax;
  for (l=0;l<nbrmax;l++){
    logpo[l]=exp(logpo[l]-maxlogpost);
  }
  
  double sumpost=nbrmax*mean(nbrmax,logpo);
  for (l=0;l<nbrmax;l++){
    logpo[l]=logpo[l]/sumpost;
    postgam[l]=logpo[l];
  }
  int ll=0;
  
  for (t=0;t<nbrmax;t++){
    int rm=0;
    for (m=0;m<Np;m++){
      for (j=0;j<P[m];j++){
        for (l=0;l<r;l++){
          Gam[m][j][l]=UniqModel[highmodelidx[t]][rm];
          rm++;
        }
      }
    }
    
    for (m=0;m<Np;m++){
      if (IndVar[m]!=1){
        for (l=0;l<r;l++){
          rhoest[m][l]=UniqModel[highmodelidx[t]][rm];
          rm++;
        }
      } else {
        for (l=0;l<r;l++) {
          rhoest[m][l]=Gam[m][0][l];
        }
      }
    }

    for (m=0;m<Np;m++){
      EstimateLoad(rr,r, n,P[m],rhoest[m],Taumean[m],A[m],meanU,X1[m],s2Mean[m],Gam[m]);
      for (l=0;l<r;l++){
        for (j=0;j<P[m];j++) {
          EstLoadMod[ll]=A[m][l][j];
          ll++;
        }
      }
    }
  }
  free(modelidx); free(logpo); free(highmodelidx);
  free_bmatrix(UniqModel,0,countmodel-1,0,dim-1);
  
  
  double thres=0.5;
  ll=0;int ls=0;
  for (m=0;m<Np;m++){
    for (l=0;l<r;l++){
      if (rhomean[m][l]>=thres) rhoest[m][l]=1; else rhoest[m][l]=0; 
      for (j=0;j<P[m];j++){
        if (rhoest[m][l]==1) {
          if (Gammean[m][j][l]>=thres) Gam[m][j][l]=1; else Gam[m][j][l]=0;
        } else {
          Gam[m][j][l]=0;
        }
      }
    }
    EstimateLoad(rr,r, n,P[m],rhoest[m],Taumean[m],A[m],meanU,X[m],s2Mean[m],Gam[m]);
    
    for (l=0;l<r;l++){
      for (j=0;j<P[m];j++){
        EstLoad[ll]=A[m][l][j];
        ll+=1;
      }
    }
  }
  
  free_dmatrix(rhomean,0,Np-1,0,r-1);
  
  for (m=0;m<Np;m++){
    for (j=0;j<P[m];j++){
      EstSig2[ls]=s2Mean[m][j];
      ls+=1;
    }
  }
  free_dmatrix(B0,0,Np-1,0,r-1);
  free_dmatrix(B0mean,0,Np-1,0,r-1);
  
  for (m=0;m<Np;m++){
    free_bmatrix(Gam[m],0,P[m]-1,0,r-1);free_dmatrix(Gammean[m],0,P[m]-1,0,r-1);
    free_bmatrix(R[m],0,r-1,0,K[m]-1);free_bmatrix(Path[m],0,P[m]-1,0,K[m]-1);
    free_dmatrix(Rmean[m],0,r-1,0,K[m]-1);free_dmatrix(Bmean[m],0,r-1,0,K[m]-1);
    free_dmatrix(B[m],0,r-1,0,K[m]-1);free_dmatrix(lambda2[m],0,r-1,0,P[m]-1);
    free_dmatrix(AcceptR[m],0,r-1,0,K[m]-1);free_dmatrix(Tau[m],0,r-1,0,P[m]-1);
    free_dmatrix(Taumean[m],0,r-1,0,P[m]-1);
    free_dmatrix(A[m],0,r-1,0,P[m]-1);free_dmatrix(X[m],0,n-1,0,P[m]-1);
    free_dmatrix(X1[m],0,n-1,0,P[m]-1);
    free(s2[m]);free(s2Mean[m]);
    free(rhoest[m]);free(Gamvs[m]);
    free(quadForm[m]);free(loggauss[m]);free(AcceptB0[m]);
    free(qv[m]);free(wg[m]);
  }
  
  free(Gam);free(Gamvs);
  free(Tau);free(Taumean);
  free(A);free(X);free(X1);free(Gammean);free(R);
  free_dmatrix(U,0,n-1,0,r-1);
  free_dmatrix(meanU,0,n-1,0,r-1);
  free(s2);free(s2Mean);
  free(rhoest);free(q);free(qv);free(wg);
  free(quadForm);free(loggauss);free(B);free(lambda2);
  gsl_rng_free (rr);free(AcceptB0);
  free(Path);free(AcceptR);free(Rmean);free(Bmean);
  
  t1 = clock() - t1;
  double  time_taken = ((double)t1)/CLOCKS_PER_SEC; // in seconds
  printf("\n\nTime taken in seconds is %f\n",time_taken);
  printf("\nTime taken in minutes is %f\n",time_taken/60);
  
  return Rcpp::List::create(
    Rcpp::Named("VarSelMeanGlobal") = VarSelMeanGlobal,
    Rcpp::Named("VarSelMean") = VarSelMean,
    Rcpp::Named("GrpSelMean") = GrpSelMean,
    Rcpp::Named("GrpEffectMean") = GrpEffectMean,
    Rcpp::Named("EstLoad") = EstLoad,
    Rcpp::Named("EstU") = EstU,
    Rcpp::Named("CompoSelMean") = CompoSelMean,
    Rcpp::Named("InterceptMean") = InterceptMean,
    Rcpp::Named("IntGrpMean") = IntGrpMean,
    Rcpp::Named("EstSig2") = EstSig2,
    Rcpp::Named("EstLoadMod") = EstLoadMod,
    Rcpp::Named("nbrmodel1") = nbrmodel1,
    Rcpp::Named("beta_samples") = beta_samples,
    Rcpp::Named("ksi_samples") = ksi_samples,
    Rcpp::Named("theta_samples") = theta_samples,
    Rcpp::Named("sigma_ksi2_samples") = sigma_ksi2_samples,
    Rcpp::Named("sigma_theta2_samples") = sigma_theta2_samples,
    Rcpp::Named("sigma2_samples") = sigma2_samples
  );
}

/*** R
# Importing libraries
require(MASS)
require(reshape2) # For melt

# Define R wrapper
cat("Started defining BIP")
BIP <- function(dataList=dataList,IndicVar=IndicVar, groupList=NULL,Method=Method,nbrcomp=4, sample=5000, burnin=1000,nbrmaxmodels=50,
                priorcompselv=c(1,1),priorcompselo=c(1,1),priorb0=c(2,2),priorb=c(1,1),priorgrpsel=c(1,1),probvarsel=0.05,
                # Gibbs Sampler Nested Arguments start here
                y, W, study_site, family, n_iter, priors) {
  
  if (sample<=burnin){
    stop("Argument burnin must be smaller than sample: the number of MCMC iterations.")
  }
  if (sample<=20){
    stop("Please specify a larger number of MCMC iterations")
  }
  
  if (is.null(nbrcomp)){
    D=which(IndicVar==0)
    mysvd=lapply(D, function(i)  svd(dataList[[i]]))
    mysumsvd=lapply(1:length(D), function(i) cumsum(mysvd[[i]]$d)/max(cumsum(mysvd[[i]]$d))*100)
    KMax=max(unlist(lapply(1:length(D), function(i) min(which(mysumsvd[[i]]>=80, arr.ind = TRUE))))) #chooses maximum from D cumulative proportions
    nbrcomp=min(KMax+1,10)
  }
  
  
  #Method="GroupInfo"
  if (Method=="BIP"){
    meth=0
  } else if (Method=="BIPnet") {
    meth=1
  } else {
    stop("You must provide either BIP or BIPnet")
  }
  
  Np=length(dataList)
  P=NULL
  n=nrow(dataList[[1]])
  P=NULL
  MeanData=list()
  SD=list()
  for (i in 1:Np){
    dataList[[i]]=as.matrix(dataList[[i]])
    P[i]=ncol(dataList[[i]])
    if ((IndicVar[i]==0)||(IndicVar[i]==2)){
      MeanData[[i]]=apply(dataList[[i]],2,mean)
      SD[[i]]=apply(dataList[[i]],2,sd)
      dataList[[i]]=scale(as.matrix(dataList[[i]]),T,T)
    }
    dataList[[i]]=t(dataList[[i]])
  }
  datasets=unlist(dataList)
  
  if (is.null(groupList) && (meth==1)){
    stop("You must provide a list of grouping information")
  } else if (is.null(groupList)){
    for (ll in 1:length(IndicVar)){
      groupList[[ll]]=matrix(1,P[ll],1)
    }
  }
  
  
  ll=1
  K=NULL
  for (i in 1:Np){
    if (IndicVar[i]!=0){
      K[i]=1
    } else {
      K[i]=ncol(groupList[[ll]])
      ll=ll+1
    }
  }
  groupList=lapply(groupList,t)
  Paths=unlist(groupList)
  
  result <- mainfunction(
    Method=as.integer(meth),n=as.integer(n),P=as.integer(P),r=as.integer(nbrcomp),Np=as.integer(Np), 
    datasets=as.double(datasets),IndVar=as.integer(IndicVar), K=as.integer(K), Paths=as.integer(Paths),maxmodel=as.integer(nbrmaxmodels),                nbrsample=as.integer(sample), burninsample=as.integer(burnin),
    CompoSelMean=as.double(rep(0,Np*nbrcomp)),VarSelMean=as.double(rep(0,nbrcomp*sum(P))),
    VarSelMeanGlobal=as.double(rep(0,sum(P))),GrpSelMean=as.double(rep(0,nbrcomp*sum(K))),
    GrpEffectMean=as.double(rep(0,nbrcomp*sum(K))),IntGrpMean=as.double(rep(0,nbrcomp*Np)),
    EstU=as.double(rep(0,n*nbrcomp)),EstSig2=as.double(rep(0,sum(P))),InterceptMean=as.double(rep(0,1)),
    EstLoadMod=as.double(rep(0,nbrmaxmodels*nbrcomp*sum(P))),EstLoad=as.double(rep(0,nbrcomp*sum(P))),
    nbrmodel=as.integer(0),postgam=rep(0,nbrmaxmodels),priorcompsel=priorcompselv,
    priorcompselo=priorcompselo,priorb0=priorb0,priorb=as.double(priorb),priorgrpsel=priorgrpsel,probvarsel=as.double(probvarsel),
    # Gibbs Sampler Nested Arguments start here
    y=y, W=W, study_site=study_site, family=family, n_iter=n_iter, beta_prior_mean=priors$beta_prior_mean, beta_prior_var=priors$beta_prior_var, 
    sigma_ksi_prior_a=priors$sigma_ksi_prior_a, sigma_ksi_prior_b=priors$sigma_ksi_prior_b, 
    sigma_theta_prior_a=priors$sigma_theta_prior_a, sigma_theta_prior_b=priors$sigma_theta_prior_b, 
    sigma_prior_a=priors$sigma_prior_a, sigma_prior_b=priors$sigma_prior_b
  )
  
  reseffect=result$EstLoadMod
  nbrmodel=result$nbrmodel1

  EstLoadModel <- lapply(1:nbrmodel, function(x) {
    lapply(1:Np, function(y) {
      NULL
    })
  })
  for (mo in 1:nbrmodel){
    for (m in 1:Np){
      x=sum(P[1:(m-1)])
      y=sum(P[1:m])
      if (m==1) {x=0}
      init=1+x*nbrcomp+(mo-1)*nbrcomp*sum(P)
      final=y*nbrcomp+(mo-1)*nbrcomp*sum(P)
      EstLoadModel[[mo]][[m]]=matrix(reseffect[init:final],nbrcomp,byrow=T);
    }
  }

  resvarsel1=result$VarSelMeanGlobal
  resvarsel2=result$VarSelMean
  resvarsel3=result$GrpSelMean
  resvarsel4=result$GrpEffectMean
  resvarsel5=result$EstLoad
  EstimateU=matrix(result$EstU,n,byrow=T)
  CompoSelMean=matrix(result$CompoSelMean,Np,byrow=T)
  IntGrpMean=matrix(result$IntGrpMean,Np,byrow=T)
  VarSelMeanGlobal=list()
  VarSelMean=list()
  GrpSelMean=list()
  GrpEffectMean=list()
  EstLoad=list()
  EstSig2=list()
  m1=m2=m3=1
  for (m in 1:Np){
    VarSelMeanGlobal[[m]]=resvarsel1[m1:(m1-1+P[m])]
    VarSelMean[[m]]=matrix(resvarsel2[m2:(m2-1+P[m]*nbrcomp)],P[m],byrow=T)
    GrpSelMean[[m]]=matrix(resvarsel3[m3:(m3-1+K[m]*nbrcomp)],nbrcomp,byrow=T)
    GrpEffectMean[[m]]=matrix(resvarsel4[m3:(m3-1+K[m]*nbrcomp)],nbrcomp,byrow=T)
    EstLoad[[m]]=matrix(resvarsel5[m2:(m2-1+P[m]*nbrcomp)],nbrcomp,byrow=T)
    EstSig2[[m]]=result$EstSig2[m1:(m1-1+P[m])]
    m1=m1+P[m]
    m2=m2+P[m]*nbrcomp
    m3=m3+K[m]*nbrcomp

  }

  return (list(EstU=EstimateU,VarSelMean=VarSelMean,VarSelMeanGlobal=VarSelMeanGlobal,CompoSelMean=CompoSelMean,
               GrpSelMean=GrpSelMean,GrpEffectMean=GrpEffectMean,IntGrpMean=IntGrpMean,EstLoad=EstLoad,
               EstLoadModel=EstLoadModel,nbrmodel=result$nbrmodel1,EstSig2=EstSig2,EstIntcp=result$InterceptMean,
               PostGam=result$postgam,IndicVar=IndicVar,nbrcomp=nbrcomp,MeanData=MeanData,SDData=SD,
               beta_samples = result$beta_samples,
               ksi_samples = result$ksi_samples,
               theta_samples = result$theta_samples,
               sigma_ksi2_samples = result$sigma_ksi2_samples,
               sigma_theta2_samples = result$sigma_theta2_samples,
              sigma2_samples = result$sigma2_samples))
}

cat("Started simulating data")
# Simulate data & estimate associated parameters
source("simulate_random_intercept_data.R")
# nu2_site_truth <- 3
# nu2_family_truth <- 2
# simulation_results <- simulate_re_data(n_sites=30, nu2=nu2_truth, seed=2)
simulation_results <- simulate_re_data_nested(n_obs = 200, 
                                              n_sites = 10, n_families_per_site = 10,
                                              n_individs_per_family = 10,
                                              nu2_site = nu2_site_truth, 
                                              nu2_family = nu2_family_truth)
dataList <- list(simulation_results$X[[1]],
                  simulation_results$X[[2]],
                  simulation_results$Y)
n_sample <- 5000
n_burnin <- 1000

# Priors
priors <- list(beta_prior_mean = c(0, 0),
               beta_prior_var = c(100, 100),
               sigma_ksi_prior_a = 2,
               sigma_ksi_prior_b = 1,
               sigma_theta_prior_a = 2,
               sigma_theta_prior_b = 1,
               sigma_prior_a = 2,
               sigma_prior_b = 1)

cat("Calling BIP")

BA <- BIP(dataList=dataList, IndicVar=c(0,0,1), Method="BIP", nbrcomp=4, sample=n_sample, burnin=n_burnin, 
          y=simulation_results$Y, W=simulation_results$W, study_site=simulation_results$study_site, family=simulation_results$family, 
          n_iter = n_sample+n_burnin, priors)


*/

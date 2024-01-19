#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int logmultigaussianT(const gsl_vector * x, const gsl_vector * y,
                      const gsl_matrix * L,
                      double * result,double *quadF,
                      gsl_vector * work){
  const size_t M = L->size1;
  const size_t N = x->size;
  size_t i;
  double quadForm=0;        /* (x)' Sigma^{-1} (x) */
  double logSqrtDetSigma=0; /* log [ sqrt(|Sigma|) ] */
  
  /* compute: work = x - mu*/ 
  for (i = 0; i < M; ++i)
  {
    double xi = gsl_vector_get(y, i);
    //double mui = gsl_vector_get(mu, i);
    gsl_vector_set(work, i, xi);
  }
  
  
  /* compute: work = L^{-1} * (x - mu) */
  gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, L, work);
  
  /* compute: quadForm = (x - mu)' Sigma^{-1} (x - mu) */
  gsl_blas_ddot(work, work, &quadForm);
  double x2=0;
  gsl_blas_ddot(x, x, &x2);
  quadForm=x2-quadForm;
  /* compute: log [ sqrt(|Sigma|) ] = sum_i log L_{ii} */
  for (i = 0; i < M; ++i)
  {
    double Lii = gsl_matrix_get(L, i, i);
    logSqrtDetSigma += log(Lii);
  }
  
  *quadF=quadForm;
  *result = -0.5*quadForm - logSqrtDetSigma - 0.5*N*log(2.0*M_PI);
  return GSL_SUCCESS;
  
}


// [[Rcpp::export]]
void findc(int n,_Bool R[n],int a,int * IDX, int *nx)
{
  int ii_data[n];
  int idx = 0;
  int  ii = 0;
  _Bool  exitg2 = 0;
  _Bool guard2=0;
  while ((exitg2 == 0) && (ii < n)) {
    guard2 = 0;
    if (R[ii]!= a) {
      idx++;
      ii_data[idx - 1] = ii;
      if (idx >= n) {
        exitg2 = 1;
      } else {
        guard2 = 1;
      }
    } else {
      guard2 = 1;
    }
    
    if (guard2 == 1) {
      ii++;
    }
  }
  
  int loop_ub=idx;
  for (idx = 0; idx < loop_ub; idx++) {
    IDX[idx] = ii_data[idx];
  }
  *nx=loop_ub;
  
}

// [[Rcpp::export]]
void logGausQuadForm(int j,int r, int n,int p,double ** Tau, double ** U,double ** X,double s2,double* quadForm,_Bool * Gam,double * loggauss){
  int i,s,s1;
  int NZ1[r];
  int nx1=0;
  findc(r,Gam,0,NZ1, &nx1);
  
  double result=0;double quadF=0;
  if (nx1>0){
    gsl_vector *work1 =gsl_vector_alloc (nx1);
    double * Sigma1=malloc(nx1*nx1*sizeof(double));
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
  
  /*gsl_ran_multivariate_gaussian_log_pdf (xj, mu, &m.matrix, result, work);*/
  *quadForm=quadF;
  *loggauss=result;
  
}

// [[Rcpp::export]]
void SamplerhoGamma(gsl_rng * rr,int r, int n,int IndVar,int p, _Bool * rho,double ** Tau, double ** U,double ** X, double* q1,double q2,double* s2,double* quadForm,_Bool** Gam,double *loggauss){
  int l,j;
  _Bool *rhonew=malloc(r*sizeof(_Bool));
  _Bool Gamnew[p][r];
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
      //printf("RAT=%.2lf",rat);
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
        rhonew[l]=1;
        double logpostnew=0;
        double logpostold=0;
        double * quadForm1=malloc(p*sizeof(double));
        double * loggausnew1=malloc(p*sizeof(double));
        double quad1=0; double quad2=0;
        double loggaussold=0; double loggaussnew=0;
        for (j=0;j<p;j++){
          double logqj=0;double logmqj=0;
          logpostold+=loggauss[j];
          Gam[j][l]=0;
          logGausQuadForm(j,r, n,p, Tau,  U,X,s2[j],&quad1,Gam[j],&loggaussold);
          Gamnew[j][l]=1;
          logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quad2,Gamnew[j],&loggaussnew);
          double rat=-loggaussold+loggaussnew+log(q1[l])-log(1-q1[l]); // log(P_lj / Q_lj)
          //printf("%lf ",rat);
          double x1=loggaussnew+log(q1[l]); // log G_1
          double x2=loggaussold+log(1-q1[l]); // log G_0
          double maxq=MAX(x1,x2);
          double den=maxq+log(exp(x1-maxq)+exp(x2-maxq)); // log(G_0 + G_1)
          logqj=x1-den; // log P_lj
          logmqj=x2-den; // log Q_lj = 1 - log P_lj
          
          double uni=gsl_ran_flat (rr, 0, 1);
          if ((log(uni/(1-uni))<rat)){
            Gam[j][l]=1;Gamnew[j][l]=1;
            logpostnew+=loggaussnew+log(q1[l]);
            loggausnew1[j]=loggaussnew;
            logq+=logqj; // log proposal difference
            quadForm1[j]=quad2;
          } else {
            Gam[j][l]=0;Gamnew[j][l]=0;
            logq+=logmqj;
            logpostnew+=loggaussold+log(1-q1[l]);
            quadForm1[j]=quad1;
            loggausnew1[j]=loggaussold;
          }
        }
        logpostnew+=log(q2);
        logpostold+=log(1-q2);
        double un=gsl_ran_flat (rr, 0, 1);
        //printf("%lf ",logpostnew-logpostold-logq);
        double rat1=logpostnew-logpostold-logq; // log acceptance ratio
        if (log(un)<rat1){
          rho[l]=1;
          for (j=0;j<p;j++){
            quadForm[j]=quadForm1[j];
            loggauss[j]=loggausnew1[j];
          }
        } else {
          rho[l]=0;
          for (j=0;j<p;j++) Gam[j][l]=Gamnew[j][l]=0;
        }
        rhonew[l]=rho[l];
        free(quadForm1);free(loggausnew1);
      } else {
        rhonew[l]=0;
        double logpostnew=0;
        double logpostold=0;
        double * quadForm1=malloc(p*sizeof(double));
        double * loggausnew1=malloc(p*sizeof(double));
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
            double maxq=MAX(x1,x2);
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
        double * quadForm1=malloc(p*sizeof(double));
        double * loggausnew1=malloc(p*sizeof(double));
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


/*** R

*/

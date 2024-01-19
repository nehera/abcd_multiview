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
//#include "header.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

//#include <utils.h>

//// Utils

double mean(int n,double * x){
  int i;
  double me=0;
  for (i=0;i<n;i++)
    me+=x[i];
  return me/n;
}
_Bool compar(int n,_Bool *u,_Bool *v){
  int i;
  for (i=0;i<n;i++){
    if (u[i]!=v[i])
      return 0;
  }
  return 1;
}

_Bool **  UniqueModel(int nbrsample, int p, _Bool ** rhosample,int * modelidx,int *countd1){
  
  int i;int j=0;
  int countd=0;
  modelidx[countd]=0;
  for (i=0;i<nbrsample;i++){
    for (j=0;j<countd;j++){
      _Bool resultX =compar(p,rhosample[i], rhosample[modelidx[j]]);
      if (resultX==1)
        break;
    }
    if (j==countd){
      modelidx[countd]=i;
      countd++;
    }
  }
  *countd1=countd;
  
  
  _Bool ** UniqModel=malloc(countd*sizeof(_Bool*));
  for (i=0;i<countd;i++){
    UniqModel[i]=malloc(p*sizeof(_Bool));
    for (j=0;j<p;j++){
      UniqModel[i][j]=rhosample[modelidx[i]][j];
    }
  }
  printf("\n\n\nNumber of models is= %d",countd);
  printf("\n");
  return UniqModel;
}


void save1d(char *filename,int n,double *data)
{
  FILE *qfile;
  qfile=fopen(filename,"wb");
  if(qfile==NULL){
    printf("\nUnable to open file.\n");
    exit(1);
  } else {
    int i;
    for (i = 0; i < n; i++) {
      fprintf(qfile, "%f ", data[i]);
    }
    fclose(qfile);
  }
}


void save2d(char *filename,int n,int p,double ** data)
{
  FILE *qfile;
  qfile=fopen(filename,"wb");
  if(qfile==NULL){
    printf("\nUnable to open file.\n");
    exit(1);
  } else {
    int i, j;
    for (i = 0; i < n; i++) {
      for (j = 0; j < p; j++) {
        fprintf(qfile, "%lf ", data[i][j]);
      }
      fprintf(qfile,"\n");
    }
    fclose(qfile);
  }
}

int multivariate_gaussian (const gsl_rng * r,
                           const gsl_vector * mu,
                           const gsl_matrix * L,
                           gsl_vector * result)
{
  /*
   * L     matrix resulting from the Cholesky decomposition of
   *  *      the inverse of   variance-covariance matrix Sigma^-1 = L L^T (dimension d x d)
   */
  const size_t M = L->size1;
  const size_t N = L->size2;
  
  if (M != N)
  {
    GSL_ERROR("requires square matrix", GSL_ENOTSQR);
  }
  else if (mu->size != M)
  {
    GSL_ERROR("incompatible dimension of mean vector with variance-covariance matrix", GSL_EBADLEN);
  }
  else if (result->size != M)
  {
    GSL_ERROR("incompatible dimension of result vector", GSL_EBADLEN);
  }
  else
  {
    size_t i;
    
    for (i = 0; i < M; ++i)
      gsl_vector_set(result, i, gsl_ran_ugaussian(r));
    
    gsl_blas_dtrsv(CblasLower, CblasTrans, CblasNonUnit, L, result);     
    //gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, L, result);
    gsl_vector_add(result, mu);
    
    return GSL_SUCCESS;
  }
}





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




int logmultigaussian(const gsl_vector * x,
                     const gsl_vector * mu, 
                     const gsl_matrix * L,
                     double * result,double *quadF,
                     gsl_vector * work){
  const size_t M = L->size1;
  const size_t N = L->size2;
  if (M != N)
  {
    GSL_ERROR("requires square matrix", GSL_ENOTSQR);
  }
  else if (mu->size != M)
  {
    GSL_ERROR("incompatible dimension of mean vector with variance-covariance matrix", GSL_EBADLEN);
  }
  else if (x->size != M)
  {
    GSL_ERROR("incompatible dimension of quantile vector", GSL_EBADLEN);
  }
  else if (work->size != M)
  {
    GSL_ERROR("incompatible dimension of work vector", GSL_EBADLEN);
  }
  else
  {
    size_t i;
    double quadForm=0;        /* (x - mu)' Sigma^{-1} (x - mu) */
  double logSqrtDetSigma=0; /* log [ sqrt(|Sigma|) ] */
  
  /* compute: work = x - mu */
  for (i = 0; i < M; ++i)
  {
    double xi = gsl_vector_get(x, i);
    double mui = gsl_vector_get(mu, i);
    gsl_vector_set(work, i, xi - mui);
  }
  
  /* compute: work = L^{-1} * (x - mu) */
  gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, L, work);
  
  /* compute: quadForm = (x - mu)' Sigma^{-1} (x - mu) */
  gsl_blas_ddot(work, work, &quadForm);
  
  /* compute: log [ sqrt(|Sigma|) ] = sum_i log L_{ii} */
  for (i = 0; i < M; ++i)
  {
    double Lii = gsl_matrix_get(L, i, i);
    logSqrtDetSigma += log(Lii);
  }
  *quadF=quadForm;
  *result = -0.5*quadForm - logSqrtDetSigma - 0.5*M*log(2.0*M_PI);
  return GSL_SUCCESS;
  }
}

double inverseGaussian(gsl_rng * r, double mu, double lambda) {
  double v=gsl_ran_gaussian (r, 1);  // sample from a normal distribution with a mean of 0 and 1 standard deviation
  double y = v*v;
  double x = mu + (mu*mu*y)/(2*lambda) - (mu/(2*lambda)) * sqrt(4*mu*lambda*y + mu*mu*y*y);
  double test=gsl_ran_flat (r, 0, 1);    // sample from a uniform distribution between 0 and 1
  if (test <= (mu)/(mu + x))
    return x;
  else
    return (mu*mu)/x;
}
void readBoolVector(char *filename, int nRows,  _Bool * data )
{
  
  FILE *fp = fopen (filename, "r");
  if (fp==NULL)
  {
    printf("We can't open the file (%s).\n", filename);
    exit(1);
  }
  else
  { int iR;
    for (iR = 0; iR < nRows; ++iR )
    {
      int x;
      fscanf(fp, "%d" , &x );
      _Bool bb=(x!=0);
      data[iR]=bb;
    }
    
    fclose(fp);
  }
}


void sort(int n,double *x,int *idx)
{
  int i,j;
  double a;
  int id;
  for (i = 0; i < n; i++)
    idx[i]=i;
  for (i = 0; i < n; ++i)
  {
    for (j = i + 1; j < n; ++j)
    {
      if (x[i] <= x[j])
      {
        a =  x[i];
        id=idx[i];
        idx[i]=idx[j];
        x[i] = x[j];
        idx[j]=id;
        x[j] = a;
      }
    }
  }
  
}



double auc(int n, double * esti,_Bool class[n]){
  double fpr[n+2],tpr[n+2];
  double auc1=0;
  int P=0;//P=positive instances
  int i,j;
  double esti1[n];
  for (i=0;i<n;i++){
    esti1[i]=esti[i];
    if (class[i]==1)
      P+=1;
  }
  int idx[n];
  sort(n,esti1,idx);
  fpr[n+1]=1;tpr[n+1]=1;
  fpr[0]=0;tpr[0]=0;
  for (i=n;i>=1;--i){
    double af=0;double at=0;
    for (j=0;j<n;j++){
      if (esti[j]>esti1[i-1]){
        if (class[j]==0){
          af+=1;
        }
        else {
          at+=1;
        } } }
    tpr[i]=at/P;
    fpr[i]=af/(n-P);
    auc1+=(fpr[i+1]-fpr[i])*(tpr[i+1]+tpr[i]);
  }
  auc1+=(fpr[1]-fpr[0])*(tpr[1]+tpr[0]);
  auc1=0.5*(auc1);
  return auc1;
}




void readBoolArray(char *filename, int nRows, int nCols, _Bool ** data )
{
  
  FILE *fp = fopen (filename, "r");
  if (fp==NULL)
  {
    printf("We can't open the file (%s).\n", filename);
    exit(1);
  }
  else
  { int iR,iC;
    for (iR = 0; iR < nRows; ++iR )
    {
      for (iC = 0; iC < nCols; ++iC )
      {
        int x;
        fscanf(fp, "%d" , &x );
        _Bool bb=(x!=0);
        data[iR][iC]=bb;
      }
    }
    
    fclose(fp);
  }
}



void readDoubleArray(char *filename, int nRows, int nCols, double ** data )
{
  
  FILE *fp = fopen (filename, "r");
  if (fp==NULL)
  {
    printf("We can't open the file (%s).\n", filename);
    exit(1);
  }
  else
  { int iR,iC;
    for (iR = 0; iR < nRows; ++iR )
    {
      for (iC = 0; iC < nCols; ++iC )
      {
        fscanf(fp, "%lf" , &data[iR][iC] );
      }
    }
    fclose(fp);
  }
}
double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  double **m;
  
  m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++)
  {
    m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
    if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}

void nrerror(char error_text[])
{
  printf("Utils run-time error...\n");
  printf("%s\n",error_text);
  printf("...now exiting to system...\n");
  exit(1);
}
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
  int i;
  
  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  
  free((char*) (m+nrl));
}
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

_Bool **bmatrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  _Bool **m;
  
  m=(_Bool **) malloc( (nrh-nrl+1)*sizeof(_Bool*));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++)
  {
    m[i]=(_Bool *) malloc((nch-ncl+1)*sizeof(_Bool));
    if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}
void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch)
{
  int i;
  
  for(i=nrh;i>=nrl;i--) free((_Bool*) (m[i]+ncl));
  
  free((_Bool*) (m+nrl));
}

//// My Functions

void SampleIntercept(gsl_rng * rr,int n, int r, double * intercept, double* sigma2, double sigma20,double ** U, double ** A, double **y){
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

void logPostGam(double *logpo,int * IndVar,int Np,int r, int n,int* P,double *** Tau, double ** U,double *** X,double **s2,_Bool ** rho,_Bool *** Gam,double** qv,double* q){
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
    //if (m<Np-1){
    if (IndVar[m]==0){
      for (l=0;l<r;l++){
        logpost+=rho[m][l]*log(q[m])+(1-rho[m][l])*log(1-q[m]);
      }
    }
  }
  *logpo=logpost;
}

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

void sampleGam(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** X,double *s2,double* quadForm,_Bool ** Gam,double * loggauss,double q){
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
      //for (l=0;l<r;l++) Gam[j][l]=rho[l]*Gam[j][l];
      _Bool Gamnew2[r];
      for (l=0;l<r;l++) Gamnew2[l]=0;
      double logpostold=0;double logpostnew=0;
      _Bool Gamold1[nx1];
      for (s=0;s<nx1;s++){
        Gamold1[s]=Gam[j][NZ1[s]];
      }
      _Bool Gamnew1[nx1];
      proposal(nx1,Gamold1,Gamnew1,phi, rr);
      for (s=0;s<nx1;s++){
        Gamnew2[NZ1[s]]=Gamnew1[s];
        logpostold+=Gamold1[s]*log(q)+(1-Gamold1[s])*log(1-q);
        logpostnew+=Gamnew1[s]*log(q)+(1-Gamnew1[s])*log(1-q);
      }
      
      /*Logpost new*/
      double loggaussnew=0;double quadForm1=0;
      logGausQuadForm(j,r, n,p,Tau,  U,X,s2[j],&quadForm1,Gamnew2,&loggaussnew);
      //void logGausQuadForm(int j,int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** X,double s2,double* quadForm,_Bool * Gam,double * loggauss)
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

double  logpost(int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** X, double q,double* s2,double* quadForm,_Bool ** Gam,double * loggauss){
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


void rho1(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau, double ** U,double ** X, double q,double* s2,double* quadForm,_Bool** Gam,double *loggauss){
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
  double * quadForm1=malloc(p*sizeof(double));
  double * loggauss1=malloc(p*sizeof(double));
  /*double logpostold=logpost(r, n,p, rho,Tau,U,X, q,s2,quadForm1,Gam);*/
  double phi=0.5;
  _Bool rhonew[r];
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
        // Propose to turn the component on
        rhonew[l]=1;
        double logpostnew=0;
        double logpostold=0;
        double * quadForm1=malloc(p*sizeof(double));
        double * loggausnew1=malloc(p*sizeof(double));
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
          double rat=-loggaussold+loggaussnew+log(q1[l])-log(1-q1[l]); // log(P_lj / Q_lj) i.e. log(G_1)-log(G_0)
          //printf("%lf ",rat);
          double x1=loggaussnew+log(q1[l]); // log G_1 // calculated redundant to line 741
          double x2=loggaussold+log(1-q1[l]); // log G_0 // calculated redundant to line 741
          double maxq=MAX(x1,x2);
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
            logpostnew+=loggaussold+log(1-q1[l]); // loggauss doesn't change from the last iteration when we are proposing the component be on and are keeping the feature off. 
            quadForm1[j]=quad1;
            loggausnew1[j]=loggaussold;
          }
        }
        // Add the log probability that the component is on
        logpostnew+=log(q2);
        logpostold+=log(1-q2);
        double un=gsl_ran_flat (rr, 0, 1);
        //printf("%lf ",logpostnew-logpostold-logq);
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


void sigma2(gsl_rng * rr, int p,int n,double a0, double b0, double* quadForm,double * s2){
  int j;
  for (j=0;j<p;j++){
    double inv=1/(b0+0.5*s2[j]*quadForm[j]);
    s2[j]=1/gsl_ran_gamma (rr, a0+n/2.0, inv);
  }
}

void EffectZero(int l,gsl_rng * rr,int K,int p,_Bool rho, _Bool * R,double* B,double *B0,_Bool ** Path,double alphab0, double betab0,double alpha, double * lambda2,double *AcceptB0,_Bool ** Gam){
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
  /*} 
   else {
   *B0=gsl_ran_gamma (rr, p*alpha+alphab0, 1/(kl+betab0));
   *AcceptB0+=1;
   //printf("Bo=%lf",kl+betab0);
   }
   */
  
}


void TauLambda(int l,gsl_rng * rr,int K,int p,double *A,double* B,double B0,double * Tau,_Bool ** Path,double alpha, double * lambda2,double *s2,_Bool ** Gam){
  /* K is the number of pathways
   *  * R is the binary vector for pathway effect
   *   * Path is the grouping information; it is a binary matrix
   *    * B is the pathway effect
   *     * alphab and betab are priors for B effects
   *      * l is the component index
   *       */
  int j,k1;
  for (j=0;j<p;j++){
    //if (fabs(A[j])>0){
    if ((Gam[j][l]==1)&& (fabs(A[j])>0)){
      double mu=sqrt(2*lambda2[j]*s2[j]/pow(A[j],2.0));
      Tau[j]=1.0/inverseGaussian(rr,  mu, 2*lambda2[j]);
      if (Tau[j]<0)
        //printf("TauT=%lf \n",Tau[j]);
        if (((Tau[j]-Tau[j])!=0)|(Tau[j]<0)){
          Tau[j]=1/mu+1/(2*lambda2[j]);
          //printf("J=%d, MU=%lf, lamb=%lf, S2=%lf, Beta=%lf B0=%lf\n", j,mu,lambda2[j],s2[j],A[j],B0);
        }
        double kb=0;
        for (k1=0;k1<K;k1++){
          kb+=Path[j][k1]*B[k1];
        }
        lambda2[j]=gsl_ran_gamma (rr, alpha+1, 1/(B0+kb+Tau[j]));
    } else {
      Tau[j]=gsl_ran_exponential (rr, 1.0/lambda2[j]);
      //printf("Tau=%lf ",Tau[j]);
      //lambda2[j]=gsl_ran_gamma (rr, alpha+1, 1/(Tau[j]));
      lambda2[j]=gsl_ran_gamma (rr, alpha, 1/B0);
    }
  }
}


void GroupEffect(int l,gsl_rng * rr,_Bool rho,int K,int p, _Bool * R,double *A,double* B,double B0,double * Tau,_Bool ** Path,double alphab, double betab,double w,double alpha, double * lambda2,double *s2,double * AcceptR,_Bool **Gam){
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
      _Bool Rknew=0;double Bknew=0;
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
    }// close for k
  } // End for else if (rho==0)
}//close for the function



void LoadAOther(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau,double ** A, double ** U,double ** X,double* s2,_Bool ** Gam){
  int s,s1,j,i;
  for (j=0;j<p;j++){
    for (s=0;s<r;s++){
      A[s][j]=0;
    }
    int NZ1[r];
    int nx1=0;
    findc(r,Gam[j],0,NZ1, &nx1);
    if (nx1>0){
      double * SigmaInv=malloc(nx1*nx1*sizeof(double));
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




void EstimateLoad(gsl_rng * rr,int r, int n,int p, _Bool * rho,double ** Tau,double ** A, double ** U,double ** X,double* s2,_Bool ** Gam){
  int s,s1,j,i;
  for (j=0;j<p;j++){
    for (s=0;s<r;s++){
      A[s][j]=0;
    }
    int NZ1[r];
    int nx1=0;
    findc(r,Gam[j],0,NZ1, &nx1);
    if (nx1>0){
      double * SigmaInv=malloc(nx1*nx1*sizeof(double));
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
      //gsl_vector *result=gsl_vector_alloc (nx1);
      //multivariate_gaussian (rr, mu, &m.matrix,result);
      for (s=0; s<nx1;s++){
        A[NZ1[s]][j]=gsl_vector_get(mu,s);
      }
      free(SigmaInv); gsl_vector_free (mu);
      //gsl_vector_free (result);
    }
  }
}



void proposal(int n,_Bool R[n],_Bool R1[n],float phi, gsl_rng * r) {
  int i;
  for (i=0; i<n;i++) R1[i]=R[i];
  
  int n1=0;//number of 1
  int n0=0;//number of zeros
  int v1[n];
  findc(n,R,1,v1,&n0);//find indices different from 1 i.e ==0;
  int v2[n-n0];
  findc(n,R,0,v2,&n1);//find indices different of zeros i.e ==1
  double u=gsl_ran_flat (r, 0, 1);
  //printf("No==%d, N1==%d\n",n0,n1);
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

void SampleUU(gsl_rng * rr,int r, int n,int Np,int * P,double *** A, double ** U,double *** X,double** s2){
  
  //void SampleUU(gsl_rng * rr,int r, int n,int p0,int p1,int p2,double *** A, double ** U,double *** X,double** s2){	
  
  int s,s1,j,i;
  int sumMark=0;
  for (i=0;i<Np;i++) sumMark+=P[i]; 
  
  double * SigmaInv=malloc(r*r*sizeof(double));
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
  //gsl_linalg_cholesky_decomp (&m.matrix);
  TRY
  {
    //printf("\nError on the sampling of U.");
    gsl_linalg_cholesky_decomp (&m.matrix);
    //THROW;
  }
  CATCH
  {
    printf("\nError on the sampling of U.");
    //exit(1);
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
          if ((k<=j) && (j<k+P[l])) Ax[s]+=X[l][i][j-k]*A[l][s][j-k]/s2[l][j-k];//a+=A[i][s][j-k]*A[i][s1][j-k]/s2[i][j-k];
          k+=P[l];
        }
      }
      //for (j=0; j<p0+p1+p2;j++){
      //if (j<p0)
      //Ax[s]+=X[0][i][j]*A[0][s][j]/s2[0][j];
      //else if (j<p0+p1)
      //Ax[s]+=X[1][i][j-p0]*A[1][s][j-p0]/s2[1][j-p0];
      //else Ax[s]+=X[2][i][j-p1-p0]*A[2][s][j-p1-p0]/s2[2][j-p1-p0];
      //}
    }
    gsl_vector_view b= gsl_vector_view_array (Ax, r);
    //printf("Size 1 is %zu ",(&m.matrix)->size2); printf("Size 2 is %zu ",mu->size);
    gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
    gsl_vector *result=gsl_vector_alloc (r);
    multivariate_gaussian (rr, mu, &m.matrix,result);
    for (s=0; s<r;s++){
      U[i][s]=gsl_vector_get(result,s);
      //if (i<20)
      //printf("%.2lf ",U[i][s]);
    }
    gsl_vector_free (mu);gsl_vector_free (result);
  }
  free(SigmaInv);
}


void SampleU(gsl_rng * rr,int r, int n,int p0,int p1,int p2,double *** A, double ** U,double *** X,double** s2){
  int s,s1,j,i;
  double * SigmaInv=malloc(r*r*sizeof(double));
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
  //gsl_linalg_cholesky_decomp (&m.matrix);
  TRY
  {
    //printf("\nError on the sampling of U.");
    gsl_linalg_cholesky_decomp (&m.matrix);
    //THROW;
  }
  CATCH
  {
    printf("\nError on the sampling of U.");
    //exit(1);
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
    //printf("Size 1 is %zu ",(&m.matrix)->size2); printf("Size 2 is %zu ",mu->size);
    gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
    gsl_vector *result=gsl_vector_alloc (r);
    multivariate_gaussian (rr, mu, &m.matrix,result);
    for (s=0; s<r;s++){
      U[i][s]=gsl_vector_get(result,s);
      //if (i<20)
      //printf("%.2lf ",U[i][s]);
    }
    gsl_vector_free (mu);gsl_vector_free (result);
  }
  free(SigmaInv);
}

//// BIP

// [[Rcpp::export]]
void mainfunction(int *Method1,int * n1,int *P,int *r1,int *Np1,double *datasets,int * IndVar,int *K,int * Paths,int* maxmodel1, int *nbrsample1,int *burninsample1,double *CompoSelMean,double *VarSelMean,double * VarSelMeanGlobal,double * GrpSelMean,double *GrpEffectMean,double * IntGrpMean,double * EstU, double * EstSig2,double *InterceptMean,double *EstLoadMod,double *EstLoad,int *nbrmodel1,double *postgam,double* priorcompsel,double * priorcompselo,double* priorb0,double* priorb,double *priorgrpsel,double *probvarsel){
  
  setvbuf(stdout, NULL, _IONBF, 0);
  int  Method=Method1[0];
  if (Method==1)
    printf("\nThe Method is GroupInfo\n");
  else if (Method==0) {
    printf("The Method is NoGroupInfo\n");
  }
  int i,l,j,k;
  int m;
  clock_t t1 = clock();
  int Np=Np1[0]; //Nber of platforms
  int nbrsample=nbrsample1[0];
  printf("Number of MCMC samples after burn-in is %d\n",nbrsample);
  int burninsample=burninsample1[0];
  printf("Number of burn-in is %d\n",burninsample);
  int n=n1[0];
  printf("Number of samples is %d\n",n);
  //int p0=P[0];
  for (m=0;m<Np;m++){
    printf("Number of markers in platform %d is %d\n",m,P[m]);
  }
  int r=r1[0];
  printf("Number of components is %d\n",r);
  //int P[Np];P[0]=p0;P[1]=p1;;P[2]=p2;
  double *** X=malloc(Np*sizeof(double **));
  double *** X1=malloc(Np*sizeof(double **));
  //double dat=0;
  k=0;
  for (m=0;m<Np;m++){
    X[m]=dmatrix(0,n-1,0,P[m]-1);X1[m]=dmatrix(0,n-1,0,P[m]-1);
    for (i=0;i<n;i++){
      for (j=0;j<P[m];j++){
        X[m][i][j]=X1[m][i][j]=datasets[k];k++;
      }}
  }
  //int K0=K[0];int K1=K[1];int K2=K[2];
  /*
   for (m=0;m<Np;m++){
   if (IndVar[m]!=1)
   printf("Number of groups for platform %d is %d\n",m,K[m]);
   else printf("Number of groups for outcome is 1");
   }
   */
  _Bool *** Path=malloc(Np*sizeof(_Bool **));
  int m1=0;//_Bool pp=0;
  int kk=0;
  for (m=0;m<Np;m++){
    Path[m]=bmatrix(0,P[m]-1,0,K[m]-1);
    if (IndVar[m]==1){
      Path[m][0][0]=1;
    } else if (IndVar[m]==0) {
      for (j=0;j<P[m];j++){
        for (k=0;k<K[m];k++){
          Path[m][j][k]=Paths[kk];kk++;
          //if (m==0) pp=path0[j*K[m]+k]; else if (m==1) pp=path1[j*K[m]+k];
          //Path[m][j][k]=pp;
        }
      }
      m1+=1;
    } else if (IndVar[m]==2){//covariates
      for (j=0;j<P[m];j++){
        for (k=0;k<K[m];k++){
          Path[m][j][k]=1;}}
    } 
  }
  
  long seed=1;
  gsl_rng * rr = gsl_rng_alloc (gsl_rng_rand48);
  gsl_rng_set (rr, seed);
  double ** U=dmatrix(0,n-1,0,r-1);
  double ** meanU=dmatrix(0,n-1,0,r-1);
  for (i=0;i<n;i++){
    for (l=0;l<r;l++){
      U[i][l]=gsl_ran_ugaussian(rr);
      meanU[i][l]=0;
    }
  }
  
  double ***A=malloc(Np*sizeof(double **));
  for (m=0;m<Np;m++){
    A[m]=dmatrix(0,r-1,0,P[m]-1);
  }
  
  /* Hyperparameter*/
  double * q=malloc(Np*sizeof(double));
  double ** qv=malloc(Np*sizeof(double*));
  double ** wg=malloc(Np*sizeof(double*));
  for (m=0;m<Np;m++){
    q[m]=0.5;
    qv[m]=malloc(r*sizeof(double));//
    wg[m]=malloc(r*sizeof(double));//proba. for group selection
  }
  //double q=0.5;
  //double qv=0.5;
  //double a0=0.01;double b0=0.01;
  double a0=1;double b0=1;
  //double alphab=1;double betab=1; 
  double alphab=priorb[0];double betab=priorb[1];
  //double al=1;double bl=1;// Hyper for q
  double al=priorcompsel[0];double bl=priorcompsel[1];// Hyper for q
  //double al0=1;double bl0=1;// Hyper for q in the outcome
  double al0=priorcompselo[0];double bl0=priorcompselo[1];// Hyper for q in the outcome
  //double alv=1;double blv=1;// Hyper for qv
  //double alphab0=2; double betab0=2;
  //double alv=priorselv[0];double blv=priorselv[1];
  double alphab0=priorb0[0]; double betab0=priorb0[1];//Hyper for b0
  //double w=0.5;//priors for group selection
  //double alg=1;double blg=1;
  double alg=priorgrpsel[0];double blg=priorgrpsel[1];
  double alpha=1;
  //Initialization
  
  _Bool ** rhoest=malloc(Np*sizeof(_Bool*));
  double ** rhomean=malloc(Np*sizeof(double*));
  double ** Gamvs=malloc(Np*sizeof(double*));
  //_Bool ** GamvsK=malloc(Np*sizeof(_Bool*));
  _Bool *** R=malloc(Np*sizeof(_Bool**));
  _Bool *** Gam=malloc(Np*sizeof(_Bool**));
  double *** Gammean=malloc(Np*sizeof(double**));
  double *** AcceptR=malloc(Np*sizeof(double**));
  double *** Bmean=malloc(Np*sizeof(double**));
  double *** Rmean=malloc(Np*sizeof(double**));
  double *** B=malloc(Np*sizeof(double**));
  double *** lambda2=malloc(Np*sizeof(double**));
  double *** Tau=malloc(Np*sizeof(double**));
  double *** Taumean=malloc(Np*sizeof(double**));
  double **B0=malloc(Np*sizeof(double*));
  double **B0mean=malloc(Np*sizeof(double*));
  double **AcceptB0=malloc(Np*sizeof(double*));
  double ** quadForm=malloc(Np*sizeof(double*));
  double ** loggauss=malloc(Np*sizeof(double*));
  double ** s2=malloc(Np*sizeof(double*));
  double ** s2Mean=malloc(Np*sizeof(double*));
  int mj=0;
  for (m=0;m<Np;m++){
    Gamvs[m]=malloc(P[m]*sizeof(double));
    //GamvsK[m]=malloc(P[m]*sizeof(_Bool));
    Gam[m]=bmatrix(0,P[m]-1,0,r-1);
    Gammean[m]=dmatrix(0,P[m]-1,0,r-1);
    rhoest[m]=malloc(r*sizeof(_Bool));
    rhomean[m]=malloc(r*sizeof(double));
    R[m]=bmatrix(0,r-1,0,K[m]-1);
    AcceptR[m]=dmatrix(0,r-1,0,K[m]-1);
    Bmean[m]=dmatrix(0,r-1,0,K[m]-1);
    Rmean[m]=dmatrix(0,r-1,0,K[m]-1);
    B[m]=dmatrix(0,r-1,0,K[m]-1);
    lambda2[m]=dmatrix(0,r-1,0,P[m]-1);
    B0[m]=malloc(r*sizeof(double));
    B0mean[m]=malloc(r*sizeof(double));
    AcceptB0[m] =malloc(r*sizeof(double));
    Tau[m]=dmatrix(0,r-1,0,P[m]-1);
    Taumean[m]=dmatrix(0,r-1,0,P[m]-1);
    for (l=0;l<r;l++){
      rhomean[m][l]=0;B0mean[m][l]=0;B0[m][l]=0.1;AcceptB0[m][l]=0;
      for (j=0;j<P[m];j++) Gammean[m][j][l]=0;
      double  uni=gsl_ran_flat (rr, 0, 1);
      qv[m][l]=probvarsel[0];
      wg[m][l]=0.5;//Prior prob group selection
      if (IndVar[m]==1)//Response
        qv[m][l]=0.5;
      else if (IndVar[m]==2)
        qv[m][l]=1;
      if (uni<0.5) {
        rhoest[m][l]=1;
        for (j=0;j<P[m];j++){
          if (gsl_ran_flat (rr, 0, 1)<qv[m][l]) {Gam[m][j][l]=1;A[m][l][j]=0.01;} 
          else {Gam[m][j][l]=0;A[m][l][j]=0;}
        }
      } else {
        rhoest[m][l]=0;
        for (j=0;j<P[m];j++) {Gam[m][j][l]=0;A[m][l][j]=0;}
      }
      if (IndVar[m]==1){//Response
        Gam[m][0][l]=rhoest[m][l];
      } else if (IndVar[m]==2){//Clinical factors
        rhoest[m][l]=1;
        for (j=0;j<P[m];j++) {Gam[m][j][l]=1;A[m][l][j]=0.01;}
      }
      for (k=0;k<K[m];k++){
        Rmean[m][l][k]=0;Bmean[m][l][k]=0;AcceptR[m][l][k]=0;
        R[m][l][k]=0;B[m][l][k]=0;
        //if ((strcmp(Method,"GroupInfo")==0)&& (P[m]!=1)){
        if ((Method==1)&& (IndVar[m]==0)){
          if (rhoest[m][l]==1){
            double  ui=gsl_ran_flat (rr, 0, 1);
            if (ui<wg[m][l]){ R[m][l][k]=1;B[m][l][k]=0.1;}
          }
        }
      }
      for (j=0;j<P[m];j++){
        Tau[m][l][j]=1;Taumean[m][l][j]=0;
        lambda2[m][l][j]=1;
        if (IndVar[m]==2){Tau[m][l][j]=100;}
      }
    }
    s2[m]=malloc(P[m]*sizeof(double));
    s2Mean[m]=malloc(P[m]*sizeof(double)); 
    quadForm[m]=malloc(P[m]*sizeof(double));
    loggauss[m]=malloc(P[m]*sizeof(double));
    for (j=0;j<P[m];j++) {
      loggauss[m][j]=0;
      s2[m][j]=0.1;s2Mean[m][j]=0;
      quadForm[m][j]=Gamvs[m][j]=0;//GamvsK[m][j]=0;
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
  
  /*
   for (m=0;m<Np;m++){
   if (IndVar[m]==2){
   for (j=0;j<P[m];j++) {
   printf("Gam==%d ",Gam[m][j][0]);
   }
   }
   }
   */
  int t;
  int N=nbrsample+burninsample;
  double intercp;
  _Bool ** rhomodel=malloc(nbrsample*sizeof(_Bool*));
  for (t=0;t<nbrsample;t++){
    rhomodel[t]=malloc(dim*sizeof(_Bool));
  }
  //double AUC=0;
  for (t=0;t<N;t++){
    for (m=0;m<Np;m++){
      if (IndVar[m]==1){
        SampleIntercept(rr,n, r, &intercp, s2[m], 100.0,U, A[m], X[m]);
        for (i=0;i<n;i++){
          for (j=0;j<P[m];j++){
            X1[m][i][j]=X[m][i][j]-intercp;
          }}
      }
      double sumrho=0;
      for (l=0;l<r;l++){
        double sumeta=0;
        sumrho+=rhoest[m][l];
        for (j=0;j<P[m];j++){
          sumeta+=Gam[m][j][l];
        }
        //qv[m][l]=gsl_ran_beta (rr, alv+rhoest[m][l]*sumeta, blv+rhoest[m][l]*(P[m]-sumeta));
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
      
      sigma2(rr,  P[m],n,a0, b0,quadForm[m],s2[m]);
      LoadAOther(rr,r, n,P[m],rhoest[m],Tau[m],A[m],U,X1[m],s2[m],Gam[m]);
      
      for (l=0;l<r;l++){
        if (IndVar[m]!=2){
          TauLambda(l,rr,K[m],P[m],A[m][l],B[m][l],B0[m][l],Tau[m][l],Path[m],alpha,lambda2[m][l],s2[m],Gam[m]);
        }
        if ((Method==1)&&(IndVar[m]==0)){
          sumrho=0;
          for (k=0;k<K[m];k++) sumrho+=R[m][l][k];
          //if ((strcmp(Method,"GroupInfo")==0)&&(P[m]!=1))
          wg[m][l]=gsl_ran_beta(rr,alg+rhoest[m][l]*sumrho,blg+rhoest[m][l]*(K[m]-sumrho));
          GroupEffect(l,rr,rhoest[m][l],K[m],P[m], R[m][l],A[m][l],B[m][l],B0[m][l],Tau[m][l],Path[m],alphab, betab,wg[m][l],alpha, lambda2[m][l],s2[m],AcceptR[m][l],Gam[m]);}
        EffectZero(l,rr,K[m],P[m],rhoest[m][l], R[m][l],B[m][l],&B0[m][l],Path[m],alphab0, betab0,alpha, lambda2[m][l],&AcceptB0[m][l],Gam[m]);
      }
      
    }
    //SampleU(rr,r,n,p0,p1,p2,A,U,X1,s2);
    SampleUU(rr,r,n,Np,P,A,U,X1,s2);
    
    
    if (t>=burninsample){
      *InterceptMean+=intercp/nbrsample;
      int rm=0;
      for (m=0;m<Np;m++){
        for (j=0;j<P[m];j++){
          for (l=0;l<r;l++){
            rhomodel[t-burninsample][rm]=Gam[m][j][l];
            rm++;
          }}}
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
        }}
      
      //if (t%(N/20)==1){
      if (t%((N+5)/5)==1){
        //R_CheckUserInterrupt();
        //printf("Intercept for outcome is  %.3lf is =\n",intercp);
        /*
         for (m=0;m<Np;m++){
         //printf("q is %lf \n",q[m]);
         //printf("B0 for platform %d is =\n",m);
         for (l=0;l<r;l++){
         printf("%lf ",B0[m][l]);
         }
         printf("\n");
         }
         printf("\n");
         for (m=0;m<Np;m++){
         printf("qv for platform %d is =\n",m);
         for (l=0;l<r;l++){
         printf("%lf ",qv[m][l]);
         }
         printf("\n");
         }
         */
        printf("\n");
        /*
         for (m=0;m<Np;m++){
         if (IndVar[m]==0)//if (P[m]!=1)
         //printf("Sigma2 for platform %d is=\n",m);
         //else if  (IndVar[m]==1) printf("Sigma2 for outcome is=\n");
         //else if  (IndVar[m]==2) printf("Sigma2 for clin. factor is=\n");
         for (j=0;j<MIN(10,P[m]);j++){
         printf("%lf ",s2[m][j]);
         }
         printf("\n");
         }
         printf("\n");
         */
        /*
         if (Method==1){
         //if (strcmp(Method,"GroupInfo")==0){
         for (m=0;m<Np;m++){
         if (IndVar[m]==0)//if (P[m]!=1)
         printf("Group selection for platform %d is=\n",m);
         else if  (IndVar[m]==1)  printf("Group selection for outcome  is==\n");
         for (l=0;l<r;l++){
         for (k=0;k<K[m];k++){
         printf("%d ",R[m][l][k]);
         }
         printf("\n");
         }
         printf("\n");
         }
         }
         */
        /*
         for (m=0;m<Np;m++){
         if  (IndVar[m]==0)
         //if (P[m]!=1)
         printf("Component selection for platform %d is=\n",m);
         else if  (IndVar[m]==1) printf("Component selection for outcome\n");
         for (l=0;l<r;l++){
         printf("%d ",rhoest[m][l]);
         }
         printf("\n");
         }
         printf("\n");
         
         for (m=0;m<Np;m++){
         if  (IndVar[m]==0)
         //	if (P[m]!=1)
         printf("Variable selection for platform %d \n",m);
         else if  (IndVar[m]==1) printf("Variable selection for outcome\n");
         for (l=0;l<r;l++){
         for (j=0;j<MIN(10,P[m]);j++){
         printf("%d ",Gam[m][j][l]);
         }
         printf("\n");
         }
         printf("\n");
         }
         printf("\n");
         */
        printf("The number of iterations is  %d\n",t);
      }
    }
  }
  printf("\n");
  /*
   
   for (m=0;m<Np;m++){
   if (P[m]!=1)
   printf("AcceptR for platform %d is=\n",m);
   else printf("AcceptR for outcome is=\n");
   for (l=0;l<r;l++){
   for (k=0;k<MIN(K[m],10);k++){
   printf("%.2lf ",AcceptR[m][l][k]/N);
   }
   printf("\n");
   }
   printf("\n");
   }
         }
   for (m=0;m<Np;m++){
   if (P[m]!=1)
   printf("AcceptB0 for platform %d is=\n",m);
   else printf("AcceptB0 for outcome is=\n");
   for (l=0;l<r;l++){
   printf("%.2lf ",AcceptB0[m][l]/N);
   }
   printf("\n");
   }
   printf("\n");
   
   
   
   */
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
      }}
    for (j=0;j<P[m];j++){
      VarSelMeanGlobal[sumP+j]=Gamvs[m][j];
    }
    for (l=0;l<r;l++){
      IntGrpMean[m*r+l]=B0mean[m][l];
    }
    sumP+=P[m];sumK+=K[m];
  }
  
  
  /* Loading estimate for prediction using multiple models*/
  int u=0;
  for (i=0;i<n;i++){
    for (l=0;l<r;l++){
      EstU[u]=meanU[i][l];
      u+=1;
    }
  }
  
  int countmodel=0;int *modelidx=malloc(nbrsample*sizeof(double));
  _Bool **UniqModel=  UniqueModel(nbrsample, dim, rhomodel,modelidx,&countmodel);
  free_bmatrix(rhomodel,0,nbrsample-1,0,dim-1);
  double * logpo=malloc(countmodel*sizeof(double));
  for (t=0;t<countmodel;t++){
    int rm=0;
    int sumg=0;
    for (m=0;m<Np;m++){
      
      if (IndVar[m]==1){
        //		SampleIntercept(rr,n, r, &intercp, s2[m], 100.0,U, A[m], X[m]);
        for (i=0;i<n;i++){
          for (j=0;j<P[m];j++){
            X1[m][i][j]=X[m][i][j]-*InterceptMean;
          }}
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
        for (l=0;l<r;l++){ rhoest[m][l]=Gam[m][0][l];}
      } 
    }
    
    logPostGam(&logpo[t],IndVar,Np,r, n,P,Tau, meanU,X1,s2Mean,rhoest,Gam,qv,q);
  }
  int * highmodelidx=malloc(countmodel*sizeof(int));
  sort(countmodel,logpo,highmodelidx);
  double maxlogpost=logpo[0];
  int maxmodel=maxmodel1[0];
  int nbrmax=MIN(maxmodel,countmodel);
  *nbrmodel1=nbrmax;
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
          Gam[m][j][l]=UniqModel[highmodelidx[t]][rm];rm++;
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
        for (l=0;l<r;l++){ rhoest[m][l]=Gam[m][0][l];}
      }
    }
    //for (l=0;l<r;l++){ rhoest[Np-1][l]=Gam[Np-1][0][l];}
    for (m=0;m<Np;m++){
      EstimateLoad(rr,r, n,P[m],rhoest[m],Taumean[m],A[m],meanU,X1[m],s2Mean[m],Gam[m]);
      for (l=0;l<r;l++){
        for (j=0;j<P[m];j++){
          EstLoadMod[ll]=A[m][l][j];
          ll++;}}
    }
  }
  free(modelidx);free(logpo);free(highmodelidx);
  free_bmatrix(UniqModel,0,countmodel-1,0,dim-1);
  
  
  double thres=0.5;
  ll=0;int ls=0;
  for (m=0;m<Np;m++){
    for (l=0;l<r;l++){
      if (rhomean[m][l]>=thres) rhoest[m][l]=1; else rhoest[m][l]=0; 
      for (j=0;j<P[m];j++){
        if (rhoest[m][l]==1) {
          if (Gammean[m][j][l]>=thres) Gam[m][j][l]=1; else Gam[m][j][l]=0;
        } else {Gam[m][j][l]=0;}
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
    free(rhoest[m]);free(Gamvs[m]);//free(GamvsK[m]);
    free(quadForm[m]);free(loggauss[m]);free(AcceptB0[m]);
    free(qv[m]);free(wg[m]);
  }
  free(Gam);free(Gamvs);//free(GamvsK);
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
  // printf("\nTime taken in hours is %f\n",time_taken/3600);
  
}
# Load R Libraries
pacman::p_load(Rcpp, RcppArmadillo, tidyverse)

# Load C++ implementation of mainfunction
Rcpp::sourceCpp("src/BIP.cpp")

# Define R wrapper
BIP <- function(dataList, IndicVar, groupList=NULL, Method, 
                nbrcomp=4, sample=5000, burnin=1000, nbrmaxmodels=50,
                # BIP priors
                priorcompselv=c(1,1), priorcompselo=c(1,1), priorb0=c(2,2), 
                priorb=c(1,1), priorgrpsel=c(1,1), probvarsel=0.05,
                # BIPmixed design matrices
                Z_family = NULL, Z_site = NULL, 
                # BIPmixed priors
                sigma_ksi_prior=c(1,2), sigma_theta_prior=c(1,2), sigma_prior=c(1,2)) {
  
  if (2 %in% IndicVar) {
    covariates_included <- T
  } else {
    covariates_included <- F
  }
  
  # Assume 1 covar so mainfunction compiles no matter the Method
  n_covars <- 1
  mu_beta <- rep(0, n_covars)
  mu_prior_var <- 100
  beta_prior_var <- rep(100, n_covars)

  n <- sapply(dataList, nrow) %>% unique
  if (length(n) != 1) {
    stop("Each element of dataList must have a number of rows equal to the number of observations.")
  }
  
  if (Method != "BIPmixed") {
    # Matrices required for mainfunction compilation
    Z_family <- Z_site <- W <- matrix(1, nrow = nrow(dataList[[1]]), ncol = 1)
  }
  
  if (Method == "BIPmixed") {
    if (is.null(Z_family) | is.null(Z_site)) {
      stop("If Method is BIPmixed, you must provide design matrices Z_family and Z_site.")
    }
    if (covariates_included) {
      # Update fixed effect priors if covars present & Method is BIPmixed
      covar_index <- which(IndicVar == 2)
      n_covars <- ncol(dataList[[covar_index]])
      mu_beta <- rep(0, n_covars)
      mu_prior_var <- 100
      beta_prior_var <- rep(100, n_covars)
      # Remove covariates from dataList & IndicVar
      W <- dataList[[covar_index]]
      dataList <- dataList[-covar_index]
      IndicVar <- IndicVar[-covar_index]
    } else {
      W <- matrix(1, nrow = nrow(dataList[[1]]), ncol = 1)
    }
  }
  
  # Map families to study site (This needs to happen regardless of Method since passed to mainfunction)
  Z_family_to_site <- t(Z_family) %*% Z_site
  
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
    meth = 0
  } else if (Method=="BIPnet") {
    meth = 1
  } else if (Method=="BIPmixed") {
    meth = 2
  } else {
    stop("You must provide Method = BIP, BIPnet, or BIPmixed")
  }
  
  Np=length(IndicVar)
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
    Method = as.integer(meth), covariates_included = covariates_included, n = as.integer(n), P = as.integer(P), r = as.integer(nbrcomp), Np = as.integer(Np),
    datasets = as.double(datasets), IndVar = as.integer(IndicVar), K=as.integer(K), Paths = as.integer(Paths),
    maxmodel = as.integer(nbrmaxmodels), nbrsample = as.integer(sample), burninsample = as.integer(burnin),
    CompoSelMean = as.double(rep(0,Np*nbrcomp)), VarSelMean = as.double(rep(0,nbrcomp*sum(P))),
    VarSelMeanGlobal = as.double(rep(0,sum(P))), GrpSelMean = as.double(rep(0,nbrcomp*sum(K))),
    GrpEffectMean = as.double(rep(0,nbrcomp*sum(K))), IntGrpMean = as.double(rep(0,nbrcomp*Np)),
    EstU = as.double(rep(0,n*nbrcomp)), EstSig2 = as.double(rep(0,sum(P))), InterceptMean=as.double(rep(0,1)),
    EstLoadMod = as.double(rep(0,nbrmaxmodels*nbrcomp*sum(P))), EstLoad=as.double(rep(0,nbrcomp*sum(P))),
    nbrmodel = as.integer(0), postgam = rep(0,nbrmaxmodels), priorcompsel = priorcompselv,
    priorcompselo = priorcompselo, priorb0 = priorb0, priorb = as.double(priorb),
    priorgrpsel = priorgrpsel, probvarsel = as.double(probvarsel),
    Z_family, Z_site, Z_family_to_site,
    mu_prior_var, mu_beta, beta_prior_var, W,
    sigma_ksi_prior_a = sigma_ksi_prior[1],
    sigma_ksi_prior_b = sigma_ksi_prior[2],
    sigma_theta_prior_a = sigma_theta_prior[1],
    sigma_theta_prior_b = sigma_theta_prior[2],
    sigma_prior_a = sigma_prior[1],
    sigma_prior_b = sigma_prior[2]
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
  
  return (list(Method = as.integer(meth),
               EstU=EstimateU,VarSelMean=VarSelMean,VarSelMeanGlobal=VarSelMeanGlobal,
               CompoSelMean=CompoSelMean,GrpSelMean=GrpSelMean,GrpEffectMean=GrpEffectMean,
               IntGrpMean=IntGrpMean,EstLoad=EstLoad,EstLoadModel=EstLoadModel,
               nbrmodel=result$nbrmodel1,EstSig2=EstSig2,EstIntcp=result$intercepts$InterceptMean,
               PostGam=result$postgam,IndicVar=IndicVar,nbrcomp=nbrcomp,MeanData=MeanData,SDData=SD,
               mu_samples = result$samples$mu_samples,
               #alpha_samples = result$alpha_samples,
               beta_samples = result$samples$beta_samples,
               gamma_samples = result$samples$gamma_samples,
               ksi_samples = result$samples$ksi_samples,
               theta_samples = result$samples$theta_samples,
               sigma2_ksi_samples = result$samples$sigma2_ksi_samples,
               sigma2_theta_samples = result$samples$sigma2_theta_samples,
               sigma2_samples = result$samples$sigma2_samples,
               initial_values = result$initial_values,
               # sigma2_non_outcome_samples = result$sigma2_non_outcome_samples,
               beta_mean = result$intercepts$beta_mean, 
               theta_mean = result$intercepts$theta_mean
  )
  )
}


library(MASS)
BIPpredict=function(dataListNew=dataListNew, 
                    Result=Result, meth="BMA",
                    # TODO include Wnew in dataListNew and peel out
                    # TODO specify the design matrices to be for new individuals
                    # TODO add checks to ensure design matrices conform to desired dimensions
                    Wnew = NULL, Z_family_to_site=NULL, Z_family=NULL) {
  IndicVar=Result$IndicVar
  nbrcomp=Result$nbrcomp
  X_new=list();
  np=which(IndicVar==1)
  Np=length(IndicVar);

  print("Scaling new data...")
  j=1;
  for (i in  1:Np){
    X_new[[i]]=NULL;
    if (IndicVar[i]!=1){ ## We do not center the response variable
      X_new[[i]]=dataListNew[[j]];
      X_new[[i]]=t((t(X_new[[i]])-Result$MeanData[[i]])/Result$SDData[[i]])
      j=j+1;
    }
  }
  
  nbrmodel=Result$nbrmodel
  EstLoad=Result$EstLoadModel
  N_obs <- dataListNew %>% sapply(nrow) %>% unique()
  if (length(N_obs)!=1) { stop("dataListNew must have matrices with equal numbers of rows.") }
  Upredict=matrix(0, N_obs, nbrcomp)
  ypredict=rep(0, N_obs)
  nbrmodel1=nbrmodel
  if (meth!="BMA"){
    nbrmodel1=1;
  }
  
  print("Estimating loadings...")
  for (nb in 1:nbrmodel1){
    SelLoadPacked=NULL;SelLoadPackedFull=NULL
    SelVarXnew=NULL;SelVarXnewFull=NULL;
    Sig2=NULL;Sig2Full=NULL;
    if (meth=="BMA"){
      EstL=EstLoad[[nb]];
      pb=Result$PostGam[nb];
    } else {
      EstL=Result$EstLoad; 
      pb=1;
    }
    for (m in 1:Np){
      if (IndicVar[m]!=1){ # We exclude the response variable
        nc=nrow(EstL[[m]]);
        selvar=apply(EstL[[m]],2, function(x) length((which(x==0)))!=nc)
        SelLoadPacked=cbind(SelLoadPacked,EstL[[m]][,selvar]);
        Sig2=c(Sig2,Result$EstSig2[[m]][selvar])
        SelVarXnew=cbind(SelVarXnew,X_new[[m]][,selvar])
      }
    }
    SelLoadPackedFull=SelLoadPacked;
    Sig2Full=Sig2;
    SelVarXnewFull=SelVarXnew;
    SelLoadPackedFull=as.matrix(SelLoadPackedFull);
    Sigma2Full=solve(SelLoadPackedFull%*%diag(1/Sig2Full)%*%t(SelLoadPackedFull)+diag(nrow(SelLoadPackedFull)));
    Upredict1=t(Sigma2Full%*%SelLoadPackedFull%*%diag(1/Sig2Full)%*%t(SelVarXnewFull))
    Upredict=Upredict+Upredict1*pb;
    ypredict=ypredict+as.matrix(Upredict1%*%EstL[[np]])*pb;
  }
  
  if (Result$Method==2) {
    ypredict=as.vector(ypredict)
    Wbeta <- Wnew %*% Result$beta_mean
    for (s in 1:N_sites) {
      families_in_s <- which(Z_family_to_site[, s] != 0)
      for (f in families_in_s) {
        individuals_in_f <- which(Z_family[, f] == 1)
        for (ind in individuals_in_f) {
            # Add in theta_mean and beta_mean
            ypredict[ind] <- ypredict[ind] + Wbeta[ind] + Result$theta_mean[f]
          }
        }
      }
  } else {
    ypredict=as.vector(ypredict)+Result$EstIntcp
  }
  return (list(ypredict=as.vector(ypredict)+Result$EstIntcp,Upredtest=Upredict))
}

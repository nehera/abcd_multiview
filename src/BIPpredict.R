# Define BIPpredict function
BIPpredict <- function(dataListNew=dataListNew, 
                    Result=Result, meth="BMA",
                    Z_site=NULL) {
  
  print("Getting predictions...")
  
  IndicVar <- Result$IndicVar
  
  if (Result$Method == 2) {
    # Peel out covariates in the case that Method= BIPmixed & covariates included
    covariate_index <- which(IndicVar==2)
    if (length(covariate_index) == 1) {
      Wnew <- dataListNew[[covariate_index]]
      dataListNew <- dataListNew[-covariate_index]
      IndicVar <- IndicVar[-covariate_index]
    } else {
      Wnew <- matrix(0, nrow = nrow(Z_site), ncol = 1)
    }
  }
  
  nbrcomp <- Result$nbrcomp
  np <- which(IndicVar == 1)
  
  # Standardize using the mean and variance from the training data
  mean_list <- Result$MeanData[IndicVar!=1]
  sd_list <- Result$SDData[IndicVar!=1]
  X_new <- dataListNew
  for (m in 1:length(dataListNew)) {
    X_new_m <- dataListNew[[m]]
    for (j in 1:ncol(X_new_m)) {
      X_new_m[, j] <- ( X_new_m[, j] - mean_list[[m]][j] ) / sd_list[[m]][j]
    }
    X_new[[m]] <- X_new_m
  }
  
  # Get the number of platforms
  Np <- length(X_new)
  
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
    Sigma2Full=solve( SelLoadPackedFull %*% diag(1/Sig2Full) %*% t(SelLoadPackedFull) + diag(nrow(SelLoadPackedFull)))
    Upredict1=t(Sigma2Full%*%SelLoadPackedFull%*%diag(1/Sig2Full)%*%t(SelVarXnewFull))
    Upredict=Upredict+Upredict1*pb;
    ypredict=ypredict+as.matrix(Upredict1%*%EstL[[np]])*pb;
  }
  
  if (Result$Method==2) {
    # Adjust the predictions by random and fixed effect estimates
    ypredict <- as.vector(ypredict)
    Wbeta <- Wnew %*% Result$beta_mean
    for (i in 1:length(ypredict)) {
      
      site_i <- which(Z_site[i, ] == 1)
      ypredict[i] <- ypredict[i] + Wbeta[i] + Result$ksi_mean[site_i] # we use the site effect
      
    }
  } else {
    ypredict <- as.vector(ypredict) + Result$EstIntcp
  }
  return (list(ypredict=ypredict, Upredtest=Upredict))
}

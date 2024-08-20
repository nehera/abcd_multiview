# Load libraries
pacman::p_load(MASS)

# Define BIPpredict function
BIPpredict=function(dataListNew=dataListNew, 
                    Result=Result, meth="BMA",
                    Z_site=NULL, Z_family=NULL) {
  
  IndicVar <- Result$IndicVar
  
  if (Result$Method==2) {
    # Peel out covariates in the case that Method= BIPmixed & covariates included
    covariate_index <- which(IndicVar==2)
    if (length(covariate_index) == 1) {
      Wnew <- dataListNew[covariate_index]
      dataListNew <- dataListNew[-covariate_index]
      IndicVar <- IndicVar[-covariate_index]
    } else {
      Wnew <- matrix(0, nrow = nrow(Z_site), ncol = 1)
    }
    # Map families to study site (This needs to happen regardless of Method since passed to mainfunction)
    Z_family_to_site <- t(Z_family) %*% Z_site
    N_sites <- ncol(Z_site)
  }
  
  nbrcomp <- Result$nbrcomp
  X_new <- list()
  np <- which(IndicVar==1)
  Np <- length(IndicVar)
  
  print("Scaling new data...")
  
  # Normalize each matrix in the list
  X_new <- lapply(1:Np, function(i) {
    if (IndicVar[i] != 1) {
      # Subtract the mean and divide by the standard deviation
      (dataListNew[[i]] - matrix(Result$MeanData[[i]], nrow = nrow(dataListNew[[i]]), ncol = ncol(dataListNew[[i]]), byrow = TRUE)) /
        matrix(Result$SDData[[i]], nrow = nrow(dataListNew[[i]]), ncol = ncol(dataListNew[[i]]), byrow = TRUE)
    } else {
      # If IndicVar[i] == 1, return the matrix as is
      dataListNew[[i]]
    }
  })
  
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
    ypredict <- as.vector(ypredict) + Result$EstIntcp
  }
  return (list(ypredict=ypredict, Upredtest=Upredict))
}

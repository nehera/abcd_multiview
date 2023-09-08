n_views=2; n_obs=200; p_m=10 
prob_feature_importance=0.5; r=4
prob_component_importance=0.5
seed=1
library(tidyverse)

simulate_iid_data <- function(n_views=2, n_obs=200, p_m=10, r=4,
                              prob_feature_importance=0.5, 
                              prob_component_importance=0.5,
                              seed=1) {
  
  set.seed(seed)
  n_important_features <- floor(prob_feature_importance*p_m)
  n_important_components <- floor(prob_component_importance*r)

  simulate_A <- function(n_important_components, n_important_features) {
    index_important_components <- seq(to = n_important_components)
    index_important_features <- seq(to = n_important_features)
    # simulate_A assumes all features important across active components
    n_nonzero_a <- n_important_components * n_important_features 
    nonzero_a <- matrix(rnorm(n_nonzero_a), 
                        nrow = n_important_components, 
                        ncol = n_important_features)
    A <- matrix(0, nrow = r, ncol = p_m) 
    A[index_important_components, index_important_features] <- nonzero_a
    return(A)
  }
  
  A_list <- list()
  for (m in 1:n_views) {
    A_list[[m]] <- simulate_A(n_important_components, n_important_features)
  }
  
  U <- matrix(data = rnorm(n_obs*r), nrow = n_obs, ncol = r)
  
  E_list <- list()
  for (m in 1:n_views) {
    E_list[[m]] <- matrix(data = rnorm(n_obs*p_m), nrow = n_obs, ncol = p_m)
  }
  
  X_list <- list()
  for (m in 1:n_views) {
    X_list[[m]] <- U %*% A_list[[m]] + E_list[[m]]
  }

  a <- ifelse(1:r %in% index_important_components, 1, 0) %>%
    matrix(nrow = r)
  e <- matrix(rnorm(n_obs), nrow = n_obs)
    
  Y <- U %*% a + e
  
  simulation_results <- list(X_list=X_list, Y=Y,
                             index_important_components=seq(n_important_components),
                             index_important_features=seq(n_important_features),
                             U=U, A_list=A_list)
  
  return(simulation_results)
}


# Note, at present, this version only is functional for scenario 1, setting 1 at p != 500
simulate_data <- function(n=200, p=10, prob_feature_importance=0.4, r=4, scenario=1, setting=NULL, overlap=NULL, seed=1) {
  
  # TODO breakout scenario 1 from scenarios 2 and 3 to simplify arguments required.
  # scenario: We have three scenarios as described in the paper. Scenario 1 has 5 settings, and scenarios 2 and 3 have both two choices: overlap and no overlap
  # setting: Settings for scenario 1 only. It's a total of 4 settings
  # overlap: If yes, then there are some overlap between components, otherwise, there is no overlap. This input only works wen cenario is either 2 or 3
  # seed: Seed to generate random numbers 
  
  if ((scenario==1) && is.null(setting)){
    stop("You must specify the setting for scenario 1")
  }
  
  if (((scenario==2) || (scenario==3)) && (is.null(overlap))){
    stop("You must specify overlap for scenario 2 or 3. For the overlap argument, the possible options are 'yes' or 'no'.")
  }
  
  n_important_features <- floor(p * prob_feature_importance)
  n_unimportant_features <- p - n_important_features
  n_important_groups <- 1 # TODO remove hardcoding of n_important_groups
  n_features_in_group <- n_important_features # TODO remove hardcoding of n_features_in_group
  
  Group1=matrix(0, p, n_important_groups+1)
  m <- matrix(1, nrow = n_features_in_group)
  Group1[1:n_important_features, 1:n_important_groups] <- kronecker(diag(1, n_important_groups), m) # TODO share kronecker operation with Thierry
  Group1[(n_important_features+1):p, n_important_groups+1] <- 1
  
  Group2=Group1
  
  #In general, you can generate N random numbers in the interval (a,b)
  #with the formula r = a + (b-a).*rand(N,1).
  a1=-.5
  b1=-.3
  a2=.5
  b2=.3
  #latent component
  # TODO the division by 2 in the subsetting below means the n_important_features must be an even number... The default prob_feature_importance has been set to 0.4 to circumvent this issue for now.
  V1a=V1b=matrix(0,n_important_features/2,r) 
  V2a=V2b=matrix(0,n_important_features/2,r)
  for (j in 1:r){
    V1a[,j]=runif(n_important_features/2,a1,b1);
    V1b[,j]=runif(n_important_features/2,b2,a2);
    V11=rbind(V1a,V1b);
    # for second dataset;
    V2a[,j]=runif(n_important_features/2,a1,b1);
    V2b[,j]=runif(n_important_features/2,b2,a2);
    V22=rbind(V2a,V2b);
  }
  
  #multiply effect size for main networks by 2;
  for (j in seq(1, n_important_features, by = n_features_in_group)){
    V11[j,]=2*V11[j,];
    V22[j,]=2*V22[j,];
  }

  Sx <- diag(nrow = n_features_in_group, ncol = n_features_in_group)
  Sx[2:n_features_in_group, 1] <- 0.7
  Sx[1, 2:n_features_in_group] <- 0.7
  diagonal_matrix_off_diagonal_specified <- function(nrow, ncol, off_diagonal, diagonal=1) {
    # Create an empty matrix with NA values
    m <- matrix(NA, nrow, ncol)
    # Fill the lower triangular part with 0.5
    m[lower.tri(m)] <- off_diagonal
    # Fill the upper triangular part with 0.5
    m[upper.tri(m)] <- off_diagonal
    # Fill the diagonal with 1
    diag(m) <- 1
    # Return the matrix
    return(m)
  }
  Sx[2:n_features_in_group, 2:n_features_in_group] <- diagonal_matrix_off_diagonal_specified(n_features_in_group-1, n_features_in_group-1, 0.7^2)
  
  Sigma1=Matrix::bdiag(rep(list(Sx), n_important_groups))
  s=nrow(Sigma1)
  Sigma1=Matrix::bdiag(Sigma1,diag(p-s))
  Sigma2=Sigma1
  
  if (scenario==1) {
    
    if (setting==1){
      # there are 11 groups; groups 1-10 are signal variables, group 11 are noise variables
      # generate a p by 11 indicator matrix where there is 1 in group k if variable j is in that group.
      
      # orthonormalize;
      #A1=[orth(V11) ;zeros(p-n_important_features,4)];
      #A2=[orth(V22) ;zeros(p-n_important_features,4)];
      
      #not orthonomalize
      A1=rbind(V11 ,matrix(0,p-n_important_features,4));
      A2=rbind(V22 ,matrix(0,p-n_important_features,4)); # TODO REMOVE p and n_important_features
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      
      TrueVariables1=rep(0,p);
      TrueVariables1[abs(A1[,1])!=0]=1;
      
      TrueVariables2=rep(0,p);
      TrueVariables2[abs(A2[,1])!=0]=1;
    } else if (setting==2) {
      #%there are 11 groups; first 3 groups corresponding to the first 30
      # %variables contribute to associatin between X1 and X2, and to the outcome.
      #%The remaining groups are noise. 
      #%orthonormalize;
      #%A1=[orth(V11(1:30,:)) ;zeros(p-30,4)];
      #%A2=[orth(V22(1:30,:)) ;zeros(p-30,4)];
      
      #%not orthonomalize
      A1=rbind(V11[1:30,] ,matrix(0,p-30,4));
      A2=rbind(V22[1:30,] ,matrix(0,p-30,4));
      
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      
      TrueVariables1=rep(0,p);
      TrueVariables1[abs(A1[,1])!=0]=1;
      
      TrueVariables2=rep(0,p);
      TrueVariables2[abs(A2[,1])!=0]=1;
      
      
    } else if (setting==3){
      randvec=NULL;
      #at most 5 random varibles in each group do not contribute to A;
      for (j in 1:10){
        randvec=c(randvec,sample((10*(j-1)+1):(10*j),5,1)); 
        #randvec=[randvec;randi([10*(j-1)+1,10*j],5,1)]; 
      }
      V11[randvec,]=0;
      V22[randvec,]=0;
      
      
      A1=rbind(V11 ,matrix(0,p-nrow(V11),4));
      A2=rbind(V22 ,matrix(0,p-nrow(V22),4));
      
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      
      TrueVariables1=rep(0,p);
      TrueVariables1[abs(A1[,1])!=0]=1;
      
      TrueVariables2=rep(0,p);
      TrueVariables2[abs(A2[,1])!=0]=1;
      
      
    } else if (setting==4){
      randvec=NULL;
      for (j in 1:3){
        randvec=c(randvec,sample((10*(j-1)+1):(10*j),5,1));       
      }
      
      A1=rbind(V11[1:30,] ,matrix(0,p-30,4));
      A2=rbind(V22[1:30,] ,matrix(0,p-30,4));
      
      V11=V11[1:30,] ;
      V22=V22[1:30,];
      
      V11[randvec,]=0;
      V22[randvec,]=0;
      
      
      A1=rbind(V11 ,matrix(0,p-nrow(V11),4));
      A2=rbind(V22 ,matrix(0,p-nrow(V22),4));
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      
      TrueVariables1=rep(0,p);
      TrueVariables1[abs(A1[,1])!=0]=1;
      
      TrueVariables2=rep(0,p);
      TrueVariables2[abs(A2[,1])!=0]=1;
      
      
    } else if (setting==5) {
      #%100 important variables but 25 each loaded on first to fourth factors
      
      
      
      A1=rbind(V11 ,matrix(0,p-n_important_features,4));
      A2=rbind(V22 ,matrix(0,p-n_important_features,4));
      
      whic=matrix(0,p,4);
      whic[1:25,1]=1;
      whic[26:50,2]=1;
      whic[56:75,3]=1; 
      whic[76:100,4]=1;
      A1[whic[,1]!=1,1]=0;
      A1[whic[,2]!=1,2]=0;
      A1[whic[,3]!=1,3]=0;
      A1[whic[,4]!=1,4]=0;
      
      A2[whic[,1]!=1,1]=0;
      A2[whic[,2]!=1,2]=0;
      A2[whic[,3]!=1,3]=0;
      A2[whic[,4]!=1,4]=0;
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      
      
      TrueVariables1=rep(0,p);
      TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
      #TrueVariables1(sum(abs(A1),2)~=0)=1;
      
      TrueVariables2=rep(0,p);
      TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
    } 
  } else if (scenario==2){
    
    A1=rbind(V11 ,matrix(0,p-n_important_features,4));
    A2=rbind(V22 ,matrix(0,p-n_important_features,4)); 
    if (overlap=="yes"){
      # %overlap
      whic=matrix(0,p,4);
      whic[1:100,1]=1;
      whic[1:100,2]=1;
      whic[1:100,3]=1;
      whic[1:100,4]=1;
      A1[whic[,1]!=1,1]=0;
      A1[whic[,2]!=1,2]=0;
      A1[,3:4]=0;
      
      
      A2[,1:2]=0;
      A2[whic[,3]!=1,3]=0;
      A2[whic[,4]!=1,4]=0;
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      TrueVariables1=rep(0,p);
      TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
      TrueVariables2=rep(0,p);
      TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
    } else if (overlap=="no"){
      #%No overlap    
      whic=matrix(0,p,4);
      whic[1:50,1]=1;
      whic[51:100,2]=1;
      whic[1:50,3]=1;
      whic[51:100,4]=1;
      
      A1[whic[,1]!=1,1]=0;
      A1[whic[,2]!=1,2]=0;
      A1[,3:4]=0;
      
      A2[,1:2]=0;
      A2[whic[,3]!=1,3]=0;
      A2[whic[,4]!=1,4]=0;
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      TrueVariables1=rep(0,p);
      TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
      TrueVariables2=rep(0,p);
      TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
    }
    
  } else if (scenario==3){
    
    A1=rbind(V11 ,matrix(0,p-n_important_features,4));
    A2=rbind(V22 ,matrix(0,p-n_important_features,4)); 
    if (overlap=="no"){
      whic=matrix(0,p,4);
      whic[1:25,1]=1;
      whic[26:50,2]=1;
      whic[51:100,3]=1;
      whic[51:100,4]=1;
      A1[whic[,1]!=1,1]=0;
      A1[whic[,2]!=1,2]=0;
      A1[whic[,3]!=1,3]=0;
      A1[,4]=0;
      
      A2[whic[,1]!=1,1]=0;
      A2[whic[,2]!=1,2]=0;
      A2[whic[,3]!=1,4]=0;
      A2[,3]=0;
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
    } else if (overlap=="yes"){
      whic=matrix(0,p,4);
      whic[1:50,1]=1;
      whic[1:50,2]=1;
      whic[51:100,3]=1;
      whic[51:100,4]=1;
      A1[whic[,1]!=1,1]=0;
      A1[whic[,2]!=1,2]=0;
      A1[whic[,3]!=1,3]=0;
      A1[,4]=0;
      A2[whic[,1]!=1,1]=0;
      A2[whic[,2]!=1,2]=0;
      A2[whic[,3]!=1,4]=0;
      A2[,3]=0;
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      
    }
    
    TrueVariables1=rep(0,p);
    TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
    TrueVariables2=rep(0,p);
    TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
  }
  ### Generate data
  set.seed(seed)
  sigma2u=1
  U=MASS::mvrnorm(n,mu=rep(0,r),sigma2u*diag(r))
  A3=c(1, 0, 1, 1) # components 1, 3, and 4 affect response
  E1=MASS::mvrnorm(n, mu=rep(0,p), Sigma1);
  E2=MASS::mvrnorm(n, mu=rep(0,p), Sigma2);
  E3=MASS::mvrnorm(n, mu=1, sigma2u);
  X1=U%*%t(A1)+E1;X2=U%*%t(A2)+E2;Y=U%*%A3+E3
  list(X1=X1,X2=X2,Y=Y,TrueVar1=TrueVariables1,TrueVar2=TrueVariables2,U=U,LoadA1=A1,LoadA2=A2,Group1=Group1,Group2=Group2)
}

# Start with iid data
simulation_results <- simulate_iid_data()

data_list <- list(simulation_results$X_list[[1]], 
                  simulation_results$X_list[[2]], 
                  simulation_results$Y)
indic_var <- c(0, 0, 1) 
method <- "BIP" # method without grouping information to start
group_list <- NULL # default when grouping information not included

## Scale data
# TODO: Breakout scaling into separate step & choose whether or not to scale by the Frobenius norm to account for differing n_features
data_list <- lapply(data_list, scale)

## Get results
bip_0 <- BIPnet::BIP(dataList = data_list, IndicVar = indic_var, Method = method)

## Store simulated data and results for development/ testing
today_date <- Sys.Date()
fname_results <- paste0("data/", today_date, "_simulation_results.rds")
saveRDS(simulation_results, fname_results)
fname_data <- paste0("data/", today_date, "_simulation_data_list.rds")
saveRDS(data_list, fname_data)
fname_data <- paste0("data/", today_date, "_simulation_BIP_results.rds")
saveRDS(bip_0, fname_data)
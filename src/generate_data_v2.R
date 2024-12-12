# This is Thierry's data simulation function with the following edits from Apostolos on Oct 31, 2024:
# I kept Thierry’s scenarios and included the RE variances as function arguments. 
# I updated Thierry’s scenarios 2-3 just in case, instead of simply removing them, 
# but I’m honestly not sure if my updates in that part of the code make sense. 
# Finally, I kept r=4, because I didn’t feel like gambling with the dimensions of his matrices and vectors. 
# He defines r=4, but he also hardcodes matrix dimensions involving r.

Simulate2 <- function(
    scenario=1, setting=NULL, overlap=NULL, 
    seed=1, development=T,
    sigma2_ksi_true=0, 
    sigma2_theta_true=0,
    ksi=NULL # vector of site-effects only applicable to the test set
                      ) {
  
  # scenario: We have three scenarios as described in Thierry and Sandra's Biostatistics paper. 
  # Scenario 1 has 5 settings, and scenarios 2 and 3 have both two choices: overlap and no overlap
  # setting: Settings are for scenario 1 only. There's a total of 4 possible settings
  # overlap: If yes, then there are some overlap between components, otherwise, there is no overlap. 
  # This input only works when scenario is either 2 or 3
  # seed: Seed to generate random numbers 
  
  if ((scenario==1) && is.null(setting)){
    stop("You must specify the setting for scenario 1")
  }
  
  if (((scenario==2) || (scenario==3)) && (is.null(overlap))){
    stop("You must specify overlap for scenario 2 or 3. It's yes or no for overlap")
  }
  
  P1=500;
  P2=500;
  P3=500;
  P4=500;
  
  r=4; #%4 latent components;
  p1=P1/5;
  p2=P2/5;
  p3=P3/5;
  p4=P4/5;
  
  Group1=matrix(0,P1, 11);
  Group1[1:10,1]=1;
  Group1[11:20,2]=1;
  Group1[21:30,3]=1;
  Group1[31:40,4]=1;
  Group1[41:50,5]=1;
  Group1[51:60,6]=1;
  Group1[61:70,7]=1;
  Group1[71:80,8]=1;
  Group1[81:90,9]=1;
  Group1[91:100,10]=1;
  Group1[101:P1,11]=1;
  Group4=Group3=Group2=Group1;
  
  #In general, you can generate N random numbers in the interval (a,b)
  #with the formula r = a + (b-a).*rand(N,1).
  a1=-.5;
  b1=-.3;
  a2=.5;
  b2=.3;
  
  #latent component
  V1a=V1b=matrix(0,p1/2,r);
  V2a=V2b=matrix(0,p2/2,r);
  V3a=V3b=matrix(0,p3/2,r);
  V4a=V4b=matrix(0,p4/2,r);
  for (j in 1:r){
    # For 1st view
    V1a[,j]=runif(p1/2,a1,b1);
    V1b[,j]=runif(p1/2,b2,a2);
    V11=rbind(V1a,V1b);
    # For 2nd view
    V2a[,j]=runif(p2/2,a1,b1);
    V2b[,j]=runif(p2/2,b2,a2);
    V22=rbind(V2a,V2b);
    # For 3rd view
    V3a[,j]=runif(p3/2,a1,b1);
    V3b[,j]=runif(p3/2,b2,a2);
    V33=rbind(V3a,V3b);
    # For 4th view
    V4a[,j]=runif(p4/2,a1,b1);
    V4b[,j]=runif(p4/2,b2,a2);
    V44=rbind(V4a,V4b);
  }
  
  #multiply effect size for main networks by 2;
  for (j in c(1,11, 21, 31, 41, 51, 61, 71, 81, 91)){
    V11[j,]=2*V11[j,];
    V22[j,]=2*V22[j,];
    V33[j,]=2*V33[j,];
    V44[j,]=2*V44[j,];
  }
  
  Sx=matrix(1,10,10);
  for (j in 2:10){
    Sx[1,j]=.7;
    Sx[j,1]=.7;
    for (ii in 3:10){
      Sx[ii,j]=.7^2;
      Sx[j,ii]=.7^2;
      Sx[ii,ii]=1;
    }
  }
  Sigma1=Matrix::bdiag(Sx,Sx,Sx,Sx,Sx,Sx,Sx,Sx,Sx,Sx);
  s=nrow(Sigma1);
  Sigma1=Matrix::bdiag(Sigma1,diag(P1-s));
  Sigma4=Sigma3=Sigma2=Sigma1;
  
  if (scenario==1) {
    
    if (setting==1){
      
      # there are 11 groups; groups 1-10 are signal variables, group 11 are noise variables;
      # generate a P1 by 11 indicator matrix where there is 1 in group k if
      # variable j is in that group;
      
      A1=rbind(V11 ,matrix(0,P1-p1,4));
      A2=rbind(V22 ,matrix(0,P2-p2,4));
      A3=rbind(V11 ,matrix(0,P3-p3,4));
      A4=rbind(V22 ,matrix(0,P4-p4,4));
      
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      A3[abs(A3)<=10^-8]=0;
      A4[abs(A4)<=10^-8]=0;
      
      TrueVariables1=rep(0,P1);
      TrueVariables1[abs(A1[,1])!=0]=1;
      
      TrueVariables2=rep(0,P2);
      TrueVariables2[abs(A2[,1])!=0]=1;
      
      TrueVariables3=rep(0,P3);
      TrueVariables3[abs(A3[,1])!=0]=1;
      
      TrueVariables4=rep(0,P4);
      TrueVariables4[abs(A4[,1])!=0]=1;
    } else if (setting==2) {
      
      # there are 11 groups; first 3 groups corresponding to the first 30
      # variables contribute to associatin between X1 and X2, and to the outcome.
      # The remaining groups are noise. 
      
      A1=rbind(V11[1:30,] ,matrix(0,P1-30,4));
      A2=rbind(V22[1:30,] ,matrix(0,P2-30,4));
      A3=rbind(V33[1:30,] ,matrix(0,P3-30,4));
      A4=rbind(V44[1:30,] ,matrix(0,P4-30,4));
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      A3[abs(A3)<=10^-8]=0;
      A4[abs(A4)<=10^-8]=0;
      
      TrueVariables1=rep(0,P1);
      TrueVariables1[abs(A1[,1])!=0]=1;
      
      TrueVariables2=rep(0,P2);
      TrueVariables2[abs(A2[,1])!=0]=1;
      
      TrueVariables3=rep(0,P3);
      TrueVariables3[abs(A3[,1])!=0]=1;
      
      TrueVariables4=rep(0,P4);
      TrueVariables4[abs(A4[,1])!=0]=1;
      
    } else if (setting==3){
      
      randvec=NULL;
      #at most 5 random varibles in each group do not contribute to A;
      
      for (j in 1:10){
        randvec=c(randvec,sample((10*(j-1)+1):(10*j),5,1)); 
      }
      V11[randvec,]=0;
      V22[randvec,]=0;
      V33[randvec,]=0;
      V44[randvec,]=0;
      
      A1=rbind(V11 ,matrix(0,P1-nrow(V11),4));
      A2=rbind(V22 ,matrix(0,P2-nrow(V22),4));
      A3=rbind(V33 ,matrix(0,P3-nrow(V33),4));
      A4=rbind(V44 ,matrix(0,P4-nrow(V44),4));      
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      A3[abs(A3)<=10^-8]=0;
      A4[abs(A4)<=10^-8]=0;
      
      TrueVariables1=rep(0,P1);
      TrueVariables1[abs(A1[,1])!=0]=1;
      
      TrueVariables2=rep(0,P2);
      TrueVariables2[abs(A2[,1])!=0]=1;
      
      TrueVariables3=rep(0,P3);
      TrueVariables3[abs(A3[,1])!=0]=1;
      
      TrueVariables4=rep(0,P4);
      TrueVariables4[abs(A4[,1])!=0]=1;
      
    } else if (setting==4){
      randvec=NULL;
      for (j in 1:3){
        randvec=c(randvec,sample((10*(j-1)+1):(10*j),5,1));       
      }
      
      A1=rbind(V11[1:30,] ,matrix(0,P1-30,4));
      A2=rbind(V22[1:30,] ,matrix(0,P2-30,4));
      A3=rbind(V33[1:30,] ,matrix(0,P3-30,4));
      A4=rbind(V44[1:30,] ,matrix(0,P4-30,4));
      
      V11=V11[1:30,];
      V22=V22[1:30,];
      V33=V33[1:30,];
      V44=V44[1:30,];
      
      V11[randvec,]=0;
      V22[randvec,]=0;
      V33[randvec,]=0;
      V44[randvec,]=0;
      
      A1=rbind(V11 ,matrix(0,P1-nrow(V11),4));
      A2=rbind(V22 ,matrix(0,P2-nrow(V22),4));
      A3=rbind(V33 ,matrix(0,P3-nrow(V33),4));
      A4=rbind(V44 ,matrix(0,P4-nrow(V44),4));
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      A3[abs(A3)<=10^-8]=0;
      A4[abs(A4)<=10^-8]=0;
      
      TrueVariables1=rep(0,P1);
      TrueVariables1[abs(A1[,1])!=0]=1;
      
      TrueVariables2=rep(0,P2);
      TrueVariables2[abs(A2[,1])!=0]=1;
      
      TrueVariables3=rep(0,P3);
      TrueVariables3[abs(A3[,1])!=0]=1;
      
      TrueVariables4=rep(0,P4);
      TrueVariables4[abs(A4[,1])!=0]=1;
      
      
    } else if (setting==5) {
      
      #%100 important variables but 25 each loaded on first to fourth factors
      
      A1=rbind(V11 ,matrix(0,P1-p1,4));
      A2=rbind(V22 ,matrix(0,P2-p2,4));
      A3=rbind(V33 ,matrix(0,P3-p3,4));
      A4=rbind(V44 ,matrix(0,P4-p4,4));
      
      whic=matrix(0,P1,4);
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
      
      A3[whic[,1]!=1,1]=0;
      A3[whic[,2]!=1,2]=0;
      A3[whic[,3]!=1,3]=0;
      A3[whic[,4]!=1,4]=0;
      
      A4[whic[,1]!=1,1]=0;
      A4[whic[,2]!=1,2]=0;
      A4[whic[,3]!=1,3]=0;
      A4[whic[,4]!=1,4]=0;
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      A3[abs(A3)<=10^-8]=0;
      A4[abs(A4)<=10^-8]=0;
      
      
      TrueVariables1=rep(0,P1);
      TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
      #TrueVariables1(sum(abs(A1),2)~=0)=1;
      
      TrueVariables2=rep(0,P2);
      TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
      
      TrueVariables3=rep(0,P3);
      TrueVariables3[apply(A3,1,function(x) sum(abs(x)))!=0]=1;
      #TrueVariables1(sum(abs(A1),2)~=0)=1;
      
      TrueVariables4=rep(0,P4);
      TrueVariables4[apply(A4,1,function(x) sum(abs(x)))!=0]=1;
    } 
  } else if (scenario==2){
    
    A1=rbind(V11 ,matrix(0,P1-p1,4));
    A2=rbind(V22 ,matrix(0,P2-p2,4)); 
    A3=rbind(V33 ,matrix(0,P3-p3,4));
    A4=rbind(V44 ,matrix(0,P4-p4,4)); 
    if (overlap=="yes"){
      # %overlap
      whic=matrix(0,P1,4);
      whic[1:100,1]=1;
      whic[1:100,2]=1;
      whic[1:100,3]=1;
      whic[1:100,4]=1;
      A1[whic[,1]!=1,1]=0;
      A1[whic[,2]!=1,2]=0;
      A1[,3:4]=0;
      
      A2[whic[,1]!=1,1]=0;
      A2[whic[,2]!=1,2]=0;
      A2[,3:4]=0;
      
      A3[,1:2]=0;
      A3[whic[,3]!=1,3]=0;
      A3[whic[,4]!=1,4]=0;
      
      A4[,1:2]=0;
      A4[whic[,3]!=1,3]=0;
      A4[whic[,4]!=1,4]=0;
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      A3[abs(A3)<=10^-8]=0;
      A4[abs(A4)<=10^-8]=0;
      
      TrueVariables1=rep(0,P1);
      TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
      TrueVariables2=rep(0,P2);
      TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
      TrueVariables3=rep(0,P3);
      TrueVariables3[apply(A3,1,function(x) sum(abs(x)))!=0]=1;
      TrueVariables4=rep(0,P4);
      TrueVariables4[apply(A4,1,function(x) sum(abs(x)))!=0]=1;
    } else if (overlap=="no"){
      #%No overlap    
      whic=matrix(0,P1,4);
      whic[1:50,1]=1;
      whic[51:100,2]=1;
      whic[1:50,3]=1;
      whic[51:100,4]=1;
      
      A1[whic[,1]!=1,1]=0;
      A1[whic[,2]!=1,2]=0;
      A1[,3:4]=0;
      
      A2[whic[,1]!=1,1]=0;
      A2[whic[,2]!=1,2]=0;
      A2[,3:4]=0;
      
      A3[,1:2]=0;
      A3[whic[,3]!=1,3]=0;
      A3[whic[,4]!=1,4]=0;
      
      A4[,1:2]=0;
      A4[whic[,3]!=1,3]=0;
      A4[whic[,4]!=1,4]=0;
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      A3[abs(A3)<=10^-8]=0;
      A4[abs(A4)<=10^-8]=0;
      TrueVariables1=rep(0,P1);
      TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
      TrueVariables2=rep(0,P2);
      TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
      TrueVariables3=rep(0,P3);
      TrueVariables3[apply(A3,1,function(x) sum(abs(x)))!=0]=1;
      TrueVariables4=rep(0,P4);
      TrueVariables4[apply(A4,1,function(x) sum(abs(x)))!=0]=1;
    }
    
  } else if (scenario==3){
    
    A1=rbind(V11 ,matrix(0,P1-p1,4));
    A2=rbind(V22 ,matrix(0,P2-p2,4));
    A3=rbind(V33 ,matrix(0,P3-p3,4));
    A4=rbind(V44 ,matrix(0,P4-p4,4));
    if (overlap=="no"){
      whic=matrix(0,P1,4);
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
      A2[whic[,3]!=1,3]=0;
      A2[,4]=0;
      
      A3[whic[,1]!=1,1]=0;
      A3[whic[,2]!=1,2]=0;
      A3[whic[,3]!=1,4]=0;
      A3[,3]=0;
      
      A4[whic[,1]!=1,1]=0;
      A4[whic[,2]!=1,2]=0;
      A4[whic[,3]!=1,4]=0;
      A4[,3]=0;
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      A3[abs(A3)<=10^-8]=0;
      A4[abs(A4)<=10^-8]=0;
    } else if (overlap=="yes"){
      whic=matrix(0,P1,4);
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
      A2[whic[,3]!=1,3]=0;
      A2[,4]=0;
      A3[whic[,1]!=1,1]=0;
      A3[whic[,2]!=1,2]=0;
      A3[whic[,3]!=1,4]=0;
      A3[,3]=0;
      A4[whic[,1]!=1,1]=0;
      A4[whic[,2]!=1,2]=0;
      A4[whic[,3]!=1,4]=0;
      A4[,3]=0;
      
      A1[abs(A1)<=10^-8]=0;
      A2[abs(A2)<=10^-8]=0;
      A3[abs(A3)<=10^-8]=0;
      A4[abs(A4)<=10^-8]=0;
      
    }
    
    TrueVariables1=rep(0,P1);
    TrueVariables1[apply(A1,1,function(x) sum(abs(x)))!=0]=1;
    TrueVariables2=rep(0,P2);
    TrueVariables2[apply(A2,1,function(x) sum(abs(x)))!=0]=1;
    TrueVariables3=rep(0,P3);
    TrueVariables3[apply(A3,1,function(x) sum(abs(x)))!=0]=1;
    TrueVariables4=rep(0,P4);
    TrueVariables4[apply(A4,1,function(x) sum(abs(x)))!=0]=1;
  }
  
  ### Generate data
  set.seed(seed);
  
  N_sites <- 20
  if (development) {
    n_families_per_site <- 20
  } else {
    n_families_per_site <- 100
  } 
  n_individs_per_family <- 2
  
  N_families <- N_sites*n_families_per_site
  N_obs <- N_sites*n_families_per_site*n_individs_per_family
  # Specify design matrix for sites
  Z_site <- kronecker(diag(N_sites), rep(1, n_families_per_site*n_individs_per_family))
  # Specify design matrix for families nested within sites
  Z_family <- kronecker(diag(N_families), rep(1, n_individs_per_family))
  
  mu <- 1
  sigma2_ksi <- sigma2_ksi_true # Site variance.
  sigma2_theta <- rep( sigma2_theta_true, N_sites) # Family:site variance (vec of len N_sites)
  # sigma2 <- 1 # Residual variance fixed
  #sigma2 <- .1 # Residual variance fixed
  
  sigma2u=1;
  U=MASS::mvrnorm(N_obs, mu=rep(0,r), sigma2u*diag(r))
  A5=c(1, 0, 1, 1) # components 1, 3, and 4 affect response
  E1=MASS::mvrnorm(N_obs, mu=rep(0,P1), Sigma1);
  E2=MASS::mvrnorm(N_obs, mu=rep(0,P2), Sigma2);
  E3=MASS::mvrnorm(N_obs, mu=rep(0,P3), Sigma3);
  E4=MASS::mvrnorm(N_obs, mu=rep(0,P4), Sigma4);
  E5=MASS::mvrnorm(N_obs, mu=0, sigma2u);
  
  X1 = U%*% t(A1) + E1
  X2 = U%*% t(A2) + E2
  X3 = U%*% t(A3) + E3
  X4 = U%*% t(A4) + E4
  Y = U%*% A5 + E5
  
  # if(sigma2_ksi_true!=0 & sigma2_theta_true!=0) {
    
    # Simulate random effects
    if (is.null(ksi)) {
      ksi <- rnorm(N_sites, mu, sd = sqrt(sigma2_ksi)) %>% matrix(ncol = 1)
    }
    
    # Sample theta_sf|ksi_s ~ N(ksi_s, sigma2_theta_s)
    theta <- matrix(0, nrow = N_sites, ncol = n_families_per_site)
    for (s in 1:N_sites) {
      theta[s, ] <- rnorm(n_families_per_site, mean = ksi[s], sd = sqrt(sigma2_theta[s]))
    }
    
    # Family effects as a single vector
    theta <- as.vector(t(theta)) %>% matrix(ncol = 1)
    
    # Combine effects
    Y <- Y + Z_family %*% theta
    
  # }
  
  list(X1=X1,X2=X2,X3=X3,X4=X4,Y=Y,Z_site=Z_site, Z_family=Z_family,
       TrueVar1=TrueVariables1,
       TrueVar2=TrueVariables2,
       TrueVar3=TrueVariables3,
       TrueVar4=TrueVariables4,
       U=U,LoadA1=A1,LoadA2=A2,LoadA3=A3,LoadA4=A4,
       Group1=Group1,Group2=Group2,Group3=Group3,Group4=Group4,
       ksi=ksi)
}

generate_train_test_sets <- function(scenario=1, setting=1, overlap=NULL, 
                                     train_seed=1, test_seed=2, 
                                     development=T,
                                     sigma2_ksi_true=0, 
                                     sigma2_theta_true=0) {
  
  train_set <- Simulate2(scenario, setting, 
                         seed = train_seed, 
                         sigma2_ksi_true = sigma2_ksi_true, 
                         sigma2_theta_true = sigma2_theta_true,
                         development = development)
  test_set <- Simulate2(scenario, setting, 
                        seed = test_seed, 
                        sigma2_ksi_true = sigma2_ksi_true, 
                        sigma2_theta_true = sigma2_theta_true,
                        development = development,
                        ksi = train_set$ksi)
  
  return(list(
    train_set=train_set,
    test_set=test_set
  ))
}

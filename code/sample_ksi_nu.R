# Load packages
library(tidyverse)
library(MASS)

# Simulate data
source("code/simulate_random_intercept_data.R")
n_obs <- 200
J <- 10
nu2 <- 3
simulation_results <- simulate_re_data(n_obs=n_obs, n_sites=J, nu2=nu2, seed=2)
# Residualize outcome
Y <- simulation_results$Y - 
  simulation_results$alpha_0 -
  simulation_results$U %*% simulation_results$alpha
# Number of observations
n <- nrow(Y)
# Generate response Y
ksi_true <- simulation_results$xi_s
Z <- simulation_results$Z # Design matrix
n_j <- colSums(Z) # N observations per site
# We transform our data
data <- data.frame(y = Y, j = apply(Z, 1, function(z) { which(z==1) })) %>% group_by(j)
data_sufficient <- data %>% summarize(ybar_j = mean(y), n_j = n())

# TODO Define Gibb's sampling functions
update_ksi <- function() {
  
}

update_mu <- function() {
  
}

update_sigma2 <- function() {
  
}

update_tau2 <- function() {
  
}

# Initialize data storage
param_names <- c(paste("ksi[", 1:J, "]", sep=""), "mu", "sigma2", "tau2")
chains <- 5 
iter <- 10000 
sims <- array(NA, c(iter, chains, J+3)) 
dimnames(sims) <- list(NULL, NULL, param_names) 

for (m in 1:chains) {
  
  ksi <- data %>% summarise(ystart = sample(y, 1)) %>% pull(ystart)
  mu <- data_sufficient$ybar_j %>% mean()
  # nu2 <- 
  # sigma2 <- 
  
  for (t in 1:iter) {
    
    for (j in 1:J)

  } 
}

# Analyze the samples
mcmc_summary <- monitor(sims)
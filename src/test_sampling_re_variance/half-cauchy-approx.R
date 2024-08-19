# Load packages
library(tidyverse)
library(MASS) # for mvnorm
library(rstan) # for monitor

# Simulate data
set.seed(42) # For reproducibility
# Number of observations
n <- 100
# Generate predictors
sigma2_beta_true <- 1
x1 <- rnorm(n, mean = 0, sd = sigma2_beta_true)
x2 <- rnorm(n, mean = 0, sd = sigma2_beta_true)
# Generate error term
sigma2_y <- 1
epsilon <- rnorm(n, mean = 0, sd = sigma2_y)
# Generate response Y
beta_true <- matrix(c(1,1,1), ncol = 1)
X <- model.matrix(~ x1 + x2) # Design matrix
Y <- X %*% beta_true + epsilon
n <- length(Y)

# Prior for beta
beta_prior_mean <- rep(0, ncol(X))
beta_prior_var <- 100

# Initial values
beta_current <- rep(0, ncol(X))
sigma2_current <- 1
n_samples <- 10000
n_burnin <- floor(n_samples)/ 5
a <- 0.1 # Lower bound of uniform prior for sigma^2
b <- 10  # Upper bound of uniform prior for sigma^2

# Storage for samples
beta_samples <- matrix(NA, n_samples, length(beta_current))
sigma2_samples <- numeric(n_samples)

# Set flag for sampling sigma2. Otherwise, it is fixed to sigma2_current
sample_sigma2 <- FALSE

# The following is from Bayesian Data Analysis 3 (BDA3) Appendix C:
schools <- read.csv("code/test_sampling_re_variance/schools.csv", header=TRUE)
J <- nrow(schools)
y <- schools$estimate
sigma <- schools$sd
# Gibbs sampling for the t model with fixed degrees of freedom: 
theta_update <- function(){ 
  theta_hat <- (mu/V + y/sigma^2)/(1/V + 1/sigma^2) 
  V_theta <- 1/(1/V + 1/sigma^2) 
  rnorm(J, theta_hat, sqrt(V_theta)) 
} 
mu_update <- function(){ 
  mu_hat <- sum(theta/V)/sum(1/V) 
  V_mu <- 1/sum(1/V) 
  rnorm(1, mu_hat, sqrt(V_mu)) 
} 
tau_update <- function(){ 
  sqrt(rgamma(1, J*nu/2+1, (nu/2)*sum(1/V))) 
} 
V_update <- function(){
  (nu*tau^2 + (theta-mu)^2)/rchisq(J,nu+1) 
} 
# Initially we fix the degrees of freedom at 4 to provide a robust analysis of the data. 
chains <- 5 
iter <- 1000 
sims <- array(NA, c(iter, chains, J+2)) 
dimnames(sims) <- list(NULL, NULL, c(paste("theta[", 1:8, "]", sep=""), "mu", "tau")) 
nu <- 4 
for (m in 1:chains){ 
  mu <- rnorm(1, mean(y), sd(y)) 
  tau <- runif(1, 0, sd(y)) 
  V<-runif(J,0,sd(y))^2 
  for (t in 1:iter){ 
    theta <- theta_update() 
    V<-V_update() 
    mu <- mu_update() 
    tau <- tau_update() 
    sims[t,m,] <- c(theta, mu, tau) 
  } 
} 
monitor(sims)

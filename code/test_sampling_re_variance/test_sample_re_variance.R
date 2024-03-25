# Load packages
library(tidyverse)
library(MASS)

# Simulate data
source("code/simulate_random_intercept_data.R")
n_obs <- 10000
n_sites <- 20
nu2 <- 3
simulation_results <- simulate_re_data(n_obs=n_obs, n_sites=n_sites, nu2=nu2, seed=2)
# Residualize outcome
Y <- simulation_results$Y - 
  simulation_results$alpha_0 -
  simulation_results$U %*% simulation_results$alpha
# Number of observations
n <- nrow(Y)
# Generate response Y
beta_true <- simulation_results$xi_s
X <- simulation_results$Z # Design matrix
colSums(X) # N observations per site

# We transform our data to load into code from BDA3 Appendix C:
# schools <- read.csv("code/test_sampling_re_variance/schools.csv", header=TRUE)
schools <- data.frame(y = Y, school = apply(X, 1, function(x) {which(x==1) })) %>%
  group_by(school) %>% summarize(estimate = mean(y), sd = sd(y))
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
param_names <- c(paste("theta[", 1:n_sites, "]", sep=""), "mu", "tau")
chains <- 5 
iter <- 10000 
sims <- array(NA, c(iter, chains, J+2)) 
dimnames(sims) <- list(NULL, NULL, param_names) 
nu <- 1 # prior t-distribution degrees of freedom
for (m in 1:chains){ 
  mu <- rnorm(1, mean(y), sd(y)) 
  tau <- runif(1, 0, sd(y)) 
  V <- runif(J,0,sd(y))^2 
  for (t in 1:iter){ 
    theta <- theta_update() 
    V <- V_update() 
    mu <- mu_update() 
    tau <- tau_update() 
    sims[t,m,] <- c(theta, mu, tau) 
  } 
}
# Analyze the samples
mcmc_summary <- monitor(sims)
# Posterior means
estimates <- mcmc_summary$mean
names(estimates) <- param_names
print(estimates)
# True parameter values 
truth <- c(beta_true, 0, 1)
names(truth) <- param_names
print(truth)
# Check if truth in 95% CIs
mcmc_summary$`2.5%` < truth & truth < mcmc_summary$`97.5%`

# Parameter-expanded Gibbs sampler:
gamma_update <- function(){ 
  gamma_hat <- (alpha*(y-mu)/sigma^2)/(1/tau^2 + alpha^2/sigma^2) 
  V_gamma <- 1/(1/tau^2 + alpha^2/sigma^2) 
  rnorm(J, gamma_hat, sqrt(V_gamma)) 
} 
alpha_update <- function(){ 
  alpha_hat <- sum(gamma*(y-mu)/sigma^2)/sum(gamma^2/sigma^2) 
  V_alpha <- 1/sum(gamma^2/sigma^2)
  rnorm(1, alpha_hat, sqrt(V_alpha)) 
} 
mu_update <- function(){ 
  mu_hat <- sum((y-alpha*gamma)/sigma^2)/sum(1/sigma^2) 
  V_mu <- 1/sum(1/sigma^2) 
  rnorm(1, mu_hat, sqrt(V_mu)) 
} 
tau_update <- function(){ 
  sqrt(sum(gamma^2)/rchisq(1,J-1)) 
}

sims <- array(NA, c(iter, chains, J+2)) 
dimnames(sims) <- list(NULL, NULL, c(paste("theta[", 1:n_sites, "]", sep=""), "mu", "tau")) 
for (m in 1:chains){ 
  alpha <- 1 
  mu <- rnorm(1, mean(y), sd(y)) 
  tau <- runif(1, 0, sd(y)) 
  for (t in 1:iter){ 
    gamma <- gamma_update() 
    alpha <- alpha_update() 
    mu <- mu_update() 
    tau <- tau_update() 
    sims[t,m,] <- c(mu + alpha*gamma, mu, abs(alpha)*tau) 
  } 
} 
# Analyze the samples
mcmc_summary <- monitor(sims)
# Posterior means
estimates <- mcmc_summary$mean
names(estimates) <- param_names
print(estimates)
# True parameter values 
truth <- c(beta_true, 0, 1)
names(truth) <- param_names
print(truth)
# Check if truth in 95% CIs
mcmc_summary$`2.5%` < truth & truth < mcmc_summary$`97.5%`

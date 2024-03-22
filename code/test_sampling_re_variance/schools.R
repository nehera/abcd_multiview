# The following is from Bayesian Data Analysis 3 (BDA3) Appendix C:
library(tidyverse)
setwd("code/test_sampling_re_variance")
schools <- read.csv("schools.csv", header=TRUE)
J <- nrow(schools)
y <- schools$estimate
sigma <- schools$sd

library("rstan")
schools_fit <- stan(file="schools.stan",
                    data=c("J","y","sigma"), iter=1000, chains=4)
print(schools_fit)
plot(schools_fit)

# The result is a list with four elements corresponding to the five quantities saved in the
# model: theta, eta, mu, tau, lp__
schools_sim <- extract(schools_fit)

# For example, we can display posterior inference for τ:
hist(schools_sim$tau)
# Or compute the posterior probability that the effect is larger in school A than in school C:
mean(schools_sim$theta[,1] > schools_sim$theta[,3])
# For example, we can simulate posterior predictive replicated data in the original 8 schools:
n_sims <- length(schools_sim$lp__)
y_rep <- array(NA, c(n_sims, J))
for (s in 1:n_sims)
  y_rep[s,] <- rnorm(J, schools_sim$theta[s,], sigma)

# The model as programmed above has nearly uniform prior distributions on the hyperparameters µθ and σθ. 
# An alternative is a half-Cauchy for σθ, which we could implement by
# taking the Stan model on page 592 and adding the line, tau ~ cauchy(0,25);.

# Using the t model
# It is straightforward to expand the hierarchical normal distribution for the coaching effects
# to a t distribution as discussed in Section 17.4, by replacing eta ~ normal(0,1); with
# eta ~ student_t(nu,0,1); and declaring nu as a parameter that takes on a value of 1 or
# greater (real<lower=1> nu;) and assigning it a prior distribution.

# C.3 Direct simulation, Gibbs, and Metropolis in R

# Direct/ Discrete approximation: 
mu_hat <- function(tau, y, sigma){ 
  sum(y/(sigma^2 + tau^2))/sum(1/(sigma^2 + tau^2)) 
} 

V_mu <- function(tau, y, sigma) { 
  1/sum(1/(tau^2 + sigma^2)) 
}

n_grid <- 2000 
tau_grid <- seq(.01, 40, length=n_grid) 
log_p_tau <- rep(NA, n_grid) 
for (i in 1:n_grid) { 
  mu <- mu_hat(tau_grid[i], y, sigma) 
  V <- V_mu(tau_grid[i], y, sigma) 
  log_p_tau[i] <- .5*log(V) - .5*sum(log(sigma^2 + tau_grid[i]^2)) -
    .5*sum((y-mu)^2/(sigma^2 + tau_grid[i]^2)) 
}

log_p_tau <- log_p_tau - max(log_p_tau) 
p_tau <- exp(log_p_tau) 
p_tau <- p_tau/sum(p_tau) 
n_sims <- 1000 
# Draw the simulations of τ from the approximate discrete distribution
tau <- sample(tau_grid, n_sims, replace=TRUE, prob=p_tau)

# The remaining steps are sampling from normal conditional distributions for μ and the θj’s:
mu <- rep(NA, n_sims) 
theta <- array(NA, c(n_sims,J)) 
for (i in 1:n_sims){ 
  mu[i] <- rnorm(1, mu_hat(tau[i],y,sigma), sqrt(V_mu(tau[i],y,sigma))) 
  theta_mean <- (mu[i]/tau[i]^2 + y/sigma^2)/(1/tau[i]^2 + 1/sigma^2) 
  theta_sd <- sqrt(1/(1/tau[i]^2 + 1/sigma^2)) 
  theta[i,] <- rnorm(J, theta_mean, theta_sd) 
}

## -- Gibbs sampler for the normal model assuming that variance sigma_j's known:
theta_update <- function() { 
  theta_hat <- (mu/tau^2 + y/sigma^2)/(1/tau^2 + 1/sigma^2) 
  V_theta <- 1/(1/tau^2 + 1/sigma^2) 
  rnorm(J, theta_hat, sqrt(V_theta)) 
  } 
mu_update <- function() { 
  rnorm(1, mean(theta), tau/sqrt(J)) 
  } 
tau_update <- function() { 
  sqrt(sum((theta-mu)^2)/rchisq(1,J-1)) 
}

chains <- 5 
iter <- 1000 
sims <- array(NA, c(iter, chains, J+2)) 
dimnames(sims) <- list(NULL, NULL, c(paste("theta[", 1:8, "]", sep=""), "mu", "tau")) 
for (m in 1:chains){ 
  mu <- rnorm(1, mean(y), sd(y)) 
  tau <- runif(1, 0, sd(y)) 
  for (t in 1:iter){ 
    theta <- theta_update() 
    mu <- mu_update() 
    tau <- tau_update() 
    sims[t,m,] <- c(theta, mu, tau) 
  } 
}

monitor(sims) # from library("rstan")

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
dimnames(sims) <- list(NULL, NULL, c(paste("theta[", 1:8, "]", sep=""), "mu", "tau")) 
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
monitor(sims)

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

# Sampling with t model with unknown degrees of freedom: 
log_post <- function(theta, V, mu, tau, nu, y, sigma){ 
  sum(dnorm(y, theta, sigma, log=TRUE)) + 
    sum(dnorm(theta, mu, sqrt(V), log=TRUE)) + 
    sum(.5*nu*log(nu/2) + nu*log(tau) -
        lgamma(nu/2) - (nu/2+1)*log(V) - .5*nu*tau^2/V) 
}

nu_update <- function(sigma_jump_nu){ 
  nu_inv_star <- rnorm(1, 1/nu, sigma_jump_nu) 
  if (nu_inv_star<=0 | nu_inv_star>1)
    p_jump <- 0 else { 
      nu_star <- 1/nu_inv_star 
      log_post_old <- log_post(theta, V, mu, tau, nu, y, sigma) 
      log_post_star <- log_post(theta, V, mu, tau, nu_star,y,sigma) 
      r<-exp(log_post_star-log_post_old) 
      nu <- ifelse(runif(1) < r, nu_star, nu) 
      p_jump <- min(r,1) } 
  return(list(nu=nu, p_jump=p_jump))
}

set.seed(1)
sigma_jump_nu <- 1 # Hyperparameter that defines the average jumping probability 
p_jump_nu <- array(NA, c(iter, chains)) 
sims <- array(NA, c(iter, chains, J+3)) 
dimnames(sims) <- list(NULL, NULL, c(paste("theta[", 1:8, "]", sep=""), "mu", "tau", "nu")) 
for (m in 1:chains){ 
  mu <- rnorm(1, mean(y), sd(y)) 
  tau <- runif(1, 0, sd(y)) 
  V<-runif(J,0,sd(y))^2 
  nu <- 1/runif(1, 0, 1) 
  for (t in 1:iter){ 
    theta <- theta_update() 
    V<-V_update() 
    mu <- mu_update() 
    tau <- tau_update() 
    temp <- nu_update(sigma_jump_nu) 
    nu <- temp$nu 
    p_jump_nu[t,m] <- temp$p_jump 
    sims[t,m,] <- c(theta, mu, tau, nu) 
    } 
  } 
print(mean(p_jump_nu)) 
monitor(sims)

# Finally, we can make the computations for the t model more efficient by applying parameter expansion.
gamma_update <- function(){ 
  gamma_hat <- (alpha*(y-mu)/sigma^2)/(1/V + alpha^2/sigma^2) 
  V_gamma <- 1/(1/V + alpha^2/sigma^2) 
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
  sqrt(rgamma(1, J*nu/2+1,(nu/2)*sum(1/V))) 
} 
V_update <- function(){ 
  (nu*tau^2 + gamma^2)/rchisq(J,nu+1) 
}

nu_update <- function(sigma_jump){ 
  nu_inv_star <- rnorm(1, 1/nu, sigma_jump) 
  if (nu_inv_star<=0 | nu_inv_star>1) 
    p_jump <- 0 
  else { 
    nu_star <- 1/nu_inv_star 
    log_post_old <- log_post(mu+alpha*gamma, alpha^2*V, mu, abs(alpha)*tau, nu, y, sigma) 
    log_post_star <- log_post(mu+alpha*gamma, alpha^2*V, mu, abs(alpha)*tau, nu_star, y, sigma) 
    r<-exp(log_post_star-log_post_old) 
    nu <- ifelse(runif(1) < r, nu_star, nu) 
    p_jump <- min(r,1) 
  } 
  return(list(nu=nu, p_jump=p_jump))
}

set.seed(1)
sigma_jump_nu <- 1 # Hyperparameter that defines the average jumping probability 
p_jump_nu <- array(NA, c(iter, chains)) 
sims <- array(NA, c(iter, chains, J+3)) 
dimnames(sims) <- list(NULL, NULL, c(paste("theta[", 1:8, "]", sep=""), "mu", "tau", "nu")) 
for (m in 1:chains){ 
  mu <- rnorm(1, mean(y), sd(y)) 
  tau <- runif(1, 0, sd(y)) 
  V<-runif(J,0,sd(y))^2 
  nu <- 1/runif(1, 0, 1) 
  gamma <- rnorm(J, 0, 1) 
  alpha <- rnorm(1, 0, 1)
  for (t in 1:iter){ 
    theta <- theta_update() 
    V<-V_update() 
    mu <- mu_update() 
    tau <- tau_update() 
    temp <- nu_update(sigma_jump_nu) 
    nu <- temp$nu 
    p_jump_nu[t,m] <- temp$p_jump 
    gamma <- gamma_update() 
    alpha <- alpha_update()
    sims[t,m,] <- c(mu+alpha*gamma, mu, abs(alpha)*tau, nu)
    
  } 
} 
print(mean(p_jump_nu)) 
monitor(sims)

# In the Metropolis updating function nu update() for the t degrees of freedom in Section C.3, 
# the log posterior density can be saved so that it does not need to be calculated twice at each step
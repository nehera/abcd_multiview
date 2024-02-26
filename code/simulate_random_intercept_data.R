library(tidyverse)

set.seed(1)

alpha_0 <- rnorm(1)
r <- 4
n <- 2000
U <-  matrix(data = rnorm(n*r), nrow = n, ncol = r)

prob_component_importance <- 0.5
n_important_components <- floor(prob_component_importance*r)
gamma_y <- rep(0, r)
gamma_y[1:n_important_components] <- 1
eta_y <- gamma_y # eta identical to gamma for the outcome view
p_m <- 1 # 1 feature in outcome view
A <- matrix(0, nrow = r, ncol = p_m) 

index_important_components <- seq(to = n_important_components)
n_nonzero_a <- n_important_components
nonzero_a <- matrix(rnorm(n_nonzero_a), nrow = n_important_components, ncol = p_m)
A[index_important_components, p_m] <- nonzero_a

n_covars <- 1 # Assume 1 clinical covariate for now
W <- matrix(rnorm(n), nrow = n, ncol = n_covars)
beta <- rnorm(n_covars) %>%
  matrix(ncol = 1) # Ensure column vector

epsilon <- rnorm(n) %>%
  matrix(ncol = 1) # Ensure column vector

nu2 <- 1
n_sites <- 22
xi_s <- rnorm(n_sites, sd = sqrt(nu2)) %>% matrix(ncol = 1) # make into column vector

# Sample random intercept design matrix
Z <- rmultinom(n, size = 1, prob = rep(1/ n_sites, n_sites)) %>% t() # convert into n x n_sites design matrix
n_per_site <- colSums(Z)

# Calculate outcome
y <- alpha_0 + Z%*%xi_s + U%*%A + W%*%beta + epsilon
y_tilde <- y - alpha_0 - U%*%A - W%*%beta # Essentially Z%*%xi_s + epsilon

est_conditional_params <- function(prior_mu, prior_var, post_var, x) {
  n <- nrow(x)
  conditional_var <- (1/ prior_var + n/ post_var)^(-1)
  conditional_mu <- conditional_var * (prior_mu/ prior_var + sum(x)/ post_var)
  return(c(conditional_mu, conditional_var))
}

# Get conditional distribution paramters for each site
conditional_params <- matrix(NA, nrow = n_sites, ncol = 2)
for (s in 1:n_sites) {
  site_index <- which(Z[, s]==1)
  x_site <- y_tilde[site_index, , drop = FALSE]
  conditional_params[s, ] <- est_conditional_params(prior_mu = 0, prior_var = 1, post_var = 1, x = x_site)
}

apply(conditional_params, 2, mean)

# TODO Visualize results


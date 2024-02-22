library(tidyverse)

set.seed(1)

alpha_0 <- rnorm(1)
r <- 4
n <- 200
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
beta <- rnorm(n_covars)

epsilon <- rnorm(n)

nu2_1 <- nu2_2 <- 1
n_sites <- 30
n_families_per_site <- 5
xi_s <- rnorm(n_sites, sd = sqrt(nu2_2))
xi_fs <- matrix(rnorm(n_sites*n_families_per_site, sd = sqrt(nu2_1)), nrow = n_sites, ncol = n_families_per_site)
for (s in 1:n_sites) {
  # Center family:site effects at corresponding site effects
  xi_fs[s, ] <- xi_fs[s, ] + xi_s[s]
}

# TODO draft random intercept design matrix simulation

# TODO calculate outcome

# TODO sample random intercepts from full conditionals

# TODO evaluate diff between sampled and true intercepts


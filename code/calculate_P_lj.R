library(tidyverse)

r <- 4
n <- 200
p_m <- 500
prob_component_selection <- 0.5
prob_feature_selection <- 0.05

initialize_gamma <- function(r=4, prior_component_selection=0.5) {
  gamma <- rbinom(n = r, size = 1, prob = prior_component_selection)
  return(gamma)
}
initialize_eta <- function(gamma, r=4, p_m=500, prior_feature_selection=0.05) {
  eta <- matrix(nrow = r, ncol = p_m)
  for (l in 1:r) {
    if (gamma[l] == 1) {
      eta[l, ] <- rbinom(n = p_m, size = 1, prior_feature_selection)
    } else {
      eta[l, ] <- rep(0, p_m)
    }
  }
  return(eta)
}

set.seed(1)
gamma <- initialize_gamma()
eta <- initialize_eta(gamma)

sigma2 <- rep(1, p_m)
tau2 <- matrix(1, nrow = r, ncol = p_m)
simulation_results <- readRDS("data/2023-07-03_simulation_results.rds")
U <- simulation_results$U # U is a fixed set of covariates

m <- 1 # Start with the 1st view
data_list <- readRDS("data/2023-07-03_simulation_data_list.rds")
x <- data_list[[m]]

calculate_mvnorm_var <- function(gamma, sigma2_j, tau2_j, U, n) {
  active_gamma <- which(gamma==1) 
  Sigma2_j <- U[, active_gamma, drop = FALSE] %*% 
    diag(tau2_j[active_gamma]) %*% 
    t(U[, active_gamma, drop = FALSE]) + diag(n) 
  mvnorm_var <- sigma2_j * Sigma2_j
  return(mvnorm_var)
}

calculate_log_dmvnorm_j <- function(gamma, sigma2_j, tau2_j, U, n, x_j) {
  mvnorm_var_j <- calculate_mvnorm_var(gamma, sigma2_j, tau2_j, U, n)
  log_dmvnorm_j <- mvtnorm::dmvnorm(x = x_j, mean = rep(0, n), sigma = mvnorm_var_j, log = TRUE)
  return(log_dmvnorm_j)
}

calculate_log_G_j <- function(gamma, eta_j, sigma2_j, tau2_j, U, n, x_j, prob_feature_selection) {
  log_dmvnorm_j <- calculate_log_dmvnorm_j(gamma, sigma2_j, tau2_j, U, n, x_j)
  active_gamma <- which(gamma == 1)
  n_active_gamma <- length(active_gamma) # TODO account for all inactive
  n_active_eta_given_gamma_1 <- eta_j[active_gamma] %>% sum()
  log_G_j <- log_dmvnorm_j + 
    n_active_eta_given_gamma_1 * log(prob_feature_selection) + 
    (n_active_gamma - n_active_eta_given_gamma_1) * log(1 - prob_feature_selection)
  return(log_G_j) 
}

# eta_j <- eta[, j]
# sigma2_j <- sigma2[j]
# tau2_j <- tau2[, j]
# x_j <- x[, j]

calculate_P_lj <- function(l, gamma, eta_j, sigma2_j, tau2_j, U, n, x_j, prob_feature_selection) {
  gamma_1 <- replace(gamma, l, 1)
  eta_1 <- replace(eta_j, l, 1)
  eta_0 <- replace(eta_j, l, 0)
  log_G_j_1 <- calculate_log_G_j(gamma_1, eta_1, sigma2_j, tau2_j, U, n, x_j, prob_feature_selection)
  log_G_j_0 <- calculate_log_G_j(gamma_1, eta_0, sigma2_j, tau2_j, U, n, x_j, prob_feature_selection)
  max_arg <- max(log_G_j_1, log_G_j_0) # Use logsumexp trick
  log_P_lj <- log_G_j_1 - (max_arg + log( exp(log_G_j_1 - max_arg) + exp(log_G_j_0 - max_arg) ) )
  P_lj <- exp(log_P_lj)
  # a <- log_G_j_1 - max_arg
  # b <- log_G_j_0 - max_arg
  # P_lj <- exp(a) / (exp(a)+exp(b))
  return(P_lj)
}

## -- Testing
# Component off and feature off across all components
l <- 1 
gamma[l]
j <- 1
eta[, j]
calculate_P_lj(l, gamma, eta[, j], sigma2[j], tau2[, j], U, n, x[, j], prob_feature_selection)

# Component off and feature on in another component
l <- 1 
gamma[l]
j <- 6
eta[, j]
calculate_P_lj(l, gamma, eta[, j], sigma2[j], tau2[, j], U, n, x[, j], prob_feature_selection)

# Component on and feature off across all components
l <- 3
gamma[l]
j <- 1
eta[, j]
calculate_P_lj(l, gamma, eta[, j], sigma2[j], tau2[, j], U, n, x[, j], prob_feature_selection)

# Component on and feature on in another component
l <- 3
gamma[l]
j <- 147
eta[, j]
calculate_P_lj(l, gamma, eta[, j], sigma2[j], tau2[, j], U, n, x[, j], prob_feature_selection)

## Note, only two components are on at set.seed(1) so the symmetry might cause an issue?
n_starting_points <- 30
Gamma <- matrix(nrow = r, ncol = n_starting_points)
Eta <- array(dim = c(r, p_m, n_starting_points))
for (s in 1:n_starting_points) {
  set.seed(s) 
  Gamma[,s] <- initialize_gamma()
  Eta[,,s] <- initialize_eta(Gamma[,s])
}
which(colSums(Gamma)==3)
s <- 14
gamma <- Gamma[,s]
eta <- Eta[,,s]
l <- 1
j <- 1
calculate_P_lj(l, gamma, eta[, j], sigma2[j], tau2[, j], U, n, x[, j], prob_feature_selection)

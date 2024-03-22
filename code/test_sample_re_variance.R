# Load packages
library(tidyverse)
library(MASS)

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
X <- model.matrix(~ x1 + x2, data=data) # Design matrix
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

# MCMC loop
for(i in 1:n_samples) {
  
  # Gibbs step for beta
  Sigma_beta_Y_inv <- solve(t(X) %*% X / sigma2_current + diag(1/beta_prior_var, ncol(X)))
  mu_beta_Y <- Sigma_beta_Y_inv %*% (t(X) %*% Y / sigma2_current)
  beta_current <- mvrnorm(1, mu = mu_beta_Y, Sigma = solve(Sigma_beta_Y_inv)) # Requires MASS package for mvrnorm
  
  if (sample_sigma2==TRUE) {
    # Metropolis step for theta (transformed sigma2)
    theta_current <- log(sigma2_current)
    theta_proposal <- rnorm(1, mean = theta_current, sd = 0.1) # Proposal in log-space
    sigma2_proposal <- exp(theta_proposal) # Transform back to original space
    # Calculate likelihood ratio, adjusting for Jacobian
    likelihood_ratio <- exp(-0.5 / (sigma2_proposal+sigma2_y) * sum((Y - X %*% beta_current)^2) + 
                              0.5 / (sigma2_current+sigma2_y) * sum((Y - X %*% beta_current)^2))
    # Jacobian adjustment
    jacobian_adjustment <- exp(theta_proposal - theta_current)
    # Calculate acceptance ratio, including Jacobian adjustment
    alpha <- min(1, likelihood_ratio * jacobian_adjustment)
    # Accept or reject
    if(runif(1) < alpha) {
      sigma2_current <- sigma2_proposal
    }
  }
  
  # Store samples
  beta_samples[i, ] <- beta_current
  sigma2_samples[i] <- sigma2_current
}

# After MCMC, analyze the samples

# Create a dataframe from the matrix for easier manipulation
beta_df <- as.data.frame(beta_samples)
names(beta_df) <- c("Beta_0", "Beta_1", "Beta_2")

# Add the iteration number as a column
beta_df$Iteration <- 1:nrow(beta_df)

# Convert from wide to long format for ggplot
beta_long <- pivot_longer(beta_df, cols = c("Beta_0", "Beta_1", "Beta_2"), names_to = "Parameter", values_to = "Value")

# Traceplot betas
ggplot(beta_long, aes(x = Iteration, y = Value, color = Parameter)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Traceplot of Beta Samples", x = "Iteration", y = "Beta Value") +
  theme(legend.title = element_blank())

# Summary statistics for betas
beta_long %>% 
  group_by(Parameter) %>% 
  filter(Iteration>n_burnin) %>% 
  summarise(mean(Value))

# Assuming sigma2_samples contains your MCMC output for sigma^2
sigma2_df <- data.frame(Iteration = 1:length(sigma2_samples),
                        Sigma2 = sigma2_samples,
                        Log_Sigma2 = log(sigma2_samples))

# Convert from wide to long format for ggplot
sigma2_long <- pivot_longer(sigma2_df, cols = c("Sigma2", "Log_Sigma2"), names_to = "Parameter", values_to = "Value")

# Plot for Sigma2
ggplot(sigma2_long[sigma2_long$Parameter == "Sigma2", ], aes(x = Iteration, y = Value)) +
  geom_line(color = "blue") +
  theme_minimal() +
  labs(title = "Traceplot of Sigma^2 Samples", x = "Iteration", y = "Sigma^2 Value")

# Plot for Log(Sigma2)
ggplot(sigma2_long[sigma2_long$Parameter == "Log_Sigma2", ], aes(x = Iteration, y = Value)) +
  geom_line(color = "red") +
  theme_minimal() +
  labs(title = "Traceplot of Log(Sigma^2) Samples", x = "Iteration", y = "Log(Sigma^2) Value")

# Summary statistics for sigma2
sigma2_long %>% 
  group_by(Parameter) %>% 
  filter(Iteration>n_burnin) %>% 
  summarise(mean(Value))

library(data.table)
library(tidyverse)

# Read the text files into data tables
eta_chain_data <- fread("code/Eta_chain.txt", header = FALSE)
gamma_chain_data <- fread("code/gamma_chain.txt", header = FALSE)

# Reshape the data into matrices or arrays
r <- 4  # Replace with the actual number of rows
p_m <- 10  # Replace with the actual number of columns
n_burnin <- 1000 # Replace with the actual number of burnin iterations
n_iterations <- 5000  # Replace with the actual number of iterations

# Reshape gamma_chain data into a matrix
gamma_chain <- matrix(gamma_chain_data$V1, nrow = r, ncol = n_iterations)

# Reshape Eta_chain data into a 3D array (cube)
eta_chain <- array(eta_chain_data$V1, dim = c(r, p_m, n_iterations))

## -- Analyze Component Selection
print("Component Selection Mean:")
gamma_chain[, (n_burnin+1):n_iterations] %>% apply(MARGIN = 1, FUN = mean)

print("Variable Selection Mean:")
eta_chain[,,(n_burnin+1):n_iterations] %>% apply(MARGIN = c(1,2), FUN = mean)

get_gamma_df <- function(gamma_chain) {
  n_iterations <- ncol(gamma_chain)
  gamma_df <- gamma_chain[, 1:n_iterations] %>% 
    apply(MARGIN = 1, FUN = cummean) %>% 
    as.data.frame() %>% mutate(iteration = 1:n_iterations) %>% 
    gather(key = "gamma", value = "MPP", -iteration) %>%
    mutate(gamma = gamma %>% as.factor()) 
  return(gamma_df)
}

gamma_df <- get_gamma_df(gamma_chain)

# TODO use latex in the title and labels e.g. latex2exp::TeX("$\\alpha$")
gamma_plot <- ggplot(gamma_df, 
                     aes(x = iteration, y = MPP, color = gamma)) + 
  geom_line() + geom_vline(xintercept = n_burnin, 
                           linetype = "dashed", color = "red") +
  labs(x = "iteration", y = "MPP", 
       title = "Trace plot for gamma")

print(gamma_plot)

## -- Analyze Feature Selection
features_of_interest <- c(1, 2, 6, 7)
feature_names <- paste("j", features_of_interest, sep = "_")
# TODO make trace plot function to reduce copy paste going forward
eta_df <- eta_chain[1, features_of_interest, 1:n_iterations] %>%
  apply(MARGIN = 1, FUN = cummean) %>% 
  as.data.frame() %>% rename_at(vars(names(.)), ~ feature_names) %>% mutate(iteration = 1:n_iterations) %>% 
  gather(key = "eta", value = "MPP", -iteration) %>%
  mutate(eta = eta %>% as.factor()) 

eta_1_plot <- ggplot(eta_df, 
                     aes(x = iteration, y = MPP, color = eta)) + 
  geom_line() + geom_vline(xintercept = n_burnin, 
                           linetype = "dashed", color = "red") +
  labs(x = "iteration", y = "MPP", 
       title = "Trace plot for eta_1.")

print(eta_1_plot)
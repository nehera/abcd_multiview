## First read in the arguments listed at the command line
args <- commandArgs(trailingOnly=TRUE) 
## args is now a list of character vectors 
## First check to see if arguments are passed. 
## Then cycle through each element of the list and evaluate the expressions. 
if(length(args) == 0){ 
  print("No arguments supplied.") 
  ## supply default values 
  slurm_id <- 1
} else{ 
  # set i to the first arg 
  slurm_id = args[1] 
} 
cat("Slurm ID:", slurm_id)

# Load libraries
if (("pacman" %in% installed.packages()[,"Package"]) == FALSE) { install.packages("pacman") }
pacman::p_load(tidyverse, reshape2, parallel, lme4) 
if (("BIPnet" %in% installed.packages()[,"Package"]) == FALSE) {
  pacman::p_load(devtools)
  devtools::install_github('chekouo/BIPnet')
}

# Define data analysis conditions
possible_r <- c(8, 10) # c(6, 8)
outcome_labels <- c("Internalizing Problems", "Externalizing Problems")
outcome_varnames <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r")
outcome_labels_and_varnames <- paste(outcome_labels, outcome_varnames, sep = "_and_")
normalize_response_flags <- c(FALSE, TRUE)
analysis_methods <- c("BIP") # , "BIPmixed")
analyses_df <- expand.grid(r = possible_r, outcome_label_and_varname = outcome_labels_and_varnames, 
                           normalize_response = normalize_response_flags,
                           analysis_method = analysis_methods)

# Load training data
train_list <- readRDS("data/2024-06-24_train_list.rds")

# Option to turn off subsetting
apply_subsetting <- FALSE

# Define analysis date
date_today <- Sys.Date()

# Sample ids from train_list$outcomes$src_subject_id
set.seed(123) # Setting seed for reproducibility
sampled_ids <- sample(train_list$outcomes$src_subject_id, 50)

# Function to subset data frames by sampled_ids
subset_data <- function(df, id_column) {
  df %>%
    filter(!!sym(id_column) %in% sampled_ids)
}

# Function to order data frames by a specified column
order_data <- function(df, id_column) {
  df %>%
    arrange(!!sym(id_column))
}

# Apply subsetting if needed
if (apply_subsetting) {
  train_list_subset <- train_list %>%
    map(~ subset_data(.x, "src_subject_id"))
} else {
  train_list_subset <- train_list
}

# Order each data frame in the (possibly subsetted) train_list
train_list_ordered <- train_list_subset %>%
  map(~ order_data(.x, "src_subject_id"))

# Check if all data frames are ordered by src_subject_id in the same way
check_order <- function(train_list) {
  ids <- map(train_list, ~ .x$src_subject_id)
  identical(ids[[1]], reduce(ids[-1], intersect))
}

# Ensure that all data frames are ordered by src_subject_id in the same way
all_ordered <- check_order(train_list_subset)

if (all_ordered) {
  print("All data frames are ordered by src_subject_id in the same way.")
  
  # Drop src_subject_id column from all data frames
  train_list_subset <- train_list_subset %>%
    map(~ select(.x, -src_subject_id))
  
  print("src_subject_id column has been dropped from all data frames.")
  
  # Drop design matrices matrices (views that start with Z_)
  train_list_subset <- train_list_subset[!grepl("^Z_", names(train_list_subset))]
  
  print(paste("We have dropped the design matrices."))
  
  # Identify columns with zero variance in each data frame
  zero_variance <- map_df(names(train_list_subset), ~ {
    df_name <- .x
    df <- train_list_subset[[df_name]]
    zero_var_cols <- df %>%
      select(where(~ var(.) == 0)) %>%
      colnames()
    
    tibble(
      data.type = df_name,
      variable.name = zero_var_cols
    )
  })
  
  if (nrow(zero_variance)>0) {
    # Print the zero_variance data frame to see which columns will be removed
    print("Features that have been removed for having zero variance:")
    print(zero_variance)
  } else {
    print("No features removed for having zero variance.")
  }
  
  # Remove columns with zero variance from the list of data frames
  train_list_subset <- map(train_list_subset, ~ {
    df <- .x
    zero_var_cols <- df %>%
      select(where(~ var(.) == 0)) %>%
      colnames()
    df %>% select(-all_of(zero_var_cols))
  })
  
  # Convert all data frames in train_list_subset to matrices
  train_list_matrices <- train_list_subset %>%
    map(~ as.matrix(.x))
  
} else {
  stop("Data frames are not ordered by src_subject_id in the same way.")
}

# Let's checkout the dimension of our matrices post-processing
lapply(train_list_matrices, dim)

# Function to check for NA, NaN, or Inf values in a matrix
check_na_nan_inf <- function(mat) {
  any(is.na(mat) | is.nan(mat) | is.infinite(mat))
}

# Apply the function to each matrix in train_list_matrices
check_results <- map_lgl(train_list_matrices, check_na_nan_inf)

# Print the NA check results
if (any(check_results)) {
  cat("The following matrices contain NA, NaN, or Inf values:\n")
  print(names(train_list_matrices)[check_results])
} else {
  cat("No matrices contain NA, NaN, or Inf values.\n")
}

# Do quick sensitivity analysis

# Convert each matrix to a data frame and concatenate them
df <- do.call(cbind, lapply(train_list_matrices, as.data.frame))

# Load necessary library for combining plots
library(patchwork)

# Standardizing the data with outcomes.cbcl_scr_syn_internal_r included
df_internal <- df %>%
  select(-outcomes.cbcl_scr_syn_external_r) %>%
  mutate(across(everything(), ~ scale(.) %>% as.vector()))

# Calculating the covariance matrix for internal
cov_matrix_internal <- cov(df_internal)

# Calculating the eigenvalues for internal
eigen_values_internal <- eigen(cov_matrix_internal)$values

# Prepare data for plotting internal scree plot
eigen_df_internal <- tibble(
  latent_component = 1:length(eigen_values_internal),
  eigenvalue = eigen_values_internal
)

# Plotting the scree plot for internal
scree_plot_internal <- eigen_df_internal %>%
  filter(latent_component %in% 1:20) %>%
  ggplot(aes(x = latent_component, y = eigenvalue)) +
  geom_point(shape = 22, size = 3, color = "blue", fill = "blue") +
  geom_line(color = "black", linewidth = 1) +
  labs(
    x = NULL,  # Remove x-axis label
    y = NULL,  # Remove y-axis label
    title = "Internalizing Behaviors Included"
  ) +
  theme_minimal() 

# Standardizing the data with outcomes.cbcl_scr_syn_external_r included
df_external <- df %>%
  select(-outcomes.cbcl_scr_syn_internal_r) %>%
  mutate(across(everything(), ~ scale(.) %>% as.vector()))

# Calculating the covariance matrix for external
cov_matrix_external <- cov(df_external)

# Calculating the eigenvalues for external
eigen_values_external <- eigen(cov_matrix_external)$values

# Prepare data for plotting external scree plot
eigen_df_external <- tibble(
  latent_component = 1:length(eigen_values_external),
  eigenvalue = eigen_values_external
)

# Plotting the scree plot for external
scree_plot_external <- eigen_df_external %>%
  filter(latent_component %in% 1:20) %>%
  ggplot(aes(x = latent_component, y = eigenvalue)) +
  geom_point(shape = 22, size = 3, color = "blue", fill = "blue") +
  geom_line(color = "black", linewidth = 1) +
  labs(
    x = "Number of Latent Components, r",  # Include x-axis label only here
    y = NULL,  # Remove y-axis label
    title = "Externalizing Behaviors Included"
  ) +
  theme_minimal() 

# Combine the plots
combined_plot <- scree_plot_internal / scree_plot_external +
  plot_layout(ncol = 1)

print(combined_plot)

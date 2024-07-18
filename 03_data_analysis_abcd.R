# Load necessary libraries
library(tidyverse)
library(BIPnet)

# Load training data
train_list <- readRDS("data/2024-06-24_train_list.rds")

# Sample 200 ids from train_list$outcomes$src_subject_id
set.seed(123) # Setting seed for reproducibility
sampled_ids <- sample(train_list$outcomes$src_subject_id, 200)

# Function to subset and order data frames by sampled_ids
subset_and_order <- function(df, id_column) {
  df %>%
    filter(!!sym(id_column) %in% sampled_ids) %>%
    arrange(!!sym(id_column))
}

# Subset and order each data frame in train_list
train_list_subset <- train_list %>%
  map(~ subset_and_order(.x, "src_subject_id"))

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
  
  # Convert all data frames in train_list_subset to matrices
  train_list_matrices <- train_list_subset %>%
    map(~ as.matrix(.x))
  
  # Drop the 2nd column of the outcomes matrix and rename it to y to start with internalizing problems
  train_list_matrices$y <- matrix(train_list_matrices$outcomes[, -2], ncol = 1)
  
  # Drop any matrices that start with Z_ and outcomes from the train_list_matrices object
  train_list_matrices <- train_list_matrices[!grepl("^Z_|^outcomes", names(train_list_matrices))]
  
} else {
  stop("Data frames are not ordered by src_subject_id in the same way.")
}

# Let's checkout the dimension of our matrices post-processing
lapply(train_list_matrices, dim)

# Run MCMC sampler
# Error in BIPnet::BIP(dataList = train_list_matrices, IndicVar = c(2, rep(0,  : 
# NA/NaN/Inf in foreign function call (arg 6)
# abcd_result <- BIPnet::BIP(dataList = train_list_matrices, 
#                            IndicVar = c(2, rep(0, 4), 1),
#                            groupList = NULL, Method = "BIP", 
#                            nbrcomp = 4, sample = 5000, burnin = 1000)

# Function to check for NA, NaN, or Inf values in a matrix
check_na_nan_inf <- function(mat) {
  any(is.na(mat) | is.nan(mat) | is.infinite(mat))
}

# Apply the function to each matrix in train_list_matrices
check_results <- map_lgl(train_list_matrices, check_na_nan_inf)

# Print the results
if (any(check_results)) {
  cat("The following matrices contain NA, NaN, or Inf values:\n")
  print(names(train_list_matrices)[check_results])
} else {
  cat("No matrices contain NA, NaN, or Inf values.\n")
}

# The issue is not on the data end
# On successfully running the function call below for only 3 views,
# and being able to sample with covariates included, 
# I conclude there's a bug in Thierry's code preventing running BIP on >3 views
abcd_result <- BIP(dataList = train_list_matrices[-(1:3)],
                           IndicVar = c(0,0,1),
                           groupList = NULL, Method = "BIP",
                           nbrcomp = 4, sample = 5000, burnin = 1000)

# TODO Should we normalize the outcome data as well? 
abcd_result$EstSig2[[3]]
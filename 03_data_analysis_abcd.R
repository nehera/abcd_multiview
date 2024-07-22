# Load libraries
library(tidyverse)

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
  
  # Drop the 2nd column of the outcomes matrix and make an element called y to start with internalizing problems
  train_list_subset$y <- train_list_subset$outcomes %>%
    select(cbcl_scr_syn_internal_r)
  
  # Drop design matrices matrices (views that start with Z_) and outcomes
  train_list_subset <- train_list_subset[!grepl("^Z_|^outcomes", names(train_list_subset))]
  
  print("We have subsetted the outcome view to internalizing problems and dropped design matrices.")
  
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
  
  # Print the zero_variance data frame to see which columns will be removed
  print("Features that have been removed for having zero variance.")
  print(zero_variance)
  
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

# Time taken in minutes is 1.66 for 200 observations, r=4, all views, BIP, 6K iterations
abcd_result <- BIPnet::BIP(dataList = train_list_matrices,
                           IndicVar = c(2,rep(0,4),1),
                           groupList = NULL, Method = "BIP",
                           nbrcomp = 4, sample = 5000, burnin = 1000)

# Let's checkout the BIP results

## -- Make a heatmap of U
EstU_matrix <- abcd_result$EstU

# Convert the matrix to a long format data frame
EstU_df <- as.data.frame(EstU_matrix) %>%
  rownames_to_column(var = "Individual") %>%
  pivot_longer(cols = -Individual, names_to = "LatentFactorComponent", values_to = "Value")

# Rename the Latent Factor Components to be numeric
EstU_df$LatentFactorComponent <- as.integer(gsub("V", "", EstU_df$LatentFactorComponent))

# Create the heatmap
ggplot(EstU_df, aes(x = LatentFactorComponent, y = Individual, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(
    title = "Heatmap of EstU",
    x = "Latent Factor Component",
    y = "Individual",
    fill = "Value"
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)
  )

## TODO Visualize VarSelMean, VarSelMeanGlobal, CompoSelMean

## TODO Visualize EstLoad

## TODO Make EstSig2, EstIntcp Table

## TODO Make & write-out data.frame of variables selected by thresholding VarSelMean & VarSelMeanGlobal

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

# Standardizing the outcomes and covariates
outcomes_standardized <- train_list_matrices$outcomes %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ scale(.) %>% as.vector()))

# Convert the marital status variable to a factor with labeled levels
covars_df <- train_list_matrices$covariates %>%
  as.data.frame %>%
  mutate(demo_prnt_marital_v2_bl = factor(
  demo_prnt_marital_v2_bl,
  levels = 1:6,
  labels = c("Married", "Widowed", "Divorced", "Separated", "Never_married", "Living_with_partner")
))

# One-hot encode the marital status variable
marital_status_encoded <- model.matrix(~ demo_prnt_marital_v2_bl - 1, data = covars_df)

# Convert the result back to a data frame
marital_status_encoded <- as.data.frame(marital_status_encoded)

# Inspect the one-hot encoded data
head(marital_status_encoded)

covars_standardized <- covars_df %>%
  select(-demo_prnt_marital_v2_bl) %>%
  cbind(marital_status_encoded) %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ scale(.) %>% as.vector()))

# Create site membership vector
Z_site <- train_list$Z_site[, -1] # Exclude subject id
site_membership <- apply(Z_site, 1, function(x) which(x == 1)) %>%
  as.factor()

# Create family membership vector
Z_family <- train_list$Z_family[, -1] # Exclude subject id
family_membership <- apply(Z_family, 1, function(x) which(x == 1)) %>%
  as.factor()

# Fit the model for internal outcome (1st column)
df_internal <- data.frame(
  y = outcomes_standardized[, 1],
  covars_standardized,
  site = site_membership,
  family = family_membership
)

model_internal <- lmer(y ~ . + (1 | site/family), data = df_internal)
# fixed-effect model matrix is rank deficient so dropping 21 columns / coefficients
# Warning messages:
#   1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                     unable to evaluate scaled gradient
#                   2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                     Model failed to converge: degenerate  Hessian with 2 negative eigenvalues
summary_internal <- summary(model_internal)
# Despite getting a convergence error/ issue in the fixed effects, on normalizing the outcome, I got the following var estimates:
# Random effects:
#   Groups      Name        Variance Std.Dev.
# family:site (Intercept) 0.75448  0.8686  
# site        (Intercept) 0.03206  0.1791  
# Residual                0.94222  0.9707 
print(summary_internal)

# Fit the model for external outcome (2nd column)
df_external <- data.frame(
  y = outcomes_standardized[, 2],
  covars_standardized,
  site = site_membership,
  family = family_membership
)

model_external <- lmer(y ~ . + (1 | site/family), data = df_external)
# I also get fixed-effect model matrix is rank deficient so dropping 21 columns / coefficients when fitting the model to external data, 
# I'm wondering which fixed effects are of concern. 
summary_external <- summary(model_external)

# Groups      Name        Variance  Std.Dev.
# family:site (Intercept) 17.791193 4.2180  
# site        (Intercept)  0.009101 0.0954  
# Residual                 0.911504 0.9547 

print(summary_external)


# Load necessary libraries
library(ggplot2)
library(car)

# Extract the fixed-effect design matrix from the model
X_fixed <- covars_standardized

# Calculate the correlation matrix of the design matrix
cor_matrix <- cor(X_fixed)

# Convert the correlation matrix to a format suitable for ggplot2
cor_data <- as.data.frame(as.table(cor_matrix))

# Plot the correlation matrix
cor_plot <- ggplot(cor_data, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", title = paste("Covariate Correlation Matrix"))

print(cor_plot)

# Function to extract fixed-effect design matrix, plot correlation matrix, and calculate VIF
analyze_model <- function(model, title_prefix) {
  # Calculate VIF for the fixed effects in the model
  vif_values <- vif(model)
  cat(paste(title_prefix, "VIF values:\n"))
  print(vif_values)
}

# Analyze model_internal
analyze_model(model_internal, "Internal Outcome")

# Analyze model_external
analyze_model(model_external, "External Outcome")

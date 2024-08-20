## First read in the arguments listed at the command line
args <- commandArgs(trailingOnly=TRUE) 
## args is now a list of character vectors 
## First check to see if arguments are passed. 
## Then cycle through each element of the list and evaluate the expressions. 
if(length(args) == 0){ 
  print("No arguments supplied.") 
  ## supply default values 
  slurm_id <- 9
  } else{ 
    # set i to the first arg 
    slurm_id = args[1] 
} 
cat("Slurm ID:", slurm_id)

# Load libraries
if (("pacman" %in% installed.packages()[,"Package"]) == FALSE) { install.packages("pacman") }
pacman::p_load(tidyverse, reshape2, parallel) 
# if (("BIPnet" %in% installed.packages()[,"Package"]) == FALSE) {
#   pacman::p_load(devtools)
#   devtools::install_github('chekouo/BIPnet')
# }

# Source methods
source("src/BIP.R") # Includes BIPmixed implementation
source("src/BIPpredict.R")

# Define data analysis conditions
possible_r <- c(5, 6)
outcome_labels <- c("Internalizing Problems", "Externalizing Problems")
outcome_varnames <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r")
outcome_labels_and_varnames <- paste(outcome_labels, outcome_varnames, sep = "_and_")
analysis_methods <- c("BIP", "BIPmixed")
ELA_contin_flags <- c(FALSE, TRUE) # Flag to indicate whether or not to artificially include a continuous var in the ELA view
analyses_df <- expand.grid(r = possible_r, outcome_label_and_varname = outcome_labels_and_varnames, 
                           analysis_method = analysis_methods, ELA_contin_flag = ELA_contin_flags) %>%
  filter(!(analysis_method == "BIP" & r == 5)) # There are 6 views when applying BIP, and it's important to not have n components < n views. 

# Do analysis_1 by default
analysis_conditions <- analyses_df[slurm_id, ]

# Extract analysis conditions
r <- analysis_conditions$r
analysis_conditions <- analysis_conditions %>%
  mutate(outcome_label_and_varname = as.character(outcome_label_and_varname)) %>%
  separate(outcome_label_and_varname, into = c("outcome_label", "varname"), sep = "_and_")
outcome_label <- analysis_conditions$outcome_label
outcome_varname <- analysis_conditions$varname
normalize_response <- analysis_conditions$normalize_response
analysis_method <- analysis_conditions$analysis_method
ELA_contin_flag <- analysis_conditions$ELA_contin_flag

# Load training data
train_list <- readRDS("data/2024-06-24_train_list.rds")

# Option to turn off subsetting
apply_subsetting <- FALSE
n_subjects <- 500

# Define analysis date
date_today <- Sys.Date()

# Sample ids from train_list$outcomes$src_subject_id
set.seed(123) # Setting seed for reproducibility
sampled_ids <- sample(train_list$outcomes$src_subject_id, n_subjects)

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
    map(~ dplyr::select(.x, -src_subject_id))
  
  print("src_subject_id column has been dropped from all data frames.")
  
  # Extract design matrices
  Z_family_train <- train_list_subset$Z_family
  Z_site_train <- train_list_subset$Z_site
  
  # Drop design matrices matrices (views that start with Z_)
  train_list_subset <- train_list_subset[!grepl("^Z_", names(train_list_subset))]
  
  print(paste("We have dropped the design matrices."))
  
  # Identify columns with zero variance in each data frame
  zero_variance <- map_df(names(train_list_subset), ~ {
    df_name <- .x
    df <- train_list_subset[[df_name]]
    zero_var_cols <- df %>%
      dplyr::select(where(~ var(.) == 0)) %>%
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
      dplyr::select(where(~ var(.) == 0)) %>%
      colnames()
    df %>% dplyr::select(-all_of(zero_var_cols))
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

# Focus on the outcome of interest & Drop the outcome that's not
outcome_index <- which( colnames(train_list_matrices$outcomes) == outcome_varname )
trainList <- train_list_matrices[!grepl("^outcomes", names(train_list_matrices))]
trainList$y <- train_list_matrices$outcomes[, outcome_index, drop = FALSE]

# Potentially add continuous variable to the ELA view
if (ELA_contin_flag) {
  trainList$ela_view <- cbind(trainList$ela_view, contin_var = rnorm(nrow(trainList$ela_view)))
  print("Continous variable added to the ELA view.")
}

# Fit BIP to ABCD Study Training Set
print("Fitting model...")

# Define the indicators where 2= covariates, 0= omics, and 1= outcome.
IndicVar <- c(2, rep(0,4), 1)
n_sample <- 5000
n_burnin <- 1000

if (analysis_method == "BIP") {
  
  # Fit models
  BIP_start_time <- Sys.time()
  model_fit <- BIP(dataList = trainList, IndicVar = IndicVar, Method = "BIP",
                    nbrcomp = r, sample = n_sample, burnin = n_burnin)
  BIP_end_time <- Sys.time()
  print("BIP required:")
  print(BIP_end_time - BIP_start_time)
  
} else if (analysis_method == "BIPmixed") {
  
  # Fit BIPmixed to the data
  BIPmixed_start_time <- Sys.time()
  model_fit <- BIP(dataList = trainList, IndicVar = IndicVar, Method = "BIPmixed",
                         nbrcomp = r, sample = n_sample, burnin = n_burnin,
                         Z_family = Z_family_train, Z_site = Z_site_train)
  BIPmixed_end_time <- Sys.time()
  print("BIPmixed required")
  print(BIPmixed_end_time - BIPmixed_start_time)
  
} else {
  stop("You must provide a valid analysis_method.")
}

# model_fit <- readRDS("models/2024-08-19_Internalizing_r_5_method_BIPmixed_model_fit.rds")

# Let's make the output directory name
# Collapse column names and values into a string
collapsed_string <- apply(dplyr::select(analysis_conditions, -outcome_label), 1, function(row) {
  paste(names(row), row, sep = "_", collapse = "-")
})

# Define the output directory string
output_dir <- file.path("data_analysis_results", collapsed_string)
# Directories for models, figures, and tables
models_dir <- file.path(output_dir, "models")
figures_dir <- file.path(output_dir, "figures")
tables_dir <- file.path(output_dir, "tables")

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
} else {
  stop("output_dir already exists.")
}

# Create the subdirectories
dir.create(models_dir)
dir.create(figures_dir)
dir.create(tables_dir)

# Define the prefix for any files generated during this analysis
file_prefix <- paste(date_today, str_split(outcome_label, " ")[[1]][1],
                               "r", r, "method", analysis_method, sep = "_")

# Store the model fit
saveRDS(model_fit, paste0(models_dir, "/", file_prefix, "_model_fit.rds"))

print("Making preliminary visualizations...")

# Create a named vector for the view labels, removing the "_view" suffix
view_labels <- c("covariates", "ELA", "sMRI_SA", "sMRI_CT", "fMRI", outcome_label)
names(view_labels) <- 1:length(view_labels)

## Visualize VarSelMean, VarSelMeanGlobal, CompoSelMean

# Extract the CompoSelMean matrix from the model_fit object
CompoSelMean <- model_fit$CompoSelMean

# Convert the matrix to a long format suitable for ggplot2
CompoSelMean_long <- melt(CompoSelMean)
colnames(CompoSelMean_long) <- c("View", "Component", "Probability")

# Create the heatmap with probabilities printed and custom axis labels
heatmap_CompoSelMean <- ggplot(CompoSelMean_long, aes(x = Component, y = View, fill = Probability)) +
  geom_tile() +
  geom_text(aes(label = round(Probability, 2)), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  scale_y_continuous(breaks = 1:length(view_labels), labels = view_labels) +
  theme_minimal() +
  labs(title = "Component Selection Probabilities",
       x = "Component",
       y = "View",
       fill = "Probability")

# Print the heatmap
print(heatmap_CompoSelMean)

# Assume libraries have been loaded and view_labels has been constructed upstream

# Extract the VarSelMean and VarSelMeanGlobal lists from the model_fit object
VarSelMean <- model_fit$VarSelMean
VarSelMeanGlobal <- model_fit$VarSelMeanGlobal

# Convert the VarSelMean list to a data frame with appropriate labels
VarSelMean_df <- bind_rows(
  lapply(1:length(VarSelMean), function(i) {
    data.frame(View = i, VarSelMean[[i]])
  })
)

# Rename the columns to indicate component ids
colnames(VarSelMean_df) <- c("View", paste("Component", 1:r, sep = "_"))

# Convert the VarSelMeanGlobal list to a data frame with a new column for Component
VarSelMeanGlobal_df <- bind_rows(
  lapply(1:length(VarSelMeanGlobal), function(i) {
    data.frame(View = i, Component = "Global", Probability = VarSelMeanGlobal[[i]])
  })
)

# Convert VarSelMean_df to long format for ggplot2
VarSelMean_long <- VarSelMean_df %>%
  pivot_longer(cols = starts_with("Component"), names_to = "Component", values_to = "Probability") %>%
  mutate(Component = gsub("Component_", "", Component))

# Combine VarSelMean_long and VarSelMeanGlobal_df into a single data frame
combined_df <- bind_rows(VarSelMean_long, VarSelMeanGlobal_df)

# Create the grouped boxplots by view using ggplot2
boxplot_VarSelMean <- ggplot(combined_df, aes(x = factor(View), y = Probability, fill = factor(Component))) +
  geom_boxplot() +
  scale_x_discrete(labels = view_labels) +
  scale_fill_discrete(name = "Component") +
  theme_minimal() +
  labs(title = "Variable Selection Probabilities by View",
       x = "View",
       y = "Probability")

# Print the boxplots
print(boxplot_VarSelMean)

## Make EstSig2 Visualization

# Assume libraries have been loaded and view_labels has been constructed upstream

# Extract the EstSig2 list from the model_fit object
EstSig2 <- model_fit$EstSig2

# Convert the EstSig2 list to a data frame with appropriate labels
EstSig2_df <- bind_rows(
  lapply(1:length(EstSig2), function(i) {
    data.frame(View = i, Value = EstSig2[[i]])
  })
)

# Separate the outcome view value
outcome_view_value <- EstSig2_df %>% filter(View == length(view_labels)) %>% pull(Value)

# Filter out the outcome view from the EstSig2_df
EstSig2_df <- EstSig2_df %>% filter(View != length(view_labels))

# Create the grouped boxplots by view using ggplot2
boxplot_EstSig2 <- ggplot(EstSig2_df, aes(x = factor(View), y = Value, fill = factor(View))) +
  geom_boxplot() +
  scale_x_discrete(labels = view_labels[-length(view_labels)]) +
  scale_fill_discrete(name = "View") +
  theme_minimal() +
  labs(title = bquote(hat(sigma)^2 ~ by ~ View ~ ";" ~ hat(sigma)^2(0) ~ "=" ~ .(round(outcome_view_value, 2))),
       x = "View",
       y = "Value")

# Print the boxplots
print(boxplot_EstSig2)

## Write out visualizations

# Create a list of ggplot objects and their respective names
plots <- list(
  heatmap_CompoSelMean = heatmap_CompoSelMean,
  boxplot_VarSelMean = boxplot_VarSelMean,
  boxplot_EstSig2 = boxplot_EstSig2
)

# Loop through each plot and save it
for (plot_name in names(plots)) {
  plot <- plots[[plot_name]]
  filename <- paste0(figures_dir, "/", file_prefix, "_", plot_name, ".png")
  ggsave(filename, plot = plot, width = 8, height = 6, units = "in")
}

## Make & write-out data.frame of variables selected by thresholding VarSelMean & VarSelMeanGlobal
print("Making lists of variables selected...")

# Set the variable selction threshold
mpp_threshold <- 0 # We can filter in post-processing of outputs

# Extract the VarSelMean and VarSelMeanGlobal lists from the model_fit object
VarSelMean <- model_fit$VarSelMean
VarSelMeanGlobal <- model_fit$VarSelMeanGlobal

# Convert the VarSelMean list to a data frame with appropriate labels
VarSelMean_df <- bind_rows(
  lapply(1:length(VarSelMean), function(i) {
    data.frame(View = view_labels[i], VarSelMean[[i]])
  })
)

# Rename the columns to indicate component ids
colnames(VarSelMean_df) <- c("View", paste("Component", 1:r, sep = "_"))

# Extract column names from each matrix
all_colnames <- unlist(lapply(trainList[!(names(trainList) %in% c("Z_family", "Z_site"))], colnames))

# Add feature label
VarSelMean_df$Feature <- all_colnames

# Convert VarSelMean_df to long format for ggplot2
VarSelMean_by_component <- VarSelMean_df %>%
  pivot_longer(cols = starts_with("Component"), names_to = "Component", values_to = "Probability") %>%
  mutate(Component = gsub("Component_", "", Component))

# Add the global probabilities
VarSelMean_globally <- bind_rows(
  lapply(1:length(VarSelMeanGlobal), function(i) {
    data.frame(View = view_labels[i], Probability = VarSelMeanGlobal[[i]])
  })
)

# Add feature label
VarSelMean_globally$Feature <- all_colnames

VarSelMean_list <- list("VarSelMean_by_component" = VarSelMean_by_component,
                        "VarSelMean_globally" = VarSelMean_globally)

# Filter variables based on the threshold
variables_selected <- lapply(VarSelMean_list, function(x) {
  x %>% filter(Probability >= mpp_threshold) %>%
    arrange(desc(Probability))
})

# Loop through the list and write each element to a separate CSV file
for (table_name in names(variables_selected)) {
  # Generate filename
  filename <- paste0(tables_dir, "/", file_prefix, "_", table_name, ".csv")
  # Write the data frame to a CSV file
  write.csv(variables_selected[[table_name]], filename, row.names = FALSE)
}

# Now let's make predictions. We evaluate and aggregate their error later on. 

# Load test data
test_list <- readRDS("data/2024-06-24_test_list.rds")

# Sample ids from test_list$outcomes$src_subject_id
set.seed(123) # Setting seed for reproducibility
sampled_ids <- sample(test_list$outcomes$src_subject_id, n_subjects)

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
  test_list_subset <- test_list %>%
    map(~ subset_data(.x, "src_subject_id"))
} else {
  test_list_subset <- test_list
}

# Order each data frame in the (possibly subsetted) test_list
test_list_ordered <- test_list_subset %>%
  map(~ order_data(.x, "src_subject_id"))

# Check if all data frames are ordered by src_subject_id in the same way
check_order <- function(test_list) {
  ids <- map(test_list, ~ .x$src_subject_id)
  identical(ids[[1]], reduce(ids[-1], intersect))
}

# Ensure that all data frames are ordered by src_subject_id in the same way
all_ordered <- check_order(test_list_subset)

if (all_ordered) {
  print("All data frames are ordered by src_subject_id in the same way.")
  
  # Drop src_subject_id column from all data frames
  test_list_subset <- test_list_subset %>%
    map(~ dplyr::select(.x, -src_subject_id))
  
  print("src_subject_id column has been dropped from all data frames.")
  
  # Extract design matrices
  Z_family_test <- test_list_subset$Z_family
  Z_site_test <- test_list_subset$Z_site
  
  # Drop design matrices matrices (views that start with Z_)
  test_list_subset <- test_list_subset[!grepl("^Z_", names(test_list_subset))]
  
  print(paste("We have dropped the design matrices."))
  
  # Identify columns with zero variance in each data frame
  zero_variance <- map_df(names(test_list_subset), ~ {
    df_name <- .x
    df <- test_list_subset[[df_name]]
    zero_var_cols <- df %>%
      dplyr::select(where(~ var(.) == 0)) %>%
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
  test_list_subset <- map(test_list_subset, ~ {
    df <- .x
    zero_var_cols <- df %>%
      dplyr::select(where(~ var(.) == 0)) %>%
      colnames()
    df %>% dplyr::select(-all_of(zero_var_cols))
  })
  
  # Convert all data frames in test_list_subset to matrices
  test_list_matrices <- test_list_subset %>%
    map(~ as.matrix(.x))
  
} else {
  stop("Data frames are not ordered by src_subject_id in the same way.")
}

# Reshape test data
testList <- test_list_matrices[!grepl("^outcomes", names(test_list_matrices))] # Exclude outcomes

# Get predictions
if (analysis_method == "BIP") {
  y_preds <- BIPpredict(dataListNew = testList, Result = model_fit, meth = "BMA")$ypredict
} else if (analysis_method == "BIPmixed") {
  y_preds <- BIPpredict(dataListNew = testList, Result = model_fit, meth = "BMA", 
                                 Z_site = Z_site_test, Z_family = Z_family_test)$ypredict
} else {
  stop("You must provide a valid analysis_method.")
}

Y_true <- pull(test_list$outcomes, outcome_varname)

# Create a data frame with true and predicted values
results_df <- data.frame(Y_true = Y_true, y_preds = y_preds)

# Calculate MSPE
mspe <- mean((Y_true - y_preds)^2)

# Add MSPE as a new column to the data frame
results_df$MSPE <- mspe

# Write the data frame to a CSV file
output_csv_path <- file.path(tables_dir, "y_true_vs_y_preds_mspe.csv")
write.csv(results_df, output_csv_path, row.names = FALSE)

# Create the scatter plot
p <- ggplot(results_df, aes(x = Y_true, y = y_preds)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + # Reference line for perfect prediction
  labs(title = "True Y vs Predicted Y", 
       x = "True Y", 
       y = "Predicted Y") +
  annotate("text", x = Inf, y = Inf, label = paste("MSPE =", round(mspe, 4)), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black") +
  theme_minimal()

# Save the plot
output_plot_path <- file.path(figures_dir, "true_y_vs_pred_y_mspe_plot.png")
ggsave(output_plot_path, plot = p, width = 8, height = 6)

# Store random effect inference, if applicable
if (analysis_method == "BIPmixed") {
  # Calculate posterior mean and credible intervals for sigma2_ksi_samples
  sigma2_ksi_samples_post_burnin <- model_fit$sigma2_ksi_samples[(n_burnin + 1):nrow(model_fit$sigma2_ksi_samples), , drop = FALSE]
  
  # Posterior mean
  posterior_mean_sigma2_ksi <- colMeans(sigma2_ksi_samples_post_burnin)
  
  # 2.5% and 97.5% credible intervals (assuming you want a 95% CI)
  credible_interval_sigma2_ksi <- apply(sigma2_ksi_samples_post_burnin, 2, quantile, probs = c(0.025, 0.975))
  
  # Calculate posterior mean and credible intervals for sigma2_theta_samples
  sigma2_theta_samples_post_burnin <- model_fit$sigma2_theta_samples[(n_burnin + 1):nrow(model_fit$sigma2_theta_samples), , drop = FALSE]
  
  # Posterior mean
  posterior_mean_sigma2_theta <- colMeans(sigma2_theta_samples_post_burnin)
  
  # 2.5% and 97.5% credible intervals (95% CI)
  credible_interval_sigma2_theta <- apply(sigma2_theta_samples_post_burnin, 2, quantile, probs = c(0.025, 0.975))
  
  # Combine the results into a data frame for easy viewing
  results_sigma2_ksi <- data.frame(
    Posterior_Mean = posterior_mean_sigma2_ksi,
    CI_Lower = credible_interval_sigma2_ksi[1, ],
    CI_Upper = credible_interval_sigma2_ksi[2, ]
  )
  
  results_sigma2_theta <- data.frame(
    Posterior_Mean = posterior_mean_sigma2_theta,
    CI_Lower = credible_interval_sigma2_theta[1, ],
    CI_Upper = credible_interval_sigma2_theta[2, ]
  )
  
  # Save these results to CSV files
  write.csv(results_sigma2_ksi, file.path(tables_dir, "sigma2_ksi_posterior_summary.csv"), row.names = FALSE)
  write.csv(results_sigma2_theta, file.path(tables_dir, "sigma2_theta_posterior_summary.csv"), row.names = FALSE)
}

print("Data analysis performed & results stored.")
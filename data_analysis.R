# User Options
transform_outcomes <- FALSE # sqrt transform outcomes
n_sample <- 10000
n_burnin <- 5000

# Ensure we are in the correct working directory
setwd("/users/4/neher015/abcd_multiview")

# Set the library path using the full file path
.libPaths("/users/4/neher015/R/x86_64-pc-linux-gnu-library/4.3")

# Ensure library path has been set
print(.libPaths())

# Load libraries
if (!("pacman" %in% installed.packages()[,"Package"])) {
  install.packages("pacman", lib = "/users/4/neher015/R/x86_64-pc-linux-gnu-library/4.3")
}

pacman::p_load(tidyverse, reshape2, parallel) 

# Source methods
source("src/utility_functions.R") # For data splitting
source("src/BIP.R") # Includes BIPmixed implementation
source("src/BIPpredict.R")

# First read in the arguments listed at the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 1) {
  slurm_array_task_id <- as.numeric(args[1])
} else {
  print("commandArgs have not been supplied.")
  # Supply default values
  slurm_array_task_id <- 10
}

output_dir_name <- file.path("data_analysis_results",
                             paste0(Sys.Date(), "_data_analysis"))

# Check if the directory exists, and if not, create it
if (!dir.exists(output_dir_name)) {
  dir.create(output_dir_name)
}

print("Slurm Array Task ID:")
print(slurm_array_task_id)
print("Main Output Directory:")
print(output_dir_name)

# Define data analysis conditions
possible_r <- c(6)
outcome_labels <- c("Internalizing Problems (R)", "Externalizing Problems (R)",
                    "Internalizing Problems (T)", "Externalizing Problems (T)",
                    "Body Mass Index (BMI)")
outcome_varnames <- c("cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r",
                      "cbcl_scr_syn_internal_t", "cbcl_scr_syn_external_t",
                      "BMI")
outcome_labels_and_varnames <- paste(outcome_labels, outcome_varnames, sep = "_and_")
analysis_methods <- c("BIP", "BIPmixed")
data_splits <- 1:20
analyses_df <- expand.grid(r = possible_r, outcome_label_and_varname = outcome_labels_and_varnames, 
                           analysis_method = analysis_methods, split_seed = data_splits,
                           check_convergence = c(T), include_covars = c(T, F))
# # We only check convergence on 1 of the splits/ method
# analyses_df <- analyses_df %>% 
#   filter( (split_seed != 1 & check_convergence==F) | ( split_seed==1 & check_convergence == T ) )
# Create subdirectory by the slurm_array_task_id
subdir_names <- file.path(output_dir_name, paste0(Sys.Date(), "_task_", 1:nrow(analyses_df)))
analyses_df$subdir_name <- subdir_names

# Save the analysis conditions
about_analyses_file <- file.path(output_dir_name, "analysis_conditions.csv")
analyses_df %>% write.csv(about_analyses_file)

# Save the R environment details to a file in the output directory
session_info_file <- file.path(output_dir_name, "session_info.txt")
sessionInfo() %>%
  capture.output() %>%
  writeLines(session_info_file)

# Extract analysis conditions
analysis_conditions <- analyses_df[slurm_array_task_id, ]
subdir_name <- subdir_names[slurm_array_task_id]
dir.create(subdir_name, showWarnings = FALSE, recursive = TRUE)
r <- analysis_conditions$r
analysis_conditions <- analysis_conditions %>%
  mutate(outcome_label_and_varname = as.character(outcome_label_and_varname)) %>%
  separate(outcome_label_and_varname, into = c("outcome_label", "varname"), sep = "_and_")
outcome_label <- analysis_conditions$outcome_label
outcome_varname <- analysis_conditions$varname
analysis_method <- analysis_conditions$analysis_method
split_seed <- analysis_conditions$split_seed
check_convergence <- analysis_conditions$check_convergence
include_covars <- analysis_conditions$include_covars

# Set seed for data splitting reproducibility purposes
set.seed(split_seed)

complete_data_list <- readRDS("data/2024-10-17_complete_data_list.csv")

if (transform_outcomes) {
  complete_data_list$outcomes[, -1] <- apply(complete_data_list$outcomes[, -1], 2, sqrt)
}

# Split data into 80:20 train:test family-wise split stratified by study site
sample_key <- complete_data_list$outcomes %>% pull(src_subject_id)
cluster_data_path <- 'data/2024-09-11_cluster_data.csv'
cluster_data <- read.csv(cluster_data_path) 

# Function to create train/test split for each site
split_train_test <- function(data, train_frac = 0.8) {
  unique_families <- unique(data$rel_family_id)
  train_families <- sample(unique_families, size = floor(train_frac * length(unique_families)))
  train_data <- data %>% filter(rel_family_id %in% train_families)
  test_data <- data %>% filter(!rel_family_id %in% train_families)
  list(train = train_data, test = test_data)
}

# Split data by site
split_by_site <- split(cluster_data, cluster_data$site_id_l)

# Apply split function to each site
split_data <- lapply(split_by_site, split_train_test)

# Combine train and test sets
train_data <- bind_rows(lapply(split_data, `[[`, "train"))
test_data <- bind_rows(lapply(split_data, `[[`, "test"))

# Verify the split
train_data_summary <- train_data %>% group_by(site_id_l) %>% summarize(n_train = n())
test_data_summary <- test_data %>% group_by(site_id_l) %>% summarize(n_test = n())
# Ensure no families are split across train and test
intersect(train_data$rel_family_id, test_data$rel_family_id)

# Print summaries
print(train_data_summary)
print(test_data_summary)

# Check overall n's and split proportion
n_train <- train_data_summary$n_train %>% sum
n_test <- test_data_summary$n_test %>% sum
cat("n_train:", n_train)
cat("n_test:", n_test)
cat("n_train + n_test:", n_train + n_test)
write.csv(data.frame(n_train=n_train, n_test=n_test), file.path(subdir_name, "train_test_sizes.csv"))

# Assuming you have train_data and test_data data frames with 'src_subject_id'
train_ids <- pull(train_data, src_subject_id)
test_ids <- pull(test_data, src_subject_id)

# Split complete_data_list into train and test lists
split_data <- function(df, train_ids, test_ids) {
  train_df <- df %>% filter(src_subject_id %in% train_ids)
  test_df <- df %>% filter(src_subject_id %in% test_ids)
  list(train = train_df, test = test_df)
}

# Assuming complete_data_list is your list of data frames
split_lists <- lapply(complete_data_list, split_data, train_ids = train_ids, test_ids = test_ids)

# Extract train and test lists
train_list <- lapply(split_lists, `[[`, "train")
test_list <- lapply(split_lists, `[[`, "test")

find_zero_variance_columns <- function(data_list, data_type = c("train", "test"), data_dir = ".", ignore_columns = c("src_subject_id")) {
  # Ensure the data_type is valid
  data_type <- match.arg(data_type)
  
  # Check if 'ela_view' is present in data_list
  if (!"ela_view" %in% names(data_list)) {
    stop("'ela_view' is not present in the data_list.")
  }
  
  df <- data_list$ela_view
  
  # Exclude specified columns from the data frame
  df_filtered <- df %>% dplyr::select(-all_of(ignore_columns))
  
  # Identify columns with zero variance
  zero_var_cols <- df_filtered %>%
    dplyr::select(where(~ var(., na.rm = TRUE) == 0)) %>%
    colnames()
  
  # Create a data frame for zero variance columns
  zero_variance_df <- tibble(
    variable.name = zero_var_cols
  )
  
  # Write the zero_variance data frame to a CSV file in the specified directory
  zero_variance_file <- file.path(data_dir, paste0("zero_variance_", data_type, ".csv"))
  write_csv(zero_variance_df, zero_variance_file)
  
  # Print a message to confirm the file has been written
  cat("Zero variance columns information has been written to:", zero_variance_file, "\n")
  
  return(zero_var_cols)
}

train_zero_variance_cols <- find_zero_variance_columns(train_list, "train", subdir_name)
test_zero_variance_cols <- find_zero_variance_columns(test_list, "test", subdir_name)
zero_variance_cols <- union(train_zero_variance_cols, test_zero_variance_cols)

# Function to remove zero variance columns from data_list
remove_zero_variance_columns <- function(data_list, zero_var_cols, ignore_columns = c("src_subject_id")) {
  if (!"ela_view" %in% names(data_list)) {
    stop("'ela_view' is not present in the data_list.")
  }
  
  df <- data_list$ela_view
  zero_var_cols <- union(zero_var_cols, ignore_columns)
  
  df_cleaned <- df %>%
    dplyr::select(-all_of(zero_var_cols))
  
  data_list$ela_view <- cbind(dplyr::select(data_list$ela_view, "src_subject_id"),
                              df_cleaned)
  
  return(data_list)
}

# Remove zero variance columns if there are any
if (length(zero_variance_cols) > 0) {
  train_list <- remove_zero_variance_columns(train_list, zero_variance_cols, ignore_columns = c("src_subject_id"))
  test_list <- remove_zero_variance_columns(test_list, zero_variance_cols, ignore_columns = c("src_subject_id"))
  
  cat("Zero variance columns have been removed from both training and testing datasets.\n")
} else {
  cat("No zero variance columns to remove.\n")
}

# Option to go into development mode where we subset the data to ensure quicker computation
dev_set <- F
n_subjects <- 1000

# Sample ids from train_list$outcomes$src_subject_id
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
if (dev_set) {
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
  Z_family_train <- as.matrix(train_list_subset$Z_family)
  Z_site_train <- as.matrix(train_list_subset$Z_site)
  
  # Drop design matrices matrices (views that start with Z_)
  train_list_subset <- train_list_subset[!grepl("^Z_", names(train_list_subset))]
  
  print(paste("We have dropped the design matrices."))
  
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

# Fit BIP to ABCD Study Training Set
print("Fitting model...")

# Define the indicators where 2= covariates, 0= omics, and 1= outcome.
if (include_covars) {
  print("Including covariates in model fitting.")
  IndicVar <- c(2, rep(0,4), 1)
} else {
  print("Excluding covariates from model fitting.")
  # Remove covars from IndicVar, trainList and testList
  IndicVar <- c(rep(0,4), 1)
  trainList <- trainList[-1]
}


# When check_covergence, we fit more than 1 chain and store for convergence checks 
# model_fit <- readRDS("~/abcd_multiview/data_analysis_results/2024-09-12_data_analysis/2024-09-12_task_3_job_23529704/model_fit.rds")
# model_fit <- readRDS("~/abcd_multiview/data_analysis_results/2024-09-12_data_analysis/2024-09-12_task_1_job_23529702/model_fit.rds")

if (analysis_method == "BIP") {

  # Fit models
  BIP_start_time <- Sys.time()
  model_fit <- BIP(dataList = trainList, IndicVar = IndicVar, Method = "BIP",
                    nbrcomp = r, sample = n_sample, burnin = n_burnin)
  BIP_end_time <- Sys.time()
  print("BIP required:")
  print(BIP_end_time - BIP_start_time)

  # if (check_convergence) {
  #   model_fit_2 <- BIP(dataList = trainList, IndicVar = IndicVar, Method = "BIP",
  #                    nbrcomp = r, sample = n_sample, burnin = n_burnin)
  # }

} else if (analysis_method == "BIPmixed") {

  # Fit BIPmixed to the data
  BIPmixed_start_time <- Sys.time()
  model_fit <- BIP(dataList = trainList, IndicVar = IndicVar, Method = "BIPmixed",
                         nbrcomp = r, sample = n_sample, burnin = n_burnin,
                         Z_family = Z_family_train, Z_site = Z_site_train)
  BIPmixed_end_time <- Sys.time()
  print("BIPmixed required")
  print(BIPmixed_end_time - BIPmixed_start_time)

  # if (check_convergence) {
  #   model_fit_2 <- BIP(dataList = trainList, IndicVar = IndicVar, Method = "BIPmixed",
  #                      nbrcomp = r, sample = n_sample, burnin = n_burnin,
  #                      Z_family = Z_family_train, Z_site = Z_site_train)
  # }

} else {
  stop("You must provide a valid analysis_method.")
}

if (check_convergence) {
  
  # Store the model fits for later evaluation
  saveRDS(model_fit, paste0(subdir_name, "/model_fit.rds"))
  # saveRDS(model_fit_2, paste0(subdir_name, "/model_fit_2.rds"))
  
}

print("Making preliminary visualizations...")

# Create a named vector for the view labels, removing the "_view" suffix
view_labels <- c("ELA", "sMRI_SA", "sMRI_CT", "fMRI", outcome_label) # Excluding "covariates"
names(view_labels) <- 1:length(view_labels)

## Visualize VarSelMean, VarSelMeanGlobal, CompoSelMean

# Extract the CompoSelMean matrix from the model_fit object
if (analysis_method=="BIPmixed") {
  CompoSelMean <- model_fit$CompoSelMean
} else if (analysis_method=="BIP") {
  # Covariate results removed from model output
  CompoSelMean <- model_fit$CompoSelMean[IndicVar!=2, ]
}


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
if (analysis_method=="BIPmixed") {
  VarSelMean <- model_fit$VarSelMean
  VarSelMeanGlobal <- model_fit$VarSelMeanGlobal
} else if (analysis_method=="BIP") {
  # Covariate results removed from model output
  VarSelMean <- model_fit$VarSelMean[IndicVar!=2]
  VarSelMeanGlobal <- model_fit$VarSelMeanGlobal[IndicVar!=2]
}

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
  filename <- paste0(subdir_name, "/", plot_name, ".png")
  ggsave(filename, plot = plot, width = 8, height = 6, units = "in")
}

## Make & write-out data.frame of variables selected by thresholding VarSelMean & VarSelMeanGlobal
print("Making lists of variables selected...")

# Set the variable selction threshold
mpp_threshold <- 0 # We can filter in post-processing of outputs

feature_names <- trainList[IndicVar!=2] %>% 
  sapply(function(x) colnames(x)) %>% unlist

# Add feature label
VarSelMean_df$Feature <- feature_names

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
VarSelMean_globally$Feature <- feature_names

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
  filename <- paste0(subdir_name, "/", table_name, ".csv")
  # Write the data frame to a CSV file
  write.csv(variables_selected[[table_name]], filename, row.names = FALSE)
}

# Now let's make predictions. We evaluate and aggregate their error later on. 

# Function to order data frames by a specified column
order_data <- function(df, id_column) {
  df %>%
    arrange(!!sym(id_column))
}


if (dev_set) {

  train_ela_features <- c("src_subject_id", trainList$ela_view %>% colnames())
  # Filter the columns of test_list$ela_view
  test_list_subset <- test_list
  test_list_subset$ela_view <- test_list$ela_view[, train_ela_features, drop = FALSE]
  
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
  Z_family_test <- as.matrix(test_list_subset$Z_family)
  Z_site_test <- as.matrix(test_list_subset$Z_site)
  
  # Drop design matrices matrices (views that start with Z_)
  test_list_subset <- test_list_subset[!grepl("^Z_", names(test_list_subset))]
  
  print(paste("We have dropped the design matrices."))
  
  # Convert all data frames in test_list_subset to matrices
  test_list_matrices <- test_list_subset %>%
    map(~ as.matrix(.x))
  
} else {
  stop("Data frames are not ordered by src_subject_id in the same way.")
}

# Reshape test data
testList <- test_list_matrices[!grepl("^outcomes", names(test_list_matrices))] # Exclude outcomes

# Define the indicators where 2= covariates, 0= omics, and 1= outcome.
if (include_covars) {
  print("Including covariates in prediction.")
} else {
  print("Excluding covariates from prediction.")
  # Remove covars from IndicVar, trainList and testList
  testList <- testList[-1]
}

# Get predictions
if (analysis_method == "BIP") {
  y_preds <- BIPpredict(dataListNew = testList, Result = model_fit, meth = "BMA")$ypredict
} else if (analysis_method == "BIPmixed") {
  y_preds <- BIPpredict(dataListNew = testList, Result = model_fit, meth="BMA", 
                                 Z_site = Z_site_test)$ypredict
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
output_csv_path <- file.path(subdir_name, "y_true_vs_y_preds_mspe.csv")
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
output_plot_path <- file.path(subdir_name, "true_y_vs_pred_y_mspe_plot.png")
ggsave(output_plot_path, plot = p, width = 8, height = 6)

# # Store random effect inference, if applicable
# if (analysis_method == "BIPmixed") {
#   # Calculate posterior mean and credible intervals for sigma2_ksi_samples
#   sigma2_ksi_samples_post_burnin <- model_fit$sigma2_ksi_samples[(n_burnin + 1):nrow(model_fit$sigma2_ksi_samples), , drop = FALSE]
#   
#   # Posterior mean
#   posterior_mean_sigma2_ksi <- colMeans(sigma2_ksi_samples_post_burnin)
#   
#   # 2.5% and 97.5% credible intervals (assuming you want a 95% CI)
#   credible_interval_sigma2_ksi <- apply(sigma2_ksi_samples_post_burnin, 2, quantile, probs = c(0.025, 0.975))
#   
#   # Calculate posterior mean and credible intervals for sigma2_theta_samples
#   sigma2_theta_samples_post_burnin <- model_fit$sigma2_theta_samples[(n_burnin + 1):nrow(model_fit$sigma2_theta_samples), , drop = FALSE]
#   
#   # Posterior mean
#   posterior_mean_sigma2_theta <- colMeans(sigma2_theta_samples_post_burnin)
#   
#   # 2.5% and 97.5% credible intervals (95% CI)
#   credible_interval_sigma2_theta <- apply(sigma2_theta_samples_post_burnin, 2, quantile, probs = c(0.025, 0.975))
#   
#   # Combine the results into a data frame for easy viewing
#   results_sigma2_ksi <- data.frame(
#     Posterior_Mean = posterior_mean_sigma2_ksi,
#     CI_Lower = credible_interval_sigma2_ksi[1, ],
#     CI_Upper = credible_interval_sigma2_ksi[2, ]
#   )
#   
#   results_sigma2_theta <- data.frame(
#     Posterior_Mean = posterior_mean_sigma2_theta,
#     CI_Lower = credible_interval_sigma2_theta[1, ],
#     CI_Upper = credible_interval_sigma2_theta[2, ]
#   )
#   
#   # Save these results to CSV files
#   write.csv(results_sigma2_ksi, file.path(subdir_name, "sigma2_ksi_posterior_summary.csv"), row.names = FALSE)
#   write.csv(results_sigma2_theta, file.path(subdir_name, "sigma2_theta_posterior_summary.csv"), row.names = FALSE)
# }

print("Data analysis performed & results stored.")

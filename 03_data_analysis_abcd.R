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
pacman::p_load(tidyverse, reshape2, parallel) 
if (("BIPnet" %in% installed.packages()[,"Package"]) == FALSE) {
  pacman::p_load(devtools)
  devtools::install_github('chekouo/BIPnet')
}

# Define data analysis conditions
possible_r <- c(6, 8)
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
  
# Do analysis_1 by default
analysis_conditions <- analyses_df[slurm_id, ]

print("Analysis started.")

# Extract analysis conditions
r <- analysis_conditions$r
analysis_conditions <- analysis_conditions %>%
  mutate(outcome_label_and_varname = as.character(outcome_label_and_varname)) %>%
  separate(outcome_label_and_varname, into = c("outcome_label", "varname"), sep = "_and_")
outcome_label <- analysis_conditions$outcome_label
outcome_varname <- analysis_conditions$varname
normalize_response <- analysis_conditions$normalize_response
analysis_method <- analysis_conditions$analysis_method

# Focus on the outcome of interest & Drop the outcome that's not
outcome_index <- which( colnames(train_list_matrices$outcomes) == outcome_varname )
dataList <- train_list_matrices[!grepl("^outcomes", names(train_list_matrices))]
dataList$y <- train_list_matrices$outcomes[, outcome_index, drop = FALSE]

# Potentially normalize the response variable
if (normalize_response) {
  mean_y <- mean(dataList$y)
  sd_y <- sd(dataList$y)
  dataList$y <- scale(dataList$y)
  print("Response normalized.")
}

# Ensure that r is at least equal to the number of views. 
if (r < length(train_list_matrices)) {
  stop("Number of latent factor components r must be at least equal to the number of views.")
}

# Fit BIP to ABCD Study Training Set
print("Fitting BIP...")

if (analysis_method == "BIP") {
  # Define the indicators where 2= covariates, 0= omics, and 1= outcome.
  IndicVar <- c(2, rep(0,4), 1)
  model_fit <- BIPnet::BIP(dataList, IndicVar, groupList = NULL, Method = "BIP",
                             nbrcomp = r, sample = 5000, burnin = 1000)
} else if (analysis_method == "BIPmixed") {
  # TODO Source BIPmixed wrapper
  # TODO Fit BIPmixed to the data
} else {
  stop("You must provide a valid analysis_method.")
}

# Directory to save the model fits
models_dir <- "models"
# Directory to save the figures
figures_dir <- "figures"
# Define tables directory
tables_dir <- "tables"
# Define the prefix for any files generated during this analysis
file_prefix <- paste(date_today, "outcome", str_split(outcome_label, " ")[[1]][1], 
                     "r", r, "NormalY", normalize_response, "method", analysis_method, sep = "_")

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
mpp_threshold <- 0.5

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
all_colnames <- unlist(lapply(dataList, colnames))

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

print("Data analysis performed & results stored.")
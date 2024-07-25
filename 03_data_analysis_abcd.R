# Load libraries
library(tidyverse)
library(reshape2)

# Load training data
train_list <- readRDS("data/2024-06-24_train_list.rds")

# Option to turn off subsetting
apply_subsetting <- FALSE

# Define the outcome label and today's date
outcome_label <- "Internalizing Problems"
outcome_varname <- "cbcl_scr_syn_internal_r"
date_today <- Sys.Date()

# Sample 200 ids from train_list$outcomes$src_subject_id
set.seed(123) # Setting seed for reproducibility
sampled_ids <- sample(train_list$outcomes$src_subject_id, 200)

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
  
  # Drop the 2nd column of the outcomes matrix and make an element called y to start with internalizing problems
  train_list_subset$y <- train_list_subset$outcomes %>%
    select(all_of(outcome_varname))
  
  # Drop design matrices matrices (views that start with Z_) and outcomes
  train_list_subset <- train_list_subset[!grepl("^Z_|^outcomes", names(train_list_subset))]
  
  print(paste("We have subsetted the outcome view to", outcome_label, "and dropped design matrices."))
  
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

r <- length(train_list_matrices) # Ensure that r is at least equal to the number of views. 
TrainIndicVar <- c(2,rep(0,4),1)

# Fit BIP to ABCD Study Training Set
abcd_result <- BIPnet::BIP(dataList = train_list_matrices,
                           IndicVar = TrainIndicVar,
                           groupList = NULL, Method = "BIP",
                           nbrcomp = r, sample = 5000, burnin = 1000)

#### ---- Explore BIP Training Results

# Create a named vector for the view labels, removing the "_view" suffix
view_names <- c("covariates", "ela_view", "SA_sMRI_view", "CT_sMRI_view", "fMRI_view", "y")
view_labels <- gsub("_view", "", view_names)
names(view_labels) <- 1:length(view_names)

## -- Make a heatmap of U
EstU_matrix <- abcd_result$EstU

# Convert the matrix to a long format data frame
EstU_df <- as.data.frame(EstU_matrix) %>%
  rownames_to_column(var = "Individual") %>%
  pivot_longer(cols = -Individual, names_to = "LatentFactorComponent", values_to = "Value")

# Rename the Latent Factor Components to be numeric
EstU_df$LatentFactorComponent <- as.integer(gsub("V", "", EstU_df$LatentFactorComponent))

# Create the heatmap
heatmap_EstU <- ggplot(EstU_df, aes(x = LatentFactorComponent, y = Individual, fill = Value)) +
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

print(heatmap_EstU)

## Visualize VarSelMean, VarSelMeanGlobal, CompoSelMean

# Extract the CompoSelMean matrix from the abcd_result object
CompoSelMean <- abcd_result$CompoSelMean

# Convert the matrix to a long format suitable for ggplot2
CompoSelMean_long <- melt(CompoSelMean)
colnames(CompoSelMean_long) <- c("View", "Component", "Probability")

# Create the heatmap with probabilities printed and custom axis labels
heatmap_CompoSelMean <- ggplot(CompoSelMean_long, aes(x = Component, y = View, fill = Probability)) +
  geom_tile() +
  geom_text(aes(label = round(Probability, 2)), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  scale_y_continuous(breaks = 1:length(view_names), labels = view_labels) +
  theme_minimal() +
  labs(title = "Component Selection Probabilities",
       x = "Component",
       y = "View",
       fill = "Probability")

# Print the heatmap
print(heatmap_CompoSelMean)

# Assume libraries have been loaded and view_labels has been constructed upstream

# Extract the VarSelMean and VarSelMeanGlobal lists from the abcd_result object
VarSelMean <- abcd_result$VarSelMean
VarSelMeanGlobal <- abcd_result$VarSelMeanGlobal

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

# Extract the EstSig2 list from the abcd_result object
EstSig2 <- abcd_result$EstSig2

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
  heatmap_EstU = heatmap_EstU,
  heatmap_CompoSelMean = heatmap_CompoSelMean,
  boxplot_VarSelMean = boxplot_VarSelMean,
  boxplot_EstSig2 = boxplot_EstSig2
)

# Directory to save the figures
figures_dir <- "figures"

# Loop through each plot and save it
for (plot_name in names(plots)) {
  plot <- plots[[plot_name]]
  filename <- paste0(figures_dir, "/", date_today, "_", outcome_label, "_", plot_name, ".png")
  ggsave(filename, plot = plot, width = 8, height = 6, units = "in")
}

## Make & write-out data.frame of variables selected by thresholding VarSelMean & VarSelMeanGlobal
# Assume libraries have been loaded and view_labels has been constructed upstream

# Set the variable selction threshold
mpp_threshold <- 0.5

# Extract the VarSelMean and VarSelMeanGlobal lists from the abcd_result object
VarSelMean <- abcd_result$VarSelMean
VarSelMeanGlobal <- abcd_result$VarSelMeanGlobal

# Convert the VarSelMean list to a data frame with appropriate labels
VarSelMean_df <- bind_rows(
  lapply(1:length(VarSelMean), function(i) {
    data.frame(View = view_labels[i], VarSelMean[[i]])
  })
)

# Rename the columns to indicate component ids
colnames(VarSelMean_df) <- c("View", paste("Component", 1:r, sep = "_"))

# Extract column names from each matrix in train_list_matrices
all_colnames <- unlist(lapply(train_list_matrices, colnames))

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
for (name in names(variables_selected)) {
  # Generate filename
  filename <- paste0("tables/", date_today, "_", outcome_label, "_", name, ".csv")
  # Write the data frame to a CSV file
  write.csv(variables_selected[[name]], filename, row.names = FALSE)
  # Print the data frame to the console for verification
  print(variables_selected[[name]])
}

## TODO Visualize EstLoad
# EstLoad is the Estimated Loading obtained using a threshold on the MPP of inclusions for each feature
# Threshold is assumed to be 0.5 by default

## TODO Test prediction & write-out performance statistics

## TODO Generalize analysis to both outcomes
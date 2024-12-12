# setwd("/users/4/neher015/abcd_multiview")

# User Options
base_path <- "~/abcd_multiview/data_analysis_results/2024-09-18_data_analysis"
threshold <- 0.5 # Feature selection threshold
n_samples <- 10000 # Total number of MCMC samples
n_burnin <- 5000 # Number of burnin samples

# Load libraries
library(tidyverse)
library(reshape2) # For melt
library(xtable)
library(rstan) # For convergence checks
library(ggalluvial) # For sankey_plot generation
library(RColorBrewer) # For color palettes
library(grid) # For combining sankey_plot and boxplot .png files 
library(gridExtra) # For combining sankey_plot and boxplot .png files 
library(png) # For reading in the component selection plot

# Define views
views <- c("ELA", "fMRI", "sMRI_CT", "sMRI_SA")

analysis_conditions <- file.path(base_path, "analysis_conditions.csv") %>% read.csv()

# Get all task folder paths (excluding the base folder itself)
task_folders <- list.dirs(path = base_path, full.names = TRUE, recursive = FALSE)

# Initialize a data frame to store the aggregated results
global_feature_results <- data.frame(
  Method = character(),
  Outcome = character(),
  View = character(),
  Split = integer(),
  Feature = character(),
  Probability = numeric(),
  stringsAsFactors = FALSE
)

component_feature_results <- data.frame(
  Method = character(),
  Outcome = character(),
  View = character(),
  Split = integer(),
  Component = integer(),
  Feature = character(),
  Probability = numeric(),
  stringsAsFactors = FALSE
)

component_selection_results <- data.frame(
  Method = character(),
  Outcome = character(),
  View = character(),
  Split = integer(),
  Component = integer(),
  Probability = numeric(),
  stringsAsFactors = FALSE
)

prediction_results <- data.frame(
  Method = character(),
  Outcome = character(),
  Split = integer(),
  Y_true = numeric(),
  y_preds = numeric(),
  MSPE = numeric(),
  stringsAsFactors = FALSE
)

# Convergence check results
convergence_results <- data.frame(
  Method = character(),
  Outcome = character(),
  RHat = numeric(),
  stringsAsFactors = FALSE
)

train_test_sizes <- data.frame(Seed = integer(),
                               n_train = integer(),
                               n_test = integer())

# Iterate over each task folder to process the data
for (i in 1:length(task_folders)) {
  
  # i <- 23
  
  task_folder <- task_folders[i]
  
  # Extract the number between "task" and "job"
  task_number <- str_extract(task_folder, "(?<=task_)\\d+(?=_job)") %>%
    as.numeric()
  
  # Extract the method and outcome info from the folder name (if it's encoded in the folder name)
  method <- analysis_conditions$analysis_method[task_number]  # Replace with appropriate extraction logic if needed
  outcome <- analysis_conditions$outcome_label_and_varname[task_number]  # Replace with appropriate extraction logic if needed
  split <- analysis_conditions$split_seed[task_number]
  
  # Convergence check 
  convergence_check <- analysis_conditions$check_convergence[task_number]
  if (convergence_check & method == "BIPmixed") {
    
    model_fit_1 <- readRDS(file.path(task_folder, "model_fit.rds"))
    model_fit_2 <- readRDS(file.path(task_folder, "model_fit_2.rds"))
    
    # Combine the MCMC samples for both chains
    sigma2_o_samples <- rbind(
      model_fit_1$sigma2_samples,
      model_fit_2$sigma2_samples
    )
    
    sigma2_ksi_samples <- rbind(
      model_fit_1$sigma2_ksi_samples,
      model_fit_2$sigma2_ksi_samples
    )
    
    sigma2_theta_samples <- rbind(
      model_fit_1$sigma2_theta_samples,
      model_fit_2$sigma2_theta_samples
    )
    
    # Combine all parameters into a single matrix
    combined_samples <- cbind(sigma2_o_samples, sigma2_ksi_samples, sigma2_theta_samples)
    
    # Reshape the matrix into an array that is expected by rstan::monitor (chains, iterations, parameters)
    combined_samples_array <- array(
      data = combined_samples, 
      dim = c(n_samples, 2, 24)  # 2 chains, n_samples iterations, and 24 parameters
    )
    
    # Create the parameter names
    param_names <- c("sigma2_o", "sigma2_ksi", paste0("sigma2_theta", 1:22))
    
    # Assign parameter names to the third dimension of the array (the parameters dimension)
    dimnames(combined_samples_array) <- list(NULL, NULL, param_names)
    
    # Use rstan::monitor to estimate Gelman-Rubin stats and credible intervals
    monitor_results <- monitor(
      combined_samples_array, 
      warmup = n_burnin,  # Specify the burn-in period
      print = TRUE
    )
    
  }
  
  train_test_size_i <- file.path(task_folder, "train_test_sizes.csv") %>%
    read.csv()
  train_test_sizes <- rbind(train_test_sizes, data.frame(Seed = split, train_test_size_i))
  
  # Read the MSPE from y_true_vs_y_preds_mspe.csv
  pred_file_i <- file.path(task_folder, "y_true_vs_y_preds_mspe.csv")
  pred_data_i <- read.csv(pred_file_i)
  # Extract the MSPE value
  prediction_results <- rbind(prediction_results, data.frame(Method = method,
                                                             Outcome = outcome, 
                                                             Split = split,
                                                             pred_data_i))

  # Read the globally selected features file (VarSelMean_globally.csv)
  global_file <- file.path(task_folder, "VarSelMean_globally.csv")
  global_data <- read.csv(global_file) %>% filter(Probability > threshold)
  global_feature_results <- rbind(global_feature_results, data.frame(Method = method,
                                                     Outcome = outcome,
                                                     View = global_data$View,
                                                     Split = split,
                                                     Feature = global_data$Feature,
                                                     Probability = global_data$Probability))
  
  component_file <- file.path(task_folder, "VarSelMean_by_component.csv")
  component_data <- read.csv(component_file) %>% filter(Probability > threshold)
  component_feature_results <- rbind(component_feature_results, data.frame(Method = method,
                                                     Outcome = outcome,
                                                     View = component_data$View,
                                                     Split = split,
                                                     Component = component_data$Component,
                                                     Feature = component_data$Feature,
                                                     Probability = component_data$Probability))
  
}

# Calculate avg n_train/ n_test per split for consort diagram
train_test_sizes %>% dplyr::distinct() %>% summarise(n_train_avg = mean(n_train), n_test_avg = mean(n_test))

# Calculate the avg and sd number of distinct features selected over 20 splits
ela_avg_feature_selection <- global_feature_results %>%
  group_by(Method, Outcome, View, Split) %>%  # Group by Split to count distinct features per split
  summarise(NumSelectedFeatures = n_distinct(Feature), .groups = 'drop') %>%
  filter(View == "ELA") %>%  # Filter for ELA view
  group_by(Method, Outcome) %>%
  summarise(
    AvgNumSelectedFeatures = mean(NumSelectedFeatures),
    SdNumSelectedFeatures = sd(NumSelectedFeatures),
    .groups = 'drop'
  )

# Filter features selected in at least 12 splits
ela_features_selected_12_or_more <- global_feature_results %>%
  filter(View == "ELA") %>%  # Filter for ELA view
  group_by(Method, Outcome, View, Feature) %>%  # Group by the relevant variables
  summarise(SplitCount = n_distinct(Split), .groups = 'drop') %>%  # Count the number of splits each feature appears in
  filter(SplitCount >= 12)  # Keep only those features selected in at least 12 splits

# Count the number of features selected in at least 12 splits for each Method and Outcome
ela_feature_count_summary <- ela_features_selected_12_or_more %>%
  group_by(Method, Outcome) %>%
  summarise(FeatureCount = n_distinct(Feature), .groups = 'drop')

# Combine ela feature selection summary into 1 table
ela_feature_selection_summary_combined <- left_join(ela_avg_feature_selection, ela_feature_count_summary)

# Get prediction result summary
prediction_summary <- prediction_results %>%
  group_by(Method, Outcome) %>%
  summarise(mspe_avg = mean(MSPE),
            mspe_se = sd(MSPE))

# Combine data analysis results data raw inputs
data_analysis_results_top <- left_join(ela_feature_selection_summary_combined, prediction_summary)

# Modify the Outcome column and format the table
data_analysis_results_top <- data_analysis_results_top %>%
  # Extract the part of Outcome before "_and_"
  mutate(Outcome = sub("_and_.*", "", Outcome),
         Outcome = str_replace(Outcome, " Problems", "")  # Remove " Problems"
  ) %>%
  # Round the numerical values to 3 decimal places and format Mean (SD) columns
  mutate(
    AvgNumSelectedFeatures = round(AvgNumSelectedFeatures, 2),
    SdNumSelectedFeatures = round(SdNumSelectedFeatures, 3),
    mspe_avg = round(mspe_avg, 3),
    mspe_se = round(mspe_se, 3),
    `ELA Selected` = paste0(AvgNumSelectedFeatures, " (", SdNumSelectedFeatures, ")"),
    `MSPE` = paste0(mspe_avg, " (", mspe_se, ")")
  ) %>%
  # Select the necessary columns
  dplyr::select(Outcome, Method, `ELA Selected`, FeatureCount, `MSPE`) %>%
  rename(`>= 12 Splits` = FeatureCount)

# Convert the modified table to a LaTeX table
data_analysis_latex_table_top <- xtable(data_analysis_results_top)

print(data_analysis_latex_table_top)

# Summarize important features for data analysis results table

# Define view labels
view_labels <- c("ELA", "sMRI_SA", "sMRI_CT", "fMRI")
names(view_labels) <- 1:length(view_labels)

# Rename the values in the View column according to view_labels and filter out rows where View is NA
component_feature_results <- component_feature_results %>%
  mutate(View = view_labels[as.character(View)]) %>%
  filter(!is.na(View))  # This will remove rows with NA in the View column

# Filter, count features with MPP > 0.5, and select top 5 features by MPP for each Component and View
filtered_results <- component_feature_results %>%
  filter(Method == "BIPmixed",
         Component == 3 | Component == 5,
         Split == 1,
         grepl("^Externalizing", Outcome), 
         Probability > 0.5) %>%
  group_by(Component, View) %>%
  summarise(Num_Important_Features = sum(Probability > 0.5),  # Count MPP > 0.5
            Top_5_Features = list(Feature[order(-Probability)][1:5])) %>%  # Select top 5 features by MPP
  ungroup()

# Create a separate data frame with Component, View, and top 5 features
top_5_features_df <- filtered_results %>%
  unnest(cols = c(Top_5_Features))  # Unnest the list column into individual rows

# Create a comma-separated string for the Top_5_Features column and replace underscores with escaped underscores
top_5_features_df <- top_5_features_df %>%
  group_by(Component, View, Num_Important_Features) %>%
  summarise(Top_5_Features = paste(Top_5_Features[!is.na(Top_5_Features)], collapse = ", "), .groups = 'drop') %>%
  mutate(
    View = str_replace_all(View, "_", "\\\\_"),  # Replace underscores with escaped underscores in View
    Top_5_Features = str_replace_all(Top_5_Features, "_", "\\\\_")  # Replace underscores in Top_5_Features
  ) %>%
  rename(`N Important Features` = Num_Important_Features, `Top 5 Features` = Top_5_Features)

# Replace repeated Component values with blanks
top_5_features_df <- top_5_features_df %>%
  group_by(Component) %>%
  mutate(Component = if_else(row_number() == 1, as.character(Component), ""))

data_analysis_table_bottom <- top_5_features_df %>%
  mutate(OutcomeMethod = if_else(row_number() == 1, "BIPmixed on Externalizing Problems", ""))

# Create the LaTeX table using xtable
data_analysis_latex_table_bottom <- xtable(data_analysis_table_bottom, align = c("l", "l", "l", "c", "c", "p{5cm}"))

# Print the LaTeX table with the new column
print(data_analysis_latex_table_bottom, include.rownames = FALSE, sanitize.text.function = identity)

# Generate LaTeX for the top table and remove \begin{table}, \end{table}, and \centering
latex_top <- print(xtable(data_analysis_latex_table_top), 
                   include.rownames = FALSE, 
                   print.results = FALSE)

latex_top_clean <- gsub("\\\\begin\\{table\\}\\[ht\\]|\\\\centering|\\\\end\\{table\\}|\\\\end\\{tabular\\}", "", latex_top)

# Generate LaTeX for the bottom table and remove \begin{table}, \end{table}, and \centering
latex_bottom <- print(xtable(data_analysis_latex_table_bottom, 
                             align = c("l", "l", "l", "c", "c", "p{5cm}")),
                      include.rownames = FALSE, 
                      sanitize.text.function = identity, 
                      print.results = FALSE)

latex_bottom_clean <- gsub("\\\\begin\\{table\\}\\[ht\\]|\\\\begin\\{tabular\\}|\\\\centering|\\\\end\\{table\\}|\\\\end\\{tabular\\}", "", latex_bottom)

# Combine the two tables with a single \begin{table} and \end{table}
combined_latex <- paste0(
  "\\begin{table}[ht]\n\\centering\n",
  latex_top_clean,  # Cleaned LaTeX code for the first table
  "\\hline\n",  # Insert \hline to separate tables
  latex_bottom_clean,  # Cleaned LaTeX code for the second table
  "\\end{tabular}\n\\end{table}"
)

# Print the combined LaTeX code
cat(combined_latex)

# Group by Method, Outcome, View, and Feature, and count repetitions
component_counts <- component_feature_results %>%
  # We filter to 1 BIPmixed result to use for representation
  filter(Method == "BIPmixed" & Split == 1) %>%
  group_by(Method, Outcome, View, Component) %>%
  summarise(FeatureCount = n_distinct(Feature), .groups = 'drop') %>%
  # ad hoc mutation to update view naming
  mutate(View = case_when(
    View == 1 ~ "ELA",
    View == 2 ~ "fMRI",
    View == 3 ~ "sMRI_CT",
    View == 4 ~ "sMRI_SA"
  )) %>% filter(View %in% views)

# Define the function
generate_feature_select_plots <- function(filtered_data, include_relative_contribution_title, outcome_label, 
                                          components_important_to_outcome) {
  
  # Prepare data for Sankey plot: Create source and target pairs and count the number of features for each connection
  sankey_data <- filtered_data %>%
    group_by(View, Component) %>%
    mutate(Component = as.character(Component)) # Convert Component to character for Sankey
  
  # Define a pastel color palette for the components using RColorBrewer
  custom_palette <- brewer.pal(n = 4, name = "Dark2")  
  
  # Ensure the custom_palette is named properly
  views <- unique(filtered_data$View)
  names(custom_palette) <- views
  
  # Step 2: Create the Sankey plot
  sankey_plot <- ggplot(sankey_data, 
                        aes(axis1 = View, axis2 = Component, y = FeatureCount)) +
    geom_alluvium(aes(fill = View), width = 0.2) +
    geom_stratum(width = 1/12, fill = "white", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("View", "Component"), expand = c(0.05, 0.05)) +
    scale_fill_manual(values = custom_palette) +  # Apply the custom color palette
    theme_minimal() +
    theme(
      legend.position = "none",  # Remove legend
      axis.title = element_blank(),  # Remove axis titles
      axis.text.y = element_blank(),  # Show y-axis text
      axis.ticks.y = element_blank(),  # Show y-axis ticks
      panel.grid = element_blank()  # Remove gridlines
    ) +
    labs(x = "", y = outcome_label)  # Show outcome_label on the y-axis
  
  # Step 3: Calculate the relative contribution of each view to each component
  relative_contributions <- filtered_data %>%
    group_by(Component) %>%
    mutate(RelativeContribution = FeatureCount / sum(FeatureCount)) %>%
    ungroup()
  
  # Step 4: Create the bar plot for relative contributions
  if (include_relative_contribution_title) {
    relative_contribution_title <- "View Contribution by Features Selected"
    x_lab_title <- NULL
    show_legend <- T
  } else {
    relative_contribution_title <- NULL
    x_lab_title <- "Component"
    show_legend <- F
  }
  
  # Create the bar plot with highlighted components
  relative_contributions_plot <- ggplot(relative_contributions, aes(x = factor(Component), y = RelativeContribution, fill = View)) +
    geom_bar(stat = "identity", position = "fill", color = "black", alpha = 0.5) +  # Add black outline to each bar
    scale_y_continuous(labels = scales::percent) +    # Display y-axis in percentages
    scale_fill_manual(values = custom_palette) +  # Apply the custom color palette
    labs(x = x_lab_title, y = NULL, title = relative_contribution_title, fill = "View") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),      # Remove all grid lines
      panel.background = element_rect(fill = "white", color = NA),  # Ensure a white background
      plot.background = element_rect(fill = "white", color = NA),    # Ensure the plot area background is white
      legend.position = ifelse(show_legend, "right", "none")  # Show or hide legend based on condition
    ) 
  
  # Return both plots as a list
  return(list(sankey_plot = sankey_plot, bar_plot = relative_contributions_plot))
}

# Filter for now
externalizing_data <- component_counts %>%
  filter(Method == "BIPmixed" & Outcome == "Externalizing Problems_and_cbcl_scr_syn_external_r") %>%
  dplyr::select(-Method, -Outcome)

internalizing_data <- component_counts %>%
  filter(Method == "BIPmixed" & Outcome == "Internalizing Problems_and_cbcl_scr_syn_internal_r") %>%
  dplyr::select(-Method, -Outcome)

feature_select_plots <- list()
feature_select_plots[["externalizing"]] <- generate_feature_select_plots(externalizing_data, include_relative_contribution_title = T, outcome_label = "Externalizing Problems", components_important_to_outcome = c(3, 5))
feature_select_plots[["internalizing"]] <- generate_feature_select_plots(internalizing_data, include_relative_contribution_title = F, outcome_label = "Internalizing Problems", components_important_to_outcome = NULL)

# Let's consider the variance parameter ratios

# 1. Extract slices for sigma2_o, sigma2_ksi, and sigma2_theta
sigma2_o <- combined_samples_array[, , 1]  # sigma2_o is the 1st parameter (3D array: iterations x chains x 1)
sigma2_ksi <- combined_samples_array[, , 2]  # sigma2_ksi is the 2nd parameter (3D array: iterations x chains x 1)
sigma2_theta <- combined_samples_array[, , 3:24]  # sigma2_theta corresponds to the next 22 parameters

sigma2_o_estimate <- monitor(sigma2_o, warmup = n_burnin)
sigma2_ksi_estimate <- monitor(sigma2_ksi, warmup = n_burnin)

# # 2. Calculate the ratio for the first array (sigma2_theta / sigma2_o)
# sigma2_theta_o_ratio <- sigma2_theta / array(sigma2_o, dim = c(dim(sigma2_theta)[1], dim(sigma2_theta)[2], 22))

# 3. Calculate the ratio for the second array (sigma2_theta / sigma2_ksi)
sigma2_theta_ksi_ratio <- sigma2_theta / array(sigma2_ksi, dim = c(dim(sigma2_theta)[1], dim(sigma2_theta)[2], 22))

# Use rstan::monitor to estimate Gelman-Rubin stats and credible intervals
monitor_results <- monitor(
  sigma2_theta_ksi_ratio, # Start with sigma2_theta_ksi_ratio
  warmup = n_burnin,  # Specify the burn-in period
  print = TRUE
)

# Create a new column to classify parameters into groups and assign numeric indices for plotting
monitor_results <- monitor_results %>%
  as.data.frame() %>%
  mutate(Parameter = rownames(monitor_results),  # Get parameter names from rownames
         ParameterType = case_when(
           Parameter == "sigma2_o" ~ "sigma2_o",
           Parameter == "sigma2_ksi" ~ "sigma2_ksi",
           grepl("sigma2_theta", Parameter) ~ "sigma2_theta",
           TRUE ~ "Other"  # Add a fallback for any unclassified cases
         )) %>%
  arrange(ParameterType) %>%  # Ensure parameters are grouped by type
  mutate(ParameterIndex = row_number())  # Create an index for each parameter to use in plotting

# Extract the posterior mean estimate for sigma2_ksi
sigma2_ksi_mean <- monitor_results %>%
  filter(Parameter == "sigma2_ksi") %>%
  pull(mean)

# Extract the y-axis position (index) for sigma2_ksi
sigma2_ksi_index <- monitor_results %>%
  filter(Parameter == "sigma2_ksi") %>%
  pull(ParameterIndex)

variance_forest_plot <- ggplot() +
  # Plot the points and error bars
  geom_point(data = monitor_results, aes(x = mean, y = ParameterIndex)) + # , color = ParameterType)) +
  geom_errorbarh(data = monitor_results, aes(xmin = `2.5%`, xmax = `97.5%`, y = ParameterIndex), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  # Set up axes and labels
  scale_y_continuous(breaks = NULL) +  # Remove parameter names on the y-axis
  labs(title = "Within site variance/ between site variance",
       x = "Mean and 95% CI",
       y = "Site s") +
  theme_minimal() +
  coord_flip()

# Read in component selection plot
# TODO Fix hardcoding
base_path_for_component_selection <- "~/abcd_multiview/data_analysis_results/2024-09-18_data_analysis/2024-09-18_task_4_job_23845260"
model_for_component_selection <- readRDS(file.path(base_path_for_component_selection, "model_fit.rds"))
CompoSelMean <- model_for_component_selection$CompoSelMean

# Create a named vector for the view labels, removing the "_view" suffix
view_labels <- c("ELA", "Externalizing") # Excluding "covariates", and "sMRI_SA", "sMRI_CT", "fMRI" (MPP > 0.95 for imaging for all components)
names(view_labels) <- 1:length(view_labels)

# Convert the matrix to a long format suitable for ggplot2
CompoSelMean_long <- melt(CompoSelMean) 
colnames(CompoSelMean_long) <- c("View", "Component", "Probability")
CompoSelMean_long <- CompoSelMean_long %>% 
  filter(View == 1 | View == 5) %>%
  mutate(View = ifelse(View == 1, view_labels[1], view_labels[2]) %>%
           as.factor())

# Create the heatmap with probabilities printed and custom axis labels
heatmap_CompoSelMean <- ggplot(CompoSelMean_long, aes(x = Component, y = View, fill = Probability)) +
  geom_tile() +
  geom_text(aes(label = round(Probability, 2)), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = expression(hat(gamma[l])^(m)),
       x = "Component",
       y = NULL,
       fill = "Probability")

# Define a custom layout matrix where the forest plot spans the entire last row
layout_matrix <- rbind(
  c(1, 2),  # First row: two columns for externalizing feature selection plots
  c(3, 3)
)

# Combine the plots using the custom layout
combined_plot <- grid.arrange(
  feature_select_plots[["externalizing"]]$bar_plot, 
  feature_select_plots[["externalizing"]]$sankey_plot, 
  variance_forest_plot,
  layout_matrix = layout_matrix
)

# Save the combined plot as a PNG
ggsave("figures/combined_feature_selection_plot.png", plot = combined_plot, width = 12, height = 6)

# Display the combined plot
print(combined_plot)

# We generate a caption for the combined_plot

# Extract necessary values for sigma2_ksi_estimate and sigma2_o_estimate
sigma2_ksi_mean <- round(sigma2_ksi_estimate$mean, 3)
sigma2_ksi_ci <- paste0(round(sigma2_ksi_estimate$`2.5%`, 3), ", ", round(sigma2_ksi_estimate$`97.5%`, 3))

sigma2_o_mean <- round(sigma2_o_estimate$mean, 3)
sigma2_o_ci <- paste0(round(sigma2_o_estimate$`2.5%`, 3), ", ", round(sigma2_o_estimate$`97.5%`, 3))

# Important components to outcome
components_important_to_outcome <- CompoSelMean_long %>% 
  filter(View == "Externalizing" & Probability > 0.5) %>% pull(Component)
components_text <- paste(components_important_to_outcome, collapse = ", ")

# Create the caption
combined_plot_caption <- paste0(
  "BIPmixed analysis of the ABCD Study dataset with outcome \\( \\sqrt{y} \\) raw externalizing problems. ",
  "Internalizing problems results omitted. ",
  "\\textbf{Panel A}. View contributions to latent factor components where contribution is defined as the number of important features, ",
  "those with marginal posterior probabilities \\( >0.5 \\). Views: Early Life Adversity (ELA), functional MRI (fMRI) functional connectivity, ",
  "and 2 from the structural MRI (sMRI) modality, Cortical Thickness (CT) and Surface Area (SA). ",
  "Important components related to the outcome are highlighted with a red dashed box: components ", components_text, ". ",
  "\\textbf{Panel B}. Sankey plot important feature mapping from views to latent components with a red dashed box around important components. ",
  "\\textbf{Panel C}. Within study site variances \\( \\sigma^2_{\\theta_s} \\) to between study site variance \\( \\sigma^2_\\xi \\) credible intervals, ",
  "with the dashed line indicating within and between site variance equivalence. ",
  "Posterior mean (credible interval) for outcome model residual variance \\( \\sigma^{2(0)} \\) is ", sigma2_o_mean, " (", sigma2_o_ci, "), ",
  "and for between study site variance \\( \\sigma^2_\\xi \\) is ", sigma2_ksi_mean, " (", sigma2_ksi_ci, ")."
)

# Print the caption
cat(combined_plot_caption)

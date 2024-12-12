setwd("/users/4/neher015/abcd_multiview")

library(tidyverse)
library(vctrs)
library(xtable)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
source("src/utility_functions.R")

base_path <- "simulation_study_results/2024-11-14_simulation_study"

scenarios <- file.path(base_path, "scenarios.csv") %>% read.csv() %>% rename(task_number = X) 

coverage_results <- data.frame(
  Scenario = integer(),
  Seed = integer(),
  param_type = character(),
  total = integer(),
  within_ci = integer(),
  proportion_within_ci = numeric(),
  stringsAsFactors = FALSE
)

# coverage <- file.path(task_folder, paste0("variance_param_coverage_", task_number, ".csv")) %>%
#   read.csv()
# coverage_results <- rbind(coverage_results, data.frame(Scenario = scenario_id, Seed = seed, coverage))

# Create a credible interval coverage table
# # Process all files and calculate the summary statistics
# coverage_summary <- coverage_results %>%
#   group_by(Scenario, param_type) %>%
#   summarise(
#     proportion_within_ci_mean = mean(proportion_within_ci),
#     proportion_within_ci_sd = sd(proportion_within_ci),
#     .groups = "drop"
#   ) %>% filter(!is.na(scenario_id)) %>%
#   mutate(Method = "BIPmixed")

prediction_results <- data.frame(
  Scenario = integer(),
  Seed = integer(),
  Method = character(),
  MSE = numeric(),
  Bias2 = numeric(),
  Variance = numeric(),
  Mean = numeric(), 
  Correlation = numeric(),
  stringsAsFactors = FALSE
)

prediction_granular_results <- data.frame(
  Scenario = integer(),
  True_Y = numeric(),
  Predicted_Y = numeric(),
  Method = character(),
  Family = integer(),
  Seed = integer(),
  stringsAsFactors = FALSE
) 

feature_selection_results <- data.frame(
  Scenario = integer(),
  Seed = integer(),
  FalsePosRate = numeric(),
  FalseNegRate = numeric(),
  F1measure = numeric(),
  AUC = numeric(),
  Method = character(), 
  View = integer(),
  stringsAsFactors = FALSE
)

# Process all relevant files
for (i in 1:nrow(scenarios)) {
  
  task_folder <- scenarios$subdir_name[i]
  task_number <- scenarios$task_number[i]
  scenario_id <- scenarios$scenario_id[i]
  seed <- scenarios$train_seed[i]
  
  # Check if the folder is empty
  files <- list.files(task_folder)
  
  if (length(files) == 0) {
    # Skip this iteration if the folder is empty
    print(paste("Skipping", task_folder))
    next
  }
  
  # Continue with your operations if the folder is not empty
  print(paste("Processing files in", task_folder))
  
  # Prediction
  prediction_file_path <- file.path(task_folder, paste0("prediction_performance_by_method_", task_number, ".csv")) 
  if (file.exists(prediction_file_path)) {
    prediction <- prediction_file_path  %>%
      read.csv()
    
  } else {
    # Skip this iteration if the folder is empty
    print(paste("Prediction results don't exist for:", task_folder))
    next
  }

  prediction_results <- rbind(prediction_results, data.frame(Scenario = scenario_id, Seed = seed, prediction))
  
  prediction_granular <- file.path(task_folder, paste0("prediction_data_", task_number, ".csv")) %>%
    read.csv()
  prediction_granular_results <- rbind(prediction_granular_results, data.frame(Scenario = scenario_id, Seed = seed, prediction_granular))
  
  # Construct the file path
  feature_selection_path <- file.path(task_folder, paste0("variable_selection_performance_", task_number, ".csv"))
  
  # Check if the file exists
  if (file.exists(feature_selection_path)) {
    # If the file exists, read it and append to feature_selection_results
    feature_selection_data <- read.csv(feature_selection_path)
    feature_selection_results <- rbind(feature_selection_results, data.frame(Scenario = scenario_id, Seed = seed, feature_selection_data))
  }
  
}

# Process all files and calculate the summary statistics
feature_selection_summary <- feature_selection_results %>%
  group_by(Scenario, Method, View) %>%
  summarise(across(c(FalsePosRate, FalseNegRate, F1measure, AUC), 
                   list(mean = ~mean(.), sd = ~sd(.))),
            .groups = "drop") %>% filter(!is.na(scenario_id)) %>% filter(!is.na(scenario_id))

# 1. Extract the results for the 1st view
first_view_selection <- feature_selection_summary %>%
  filter(View == 1) %>%
  rename_with(~paste0(., "_X1"), -c(Scenario, Method, View))

# 2. Calculate the averaged results across all 4 views
average_view_selection <- feature_selection_summary %>%
  group_by(Scenario, Method) %>%
  summarise(across(c(FalsePosRate_mean, FalseNegRate_mean, F1measure_mean, AUC_mean, 
                     FalsePosRate_sd, FalseNegRate_sd, F1measure_sd, AUC_sd), 
                   list(avg = ~mean(.)), .names = "{col}"), .groups = "drop")

# 3. Combine the first view and averaged results into a single data frame
combined_feature_selection_summary <- first_view_selection %>%
  left_join(average_view_selection, by = c("Scenario", "Method"))

# Create a prediction performance table

# Remove replicates with NA prediction results
filtered_prediction_results <- prediction_results %>% 
  # Remove replicates with missing MSE
  filter(!is.na(MSE)==T) %>%
  # Remove outliers
  group_by(Method, Scenario) %>%
  mutate(
    Q1 = quantile(MSE, 0.25, na.rm = TRUE),
    Q3 = quantile(MSE, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = MSE < lower_bound | MSE > upper_bound
  ) %>%
  filter(!is_outlier)

print("N tasks removed due to missing/ outlier MSE:")
print(nrow(prediction_results)-nrow(filtered_prediction_results))
print("Out of N tasks overall:")
print(nrow(prediction_results))

# Count the number of rows for each Method
prediction_results %>%
  group_by(Method, Scenario) %>%
  summarise(n = n())

# Count the number of rows for each Method
filtered_prediction_results %>%
  group_by(Method, Scenario) %>%
  summarise(n = n())

# Calculate the summary statistics
prediction_summary <- filtered_prediction_results %>%
  group_by(Scenario, Method) %>%
  summarise(
    MSE_mean = mean(MSE),
    MSE_sd = sd(MSE),
    Bias2_mean = mean(Bias2),
    Bias2_sd = sd(Bias2),
    Variance_mean = mean(Variance),
    Variance_sd = sd(Variance),
    Mean_mean = mean(Mean),
    Mean_sd = sd(Mean),
    Correlation_mean = mean(Correlation),
    Correlation_sd = sd(Correlation),
    .groups = "drop"
  )

# # Calculate quartiles of Predicted_Y and filter to smallest (1st) and largest (4th) quartiles
# quartile_summary_mspe <- prediction_granular_results %>%
#   group_by(Scenario, Seed, Method) %>%
#   mutate(Quartile = ntile(Predicted_Y, 4)) %>%  # Assign quartiles based on Predicted_Y
#   group_by(Scenario, Seed, Method, Quartile) %>%
#   reframe(MSPE = (True_Y - Predicted_Y)^2) %>%   # Calculate MSPE for given replicate
#   group_by(Scenario, Method, Quartile) %>%
#   summarise(
#     avg_mspe = mean(MSPE),                      # Average MSPE in the quartile
#     sd_mspe = sd(MSPE),                         # Standard deviation of MSPE in the quartile
#     .groups = 'drop'
#   )
# 
# # View the filtered quartile summary for MSPE
# head(quartile_summary_mspe)

# Step 1: Combine mean and sd for each metric
results_aggregated <- left_join(prediction_summary,
                                combined_feature_selection_summary) 
results_formatted <- results_aggregated %>%
  mutate(
    MSE = paste0(sprintf("%.3f", MSE_mean), " (", sprintf("%.3f", MSE_sd), ")"),
    Bias2 = paste0(sprintf("%.3f", Bias2_mean), " (", sprintf("%.3f", Bias2_sd), ")"),
    Variance = paste0(sprintf("%.3f", Variance_mean), " (", sprintf("%.3f", Variance_sd), ")"),
    Mean = paste0(sprintf("%.3f", Mean_mean), " (", sprintf("%.3f", Mean_sd), ")"),
    Correlation = paste0(sprintf("%.3f", Correlation_mean), " (", sprintf("%.3f", Correlation_sd), ")"),
    FalsePosRate = paste0(sprintf("%.3f", FalsePosRate_mean), " (", sprintf("%.3f", FalsePosRate_sd), ")"),
    FalseNegRate = paste0(sprintf("%.3f", FalseNegRate_mean), " (", sprintf("%.3f", FalseNegRate_sd), ")"),
    F1measure = paste0(sprintf("%.3f", F1measure_mean), " (", sprintf("%.3f", F1measure_sd), ")"),
    AUC = paste0(sprintf("%.3f", AUC_mean), " (", sprintf("%.3f", AUC_sd), ")"),
    # proportion_within_ci = paste0(sprintf("%.3f", proportion_within_ci_mean), " (", sprintf("%.3f", proportion_within_ci_sd), ")")
  ) %>%
  # We exclude the predicted mean
  dplyr::select(Scenario, Method, MSE, -Bias2, Variance, -Mean, -Correlation, FalsePosRate, FalseNegRate, -F1measure, AUC) # %>% # , proportion_within_ci) %>%
  # rename(`Variance Parameter Coverage` = proportion_within_ci)

# Step 2: Pivot the data to long format
results_long <- results_formatted %>%
  pivot_longer(cols = -c(Scenario, Method), names_to = "Metric", values_to = "Value")

# Check for duplicates
results_long %>%
  group_by(Metric, Scenario, Method) %>%
  summarise(n = n(), .groups = 'drop') %>%
  filter(n > 1)

# Handle duplicates by taking the first occurrence (or mean, max, etc.)
results_wide <- results_long %>%
  unite("Scenario_Method", Scenario, Method, sep = " - ") %>%  # Combine Scenario and Method
  pivot_wider(names_from = Scenario_Method, values_from = Value, values_fn = list(Value = first)) %>%
  t()

# Assuming `results_wide` is your matrix:
results_df <- as.data.frame(results_wide[-1, ]) 
colnames(results_df) <- results_wide[1, ]           # Assign the first row as column names

results_df_named <- cbind(data.frame("Scenario - Method" = rownames(results_df)), results_df)
rownames(results_df_named) <- NULL

# Step 4: Convert to LaTeX table using xtable
latex_table <- xtable(results_df_named, 
                      # caption = "Summary of model performance metrics by Scenario and Method. Metrics represent the mean (standard deviation).",
                      digits = 3)  # Ensure 3 digits are printed)

# Print LaTeX table with caption placed at the top
print(latex_table, include.rownames = FALSE, caption.placement = "top")

# Make prediction scatterplots
# Filter the data for Seed == 1
filtered_data <- prediction_granular_results %>%
  filter(Seed == 6)

# Create a named vector for the facet labels
scenario_labels <- c("1" = "Scenario 1", "2" = "Scenario 2", "3" = "Scenario 3")

# Create the scatterplots for each scenario, colored by Method
scatter_plot <- filtered_data %>% 
  filter(Method != "RandMVLearn") %>%
  mutate(Method = ifelse(Method=="2step", "PCA2Step", Method))   %>%
  mutate(Method = ifelse(Method=="CooperativeLearning", "Cooperative", Method)) %>%
  ggplot(aes(x = Predicted_Y, y = True_Y, color = Method)) +
  geom_point(alpha = 0.15) +  # Add points with transparency
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Perfect calibration line
  facet_wrap(~ Scenario, nrow = 1, labeller = as_labeller(scenario_labels)) +  # Custom facet labels
  labs(x = expression(hat(y)), y = expression(y), title = NULL) +
  theme_minimal() +
  theme(legend.position = "top")  # Move legend to the bottom for better visualization

# Obtain the third color from the default palette
default_colors <- scale_color_discrete()$palette(4)  # Get first 3 colors from the default palette
# true_y_color <- default_colors[2]  # Select the third color

# Define a named vector for colors, assigning "True Y" the third color from the default palette
cb_palette <- c(
  "BIP" = default_colors[1], 
  "BIPmixed" = default_colors[2],
  "Cooperative" = default_colors[3],
  "PCA2Step" = default_colors[4],
  "True Y" = "red")

# Create a combined data frame with True_Y and Predicted_Y by Method for plotting
combined_data <- filtered_data %>%
  mutate(Method = ifelse(Method=="2step", "PCA2Step", Method)) %>%
  mutate(Method = ifelse(Method=="CooperativeLearning", "Cooperative", Method)) %>%
  pivot_longer(cols = c(True_Y, Predicted_Y), 
               names_to = "Type", 
               values_to = "Y_Value") %>%
  mutate(Type = ifelse(Type == "True_Y", "True Y", Method)) %>%
  filter(Method != "RandMVLearn") 

# Create the combined density plot stratified by Method using the customized palette
combined_density_plot <- ggplot(combined_data, aes(x = Y_Value, color = Type)) +
  geom_density(alpha = 0.3, size = 0.5) +
  facet_wrap(~ Scenario, nrow = 1) +
  scale_color_manual(values = cb_palette, guide = guide_legend(override.aes = list(color = cb_palette))) +
  labs(x = "y (True and Predicted)", y = "Density", title = NULL) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", size = 0.5),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank())  # Hide facet labels

# Print the plot
print(combined_density_plot)

# Specify the layout matrix
layout_matrix <- rbind(
  c(1),  # First row is scatter_plot spanning 2 columns
  c(1),  # Second row is scatter_plot spanning 2 columns
  c(2)   # Third row is combined_density_plot spanning 2 columns
)

# Combine the plots with grid.arrange
combined_plot <- grid.arrange(
  scatter_plot,
  combined_density_plot,
  layout_matrix = layout_matrix
)

# Display the combined plot
print(combined_plot)

# Specify the file path to save the plot
file_path <- "figures/combined_simulations_plot.png"

# Save the combined density plot to the specified path
ggsave(filename = file_path, 
       plot = combined_plot, 
       width = 8,  # Width of the plot in inches
       height = 8,  # Height of the plot in inches
       dpi = 300)  # Resolution of the plot (dots per inch)


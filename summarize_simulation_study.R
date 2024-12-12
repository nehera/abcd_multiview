library(tidyverse)
library(dplyr)
library(purrr)
library(stringr)
library(vctrs)

simulation_study_path <- "simulation_study_results/2024-09-03_test_simulation_study"
scenarios <- file.path(simulation_study_path, "scenarios.csv") %>%
  read.csv() %>%
  # Filter to scenarios in the manuscript
  filter(scenario_id %in% c(1, 6, 8) & !is.na(scenario_id)) %>%
  # Rename the scenario_id to match the manuscript
  mutate(scenario_id = case_when(
    scenario_id == 6 ~ 2,
    scenario_id == 8 ~ 3,
    TRUE ~ scenario_id  # Keep the original value for scenario_id 1
  ))

# Function to extract task_id from file path
extract_task_id <- function(file_path) {
  match <- str_match(file_path, "_task_(\\d+)_job_")
  as.numeric(match[, 2])  # Extract the task_id from the match
}

# List all relevant files
all_files <- list.files(simulation_study_path, pattern = "variable_selection_performance_\\d+\\.csv", full.names = TRUE, recursive = TRUE)

# Extract task_id from file paths
file_info <- tibble(file_path = all_files) %>%
  mutate(task_id = map_dbl(file_path, extract_task_id))

# Join the task_id with the scenarios data to map it back to scenario_id
file_paths <- file_info %>%
  left_join(scenarios %>% mutate(task_id = row_number()), by = "task_id") %>%
  filter(!is.na(file_path))  # Remove any file without a corresponding scenario_id

# Select only the columns needed for processing
file_paths_selected <- file_paths %>%
  dplyr::select(file_path, scenario_id)

# Function to process files and add scenario_id
process_file <- function(file_path, scenario_id) {
  read.csv(file_path) %>%
    mutate(scenario_id = scenario_id)
}

# Process all files and calculate the summary statistics
summary_stats <- file_paths_selected %>%
  pmap(function(file_path, scenario_id) {
    process_file(file_path, scenario_id)
  }) %>%
  list_rbind() %>%
  group_by(scenario_id, Method, View) %>%
  summarise(across(c(FalsePosRate, FalseNegRate, F1measure, AUC), 
                   list(mean = ~mean(.), sd = ~sd(.))),
            .groups = "drop") %>% filter(!is.na(scenario_id)) %>% filter(!is.na(scenario_id))

# View the resulting summary statistics
print(summary_stats)

library(xtable)

# Assuming 'summary_stats' is the data frame from previous steps

# 1. Extract the results for the 1st view
first_view_stats <- summary_stats %>%
  filter(View == 1) %>%
  rename_with(~paste0(., "_X1"), -c(scenario_id, Method, View))

# 2. Calculate the averaged results across all 4 views
average_view_stats <- summary_stats %>%
  group_by(scenario_id, Method) %>%
  summarise(across(c(FalsePosRate_mean, FalseNegRate_mean, F1measure_mean, AUC_mean, 
                     FalsePosRate_sd, FalseNegRate_sd, F1measure_sd, AUC_sd), 
                   list(avg = ~mean(.)), .names = "{col}"), .groups = "drop") %>%
  rename_with(~paste0(., "_Xavg"), -c(scenario_id, Method))

# 3. Combine the first view and averaged results into a single data frame
combined_stats <- first_view_stats %>%
  left_join(average_view_stats, by = c("scenario_id", "Method"))

# 4. Reshape the data frame for better formatting in xtable
final_stats <- combined_stats %>%
  pivot_longer(cols = -c(scenario_id, Method, View), 
               names_to = c("Metric", "Statistic", "ViewType"),
               names_sep = "_", 
               values_to = "Value") %>%
  pivot_wider(names_from = c(ViewType, Statistic), values_from = Value) %>%
  arrange(scenario_id, Method) %>%
  dplyr::select(-View)

# View the resulting final_stats data frame
print(final_stats)

# 5. Pass the data frame to xtable to create a LaTeX table
latex_table <- xtable(final_stats)

# 6. Print the LaTeX table
print(latex_table, include.rownames = FALSE)

# Create a prediction performance table

library(dplyr)
library(purrr)
library(stringr)
library(xtable)

# Function to extract task_id from file path
extract_task_id <- function(file_path) {
  match <- str_match(file_path, "_task_(\\d+)_job_")
  as.numeric(match[, 2])  # Extract the task_id from the match
}

# List all relevant prediction performance files
all_files <- list.files(simulation_study_path, pattern = "prediction_performance_by_method_\\d+\\.csv", full.names = TRUE, recursive = TRUE)

# Extract task_id from file paths
file_info <- tibble(file_path = all_files) %>%
  mutate(task_id = map_dbl(file_path, extract_task_id))

# Join the task_id with the scenarios data to map it back to scenario_id
file_paths <- file_info %>%
  left_join(scenarios %>% mutate(task_id = row_number()), by = "task_id") %>%
  filter(!is.na(file_path))  # Remove any file without a corresponding scenario_id

# Select only the columns needed for processing
file_paths_selected <- file_paths %>%
  select(file_path, scenario_id)

# Function to process files and add scenario_id
process_file <- function(file_path, scenario_id) {
  read.csv(file_path) %>%
    mutate(scenario_id = scenario_id)
}

# Process all files and calculate the summary statistics
summary_stats <- file_paths_selected %>%
  pmap_dfr(function(file_path, scenario_id) {
    process_file(file_path, scenario_id)
  }) %>%
  group_by(scenario_id, Method) %>%
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

# Create a LaTeX table using xtable
latex_table <- xtable(summary_stats)

# Print the LaTeX table
print(latex_table, include.rownames = FALSE)


# Create a credible interval coverage table

library(dplyr)
library(purrr)
library(stringr)
library(xtable)

# Function to extract task_id from file path
extract_task_id <- function(file_path) {
  match <- str_match(file_path, "_task_(\\d+)_job_")
  as.numeric(match[, 2])  # Extract the task_id from the match
}

# List all relevant variance parameter coverage files
all_files <- list.files(simulation_study_path, pattern = "variance_param_coverage_\\d+\\.csv", full.names = TRUE, recursive = TRUE)

# Extract task_id from file paths
file_info <- tibble(file_path = all_files) %>%
  mutate(task_id = map_dbl(file_path, extract_task_id))

# Join the task_id with the scenarios data to map it back to scenario_id
file_paths <- file_info %>%
  left_join(scenarios %>% mutate(task_id = row_number()), by = "task_id") %>%
  filter(!is.na(file_path))  # Remove any file without a corresponding scenario_id

# Select only the columns needed for processing
file_paths_selected <- file_paths %>%
  select(file_path, scenario_id)

# Function to process files and add scenario_id
process_file <- function(file_path, scenario_id) {
  read.csv(file_path) %>%
    mutate(scenario_id = scenario_id)
}

# Process all files and calculate the summary statistics
summary_stats <- file_paths_selected %>%
  pmap_dfr(function(file_path, scenario_id) {
    process_file(file_path, scenario_id)
  }) %>%
  group_by(scenario_id, param_type) %>%
  summarise(
    proportion_within_ci_mean = mean(proportion_within_ci),
    proportion_within_ci_sd = sd(proportion_within_ci),
    .groups = "drop"
  ) %>% filter(!is.na(scenario_id))

# Create a LaTeX table using xtable
latex_table <- xtable(summary_stats)

# Print the LaTeX table
print(latex_table, include.rownames = FALSE)

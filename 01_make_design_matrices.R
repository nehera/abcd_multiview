library(tidyverse)

## -- User Arguments

# Set the path for raw data files
data_dir <- '/Users/aidanneher/Library/CloudStorage/Box-Box/ABCD Tabulated Data/5.1/core'
# Location of desired output directory - if NULL, will output into working directory
out_dir <- '/Users/aidanneher/Documents/GitHub/abcd_multiview/data'
# Date you used in output name - if NULL, will use output from Sys.Date() (current date)
out_date <- NULL
# Initials or other string you want in output naming - no NULL option here
out_initials <- 'AN'

## -- Make Design Matrices

# Extract clustering information
path <- file.path(data_dir, "abcd-general", "abcd_y_lt.csv")
cluster_data <- read.csv(path) %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(src_subject_id, rel_family_id, site_id_l) %>%
  arrange(site_id_l, rel_family_id)

# Create the Z_family design matrix
subject_ids <- unique(cluster_data$src_subject_id)
family_ids <- unique(cluster_data$rel_family_id)
Z_family <- matrix(0, nrow = nrow(cluster_data), ncol = length(family_ids))
rownames(Z_family) <- subject_ids
colnames(Z_family) <- family_ids

for (i in seq_along(family_ids)) {
  family_id <- family_ids[i]
  Z_family[cluster_data$rel_family_id == family_id, i] <- 1
}

# Create the Z_site design matrix
site_ids <- unique(cluster_data$site_id_l)
Z_site <- matrix(0, nrow = nrow(cluster_data), ncol = length(site_ids))
rownames(Z_site) <- subject_ids
colnames(Z_site) <- site_ids

for (i in seq_along(site_ids)) {
  site_id <- site_ids[i]
  Z_site[cluster_data$site_id_l == site_id, i] <- 1
}

## -- Verify Z_family and Z_site reflect observed n families and n sites
# Calculate sums using apply for Z_family and Z_site
family_sums <- apply(Z_family, 2, sum)
site_sums <- apply(Z_site, 2, sum)
# Convert family_sums and site_sums to data frames for comparison
family_sums_df <- data.frame(rel_family_id = names(family_sums), n_observations = family_sums)
site_sums_df <- data.frame(site_id_l = names(site_sums), n_observations = site_sums)
# Summarize by number of observations in each family
family_summary <- cluster_data %>%
  group_by(rel_family_id) %>%
  summarize(n_observations = n()) %>%
  arrange(desc(n_observations))
# Summarize by number of observations in each site
site_summary <- cluster_data %>%
  group_by(site_id_l) %>%
  summarize(n_observations = n()) %>%
  arrange(desc(n_observations))
# Compare the summaries with the sums
compare_family <- merge(family_summary, family_sums_df, by = "rel_family_id", suffixes = c("_summary", "_apply"))
compare_site <- merge(site_summary, site_sums_df, by = "site_id_l", suffixes = c("_summary", "_apply"))
# Check if the summaries match
all(compare_family$n_observations_summary == compare_family$n_observations_apply)
all(compare_site$n_observations_summary == compare_site$n_observations_apply)
# Print mismatches if any
mismatched_families <- compare_family %>%
  filter(n_observations_summary != n_observations_apply)
mismatched_sites <- compare_site %>%
  filter(n_observations_summary != n_observations_apply)
cat("n mismatched families or sites:", nrow(mismatched_families) + nrow(mismatched_sites))

## -- Map Family to Site
# Compute the transpose of Z_family
Z_family_transpose <- t(Z_family)
# Perform matrix multiplication
Z_family_to_site <- Z_family_transpose %*% Z_site
# Convert to data frame for better readability
Z_family_to_site_df <- as.data.frame(Z_family_to_site)
# Add row and column names for better interpretation
rownames(Z_family_to_site_df) <- colnames(Z_family)  # Family IDs
colnames(Z_family_to_site_df) <- colnames(Z_site)    # Site IDs
# Checkout the resulting design matrix
print(head(Z_family_to_site_df))

## -- Calculate Summary Stats
library(gridExtra)

# Calculate the number of families in each site
families_per_site <- cluster_data %>%
  group_by(site_id_l) %>%
  summarize(n_families = n_distinct(rel_family_id)) %>%
  arrange(desc(n_families))

# Calculate the number of families per site
families_per_site <- cluster_data %>%
  group_by(site_id_l) %>%
  summarize(n_families = n_distinct(rel_family_id))

# Calculate the number of individuals per family
individuals_per_family <- cluster_data %>%
  group_by(rel_family_id) %>%
  summarize(n_individuals = n())

# Calculate the number of individuals per site
individuals_per_site <- cluster_data %>%
  group_by(site_id_l) %>%
  summarize(n_individuals = n())

# Function to identify outliers
identify_outliers <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25)
  Q3 <- quantile(data[[column]], 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  outliers <- data %>%
    filter(data[[column]] < lower_bound | data[[column]] > upper_bound)
  return(outliers)
}

# Identify outliers for each dataset
outliers_families_per_site <- identify_outliers(families_per_site, "n_families")
outliers_individuals_per_family <- identify_outliers(individuals_per_family, "n_individuals")
outliers_individuals_per_site <- identify_outliers(individuals_per_site, "n_individuals")

# Create boxplots
plot_families_per_site <- ggplot(families_per_site, aes(x = 0, y = n_families)) +
  geom_boxplot(fill = "blue", color = "black", alpha = 0.7) +
  geom_text(data = outliers_families_per_site, aes(x = 0, label = n_families), 
            position = position_jitter(width = 0.2, height = 0), 
            hjust = -0.3, color = "red") +
  labs(title = "Number of Families per Site",
       y = "Number of Families",
       x = "") +
  coord_flip() +
  theme_minimal()

plot_individuals_per_site <- ggplot(individuals_per_site, aes(x = 0, y = n_individuals)) +
  geom_boxplot(fill = "orange", color = "black", alpha = 0.7) +
  geom_text(data = outliers_individuals_per_site, aes(x = 0, label = n_individuals), 
            position = position_jitter(width = 0.2, height = 0), 
            hjust = -0.3, color = "red") +
  labs(title = "Number of Individuals per Site",
       y = "Number of Individuals",
       x = "") +
  coord_flip() +
  theme_minimal()

# Calculate the counts for each number of individuals
counts <- individuals_per_family %>%
  group_by(n_individuals) %>%
  summarize(count = n())

# Create a histogram with counts labeled above each bar
plot_individuals_per_family <- ggplot(individuals_per_family, aes(x = n_individuals)) +
  geom_histogram(binwidth = 1, fill = "green", color = "black", alpha = 0.7) +
  geom_text(data = counts, aes(x = n_individuals, y = count, label = count), vjust = -0.5) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Histogram of Number of Individuals per Family",
       x = "Number of Individuals",
       y = "Count of Families") +
  theme_minimal()

# Arrange the boxplots in a grid
grid_plots <- grid.arrange(plot_individuals_per_family, plot_families_per_site, plot_individuals_per_site, ncol = 1)

# Define function to save the plot with date and initials
save_plot_with_date_initials <- function(plot, out_initials, out_dir = NULL, out_date = NULL) {
  # Use current date if out_date is NULL
  if (is.null(out_date)) { out_date <- Sys.Date() }
  
  # Construct the filename
  file_name <- sprintf("%s_%s_abcd_clustering_summary.png", out_date, out_initials)
  
  # Define the output path
  output_path <- ifelse(is.null(out_dir), file_name, file.path(out_dir, file_name))
  
  # Save the plot to a .png file
  ggsave(output_path, plot, width = 10, height = 15)
  
  # Print the output path for verification
  print(paste("Plot saved to:", output_path))
}

# Example usage:
save_plot_with_date_initials(grid_plots, out_initials = "AN", out_dir = "figures")

## -- Write Design Matrices

# Define the function to write the matrix to CSV
write_matrix_to_csv <- function(matrix_data, matrix_name, out_dir = NULL, out_initials, out_date = NULL) {
  # Use current date if out_date is NULL
  if (is.null(out_date)) { out_date <- Sys.Date() }
  
  # Construct the filename
  file_name <- sprintf("%s_%s_%s.csv", out_date, out_initials, matrix_name)
  
  # Define the output path
  output_path <- ifelse(is.null(out_dir), file_name, file.path(out_dir, file_name))
  
  # Write the matrix to CSV file
  write_csv(as.data.frame(matrix_data), output_path)
  
  # Print the output path for verification
  print(paste(matrix_name, "written to:", output_path))
}

# Write matrices to CSV files
write_matrix_to_csv(Z_family, "family_Z", out_dir, out_initials = "AN")
write_matrix_to_csv(Z_site, "site_Z", out_dir, out_initials = "AN")
write_matrix_to_csv(Z_family_to_site, "family_to_site_Z", out_dir, out_initials = "AN")
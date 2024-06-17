library(tidyverse)
library(table1)
library(htmltools)
library(gridExtra)

# Define the directory, date, processor, and file suffixes
data_dir <- "data"
figures_dir <- "figures"
data_processing_date <- "2024-06-06"
data_processor <- "AN"
file_suffix <- c("outcomes", "covariates", "ela_view", "SA_sMRI_view", "CT_sMRI_view", "fMRI_view")

# Function to read data and name list elements
read_data <- function(suffix) {
  file_path <- file.path(data_dir, paste0(data_processing_date, "_", data_processor, "_", suffix, ".csv"))
  data <- read_csv(file_path)
  return(data)
}

# Read in the data frames and create a named list
data_list <- lapply(file_suffix, read_data)
names(data_list) <- file_suffix

# Checkout what's been loaded
lapply(data_list, dim)

## -- Calculate and Save the Number of Variables per Data Type

# Calculate the number of columns excluding 'src_subject_id'
variables_per_datatype <- sapply(data_list, function(df) ncol(df) - 1)

# Print the number of variables for each data frame
print(variables_per_datatype)

# Prepare variables_per_datatype for saving
variables_per_datatype_df <- tibble(data_type = names(variables_per_datatype), n_variables = variables_per_datatype)

# Save the n_variables data to a CSV file
variables_data_path <- file.path(data_dir, paste0(data_processing_date, "_", data_processor, "_n_variables_per_datatype.csv"))
write_csv(variables_per_datatype_df, variables_data_path)

## -- Combine, Filter, & Summarize Missingness

# Calculate the number of rows with any missing data for each data frame
missing_per_datatype <- sapply(data_list, function(df) sum(!complete.cases(df)))

# Print the number of rows with missing data for each data frame
print(missing_per_datatype)

# Prepare missing_per_datatype for saving
missing_per_datatype_df <- tibble(data_type = names(missing_per_datatype), n_missing = missing_per_datatype)

# Save the n_missing data to a CSV file
missing_data_path <- file.path(data_dir, paste0(data_processing_date, "_", data_processor, "_n_missing_per_datatype.csv"))
write_csv(missing_per_datatype_df, missing_data_path)

# Calculate the number of missing values per variable for each data frame
missing_per_variable <- lapply(data_list, function(df) {
  colSums(is.na(df))
})

# Combine the results into a data frame
missing_per_variable_df <- bind_rows(missing_per_variable, .id = "data_type") %>%
  pivot_longer(-data_type, names_to = "variable", values_to = "n_missing") %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))

# Print the top 10 features with the most missing values
top_10_missing <- missing_per_variable_df %>% slice_max(n = 10, order_by = n_missing)
print(top_10_missing)

# Save the missing data per variable to a CSV file
missing_variable_path <- file.path(data_dir, paste0(data_processing_date, "_", data_processor, "_n_missing_per_variable.csv"))
write_csv(missing_per_variable_df, missing_variable_path)

# Create a list of data frames with only complete cases
complete_data_list <- lapply(data_list, function(df) df %>% drop_na())

# Find the intersection of src_subject_id across all complete data frames
common_subjects <- reduce(complete_data_list, function(df1, df2) inner_join(df1, df2, by = "src_subject_id")) %>%
  select(src_subject_id)

# filter all data frames in list to common_subjects
complete_data_list <- lapply(complete_data_list, function(df) df[df$src_subject_id %in% common_subjects$src_subject_id, ])

# Print the number of common subjects
cat("n Common Subjects Across Complete Case Outcomes, Covariates, & Views:", nrow(common_subjects), "\n")

# Save the common_subjects to a CSV file
sample_key_path <- file.path(data_dir, paste0(data_processing_date, "_", data_processor, "_sample_key.csv"))
write_csv(common_subjects, sample_key_path)

## -- Summarize Sample's Outcomes & Covariates

# Merge outcomes and covariates data frames and filter to common subjects
merged_data <- inner_join(data_list$outcomes, data_list$covariates, by = "src_subject_id") %>%
  filter(src_subject_id %in% common_subjects$src_subject_id)

# Convert all columns except 'src_subject_id' to numeric
merged_data <- merged_data %>%
  mutate_at(vars(-src_subject_id), as.numeric)

# Relevel race variables to be factors with levels "No" and "Yes"
race_vars <- colnames(merged_data)[grepl("_race$", colnames(merged_data))]
merged_data <- merged_data %>%
  mutate_at(vars(one_of(race_vars)), ~ factor(., levels = c(0, 1), labels = c("No", "Yes")))

# Relevel demo_sex_v2
merged_data <- merged_data %>%
  mutate(demo_sex_v2 = factor(demo_sex_v2, levels = c(1, 2, 3, 4), labels = c(
    "Male", "Female", "Intersex-Male", "Intersex-Female"
  )))

# Define labels
table1_labels <- c(
  "Internalizing Problems (Raw)",
  "Externalizing Problems (Raw)",
  "Sex (At Birth)",
  "Age (Months)",
  "White (Yes/No)",
  "Black (Yes/No)",
  "American Indian or Native American (Yes/No)",
  "Native Hawaiian or Pacific Islander (Yes/No)",
  "Asian (Yes/No)",
  "Other Race (Yes/No)",
  "Missing Race & Ethnicity (Yes/No)",
  # We do not consider indigenous race superset in our analysis. 
  # "Indigenous Race (Yes/No)",
  "Hispanic/Latinx (Yes/No)",
  "Total Family Income (Past 12 Months)",
  "Highest Parent Education Completed",
  "Parent Marital Status"
)

# Define factor levels and labels
merged_data <- merged_data %>%
  mutate(
    # We now treat income and highest education as continuous random variables
    # demo_comb_income_v2_bl = factor(demo_comb_income_v2_bl, levels = 1:10, labels = c(
    #   "Less than $5,000", "$5,000 through $11,999", "$12,000 through $15,999", 
    #   "$16,000 through $24,999", "$25,000 through $34,999", "$35,000 through $49,999",
    #   "$50,000 through $74,999", "$75,000 through $99,999", "$100,000 through $199,999",
    #   "$200,000 and greater")),
    # highest_demo_ed_bl = factor(highest_demo_ed_bl, levels = c(1:22), labels = c(
    #   "1st grade", "2nd grade", "3rd grade", "4th grade", "5th grade", 
    #   "6th grade", "7th grade", "8th grade", "9th grade", "10th grade", "11th grade", 
    #   "12th grade, no diploma", "High school graduate", "GED or equivalent", 
    #   "Less than 1 year of college credit/post-secondary education", "One year or more of college credit, no degree", 
    #   "Associate degree: Occupational, Technical, or Vocational", "Associate degree: Academic Program", 
    #   "Bachelor's degree (e.g., BA, AB, BS, BBA)", "Master's degree (e.g., MA, MS, MEng, MEd, MBA)", 
    #   "Professional School degree (e.g., MD, DDS, DVM, JD)", "Doctoral degree (e.g., PhD, EdD)")),
    demo_prnt_marital_v2_bl = factor(demo_prnt_marital_v2_bl, levels = c(1:6), labels = c(
      "Married", "Widowed", "Divorced", "Separated", "Never married", "Living with partner"))
  )

# Apply labels to variables
colnames(merged_data)[-1] <- table1_labels

# Create Table 1 using the table1 package
# caption <- "Table 1. Summary statistics of outcomes Internalizing Problems (Raw) and Externalizing Problems (Raw) and covariates for study sample."
table1_object <- table1(~ . , data = merged_data[, -1]) # , caption = caption)

# Save the table1 object as an HTML file
html_output_path <- file.path(figures_dir, paste0(data_processing_date, "_", data_processor, "_table1.html"))
save_html(table1_object, file = html_output_path)

# Function to save table1 object to HTML
save_html <- function(table1_object, file) {
  html <- as.character(table1_object)
  write(html, file)
}

# Date of the input - the subject ids included in the study sample
input_date <- "2024-06-06"
sample_key_path <- file.path("data", paste0(input_date, "_AN_sample_key.csv"))
sample_key <- read.csv(sample_key_path) %>% pull("src_subject_id")

## -- User Arguments

# Set the path for raw data files
data_dir <- '/Users/aidanneher/Library/CloudStorage/Box-Box/ABCD Tabulated Data/5.1/core'
# Location of desired output directory - if NULL, will output into working directory
out_dir <- '/Users/aidanneher/Documents/GitHub/abcd_multiview/data'
# Date you used in output name - if NULL, will use output from Sys.Date() (current date)
out_date <- Sys.Date()
# Initials or other string you want in output naming - no NULL option here
out_initials <- 'AN'
# Where to put plots
figures_dir <- "figures"

## -- Make Design Matrices

# Extract clustering information
path <- file.path(data_dir, "abcd-general", "abcd_y_lt.csv")
cluster_data <- read.csv(path) %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(src_subject_id, rel_family_id, site_id_l) %>%
  arrange(site_id_l, rel_family_id) %>%
  # Filter to those we are using for our analysis
  filter(src_subject_id %in% sample_key)

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

# Convert matrices to data frames
Z_site_df <- as.data.frame(Z_site)
Z_family_df <- as.data.frame(Z_family)
Z_family_to_site_df <- as.data.frame(Z_family_to_site)

# Add the first column as specified
Z_site_df$src_subject_id <- rownames(Z_site_df)
Z_family_df$src_subject_id <- rownames(Z_family_df)
Z_family_to_site_df$rel_family_id <- rownames(Z_family_to_site_df)

# Reorder columns to have the new column as the first column
Z_site_df <- Z_site_df[, c(ncol(Z_site_df), 1:(ncol(Z_site_df)-1))]
Z_family_df <- Z_family_df[, c(ncol(Z_family_df), 1:(ncol(Z_family_df)-1))]
Z_family_to_site_df <- Z_family_to_site_df[, c(ncol(Z_family_to_site_df), 1:(ncol(Z_family_to_site_df)-1))]

# Write data frames to CSV files with today's date prefix
write.csv(Z_site_df, paste0("data/", out_date, "_Z_site.csv"), row.names = FALSE)
write.csv(Z_family_df, paste0("data/", out_date, "_Z_family.csv"), row.names = FALSE)
write.csv(Z_family_to_site_df, paste0("data/", out_date, "_Z_family_to_site.csv"), row.names = FALSE)

# Append the data frames to the existing complete_data_list
complete_data_list$Z_site <- Z_site_df
complete_data_list$Z_family <- Z_family_df

## -- Summarize Sample's Observational Clustering

# Calculate the number of families in each site
families_per_site <- cluster_data %>%
  group_by(site_id_l) %>%
  summarize(n_families = n_distinct(rel_family_id)) %>%
  arrange(desc(n_families))

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

# Calculate non-outlier statistics
non_outliers_families_per_site <- families_per_site %>%
  filter(!n_families %in% outliers_families_per_site$n_families) %>%
  summarize(
    min = min(n_families),
    median = median(n_families),
    max = max(n_families)
  )

non_outliers_individuals_per_site <- individuals_per_site %>%
  filter(!n_individuals %in% outliers_individuals_per_site$n_individuals) %>%
  summarize(
    min = min(n_individuals),
    median = median(n_individuals),
    max = max(n_individuals)
  )

# Create boxplots
plot_families_per_site <- ggplot(families_per_site, aes(x = 0, y = n_families)) +
  geom_boxplot(fill = "blue", color = "black", alpha = 0.7) +
  geom_text(data = outliers_families_per_site, aes(x = 0, label = n_families), 
            position = position_jitter(width = 0.2, height = 0), 
            hjust = -0.3, color = "red") +
  geom_text(data = non_outliers_families_per_site, aes(x = 0, y = min, label = paste("Min:", min)),
            vjust = -1.5, color = "black") +
  geom_text(data = non_outliers_families_per_site, aes(x = 0, y = median, label = paste("Median:", median)),
            vjust = -1.5, color = "black") +
  geom_text(data = non_outliers_families_per_site, aes(x = 0, y = max, label = paste("Max:", max)),
            vjust = -1.5, color = "black") +
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
  geom_text(data = non_outliers_individuals_per_site, aes(x = 0, y = min, label = paste("Min:", min)),
            vjust = -1.5, color = "black") +
  geom_text(data = non_outliers_individuals_per_site, aes(x = 0, y = median, label = paste("Median:", median)),
            vjust = -1.5, color = "black") +
  geom_text(data = non_outliers_individuals_per_site, aes(x = 0, y = max, label = paste("Max:", max)),
            vjust = -1.5, color = "black") +
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
save_plot_with_date_initials(grid_plots, out_initials = "AN", out_dir = figures_dir)

# Create the table and summaries
individuals_per_family_table <- table(individuals_per_family$n_individuals)
individuals_per_site_summary <- summary(individuals_per_site$n_individuals)
families_per_site_summary <- summary(families_per_site$n_families)

# Convert them to data frames
individuals_per_family_df <- as.data.frame(individuals_per_family_table)
individuals_per_site_df <- as.data.frame(as.list(individuals_per_site_summary))
families_per_site_df <- as.data.frame(as.list(families_per_site_summary))

# Write data frames to CSV files with today's date prefix
write.csv(individuals_per_family_df, paste0("data/", out_date, "_individuals_per_family_table.csv"), row.names = FALSE)
write.csv(individuals_per_site_df, paste0("data/", out_date, "_individuals_per_site_summary.csv"), row.names = FALSE)
write.csv(families_per_site_df, paste0("data/", out_date, "_families_per_site_summary.csv"), row.names = FALSE)

# Split data into 80:20 train:test family-wise split stratified by study site

# Date of the input - the subject ids included in the study sample
input_date <- "2024-06-06"
sample_key_path <- file.path("data", paste0(input_date, "_AN_sample_key.csv"))
sample_key <- read.csv(sample_key_path) %>% pull("src_subject_id")
# Set the path for raw data files
data_dir <- '/Users/aidanneher/Library/CloudStorage/Box-Box/ABCD Tabulated Data/5.1/core'
# Location of desired output directory - if NULL, will output into working directory
out_dir <- '/Users/aidanneher/Documents/GitHub/abcd_multiview/data'
# Date you used in output name - if NULL, will use output from Sys.Date() (current date)
out_date <- Sys.Date()
# Initials or other string you want in output naming - no NULL option here
out_initials <- 'AN'
# Extract clustering information
path <- file.path(data_dir, "abcd-general", "abcd_y_lt.csv")
cluster_data <- read.csv(path) %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select(src_subject_id, rel_family_id, site_id_l) %>%
  arrange(site_id_l, rel_family_id) %>%
  # Filter to those we are using for our analysis
  filter(src_subject_id %in% sample_key)

# Load data
# cluster_data <- read.csv('path_to_cluster_data.csv')

library(tidyverse)

# Load data
# cluster_data <- read.csv('path_to_cluster_data.csv')

# Function to create train/test split for each site
split_train_test <- function(data, train_frac = 0.8) {
  set.seed(123) # For reproducibility
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

# Print summaries
print(train_data_summary)
print(test_data_summary)

# Check overall n's and split proportion
n_train <- train_data_summary$n_train %>% sum
n_test <- test_data_summary$n_test %>% sum
cat("n_train:", n_train)
cat("n_test:", n_test)
cat("n_train + n_test:", n_train+n_test)

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

# Save the split lists as RDS objects
out_date <- Sys.Date()
saveRDS(train_list, paste0("data/", out_date, "_train_list.rds"))
saveRDS(test_list, paste0("data/", out_date, "_test_list.rds"))

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

## -- Load Datasets

load_sMRI_data <- function(data_dir) {
  
  ct_path <- file.path(data_dir, "imaging", "mri_y_smr_thk_dst.csv")
  ct_data <- read.csv(ct_path) %>%
    # Remove summary statistic variables: 
    # These represent mean CT in left hemisphere, right hemisphere, & brain overall respectively
    select(-mrisdp_149, -mrisdp_150, -mrisdp_151)
  
  sa_path <- file.path(data_dir, "imaging", "mri_y_smr_area_dst.csv")
  sa_data <- read.csv(sa_path) %>%
    # Remove summary statistic variables: 
    # These represent mean SA in left hemisphere, right hemisphere, & brain overall respectively
    select(-mrisdp_451, -mrisdp_452, -mrisdp_453)
  
  # Filter to baseline
  ct_data <- ct_data %>%
    filter(eventname == "baseline_year_1_arm_1") %>%
    select(-eventname)
  
  sa_data <- sa_data %>%
    filter(eventname == "baseline_year_1_arm_1") %>%
    select(-eventname)
  
  cat("Total n sMRI Scans Available at Baseline (CT):", nrow(ct_data), "\n")
  cat("Total n sMRI Scans Available at Baseline (SA):", nrow(sa_data), "\n")
  
  # Check for missing values in each column for ct_data
  missingness_summary_ct <- ct_data %>%
    summarize_all(~ sum(is.na(.)))
  
  # Print the summary of missing values for ct_data
  cat("The following variables in CT data have some missingness:",
      (names(missingness_summary_ct)[missingness_summary_ct != 0]), "\n")
  
  cat("n Missing for the variable with the most missingness in CT data:",
      missingness_summary_ct %>% max, "\n")
  
  # Check for missing values in each column for sa_data
  missingness_summary_sa <- sa_data %>%
    summarize_all(~ sum(is.na(.)))
  
  # Print the summary of missing values for sa_data
  cat("The following variables in SA data have some missingness:",
      (names(missingness_summary_sa)[missingness_summary_sa != 0]), "\n")
  
  cat("n Missing for the variable with the most missingness in SA data:",
      missingness_summary_sa %>% max, "\n")
  
  return(list(ct_data = ct_data, sa_data = sa_data))
}

# Example usage:
# data_list <- load_sMRI_view("path_to_data_directory")
# ct_data <- data_list$ct_data
# sa_data <- data_list$sa_data

load_fMRI_view <- function(data_dir) {

  # Load Data
  table_names <- c("mri_y_rsfmr_cor_gp_gp") # Correlations [Gordon Network]
  
  table_dirs <- paste0("imaging/", table_names, ".csv")
  file_paths <- file.path(data_dir, table_dirs)
  read_filter_data <- function(file_path) {
    data <- read.csv(file_path) %>% 
      filter(eventname == "baseline_year_1_arm_1") %>%
      select(-eventname)
  }
  
  data_list <- lapply(file_paths, read_filter_data)
  
  # Calculate Summary Statistics
  cat("Table Names:", table_names, "\n")
  n_features_per_table <- data_list %>% sapply(function(df) ncol(df)-1)
  names(n_features_per_table) <- table_names
  cat("n rsfMRI Variables Available per Table:", n_features_per_table, "\n")
  n_rows_with_missing <- sapply(data_list, function(df) sum(apply(df, 1, function(row) any(is.na(row)))))
  names(n_rows_with_missing) <- table_names
  cat("n Observations with Any Missing Data per Table:", n_rows_with_missing)
  
  merged_df <- reduce(data_list, function(x, y) inner_join(x, y, by = "src_subject_id"))
  
  return(merged_df)

}

smri_data <- load_sMRI_data(data_dir)

attach(smri_data)

fmri_data <- load_fMRI_view(data_dir)

# Let's checkout the fMRI data we have loaded
# Note, pairwise correlations include self-correlations
# These self-correlations do not necessarily equal 1.  
# From Mark: "The reason that the self-correlations are not equal to 1 
# is that it's a within-network correlation. 
# Essentially, each network is made of a number of parcels. 
# They then computed all pairwise correlations between the parcels 
# within a network to get all these values. So for, e.g., auditory to 
# auditory, there are still some number of parcels there, and so their 
# correlations are not exactly equal to 1 (though they will likely be large). 
# See, for example, Figure 2 here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9589037/."

# Create a lookup table for the Gordon Network abbreviations
network_lookup <- data.frame(
  Abbreviation = c("ad", "cgc", "ca", "dt", "dla", 
                   "fo", "n", "rspltp", "sa", "smh", 
                   "smm", "vta", "vs"),
  FullName = c("Auditory Network", "Cingulo-Opercular Network", 
               "Cingulo-Parietal Network", "Default Mode Network", 
               "Dorsal Attention Network", "Fronto-Parietal Network", 
               "None Network", "Retrosplenial Temporal Network", 
               "Salience Network", "Sensorimotor Hand Network", 
               "Sensorimotor Mouth Network", "Ventral Attention Network", 
               "Visual Network")
)

# Remove the "rsfmri_c_" and "ngd_" prefixes from column names
colnames(fmri_data) <- sub("^rsfmri_c_ngd_", "", colnames(fmri_data))
colnames(fmri_data) <- gsub("_ngd_", "_", colnames(fmri_data))
# Reshape the data frame to long format
fmri_long <- fmri_data %>%
  pivot_longer(
    cols = -src_subject_id,  # Exclude the subject ID column
    names_to = c("Network1", "Network2"),
    names_sep = "_",
    values_to = "Correlation"
  )

# Print the reshaped data frame
print(fmri_long)

# Calculate the average correlation for each network pair
average_correlations <- fmri_long %>%
  group_by(Network1, Network2) %>%
  summarise(AverageCorrelation = mean(Correlation, na.rm = TRUE)) %>%
  ungroup()

# Create a matrix from the average correlations
cor_matrix <- average_correlations %>%
  pivot_wider(names_from = Network2, values_from = AverageCorrelation) %>%
  column_to_rownames(var = "Network1") %>%
  as.matrix()

# Plot the heatmap
ggplot(reshape2::melt(cor_matrix), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  coord_fixed() +
  labs(x = "Network 1", y = "Network 2", title = "Average Correlation Matrix Heatmap")

## -- Write Imaging Views

# Use current date if out_date is NULL
if (is.null(out_date)) { out_date <- Sys.Date() }

# Write ct_data to CSV
ct_file_name <- sprintf("%s_%s_CT_sMRI_view.csv", out_date, out_initials)
ct_output_path <- ifelse(is.null(out_dir), ct_file_name, file.path(out_dir, ct_file_name))
write_csv(ct_data, ct_output_path)
print(paste("CT sMRI view written to:", ct_output_path))

# Write sa_data to CSV
sa_file_name <- sprintf("%s_%s_SA_sMRI_view.csv", out_date, out_initials)
sa_output_path <- ifelse(is.null(out_dir), sa_file_name, file.path(out_dir, sa_file_name))
write_csv(sa_data, sa_output_path)
print(paste("SA sMRI view written to:", sa_output_path))

# Construct the filename
file_name <- sprintf("%s_%s_fMRI_view.csv", out_date, out_initials)
# Define the output path
output_path <- ifelse(is.null(out_dir), file_name, file.path(out_dir, file_name))
# Write the fMRI_data to CSV file
write_csv(fmri_data, output_path)
# Print the output path for verification
print(paste("fMRI view written to:", output_path))


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

load_sMRI_view <- function(data_dir) {
  
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
  
  # Merge sMRI data, filter to baseline, & remove eventname since cross-sectional
  smri_data <- merge(ct_data, sa_data) %>%
    filter(eventname == "baseline_year_1_arm_1") %>%
    select(-eventname)
  cat("Total n sMRI Scans Available at Baseline:", nrow(smri_data))
  
  return(smri_data)
  
}

# TODO load_FC_view <- function(data_dir) {

  # TODO load data

  # TODO filter to baseline
  # filtered_data <- fc_data %>%
  #   filter(eventname == "baseline_year_1_arm_1")
  # cat("n fMRI Scans Available at Baseline Pre-QC Filtering:", nrow(filtered_data))

  # TODO perform qc filtering - I'm pretty sure this is only needed for fMRI data
  # qc_data <- filtered_data %>%
  #   filter()

  # cat("n fMRI Scans Available at Baseline Post-QC Filtering:", nrow(qc_data))

# }

smri_data <- load_sMRI_view(data_dir)

# Doublecheck missingness
# Check for missing values in each column
missingness_summary <- smri_data %>%
  summarize_all(~ sum(is.na(.)))

# Print the summary of missing values
cat("The following variables have some missingness:",
    (names(missingness_summary)[missingness_summary!=0]))

cat("n Missing for the variable with the most missingness:",
    missingness_summary %>% max)

## -- Write Imaging View

# Use current date if out_date is NULL
if (is.null(out_date)) { out_date <- Sys.Date() }

# Construct the filename
file_name <- sprintf("%s_%s_sMRI_view.csv", out_date, out_initials)
# Define the output path
output_path <- ifelse(is.null(out_dir), file_name, file.path(out_dir, file_name))
# Write the smri_data to CSV file
write_csv(smri_data, output_path)
# Print the output path for verification
print(paste("sMRI view written to:", output_path))

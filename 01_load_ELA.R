library(tidyverse)
library(readxl)

## -- User Arguments

# Set the path for raw data files
data_dir <- '/Users/aidanneher/Library/CloudStorage/Box-Box/ABCD Tabulated Data/5.1/core'
# Location of desired output directory - if NULL, will output into working directory
out_dir <- '/Users/aidanneher/Documents/GitHub/abcd_multiview/data'
# Date you used in output name - if NULL, will use output from Sys.Date() (current date)
out_date <- NULL
# Initials or other string you want in output naming - no NULL option here
out_initials <- 'AN'

## -- Construct Output Path

# Use current date if out_date is NULL
if (is.null(out_date)) { out_date <- Sys.Date() }
# Construct the filename
file_name <- sprintf("view_ela_%s_%s.csv", out_initials, out_date)
# Define the output path
output_path <- ifelse(is.null(out_dir), file_name, file.path(out_dir, file_name))

## -- Construct ELA View

# Loads the ELA variables considered by Orendain et. al. 2023
# Dependency: "orendain_vars_from_supplement.xlsx"
load_orendain_ela <- function(data_dir, supp_features=TRUE) {
  
  # Physical and sexual violence
  mh_p_ksads_path <- file.path(data_dir, "mental-health/mh_p_ksads_ptsd.csv")
  mh_p_ksads <- read_csv(mh_p_ksads_path) %>%
    select(src_subject_id, eventname, ksads_ptsd_raw_761_p,
           ksads_ptsd_raw_762_p, ksads_ptsd_raw_763_p,
           ksads_ptsd_raw_767_p, ksads_ptsd_raw_768_p,
           ksads_ptsd_raw_769_p, ksads_ptsd_raw_760_p,
           ksads_ptsd_raw_764_p, ksads_ptsd_raw_765_p)
  
  # Parent psychopathology
  mh_p_fhx_path <- file.path(data_dir, "mental-health/mh_p_fhx.csv")
  mh_p_fhx <- read_csv(mh_p_fhx_path) %>%
    select(src_subject_id, eventname, famhx_4_p,
           fam_history_5_yes_no, fam_history_6_yes_no,
           fam_history_7_yes_no, fam_history_8_yes_no,
           fam_history_11_yes_no, fam_history_12_yes_no,
           fam_history_13_yes_no)
  
  # Neighborhood Threat
  ce_p_nsc_path <- file.path(data_dir, "culture-environment/ce_p_nsc.csv")
  ce_p_nsc <- read_csv(ce_p_nsc_path) %>%
    select(src_subject_id, eventname, neighborhood3r_p,
           neighborhood2r_p)
  
  # Prenatal Substance Exposure
  ph_p_dhx_path <- file.path(data_dir, "physical-health/ph_p_dhx.csv")
  ph_p_dhx <- read_csv(ph_p_dhx_path) %>%
    select(src_subject_id, eventname,
           devhx_9_tobacco, devhx_9_alcohol,
           devhx_9_marijuana, devhx_9_coc_crack,
           devhx_9_her_morph, devhx_9_oxycont)
  
  # Scarcity
  abcd_p_demo_path <- file.path(data_dir, "abcd-general/abcd_p_demo.csv")
  abcd_p_demo <- read_csv(abcd_p_demo_path) %>%
    select(src_subject_id, eventname, 
           demo_fam_exp1_v2, demo_fam_exp5_v2)
  
  # Household Dysfunction
  ce_y_fes_path <- file.path(data_dir, "culture-environment/ce_y_fes.csv")
  ce_y_fes <- read_csv(ce_y_fes_path) %>%
    select(src_subject_id, eventname, fes_youth_q6, 
           fes_youth_q1, fes_youth_q5)
  
  # Merging all data frames
  merged_data <- list(mh_p_ksads, mh_p_fhx, ce_p_nsc, ph_p_dhx, abcd_p_demo, ce_y_fes) %>%
    reduce(full_join, by = c("src_subject_id", "eventname"))
  
  load_supp_ela <- function() {
    
    # Read the Excel file to get the supp_lookup table
    supp_lookup <- read_xlsx("orendain_vars_from_supplement.xlsx")
    
    # Get the unique file paths
    supp_paths <- file.path(data_dir,
                            unique(supp_lookup$table_file_path))
    
    # Function to read CSV file and select the specified columns
    read_and_select <- function(file_path, var_names) {
      # Read the CSV file
      data <- read_csv(paste0(file_path, ".csv"))
      # Select the specified columns
      selected_data <- data %>% select(all_of(var_names))
      return(selected_data)
    }
    
    # Apply the function to each file path with corresponding var_names
    selected_tables <- lapply(supp_paths, function(path) {
      var_names <- supp_lookup %>%
        filter(table_name == basename(path)) %>%
        pull(var_name)
      var_names <- c("src_subject_id", "eventname", var_names)
      read_and_select(path, var_names)
    })
    
    # Merge the tables in the list by "src_subject_id" and "eventname"
    merged_data <- reduce(selected_tables, function(x, y) {
      full_join(x, y, by = c("src_subject_id", "eventname"))
    })
  }
  
  if (supp_features==TRUE) {
    supp_data <- load_supp_ela()
    merged_data <- merge(merged_data, supp_data)
  }
  
  return(merged_data %>%
           filter(eventname == "baseline_year_1_arm_1"))
}

# Loads the ELA variables considered by Brieant et. al. 2023
load_brieant_ela <- function(data_dir) {
  # Read in family variables
  rel <- read.csv(paste(data_dir, 'abcd-general/abcd_y_lt.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, rel_family_id)
  
  # Read in family substance use summary scores
  fhx <- read.csv(paste(data_dir, 'mental-health/mh_p_fhx.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, famhx_ss_fath_prob_alc_p, famhx_ss_moth_prob_alc_p, famhx_ss_fath_prob_dg_p, famhx_ss_moth_prob_dg_p)
  
  # Read in parent demographics
  pdemo <- read.csv(paste(data_dir, 'abcd-general/abcd_p_demo.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, demo_prim, demo_prnt_marital_v2, demo_prnt_ed_v2, demo_prtnr_ed_v2, demo_comb_income_v2,
           demo_fam_exp1_v2, demo_fam_exp2_v2, demo_fam_exp3_v2, demo_fam_exp4_v2, demo_fam_exp5_v2, 
           demo_fam_exp6_v2, demo_fam_exp7_v2)
  
  # Read in CRPBI
  crpbi <- read.csv(paste(data_dir, 'culture-environment/ce_y_crpbi.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, crpbi_parent1_y, crpbi_caregiver12_y, crpbi_parent2_y, crpbi_caregiver13_y,
           crpbi_parent3_y, crpbi_caregiver14_y, crpbi_parent4_y, crpbi_caregiver15_y, crpbi_parent5_y, 
           crpbi_caregiver16_y)
  
  # Read in parent report family environment scale
  fes02 <- read.csv(paste(data_dir, 'culture-environment/ce_p_fes.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, fam_enviro1_p, fam_enviro2r_p, fam_enviro3_p, fam_enviro4r_p, fam_enviro5_p,
           fam_enviro6_p, fam_enviro7r_p, fam_enviro8_p, fam_enviro9r_p)
  
  # Read in youth report family environment scale
  fes01 <- read.csv(paste(data_dir, 'culture-environment/ce_y_fes.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, fes_youth_q1, fes_youth_q2, fes_youth_q3, fes_youth_q4, fes_youth_q5, fes_youth_q6,
           fes_youth_q7, fes_youth_q8, fes_youth_q9)
  
  # Read in ksads trauma, parent interview
  ptsd <- read.csv(paste(data_dir, 'mental-health/mh_p_ksads_ptsd.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, ksads_ptsd_raw_754_p, ksads_ptsd_raw_755_p, ksads_ptsd_raw_756_p, ksads_ptsd_raw_757_p,
           ksads_ptsd_raw_758_p, ksads_ptsd_raw_759_p, ksads_ptsd_raw_760_p, ksads_ptsd_raw_761_p,
           ksads_ptsd_raw_762_p, ksads_ptsd_raw_763_p, ksads_ptsd_raw_764_p, ksads_ptsd_raw_765_p,
           ksads_ptsd_raw_766_p, ksads_ptsd_raw_767_p, ksads_ptsd_raw_768_p, ksads_ptsd_raw_769_p,
           ksads_ptsd_raw_770_p)
  
  # Read in parental monitoring
  pmq <- read.csv(paste(data_dir, 'culture-environment/ce_y_pm.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, parent_monitor_q1_y, parent_monitor_q2_y, parent_monitor_q3_y, parent_monitor_q4_y,
           parent_monitor_q5_y)
  
  # Read in neighborhood safety and crime, parents
  pnscss <- read.csv(paste(data_dir, 'culture-environment/ce_p_nsc.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, nsc_p_ss_mean_3_items)
  
  # Read in neighborhood safety and crime, youth
  ynsc <- read.csv(paste(data_dir, 'culture-environment/ce_y_nsc.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, neighborhood_crime_y)
  
  # Read in ASR (parent psychopathology)
  asr <- read.csv(paste(data_dir, 'mental-health/mh_p_asr.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, asr_scr_anxdisord_r, asr_scr_somaticpr_r, asr_scr_depress_r, asr_scr_avoidant_r,
           asr_scr_adhd_r, asr_scr_antisocial_r, asr_scr_inattention_r, asr_scr_hyperactive_r)
  
  # Read in ADI data file
  ADI <- read.csv(paste(data_dir, 'linked-external-data/led_l_adi.csv', sep='/')) %>%
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(src_subject_id, reshist_addr1_adi_wsum)
  
  # Merge all data frames
  merged_data <- full_join(fhx, asr) %>%
    full_join(pdemo) %>%
    full_join(crpbi) %>%
    full_join(fes01) %>%
    full_join(fes02) %>%
    full_join(ptsd) %>%
    full_join(pmq) %>%
    full_join(pnscss) %>%
    full_join(ynsc) %>%
    full_join(rel) %>%
    full_join(ADI)
  
  return(merged_data)
}

orendain_data <- load_orendain_ela(data_dir)
print(paste("n Columns of Orendain Data:", ncol(orendain_data)))
brieant_data <- load_brieant_ela(data_dir)
print(paste("n Columns of Brieant Data:", ncol(brieant_data)))
combined_data <- merge(orendain_data, brieant_data)
print(paste("n Columns of Combined Data:", ncol(combined_data)))

# Calculate the proportion of non-zero values per column, excluding specific columns
non_zero_data <- combined_data %>%
  select(-src_subject_id, -eventname) %>%  # Exclude these columns from the analysis
  summarise(across(everything(), ~mean(. != 0, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "prop_non_zero")
# Filter the tibble to include only those with prop_non_zero less than 0.05%
vars_to_include <- non_zero_data %>%
  filter(prop_non_zero >= 0.0005) %>%
  pull(variable)
endorsed_data <- combined_data %>%
  select(all_of(c("src_subject_id", "eventname", vars_to_include)))
print(paste("n Columns of Well-Endorsed Data:", ncol(endorsed_data)))

# Remove redundant columns
no_covariate_data <- endorsed_data %>%
  # Remove covariates, columns that start with "demo_"
  select(-starts_with("demo_")) %>% 
  # Clustering var will be represented in design matrices
  select(-rel_family_id) %>% 
  # nsc_p_ss_mean_3_items: Neighborhood Safety Protocol: Mean of Parent Report, 
  # (neighborhood1r_p + neighborhood2r_p + neighborhood3r_p)/3;
  # We remove the responses to the individual represented in this summary statistics
  select(-neighborhood2r_p, -neighborhood3r_p) %>%
  select(-eventname) # We are only taking data from baseline
print(paste("n Columns without Covariates:", ncol(no_covariate_data)))

# Recode values in all columns except 'src_subject_id' and those that start with "asr_scr_"
recoded_data <- no_covariate_data %>%
  mutate(across(
    !matches("^asr_scr_|^src_subject_id$"), 
    ~ na_if(., 999)
  )) %>%
  mutate(across(
    !matches("^asr_scr_|^src_subject_id$"), 
    ~ na_if(., 7)
  )) %>%
  mutate(across(
    !matches("^asr_scr_|^src_subject_id$"), 
    ~ na_if(., 777)
  ))

## -- Write ELA View

# Write the recoded_data to the CSV file
write_csv(recoded_data, output_path)

# Print the output path for verification
print(paste("ELA view written to:", output_path))

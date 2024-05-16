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
# For Ellery's code, we must specify what time points we want. 
# Put the string provided for all desired timepoints into a vector e.g.  c('timepoint1_name','timepoint2_name') etc.
# We stick to baseline for this analysis
timepoint_list <- c("baseline_year_1_arm_1") 


## -- Load Data Files

load_datasets <- function(data_dir) {
  # Standard files for most studies
  demog <- read_csv(file.path(data_dir, "abcd-general/abcd_p_demo.csv"))
  demog_y <- read_csv(file.path(data_dir, "mental-health/mh_y_ksads_bg.csv"))
  puberty <- read_csv(file.path(data_dir, "physical-health/ph_p_pds.csv"))
  study_covars <- read_csv(file.path(data_dir, "abcd-general/abcd_y_lt.csv"))
  sib_twin <- read_csv(file.path(data_dir, "genetics/gen_y_pihat.csv"))
  
  # MRI standard files
  mri <- read_csv(file.path(data_dir, "imaging/mri_y_adm_info.csv"))
  qc <- read_csv(file.path(data_dir, "imaging/mri_y_qc_incl.csv"))
  scan_qtns <- read_csv(file.path(data_dir, "imaging/mri_y_adm_qtn.csv"))
  
  # Standard files for clinical data
  cbcl <- read_csv(file.path(data_dir, "mental-health/mh_p_cbcl.csv"))
  bpm_y <- read_csv(file.path(data_dir, "mental-health/mh_y_bpm.csv"))
  ksad_p <- read_csv(file.path(data_dir, "mental-health/mh_p_ksads_ss.csv"))
  ksad_y <- read_csv(file.path(data_dir, "mental-health/mh_y_ksads_ss.csv"))
  
  # Combine all datasets into a list (if needed for further processing or return)
  datasets <- list(demog = demog, demog_y = demog_y, puberty = puberty,
                   study_covars = study_covars, sib_twin = sib_twin, mri = mri,
                   qc = qc, scan_qtns = scan_qtns, cbcl = cbcl, bpm_y = bpm_y,
                   ksad_p = ksad_p, ksad_y = ksad_y)
  
  return(datasets)
}

datasets <- load_datasets(data_dir)

attach(datasets)

# KSAD INITIAL CLEANING
# Why? During the data collection, a new version of the KSAD was used (KSAD2), so each KSAD variable has a corresponding KSAD2 variable which contains the data since the new version was adopted. The code below combines these variables into one.

ksad_new_df <- ksad_p %>% # initialize new data frame
  select(src_subject_id, eventname)

clean_ksads <- function(ksad_new_df, ksad, ksad2, p_or_t){ 
  # goal: combine two ksad variables into one new var
  # ksad_new_df = a new data frame created above to hold the new vars, 
  # ksad = 1st ksad variable (old ksad), 
  # ksad2 = ksad2 var (same content/question as ksad1), 
  # p_or_t = whether the variable is a p or t ksad 
  # output: a dataframe containing the new variable
  assign("ksad", paste0(ksad, sep = "_", p_or_t)) # create ksad variable from inputs and assign it the name ksad
  assign("ksad2", paste0(ksad2, sep = "_", p_or_t)) # create ksad2 variable from inputs and assign it the name ksad2
  assign("ksad_new", paste0(ksad, "_new")) # create new variable name from inputs and assign it the name ksad_new
  assign("ksad_df", paste("ksad", if_else(p_or_t == "p", "p", "y"), sep = "_")) # create data frame name from inputs and name it ksad_df
  my_cols <- c("src_subject_id", "eventname", ksad, ksad2) # create vector of necessary column names
  intermediary <- get(ksad_df) %>%
    select(all_of(my_cols)) %>% # select cols (need to do it this way because I'm referring to the cols as strings)
    mutate(!!ksad_new := if_else(is.na(get(ksad)) == T, get(ksad2), get(ksad))) # create new var
  ksad_new_df <- right_join(intermediary, ksad_new_df, by = c("src_subject_id", "eventname")) # join new var to existing data set
  ksad_new_df <- ksad_new_df %>%
    select(src_subject_id, eventname, ends_with("_new")) # select only necessary vars
  return(ksad_new_df)
}


# call the function with each pair of ksad and ksad2 variables (make sure to assign the output to the same name so you create one data set with all the new variables)
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_946", "ksads2_23_906", "p") # suicidal ideation
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_947", "ksads2_23_907", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_948", "ksads2_23_908", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_949", "ksads2_23_909", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_950", "ksads2_23_910", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_951", "ksads2_23_911", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_957", "ksads2_23_917", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_958", "ksads2_23_918", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_959", "ksads2_23_919", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_960", "ksads2_23_920", "p") 
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_961", "ksads2_23_921", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_954", "ksads2_23_914", "p") # suicidal attempts, prep
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_965", "ksads2_23_925", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_962", "ksads2_23_922", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_966", "ksads2_23_926", "p") # NO suicidal ideation/behaviors
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_955", "ksads2_23_915", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_143", "ksads2_23_134", "p") # non-suicidal self injury
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_144", "ksads2_23_135", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_956", "ksads2_23_916", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_945", "ksads2_23_905", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_840",	"ksads2_1_790", "p") # depression
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_841",	"ksads2_1_791", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_842",  "ksads2_1_792", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_843",	"ksads2_1_793", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_844",	"ksads2_1_794", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_845",	"ksads2_1_795", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_846",	"ksads2_1_796", "p")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_847",	"ksads_1_847", "p") # this one does not have a corresponding ksad2 (I put it here so it would still appear in the df)
# repeat function calls with t instead of p
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_946", "ksads2_23_906", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_947", "ksads2_23_907", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_948", "ksads2_23_908", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_949", "ksads2_23_909", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_950", "ksads2_23_910", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_957", "ksads2_23_917", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_958", "ksads2_23_918", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_959", "ksads2_23_919", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_960", "ksads2_23_920", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_961", "ksads2_23_921", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_954", "ksads2_23_914", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_965", "ksads2_23_925", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_966", "ksads2_23_926", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_955", "ksads2_23_915", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_951", "ksads2_23_911", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_962", "ksads2_23_922", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_143", "ksads2_23_134", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_144", "ksads2_23_135", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_956", "ksads2_23_916", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_945", "ksads2_23_905", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_963", "ksads2_23_923", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_964", "ksads2_23_924", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_953", "ksads2_23_913", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_23_952", "ksads2_23_912", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_840",	"ksads2_1_790", "t") 
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_841",	"ksads2_1_791", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_842",  "ksads2_1_792", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_843",	"ksads2_1_793", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_844",	"ksads2_1_794", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_845",	"ksads2_1_795", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_846",	"ksads2_1_796", "t")
ksad_new_df <- clean_ksads(ksad_new_df, "ksads_1_847",	"ksads_1_847", "t") # this one does not have a corresponding ksad2 (I put it here so it would still appear in the df)

# Demographic Cleaning Pre-Merge
# COMBINE BASELINE VARIABLES WITH LONGITUDINAL VARIABLES
# how it combines: if the observation occurred in year1/arm1 (the baseline) the combination variable takes the value of the baseline variable, if not the combo variable takes the value of the longitudinal variable
var_vect <- c("demo_brthdat_v2", "demo_gender_id_v2", "demo_prnt_ed_v2", "demo_prtnr_ed_v2" , "demo_prnt_marital_v2", "demo_comb_income_v2", "demo_roster_v2", "demo_fam_exp1_v2", "demo_fam_exp2_v2", "demo_fam_exp3_v2", "demo_fam_exp4_v2", "demo_fam_exp5_v2", "demo_fam_exp6_v2", "demo_fam_exp7_v2") # input any baseline variable name you wish to combine with long. var in this vector
combine.vars <- function(var){ # takes baseline variable name (a string) and makes a new variable using the corresponding long. var, outputs demog dataset with new combo variables (combo variables will be named the baseline variable name with "_comb" at the end)
  assign("var_comb", paste0(var, "_comb"))
  assign("var_l", paste0(var, "_l"))
  demog <- demog %>%
    mutate(!!var_comb := if_else(eventname == "baseline_year_1_arm_1", get(var), get(var_l))) # if observation occurred in year 1, arm 1 use bl var if not use long var
  return(demog)
}

for (i in 1:length(var_vect)) { # loops through vector of all variables to be combined and applies the combine.vars function
  demog <- combine.vars(var_vect[i])}

## -- Select Relevant Variables
demog_bl <- demog[demog$eventname == "baseline_year_1_arm_1",]

demog_bl <- demog_bl %>% 
  select(src_subject_id, demo_sex_v2, demo_race_a_p___10,
         demo_race_a_p___11,demo_race_a_p___12, demo_race_a_p___13, demo_race_a_p___14,
         demo_race_a_p___15, demo_race_a_p___16, demo_race_a_p___17, demo_race_a_p___18,
         demo_race_a_p___19, demo_race_a_p___20, demo_race_a_p___21, demo_race_a_p___22,
         demo_race_a_p___23, demo_race_a_p___24, demo_race_a_p___25,demo_race_a_p___77, 
         demo_race_a_p___99, demo_ethn_v2, demo_prnt_marital_v2, demo_prnt_ed_v2, demo_prtnr_ed_v2, demo_comb_income_v2)

demog_bl <- dplyr::rename(demog_bl, demo_prnt_marital_v2_bl = demo_prnt_marital_v2)
demog_bl <- dplyr::rename(demog_bl, demo_comb_income_v2_bl = demo_comb_income_v2)
demog_bl <- dplyr::rename(demog_bl, demo_prnt_ed_v2_bl = demo_prnt_ed_v2)
demog_bl <- dplyr::rename(demog_bl, demo_prtnr_ed_v2_bl = demo_prtnr_ed_v2)

demog <- demog %>%
  select(src_subject_id, eventname, demo_brthdat_v2_comb, demo_gender_id_v2_comb, 
         demo_prnt_ed_v2_comb, demo_prtnr_ed_v2_comb, demo_prnt_marital_v2_comb, demo_comb_income_v2_comb,
         demo_roster_v2_comb, demo_fam_exp1_v2_comb, demo_fam_exp2_v2_comb, demo_fam_exp3_v2_comb, 
         demo_fam_exp4_v2_comb, demo_fam_exp5_v2_comb, demo_fam_exp6_v2_comb, demo_fam_exp7_v2_comb,
         acs_raked_propensity_score)

demog_y <- demog_y %>%
  select(src_subject_id, eventname, kbi_gender, kbi_y_trans_id, kbi_y_sex_orient)

puberty <- puberty %>%
  select(src_subject_id, eventname, pds_p_ss_male_category_2, pds_p_ss_female_category_2)

fam <- study_covars[study_covars$eventname == "baseline_year_1_arm_1",]

fam <- fam %>%
  select(src_subject_id, rel_family_id)

study_covars <- study_covars %>%
  select(src_subject_id, eventname, site_id_l, interview_age)

sib_twin <- sib_twin %>% 
  select(src_subject_id, rel_relationship, rel_group_id)

mri <- mri %>%
  select(src_subject_id, eventname, mri_info_manufacturer)

qc <- qc %>%
  select(src_subject_id, eventname, imgincl_rsfmri_include)

cbcl <- cbcl %>%
  select(src_subject_id, eventname, cbcl_scr_syn_internal_r, cbcl_scr_syn_external_r, cbcl_scr_syn_totprob_r,
         cbcl_scr_dsm5_depress_r, cbcl_scr_dsm5_anxdisord_r, cbcl_scr_dsm5_adhd_r,
         cbcl_scr_syn_internal_t, cbcl_scr_syn_external_t, cbcl_scr_syn_totprob_t,
         cbcl_scr_dsm5_depress_t, cbcl_scr_dsm5_anxdisord_t, cbcl_scr_dsm5_adhd_t)

bpm_y <- bpm_y %>% 
  select(src_subject_id, eventname, bpm_y_scr_attention_r, bpm_y_scr_attention_t, bpm_y_scr_internal_r, 
         bpm_y_scr_internal_t, bpm_y_scr_external_r, bpm_y_scr_external_t, bpm_y_scr_totalprob_r, 
         bpm_y_scr_totalprob_t)

## -- Merge into single file and replace missing data codes with NA (Can change this if you do not want that)
files <- list(demog, demog_y, puberty, study_covars, 
              mri, qc, scan_qtns, cbcl, bpm_y, ksad_new_df)
abcd_data_0 <- files %>% reduce(full_join, by = c("src_subject_id", "eventname"))

files_2 <- list(demog_bl, fam, sib_twin, abcd_data_0)
abcd_data <- files_2 %>% reduce(full_join, by = "src_subject_id")

abcd_data <- abcd_data %>%
  mutate(across(where(is.numeric), ~na_if(.,777))) %>%
  mutate(across(where(is.numeric), ~na_if(.,999))) %>%
  mutate(across(where(is.numeric), ~na_if(.,555))) %>%
  mutate(across(where(is.numeric), ~na_if(.,888))) %>%
  mutate(across(where(is.character), ~na_if(.,"777"))) %>%
  mutate(across(where(is.character), ~na_if(.,"999"))) %>%
  mutate(across(where(is.character), ~na_if(.,"555"))) %>%
  mutate(across(where(is.character), ~na_if(.,"888"))) %>%
  mutate(across(where(is.character), ~na_if(.,""))) 

# Clean and Recode Age
abcd_data <- abcd_data %>%
  mutate(demo_brthdat_v2_comb_clean = if_else(demo_brthdat_v2_comb > 21, 
                                              demo_brthdat_v2_comb/12, demo_brthdat_v2_comb), # convert months to years
         demo_brthdat_v2_comb_clean = if_else(demo_brthdat_v2_comb_clean < 8, NA, demo_brthdat_v2_comb_clean), # younger than 8 --> NA
         demo_brthdat_v2_comb_clean = trunc(demo_brthdat_v2_comb_clean)) %>% # remove decimals 
  mutate(interview_age_b = interview_age / 12) # convert months to years


# Recode Race Variables
abcd_data <- abcd_data %>%
  mutate(demo_ethn_v2 = abs(demo_ethn_v2 - 2), # change "2" to 0 to match other vars 
         White_race = demo_race_a_p___10,
         Black_race = demo_race_a_p___11,
         AIAN_race = if_else(rowSums(select(., num_range("demo_race_a_p___", 12:13)))
                             >=1, 1, 0),
         NHPI_race = if_else(rowSums(select(., num_range("demo_race_a_p___", 14:17)))
                             >=1, 1, 0),
         Asian_race = if_else(rowSums(select(., num_range("demo_race_a_p___", 18:24)))
                              >=1, 1, 0),
         Other_race = demo_race_a_p___25,
         Missing_race = if_else(demo_race_a_p___99 == 1, 1, demo_race_a_p___77), 
         Missing_race = if_else(White_race == 1 | # if participant did not endorse any of these, mark as missing
                                  Black_race == 1 | 
                                  AIAN_race == 1 |
                                  NHPI_race == 1 |
                                  Asian_race == 1 | 
                                  Other_race == 1 |
                                  demo_ethn_v2 == 1, 0, 1),
         Missing_race = if_else(is.na(Missing_race) == T, 1, Missing_race)) %>% # mark all NAs in missing race as 1
  mutate(Indigenous_race = if_else(AIAN_race == 1 | NHPI_race == 1, 1, 0)) 

# Recode Puberty Variable
abcd_data <- abcd_data %>%
  mutate(pubertal_status = if_else(demo_sex_v2 == '1' | demo_sex_v2 == '3',
                                   pds_p_ss_male_category_2, # Set to this if the condition above is TRUE
                                   pds_p_ss_female_category_2)) # Otherwise set to this 

# Recode Parent-Related Demographics
abcd_data <- abcd_data %>%
  mutate(across(c(demo_prnt_ed_v2_comb, demo_prtnr_ed_v2_comb, 
                  demo_prnt_ed_v2_bl, demo_prtnr_ed_v2_bl), ~as.integer(.x))) %>%
  mutate(highest_demo_ed_comb = case_when(
    is.na(demo_prtnr_ed_v2_comb) == T ~ 
      demo_prnt_ed_v2_comb,
    demo_prnt_ed_v2_comb > demo_prtnr_ed_v2_comb ~ 
      demo_prnt_ed_v2_comb,
    demo_prnt_ed_v2_comb <= demo_prtnr_ed_v2_comb ~ 
      demo_prtnr_ed_v2_comb),
    highest_demo_ed_bl = case_when(
      is.na(demo_prtnr_ed_v2_bl) == T ~ demo_prnt_ed_v2_bl,
      demo_prnt_ed_v2_bl > demo_prtnr_ed_v2_bl ~ 
        demo_prnt_ed_v2_bl,
      demo_prnt_ed_v2_bl <= demo_prtnr_ed_v2_bl ~ 
        demo_prtnr_ed_v2_bl))

# Remove raw variables
abcd_data <- abcd_data %>%
  select(-num_range("demo_race_a_p___", 12:25), # race vars
         -num_range("ksads_23_", 946:950, "_t_new"), # SI, SA, NSSI vars
         -num_range("ksads_23_", 952:954, "_t_new"),
         -num_range("ksads_23_", 957:961, "_t_new"),
         -num_range("ksads_23_", 963:965, "_t_new"),
         -num_range("ksads_23_", 946:950, "_p_new"),
         -ksads_23_954_p_new,
         -num_range("ksads_23_", 957:961, "_p_new"), 
         -num_range("ksads_1_", 840:847, "_t_new"), # depression vars
         -num_range("ksads_1_", 840:847, "_p_new"))

# Remove the time points you do not want
abcd_data.selected_time <- abcd_data %>%
  filter(eventname %in% timepoint_list)

## -- Extract Covariates and Outcomes from Ellery's Processing Steps

covariates <- abcd_data.selected_time %>%
  select(src_subject_id, 
         interview_age, 
         ends_with("_race"), # Select all columns that end with "_race"
         demo_comb_income_v2_bl,
         highest_demo_ed_bl,
         demo_prnt_marital_v2_bl)
         # TODO consider including religion as covar with varname demo_relig_v2_l

outcomes <- abcd_data.selected_time %>%
  # Brient et al. 2023 use t scores (Not the raw r scores)
  select("src_subject_id", "cbcl_scr_syn_internal_t", "cbcl_scr_syn_external_t")

## -- Write out a covariates and outcome variables separately

# Creates name based on date and initials/string passed in by user. If no date given, use the current date
if(is.null(out_date)){
  out_date <- Sys.Date()
}
# If no initials/string is given, throw an error that user must input one
if(is.null(out_initials)){
  stop("No input given for 'out_initials'. 
       User must provide string input to create output file name.")
}

# Construct the filename
file_name <- sprintf("%s_%s_covariates.csv", out_date, out_initials)
# Define the output path
output_path <- ifelse(is.null(out_dir), file_name, file.path(out_dir, file_name))
write_csv(covariates, output_path)

# Construct the filename
file_name <- sprintf("%s_%s_outcomes.csv", out_date, out_initials)
# Define the output path
output_path <- ifelse(is.null(out_dir), file_name, file.path(out_dir, file_name))
write_csv(outcomes, output_path)

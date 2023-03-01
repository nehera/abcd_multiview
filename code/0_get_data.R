# Author: Aidan Neher
# Initialized: January 10, 2023
# Description: This script gets ABCD Study data from Box.
library(tidyverse)
get_box_file <- function(file_name, sep = "\t", skip_description = TRUE,
                         box_path = "~/Library/CloudStorage/Box-Box",
                         release_path = "ABCD STUDY/ABCDStudy 4.0 Release/Raw data files") {
  filepath <- paste(box_path, release_path, file_name, sep = "/")
  data <- read.csv(filepath, header = TRUE, sep = sep)
  if (skip_description == TRUE) { data <- data[-1,] }
  return(data)
}
usual_mutate <- function(data) {
  data %>%
    mutate(sex = factor(sex)) %>%
    mutate(interview_age = as.integer(interview_age)) %>%
    mutate(eventname = factor(eventname, 
                              levels = c("baseline_year_1_arm_1", "2_year_follow_up_y_arm_1")))
}
usual_vars = c("subjectkey", "sex", "interview_age", "eventname")

## -- Atlas Data
atlas_data <- ggseg::read_atlas_files(subjects_dir = "/Applications/freesurfer/7.2.0/subjects/", atlas="aparc.stats")

## -- Structural MRI Data
smri_data <- get_box_file("abcd_smrip10201.txt") %>%
  select(c(all_of(usual_vars), starts_with("smri"))) %>%
  select(-c(smri_visitid, starts_with("smri_vol_"), contains("_cf"))) %>% 
  usual_mutate() %>% 
  mutate_at(vars(starts_with("smri_")), as.numeric) %>%
  group_by(subjectkey)
summary_var_index <- str_detect(colnames(smri_data), "mean|total")
smri_data <- smri_data[,!summary_var_index]
# remove sulcal depth
smri_data <- smri_data %>%
  select(-contains("_sulc_"))

## -- Demographic Data
# Each race/ethnicity included as a separate variable where individuals are allowed to be coded as 1 across more than one variable. 
# We have been interpreting these in models as the effect for Black versus non-Black youth, as an example. 
demo_data <- get_box_file("pdem02.txt") %>%
  select(c(usual_vars, demo_comb_income_v2,
           demo_race_a_p___10:demo_race_a_p___25,demo_ethn_v2,
           demo_prnt_marital_v2,demo_prnt_ed_v2,demo_prtnr_ed_v2)) %>%
  usual_mutate() %>%
  mutate_at(vars(starts_with("demo_race")), as.numeric) %>%
  mutate(white = if_else(demo_race_a_p___10 == 1,1,0),
         black = if_else(demo_race_a_p___11 == 1,1,0),
         native = if_else(demo_race_a_p___12 == 1 | demo_race_a_p___13 == 1, 1, 0),
         pacific_islander = if_else(demo_race_a_p___14 == 1 | demo_race_a_p___15 == 1 |
                                          demo_race_a_p___16 == 1 | demo_race_a_p___17 == 1 ,1,0),
         asian = if_else(demo_race_a_p___18 == 1 | demo_race_a_p___19 == 1 |
                                demo_race_a_p___20 == 1 | demo_race_a_p___21 == 1 |
                                demo_race_a_p___22 == 1 | demo_race_a_p___23 == 1 |
                                demo_race_a_p___24 == 1, 1, 0),
         other = if_else(demo_race_a_p___25 == 1, 1, 0)) %>%
  select(-starts_with("demo_race")) 

race_index <- which(colnames(demo_data)=="white"):ncol(demo_data)
race_names <- colnames(demo_data)[race_index]
colnames(demo_data)[race_index] <- paste("race", race_names, sep = "_")

# It is important that you make sure everyone in your model has race/ethnicity coding otherwise they will all be compared to some subset with missing data (e.g., like dummy coding). 
demo_data$race_response_present <- rowSums(demo_data[,race_index], na.rm=TRUE)
demo_data$race_response_present %>% table # 171 non-responses
demo_data <- demo_data %>% 
  filter(race_response_present>0) %>%
  select(-race_response_present)

# It's also important to remove non-answers to other questions
demo_filter_indicator <- demo_data %>%
  select(starts_with("demo_")) %>%
  mutate_all(as.numeric) %>%
  apply(1, function(r) any(r %in% c(777, 999)))
demo_data <- demo_data[!demo_filter_indicator,] %>%
  mutate_at(vars(starts_with("demo_")), as.numeric) 

# We need to take the highest level of parental education
# and format demo vars as factors
demo_prnt_highest_ed <- demo_data %>%
  select(contains("_ed_")) %>%
  mutate(demo_prnt_highest_ed = if_else(is.na(demo_prtnr_ed_v2), demo_prnt_ed_v2, 
                                   if_else(demo_prnt_ed_v2 >= demo_prtnr_ed_v2,
                                           demo_prnt_ed_v2, demo_prtnr_ed_v2))) %>%
  pull(demo_prnt_highest_ed)

demo_data <- demo_data %>%
  mutate(demo_prnt_highest_ed = demo_prnt_highest_ed) %>%
  select(-contains("_ed_")) %>%
  mutate_at(vars(starts_with("demo_")), as.factor) 

# Lastly, we need to get the abcd_site and the family_id to include as random effects
deap_covars <- read.csv('data/DEAP_covariates_release4.csv', header=T, sep=',') %>%
  rename(subjectkey = src_subject_id) %>%
  select(subjectkey, abcd_site, rel_family_id) %>%
  mutate_at(vars(c(abcd_site, rel_family_id)), as.factor) %>%
  group_by(subjectkey) %>%
  distinct() %>%
  filter(row_number(subjectkey) == 1) # TODO: Determine what to do with subs with multiple sites/ families.

demo_data <- left_join(demo_data, deap_covars)

## -- Clinical Outcomes
## Suicidal Ideation Outcome (Binary)
outcome_vars = c("ksads_23_946_t","ksads_23_947_t","ksads_23_948_t","ksads_23_949_t","ksads_23_950_t",
                 "ksads_23_957_t","ksads_23_958_t","ksads_23_959_t","ksads_23_960_t","ksads_23_961_t")
si_data <- get_box_file("abcd_ksad501.txt") %>%
  select(usual_vars, all_of(outcome_vars)) %>%
  usual_mutate() %>%
  mutate_at(outcome_vars, as.numeric) %>%
  mutate(outcome_si = rowSums(select(., starts_with("ksads")))) %>%
  mutate(outcome_si = ifelse(outcome_si >= 1, 1, 0)) %>%
  select(-starts_with("ksads"))
## Internalizing Symptoms Outcome (Continuous Severity Score)
internalizing_data <- get_box_file("abcd_cbcls01.txt") %>%
  select(usual_vars, "cbcl_scr_syn_internal_r") %>%
  usual_mutate() %>%
  rename(outcome_internalizing_score = cbcl_scr_syn_internal_r) %>%
  mutate(outcome_internalizing_score = as.numeric(outcome_internalizing_score))

## -- Understand Views
data <- list(view_demo = demo_data,
             view_smri = smri_data,
             outcome_si = si_data,
             outcome_internalizing = internalizing_data)

lapply(data, dim)

lapply(data, function(x){table(x$eventname)})

tidy_data <- data$view_smri %>%
  left_join(select(data$view_demo, -c("sex", "interview_age", "eventname")), by = "subjectkey") %>% 
  left_join(data$outcome_si) %>%
  left_join(data$outcome_internalizing) %>%
  group_by(subjectkey, eventname)

table(tidy_data$eventname)

summary(tidy_data) 

## -- Save Tidy Data
saveRDS(tidy_data, file = paste0("data/", Sys.Date(), "-tidy_data.RDS"))

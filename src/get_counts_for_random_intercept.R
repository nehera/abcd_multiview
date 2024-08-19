# Author: Aidan Neher
# Initialized: Feb 26, 2024
# Description: 
library(tidyverse)

df <- read.csv("data/abcd_y_lt.csv") 

unique_counts <- df %>%
  group_by(rel_family_id, site_id_l, eventname) %>%
  summarise(unique_subjects = n_distinct(src_subject_id), .groups = 'drop')

# View the result
print(unique_counts)

# Filter to baseline
bl_counts <- unique_counts %>%
  filter(eventname == "baseline_year_1_arm_1")

# Viz results
ggplot(bl_counts, aes(x = unique_subjects)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Unique Subjects per Family/Site at Baseline",
       x = "Number of Unique Subjects",
       y = "Frequency")

bl_counts %>% pull(unique_subjects) %>% tabulate()

# Consider sites
site_counts <- bl_counts %>% 
  ungroup() %>% 
  group_by(site_id_l) %>% 
  summarise(unique_subjects = sum(unique_subjects), .groups = 'drop')

# Viz results
ggplot(site_counts, aes(x = unique_subjects)) +
  geom_boxplot(fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Boxplot of Unique Subjects per Site at Baseline",
       x = "Number of Unique Subjects",
       y = "Frequency")

# Minimum observations/ site at baseline
print(min(site_counts$unique_subjects))


# get_box_file <- function(file_name, sep = "\t", skip_description = TRUE,
#                          box_path = "~/Library/CloudStorage/Box-Box",
#                          release_path = "ABCD Tabulated Data/5.1/core/abcd-general") {
#   filepath <- paste(box_path, release_path, file_name, sep = "/")
#   data <- read.csv(filepath, header = TRUE, sep = sep)
#   if (skip_description == TRUE) { data <- data[-1,] }
#   return(data)
# }
# 
# usual_mutate <- function(data) {
#   data %>%
#     mutate(sex = factor(sex)) %>%
#     mutate(interview_age = as.integer(interview_age)) %>%
#     mutate(eventname = factor(eventname, 
#                               levels = c("baseline_year_1_arm_1", "2_year_follow_up_y_arm_1")))
# }
# 
# usual_vars = c("subjectkey", "sex", "interview_age", "eventname")

## -- Demographic Data
# Each race/ethnicity included as a separate variable where individuals are allowed to be coded as 1 across more than one variable. 
# We have been interpreting these in models as the effect for Black versus non-Black youth, as an example. 
# demo_data <- get_box_file("p_demo.csv")


# %>%
#   select(c(usual_vars, demo_comb_income_v2,
#            demo_race_a_p___10:demo_race_a_p___25,demo_ethn_v2,
#            demo_prnt_marital_v2,demo_prnt_ed_v2,demo_prtnr_ed_v2)) %>%
#   usual_mutate() %>%
#   mutate_at(vars(starts_with("demo_race")), as.numeric) %>%
#   mutate(white = if_else(demo_race_a_p___10 == 1,1,0),
#          black = if_else(demo_race_a_p___11 == 1,1,0),
#          native = if_else(demo_race_a_p___12 == 1 | demo_race_a_p___13 == 1, 1, 0),
#          pacific_islander = if_else(demo_race_a_p___14 == 1 | demo_race_a_p___15 == 1 |
#                                       demo_race_a_p___16 == 1 | demo_race_a_p___17 == 1 ,1,0),
#          asian = if_else(demo_race_a_p___18 == 1 | demo_race_a_p___19 == 1 |
#                            demo_race_a_p___20 == 1 | demo_race_a_p___21 == 1 |
#                            demo_race_a_p___22 == 1 | demo_race_a_p___23 == 1 |
#                            demo_race_a_p___24 == 1, 1, 0),
#          other = if_else(demo_race_a_p___25 == 1, 1, 0)) %>%
#   select(-starts_with("demo_race")) 

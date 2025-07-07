#############################################################################################
### Summary:
### 1. Data cleaning and correct labelling of sex, ethnicity and case/control status
### 2. Definition of adverse early life experiences (experienced before/including the age of 24)
#############################################################################################


library(dplyr)
library(tidyr)

env_data <- read.csv("/SAN/ugi/ukhls/Paul_MS_proj/DIVERGE-PaulBrandesProject_DATA_2025-05-28_1315.csv")



### Exclude all text entries and participants who are neither listed as cases or control, as well as clean up important variables ---------------------------------------------------------------

phenotype <- env_data %>%
	filter(subject_type != 3) %>% # Filter out tentative control
	filter(!grepl("est",subject_id)) %>% # Filter out test cases
	filter(!is.na(screener_age)) %>%
	mutate(
		subject_type_logical = ifelse(subject_type == 1, 1, 0),
		screener_sex = factor(screener_sex, levels = c(1, 2), labels = c("Male", "Female")), 
		ethnic_self_provincial = factor(ethnic_self_provincial, levels = c(1, 2, 3, 4, 5, 6), labels = c("KPK", "Punjab", "Sindh", "Balochistan", "GilgitBaltistan", "AzadKashmir")),
		social_support = osss_number_close_people + osss_help_from_others + osss_interest_from_others
   	)



### Define early life experiences of Stressful Life Events and Other Life Events ---------------------------------------------------------------

event_vars <- c("te_physical_assault___1", 
		"te_sexual_assault___1", 
		"te_other_unwanted_sex___1", 
		"te_natural_disaster___1", 
		"te_war_zone___1", 
		"le_courts_issues", 
		"le_forced_leave_home", 
		"le_domestic_issues", 
		"le_finance_problem", 
		"le_death_parents",
		"le_conflicts",
		"le_residence_problems",
		"le_teenage_pregnant",
		"le_chronic_illness",
		"le_parents_separated",
		"le_death_sibling_friend",
		"le_hospital_fam_members",
		"te_captivity___1",
		"te_other_serious_accident___1",
		"te_transport_accident___1",
		"te_bomb___1",
		"te_injury_smb_else___1"
	       )
age_vars <- c("te_physical_assault_age1", 
	      "te_sexual_assault_age1", 
	      "te_other_unwanted_sex_age1", 
	      "te_natural_disaster_age1", 
	      "te_war_zone_age1", 
	      "le_courts_issues_age", 
	      "le_forced_leave_age", 
	      "le_domestic_issues_age", 
	      "le_finance_problem_age", 
	      "le_death_parents_age",
	      "le_conflicts_age",
	      "le_residence_problems_age",
	      "le_teenage_pregnant_age",
	      "le_chronic_illness_age",
	      "le_parents_separated_age",
	      "le_death_sibling_friend_age",
	      "le_hospital_fam_members_age",
	      "te_captivity_age1",
	      "te_other_serious_accident_age1",
	      "te_transport_accident_age1",
	      "te_bomb_age1",
	      "te_injury_smb_else_age1"
	     )


for (i in seq_along(event_vars)) {
  event <- event_vars[i]
  age <- age_vars[i]
  # Create clean variable name
  new_col <- paste0("early_", gsub("te_|___1|le_", "", event))
  # Logical condition used for both case and control
  is_event_early <- phenotype[[event]] == 1 & 
                    phenotype[[age]] <= 24 & 
                    phenotype[[age]] < 200
  # Apply condition differently depending on case/control
  phenotype[[new_col]] <- ifelse(
    phenotype$subject_type_logical == 1,
    as.integer(is_event_early & phenotype[[age]] <= phenotype$age_onset),
    as.integer(is_event_early)
  )
}



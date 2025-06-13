library(dplyr)
library(survival)
library(MatchIt)


### Manual matching
## Unusually we have more cases than controls, which makes matching without losing a significant number of cases a challenge. Therefore 2 cases are matched to 1 control here:

reverse_match_caliper <- function(data, case_var = "subject_type_logical", match_var = "screener_age", k = 2, caliper = 5) {
  # Split dataset into cases and controls
  controls <- data[data[[case_var]] == 0, ]
  cases <- data[data[[case_var]] == 1, ]
  matched_pairs <- list()
  
  set.seed(123)  # Important e.g. when multiple cases have the exact same distance from a control. Makes it reproducible
  
  for (i in 1:nrow(controls)) {
    # Calculate age differences
    age_diff <- abs(cases[[match_var]] - controls[i, match_var])
    
    # Find cases within the caliper (≤ 5 years difference)
    eligible_cases <- which(age_diff <= caliper)
    
    # Skip if no cases meet the caliper
    if (length(eligible_cases) == 0) next
    
    # Select the k closest cases (within caliper)
    closest_cases <- eligible_cases[order(age_diff[eligible_cases])][1:min(k, length(eligible_cases))]
    
    # Store matches
    matched_pairs[[i]] <- rbind(
      cbind(controls[i, ], pair_id = i, role = "control"),
      cbind(cases[closest_cases, ], pair_id = i, role = "case")
    )
    
    # Remove matched cases
    cases <- cases[-closest_cases, ]
  }
  
  matched_df <- do.call(rbind, matched_pairs)
  return(matched_df)
}

# Run matching
matched_data_caliper <- reverse_match_caliper(phenotype, k = 2, caliper = 5)

# Run conditional logistic regression
clogit(subject_type_logical ~ early_loss + strata(pair_id), data = matched_data_caliper)



### Alternative 1:1 matching using MatchIt
#matchit_result <- matchit(
#  subject_type_logical ~ screener_age,
#  data = phenotype,
#  method = "nearest",
#  distance = "glm"
#  )
#
#matched_data <- match.data(matchit_result)



### Confirm number of exclusions
## Number of excluded cases and controls
n_original_cases <- sum(phenotype$subject_type_logical == 1)
n_matched_cases <- sum(matched_data_caliper$role == "case")

n_original_controls <- sum(phenotype$subject_type_logical == 0)
n_matched_controls <- sum(matched_data_caliper$role == "control")

cat("Retained", n_matched_cases, "of", n_original_cases, "cases (", round(n_matched_cases / n_original_cases * 100, 1), "%) and ", n_matched_controls, "of", n_original_controls, "controls (", round(n_matched_controls / n_original_controls * 100, 1), "%)\n")

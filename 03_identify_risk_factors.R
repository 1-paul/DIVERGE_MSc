#############################################################################################
### Summary:
### 1. Define list of environmental risk factors
### 2. Run regressions for each one:
###   2.1. Without any covariates
###   2.2. With age, sex, and ethnicity as covariates
###   2.3. After age matching
#############################################################################################


library(survival)
library(dplyr)


# Define risk factors
predictors <- c(
  "early_courts_issues",
  "early_forced_leave_home",
  "early_domestic_issues",
  "early_finance_problem",
  "le_refugee_status_ever",
  "early_loss", 
  "early_physical_assault",
  "early_sexual_assault",
  "early_other_unwanted_sex",
  "early_natural_disaster",
  "early_war_zone"
)



### Logistic regression without covariates -----------------------------------------------------------------------------------------------------------------------------------------
# Establish storage variables
pvals <- numeric(length(predictors))
estimates <- numeric(length(predictors))
names(pvals) <- names(estimates) <- predictors

# Run logistic regressions and extract p-values and estimates
for (i in seq_along(predictors)) {
  formula <- as.formula(paste("subject_type_logical ~", predictors[i]))
  model <- glm(formula, data = phenotype, family = binomial)
  smry <- summary(model)
  
  # Extract estimate and p-value for the predictor
  pvals[i] <- coef(smry)[predictors[i], "Pr(>|z|)"]
  estimates[i] <- exp(coef(smry)[predictors[i], "Estimate"])
}

# Apply multiple testing correction
pvals_corrected <- p.adjust(pvals, method = "bonferroni")

# Combine into a data frame
results_df <- data.frame(
  Predictor = predictors,
  Estimate = estimates,
  P_Value = pvals,
  Adjusted_P = pvals_corrected
)

print(results_df)



### Logistic regression with age as a covariate -----------------------------------------------------------------------------------------------------------------------------------------
# Establish storage variables
pvals <- numeric(length(predictors))
estimates <- numeric(length(predictors))
names(pvals) <- names(estimates) <- predictors

# Run logistic regressions with covariates and extract p-values and estimates
for (i in seq_along(predictors)) {
  formula <- as.formula(paste("subject_type_logical ~", predictors[i], "+ screener_age + screener_sex + ethnic_self_provincial"))
  model <- glm(formula, data = phenotype, family = binomial)
  smry <- summary(model)
  
  # Extract p-value
  pvals[i] <- coef(smry)[predictors[i], "Pr(>|z|)"]
  # Extract odds
  estimates[i] <- exp(coef(smry)[predictors[i], "Estimate"])
}

# Apply multiple testing correction
pvals_corrected <- p.adjust(pvals, method = "bonferroni")

# Combine into a data frame
results_df <- data.frame(
  Predictor = predictors,
  Estimate = estimates,
  P_Value = pvals,
  Adjusted_P = pvals_corrected
)

print(results_df)



### Conditional logistic regressions for age-matched case/controls after running "02_age_matching.R" ----------------------------------------------------------------------------------------------
# Establish storage variables
pvals <- numeric(length(predictors))
estimates <- numeric(length(predictors))
names(pvals) <- names(estimates) <- predictors

# Run conditional logistic regressions and extract p-values and estimates
for (i in seq_along(predictors)) {
  formula <- as.formula(paste("subject_type_logical ~", predictors[i], "+ strata(pair_id)"))
  model <- clogit(formula, data = matched_data_caliper)
  smry <- summary(model)
  
  # Extract p-value
  pvals[i] <- coef(smry)[predictors[i], "Pr(>|z|)"]
  # Extract odds
  estimates[i] <- coef(smry)[predictors[i], "exp(coef)"]
}

# Apply multiple testing correction
pvals_corrected <- p.adjust(pvals, method = "bonferroni")

# Combine into a data frame
results_df <- data.frame(
  Predictor = predictors,
  Odds = estimates,
  P_Value = pvals,
  Adjusted_P = pvals_corrected
)

print(results_df)


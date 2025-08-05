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
  'early_physical_assault',
  'early_sexual_assault2',
  'early_natural_disaster',
  'early_war_zone',
  'early_courts_issues',
  'early_forced_leave_home',
  'early_domestic_issues',
  'early_finance_problem',
  'early_death_parents',
  'early_conflicts',
  'early_residence_problems',
  'early_teenage_pregnant',
  'early_parents_separated',
  'early_death_sibling_friend',
  'early_hospital_fam_members',
  'early_captivity',
  'early_other_serious_accident',
  'early_transport_accident',
  'early_bomb',
  'early_injury_smb_else'
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



### Conditional logistic regressions for age-matched case/controls after running "02_age_matching.R" 
# Establish storage variables
pvals <- numeric(length(predictors))
estimates <- numeric(length(predictors))
lower_ci <- numeric(length(predictors))
upper_ci <- numeric(length(predictors))
names(pvals) <- names(estimates) <- names(lower_ci) <- names(upper_ci) <- predictors

# Run conditional logistic regressions and extract p-values, estimates, and CIs
for (i in seq_along(predictors)) {
  formula <- as.formula(paste("subject_type_logical ~", predictors[i], "+ strata(subclass)"))
  model <- clogit(formula, data = matched_data)
  smry <- summary(model)
  
  # Extract p-value
  pvals[i] <- coef(smry)[predictors[i], "Pr(>|z|)"]
  # Extract odds ratio
  estimates[i] <- coef(smry)[predictors[i], "exp(coef)"]
  # Extract 95% CI
  lower_ci[i] <- smry$conf.int[predictors[i], "lower .95"]
  upper_ci[i] <- smry$conf.int[predictors[i], "upper .95"]
}

# Apply multiple testing correction
pvals_corrected <- p.adjust(pvals, method = "bonferroni")

# Combine into a data frame
results_df <- data.frame(
  Predictor = predictors,
  Odds_Ratio = estimates,
  CI_lower = lower_ci,
  CI_upper = upper_ci,
  P_Value = pvals,
  Adjusted_P = pvals_corrected
)

print(results_df)

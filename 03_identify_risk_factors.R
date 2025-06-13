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



### Log regression
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



### Log regression with age
pvals <- numeric(length(predictors))
estimates <- numeric(length(predictors))
names(pvals) <- names(estimates) <- predictors

# Run logistic regressions and extract p-values and estimates
for (i in seq_along(predictors)) {
  formula <- as.formula(paste("subject_type_logical ~", predictors[i], "+ screener_age"))
  model <- glm(formula, data = phenotype, family = binomial)
  smry <- summary(model)
  
  # Extract p-value for the predictor
  pvals[i] <- coef(smry)[predictors[i], "Pr(>|z|)"]
  # get the odds
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

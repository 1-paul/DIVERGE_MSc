library(MASS)
library(dplyr)

# Filter cases
cases <- phenotype %>%
  filter(subject_type == 1)

# Define symptoms
symptoms <- c(
  "dip_op37_dysphoria",
  "dip_op39_anhedonia",
  "dip_op25_energy",
  "dip_op24_slowed",
  "dip_op41_concentration",
  "dip_op42_reproach",
  "dip_op43_suicidality",
  "dip_op48_decr_appetite",
  "dip_op49_weightloss",
  "dip_op50_incr_appetite",
  "dip_op51_weightgain",
  "dip_op44_insomn_initial",
  "dip_op46_insomn_terminal",
  "dip_op47_excess_sleep",
  "dip_op90_course"
)

# Define predictor
predictor <- "early_domestic_issues"

# Initialize results list
results_list <- list()

# Loop over each symptom
for (symptom in symptoms) {
  # Create formula
  formula <- as.formula(paste0("as.ordered(", symptom, ") ~ ", predictor, "+ screener_age"))
  # Fit ordinal logistic regression model
  model <- polr(formula, data = cases, Hess = TRUE)
  smry <- summary(model)
  # Extract coefficient table
  coef_table <- coef(smry)
  # Check if predictor is in the model (defensive programming)
  if (predictor %in% rownames(coef_table)) {
    estimate <- exp(coef_table[predictor, "Value"])  # Odds ratio (exponentiated coefficient)
    t_value <- coef_table[predictor, "t value"]
    p_value <- 2 * (1 - pnorm(abs(t_value)))        # Wald test p-value
  } else {
    estimate <- NA
    p_value <- NA
  }
  # Store result
  results_list[[length(results_list) + 1]] <- data.frame(
    Symptom = symptom,
    Predictor = predictor,
    Estimate = estimate,
    P_Value = p_value
  )
}

# Combine all results into a dataframe
results_df <- do.call(rbind, results_list)
# Adjust p-values for multiple testing (Bonferroni)
results_df$Adjusted_P <- p.adjust(results_df$P_Value, method = "bonferroni")
# Print results
print(results_df)



### Plot Age of Onset ~ Early domestic issues
cases <- cases %>%
	filter(!is.na(early_domestic_issues))

ggplot(cases, aes(x = age_onset, fill = early_domestic_issues)) +
	geom_density(alpha = 0.5, position = "identity") +
	labs(
	    title = "Age of Onset Distribution by Domestic Issues Status",
	    y = "Density (relative distribution)", 
	    fill = "Early Domestic Issues"
	) +
	geom_vline(
	    data = onset_means, 
	    aes(xintercept = grp.mean, color = early_domestic_issues), 
	    linetype = "dashed", 
	    linewidth = 1  # Use linewidth, not size
	) +
	scale_fill_manual(
		values = c("No Early Domestic Issues" = "#4dac26", "Early Domestic Issues" = "#d01c8b")
	) +
	scale_color_manual(
		values = c("No Early Domestic Issues" = "#4dac26", "Early Domestic Issues" = "#d01c8b")
	) +
	theme_minimal()

#############################################################################################
### Summary:
### 1. Definition of cases, symptoms, and risk factor
#############################################################################################


library(MASS)
library(dplyr)
library(ggplot2)


### Define symptoms -----------------------------------------------------------------------------
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
  "dip_op90_course",
  "sbq_suicide_idea_behav",
  "suicide_attempts_times",
  "sbq_suicide_thought",
  "sbq_suicide_told_smb"
)


### Adversity score with found consistent risk factors in 03_identify_risk_factors.R
phenotype <- phenotype %>%
	mutate(adversity_score = rowSums(across(c(
		early_physical_assault, 
		early_sexual_assault2,
		early_courts_issues,
		early_domestic_issues,
		early_finance_problem, 
		early_conflicts)), na.rm = TRUE)
	      )



# Define risk factor
#predictor <- "early_domestic_issues"
predictor <- "adversity_score"



cases <- phenotype %>%
  filter(subject_type == 1)


### Run ordinal logistic regressions -----------------------------------------------------------------------------------------------
## An ordinal logistic regressions is necessary since all the indicators of severity are ordered categorical variable e.g. "Low," "Medium," "High" 

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
	# Check if predictor is in the model
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



### Regression for Age of Onset -------------------------------------------------------------------------------------
# Filter cases with invalid ages of onset
cases <- cases %>%
	filter(age_onset < 200)

formula <- as.formula(paste("age_onset ~", predictor, "+ screener_age"))
model <- lm(formula, data = cases)

summary(model)



### Boxplot --------------------------------------------------------------------------------------------------------
# First, create the 6 categories
phenotype <- phenotype %>%
  mutate(num_of_adversities = case_when(
    adversity_score == 0 ~ "0",
    adversity_score == 1 ~ "1",
    adversity_score == 2 ~ "2",
    adversity_score == 3 ~ "3",
    adversity_score == 4 ~ "4",
    adversity_score >= 5 ~ "5+"
  )) %>%
  mutate(num_of_adversities = factor(num_of_adversities, 
                                    levels = c("0", "1", "2", "3", "4", "5+")))

cases <- phenotype %>%
  filter(subject_type == 1)

# Filter out NA values
cases <- cases %>%
	filter(!is.na(num_of_adversities)) %>%
	filter(age_onset < 200)

# Create a facetted boxplot
ggplot(cases, aes(x = num_of_adversities, y = age_onset, fill = num_of_adversities)) +
	geom_boxplot() +
	labs(
    		title = "Age of onset by Number of experienced Early Aversities",
    		x = "Number of Adversities",
    		y = "Age of Onset"
  	) +
  	scale_fill_manual(
	    values = c("0" = "#ffffb2", 
	               "1" = "#fed976", 
	               "2" = "#feb24c", 
	               "3" = "#fd8d3c", 
	               "4" = "#f03b20", 
	               "5+" = "#bd0026"),
	    name = "Number of Adversities"
  	) +
	theme_bw() +
	theme(legend.position = "none") +  # Remove legend since x-axis shows the categories
	stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black")  # Add mean markers



### Plot Age of Onset ~ Risk Factor(e.g. early domestic issues) ---------------------------------------------------------------------------------

## First define categories of exposure:
phenotype <- phenotype %>%
	mutate(num_of_adversities = case_when(
		adversity_score == 0 ~ "0",
    		adversity_score >= 1 & adversity_score <= 2 ~ "1-2",
    		adversity_score >= 3 ~ "3+")) %>%
	mutate(num_of_adversities = factor(num_of_adversities, levels = c("0", "1-2", "3+")))


# Filter out cases with missing risk factor
cases <- cases %>%
	filter(!is.na(num_of_adversities))


onset_means <- cases %>%
	group_by(num_of_adversities) %>% 
	summarize(grp.mean = mean(age_onset, na.rm = TRUE))

## Then plot:
ggplot(cases, aes(x = age_onset, fill = factor(num_of_adversities))) +
  geom_density(alpha = 0.5, position = "identity") +
  labs(
    title = "Age of onset distribution by number of experienced early adversities",
    x = "Age of Onset",
    y = "Percentage of Total Cases", 
    fill = "Number of Adversities"
  ) +
  geom_vline(
    data = onset_means, 
    aes(xintercept = grp.mean, color = factor(num_of_adversities)), 
    linetype = "dashed", 
    linewidth = 1
  ) +
  scale_fill_manual(
    values = c("0" = "#ffeda0", "1-2" = "#feb24c", "3+" = "#f03b20"),
    labels = c("0", "1-2", "3+")
  ) +
  scale_color_manual(
    values = c("0" = "#ffeda0", "1-2" = "#feb24c", "3+" = "#f03b20"),
    labels = c("0", "1-2", "3+")
  ) +
  theme_bw()


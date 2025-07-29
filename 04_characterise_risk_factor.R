#############################################################################################
### Summary:
### 1. Definition of cases, symptoms, and risk factor
#############################################################################################


library(MASS)
library(dplyr)
library(ggplot2)
library(patchwork)


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


### Adversity score with found consistent risk factors in 03_identify_risk_factors.R -----------------------------------------------------------------------------------------------
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



### Plot adversity score (NEEDS TO BE REWORKED) -----------------------------------------------------------------------------------------------
plot_data <- phenotype %>%
	mutate(Status = factor(subject_type_logical, levels = c(0, 1), labels = c("Control", "Case"))) %>%
  	count(Status, adversity_score) %>%
  	group_by(Status) %>%
  	mutate(prop = n / sum(n)) %>%
	ungroup() %>%  # Ungroup just in case it's grouped
	mutate(Status = factor(Status, levels = c("Control", "Case"))) %>%
	complete(adversity_score, Status, fill = list(prop = 0))


# p1: Top plot
p1 <- ggplot(plot_data, aes(x = adversity_score, y = prop, fill = Status)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8, alpha = 0.7) +
  scale_fill_manual(values = c("Control" = "#008837", "Case" = "#7b3294")) +
  scale_x_continuous(breaks = 0:6, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.8, 1)) +
  theme_bw(base_family = "CMU Serif") +
  labs(title = "Proportion per Adversity Score", x = NULL, y = NULL, fill = NULL) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7)
  )

# p2: Bottom plot
p2 <- ggplot(plot_data, aes(x = adversity_score, y = prop, fill = Status)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8, alpha = 0.7) +
  geom_hline(yintercept = 0.2, linewidth = 1.4, color = "black") +
  scale_fill_manual(values = c("Control" = "#008837", "Case" = "#7b3294")) +
  scale_x_continuous(breaks = 0:6, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(x = "Adversity Score", y = NULL) +
  theme_bw(base_family = "CMU Serif") +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.7),
    axis.line.y.left = element_line(color = "black", linewidth = 0.7),
    axis.line.x.top = element_line(color = "black", linewidth = 0.7)
  )

# Shared y-axis label
y_axis <- ggplot() +
  theme_void() +
  annotate(
    "text",
    x = 0.9, y = 0.5,
    label = "Proportion of Total Participants with Adversity Score",
    angle = 90, size = 5,
    fontface = "bold", family = "CMU Serif"
  )

# Combine using patchwork
combined <- (y_axis + plot_spacer() + (p1 / p2 + plot_layout(heights = c(1, 1.05)))) +
  plot_layout(ncol = 3, widths = c(0.05, 0.01, 1)) &
  theme(plot.margin = margin(t = 4, r = 8, b = 4, l = 8, unit = "pt"))

print(combined)



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


ggplot(cases, aes(x = num_of_adversities, y = age_onset, fill = num_of_adversities)) +
  geom_boxplot() +
  scale_fill_manual(
    values = c("0" = "#ffffb2", 
               "1" = "#fed976", 
               "2" = "#feb24c", 
               "3" = "#fd8d3c", 
               "4" = "#f03b20", 
               "5+" = "#bd0026"),
    name = "Number of Adversities"
  ) +
  labs(
    x = "Number of Adversities",
    y = "Age of Onset"
  ) +
  theme_bw(base_family = "CMU Serif") +
  theme(
    plot.title = element_blank(),  # Remove title
    legend.position = "none",      # Remove legend
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.border = element_blank(),                         # Remove full panel border
    axis.line = element_line(color = "black", linewidth = 0.8), # Show left and bottom axis lines
    axis.line.y.right = element_blank(),
    axis.line.x.top = element_blank()
  )



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


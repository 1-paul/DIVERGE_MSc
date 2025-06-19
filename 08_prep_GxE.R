library(dplyr)

# Filter phenotype
early_domestic_issues_df <- phenotype %>%
	select(subject_id, early_domestic_issues)

# Load covariates file
covariates <- read.table(
	"/cluster/project2/DIVERGE/20250605_GWAS/covariates.txt",
	header = TRUE,
	sep = "\t",          
	comment.char = ""    # Disable skipping lines starting with #
)

# Merge by ID, just adding early_domestic_issues
combined_df <- merge(covariates, early_domestic_issues_df, by.x = "IID", , by.y = "subject_id")


# Output new covariates file
write.table(combined_df, 
           "covariates2.txt", 
           sep = "\t", 
           quote = FALSE, 
           row.names = FALSE,
           col.names = TRUE)

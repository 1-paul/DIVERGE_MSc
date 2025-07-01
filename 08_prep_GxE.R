library(dplyr)


# Filter phenotype
early_domestic_issues_df <- phenotype %>%
    select(subject_id, early_domestic_issues, adversity_score)


# Load covariates file
covariates <- read.table(
    "/cluster/project2/DIVERGE/20250620_GWAS/GWAS/covariates.txt",
    header = TRUE,
    sep = "\t",          
    comment.char = ""    # Disable skipping lines starting with #
)


# Merge by ID, just adding early_domestic_issues
combined_df <- merge(covariates, early_domestic_issues_df, by.x = "IID", by.y = "subject_id")


# Load freeze2 IDs file
freeze2_ids <- read.table(
    "/cluster/project2/DIVERGE/munim_workspace/QC_pipeline_20250514/09_batch_effects/freeze2_ids.txt",
    header = FALSE,
    col.names = c("FID", "IID")  # Naming columns for clarity
)


# Create freeze variable (2 for freeze 2 individuals, 1 for freeze 1 individuals)
combined_df <- combined_df %>%
  mutate(freeze = if_else(IID %in% freeze2_ids$IID, 2, 1))


# Rename columns and remove unwanted columns
combined_df <- combined_df %>%
  rename("#FID" = FID) %>%  # Rename X.FID
  select("#FID", IID, everything())


# Output new covariates file
write.table(combined_df, 
           "covariates2.txt", 
           sep = "\t", 
           quote = FALSE, 
           row.names = FALSE,
           col.names = TRUE)

#############################################################################################
### Summary:
### 1. Merge phenotype and genetic data
#############################################################################################

library(dplyr)

covariates <- read.table(
    "/cluster/project2/DIVERGE/20250620_GWAS/GWAS/covariates.txt",
    header = TRUE,
    sep = "\t",          
    comment.char = ""
)

freeze2_ids <- read.table(
    "/cluster/project2/DIVERGE/munim_workspace/QC_pipeline_20250514/09_batch_effects/freeze2_ids.txt",
    header = FALSE,
    col.names = c("FID", "IID") 
)



### Merge relevant variables of the phenotype data with the genetic data ##############################################################################################
relevant_pheno <- phenotype %>%
    select(subject_id, early_domestic_issues, adversity_score)

combined_df <- merge(covariates, relevant_pheno, by.x = "IID", by.y = "subject_id")

combined_df <- combined_df %>%
    mutate(freeze = if_else(IID %in% freeze2_ids$IID, 2, 1)) %>% # Create batch variable (1 for freeze 1 individuals, 2 for freeze 2 individuals)
    rename("#FID" = FID) %>% 
    select("#FID", IID, everything())

# Output new covariates file
write.table(
    combined_df, 
    "covariates2.txt", 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE,
    col.names = TRUE
)

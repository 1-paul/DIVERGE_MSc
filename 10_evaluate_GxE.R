#############################################################################################
### Summary:
### 1. Format all the GxE result files correctly
### 2. Get the common SNPs with p-values above suggestive signficance
### 3. Plot p-Value of SNP vs beta of the interaction
#############################################################################################



library(dplyr)
library(ggplot2)
library(tidyr)
library(qqman)

# Define which files to use
file_gxe_results_snps <- "/home/pbrandes/20250701_GxE/gxe_results_snp_edi_snpxedi_results_only.txt"
file_gxe_results_all <- "/cluster/project2/DIVERGE/20250701_GxE/00_gwas_results.PHENO1.glm.logistic"
file_frq_data <- "/cluster/project2/DIVERGE/20250620_GWAS/QC/00_plink_files/02_call_rate_95g_95m.frq"
file_pc_results <- "/cluster/project2/DIVERGE/20250701_GxE/covariates.txt"



### Correctly format GxE results -----------------------------------------------------------------------------
gxe_results_snps <- read.table(file_gxe_results_snps, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

colnames(gxe_results_snps) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "LOG_OR", "LOG_OR_SE", "Z_STAT", "P")

gxe_results_snps <- gxe_results_snps %>% 
	mutate(CHROM = case_when(CHROM == "Y" ~ "23",CHROM == "XY" ~ "24",CHROM == "MT" ~ "25",TRUE ~ as.character(CHROM))) %>%
	filter(CHROM != 0) %>%
	mutate(CHROM = as.numeric(CHROM)) %>%
	filter(CHROM <= 22) %>% # Autosomes only
	filter(!is.na(P))



### Get minor allele frequencies of significant variants --------------------------------------------------
frq_data <- read.table(file_frq_data, header = TRUE, stringsAsFactors = FALSE)

# Filter out rare variants 
common_variants <- frq_data %>% 
	filter(MAF >= 0.05)  # Filter out MAF < 5%
	
gxe_results_snps <- gxe_results_snps %>%
	filter(ID %in% common_variants$SNP)



#### Get significant SNPs with TEST == "ADD" and P < 0.00001 --------------------------------------------------
significant_snps <- gxe_results_snps %>%
	filter(TEST == "ADD", P < 0.00001) %>%
	pull(ID)  # Extract SNP IDs

# Filter all rows (ADD, risk factor (e.g. early domestic issues), ADD x risk factor) for those SNPs
significant_results <- gxe_results_snps %>%
	filter(ID %in% significant_snps)  # Keeps all TEST types for those SNPs



### Widen to one column per SNP ---------------------------------------------------------------------------------------------------------
wider_df <- gxe_results_snps %>%
	# Will have to adapt names here, when changing risk factors
	mutate(
		suffix = case_when(
      		TEST == "ADD" ~ "_snp",
      		TEST == "adversity_score" ~ "_env",
      		TEST == "ADDxadversity_score" ~ "_gxe"
    		)
	) %>%
 	pivot_wider(
    		id_cols = c(CHROM, POS, ID, REF, ALT, A1, OBS_CT),
    		names_from = suffix,
    		values_from = c(LOG_OR, LOG_OR_SE, Z_STAT, P),
    		names_glue = "{.value}{suffix}"
	) %>%
  	mutate(log10_P_snp = log10(P_snp))



### plot beta of p-value vs beta of the interactions ---------------------------------------------------------------------------------------------------------

# png("0807_beta_vs_p_1.png", width=1200, height=1200)
ggplot(wider_df, aes(x = LOG_OR_gxe, y = log10_P_snp)) +
	geom_point(alpha = 0.7, size = 1.2) +
  
 	# Reference lines
	geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.5) +
	#geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5) + # needs to be change so there is up to 3x the sd
 	geom_hline(yintercept = log10(1e-5), linetype = "dashed", color = "blue", linewidth = 0.5) +  # y = -4
  
	# Custom y-axis (reverse breaks to put p=1 at the bottom)
	scale_y_reverse(
		breaks = log10(c(1, 0.1, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7)),  # Breaks at log10(p)
		labels = c("1", "1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7")  # Label with raw p-values
	) +
  
	# Axis labels
	labs(
		x = "Log Odds Ratio (GxE Interaction)",
		y = "P-value (log10 scale)",
		title = "GxE Interaction vs. SNP Main Effects"
	) +
  
	theme_bw() +
	theme(
		panel.grid.minor = element_blank(),
		plot.title = element_text(hjust = 0.5)
  	)

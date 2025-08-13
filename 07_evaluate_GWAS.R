#############################################################################################
### Summary:
### 1. Clean and filter GWAS results
### 2. Make manhattan and qq plots
### 3. Assess effect of principle components
#############################################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(qqman)
library(extrafont)
# Only once: font_import(paths = "~", pattern = "cmunrm", prompt = FALSE)
loadfonts(device = "pdf")
theme(text = element_text(family = "CMU Serif"))


gwas_results_snps <- read.table("/cluster/project2/DIVERGE/20250620_GWAS/GWAS/00_gwas_results_snp_results_only.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(gwas_results_snps) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "OR", "LOGOR_SE", "Z_STAT", "P")
gwas_results_all <- read.table("/cluster/project2/DIVERGE/20250620_GWAS/GWAS/00_gwas_results.PHENO1.glm.logistic")
colnames(gwas_results_all) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "OR", "LOG(OR)_SE", "Z_STAT", "P")
frq_data <- read.table("/cluster/project2/DIVERGE/20250620_GWAS/QC/00_plink_files/08_hardy_final.frq", header = TRUE, stringsAsFactors = FALSE)
pc_results_all <- read.table("/cluster/project2/DIVERGE/20250620_GWAS/GWAS/covariates.txt")
colnames(pc_results_all) <- c("FID", "IID", "SEX", "PHENO1", "FID2", "IID2", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")



### Clean GWAS results #############################################################################################
gwas_results_snps <- gwas_results_snps %>% 
	mutate(
		CHROM = case_when(CHROM == "Y" ~ "23",CHROM == "XY" ~ "24",CHROM == "MT" ~ "25",TRUE ~ as.character(CHROM)),
		CHROM = as.numeric(CHROM),
		log_or = log(OR),
		CI_lower = exp(log_or - 1.96 * LOGOR_SE),
		CI_upper = exp(log_or + 1.96 * LOGOR_SE)
	) %>%
	filter(CHROM != 0 & CHROM <= 22 & !is.na(P))



### Get significiant hits with MAF > 5% #############################################################################################
common_variants <- frq_data %>% 
	filter(MAF >= 0.05)

significant_common_hits <- gwas_results_snps %>%
	filter(ID %in% common_variants$SNP) %>%
	filter(P < 0.00001)	



### Manhattan (qqman) #############################################################################################
# For HPC env.: png("manhattan_plot.png", width=1200, height=600) 
manhattan(gwas_results_snps, chr="CHROM", bp="POS", snp="ID", p="P", col = c("#5768f1", "#9ea7f7"), ylim = c(0, 8))
# For HPC env.: dev.off()



### QQ plot (ggplot) #############################################################################################
pvals <- gwas_results_snps %>%
	pull(P) %>%        # Extract as vector
	na.omit() %>%
	sort()

n <- length(pvals)

qq_df <- data.frame(
	expected = -log10(ppoints(n)),
	observed = -log10(pvals),
	lower = -log10(qbeta(0.025, 1:n, n:1)),
	upper = -log10(qbeta(0.975, 1:n, n:1))
)

ggplot(qq_df, aes(x = expected, y = observed)) +
	geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray90", alpha = 0.7) +
	geom_point(size = 1, alpha = 0.7) +
	geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
	labs(
		x = "Expected -log10(p)",
		y = "Observed -log10(p)"
	) +
	theme_bw(base_family = "CMU Serif") +
	theme(
		plot.title = element_text(size = 20, face = "bold"),
		axis.title = element_text(size = 16),
		axis.text = element_text(size = 14),
		panel.border = element_blank(),
		axis.line = element_line(color = "black", linewidth = 0.7)
	)



### Evaluate principle components #############################################################################################
## Get the Mean, Median, Max, and Min
pc_summary <- gwas_results_all %>%
	mutate(CHROM = as.numeric(CHROM)) %>%
	filter(
		!is.na(CHROM),
		!is.na(P),
		P > 0,
		P <= 1,
		!is.infinite(P)
	) %>%
	group_by(TEST) %>% 
	summarize(
		mean = mean(P),
		median = median(P),
		maximum = max(P),
		minimum = min(P)
	)

## Null model
combined_df <- merge(phenotype, pc_results_all, by.x = "subject_id", , by.y = "IID")

model <- glm(subject_type_logical ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, data = combined_df, family = binomial)

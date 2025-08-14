#############################################################################################
### Summary:
### 1. Clean and filter all GxE result
### 2. Get the common SNPs with p-values above suggestive signficance
### 3. Plot p-Value of SNP vs beta of the interaction
### 4. Find SNPs reported in other papers
### 5. Plot p-Value of SNP vs p-Value of the interaction
#############################################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(qqman)
library(extrafont)
#font_import(paths = "~", pattern = "cmunrm", prompt = FALSE)
loadfonts(device = "pdf")
theme(text = element_text(family = "CMU Serif"))


# Change if script is to be adapted for different risk factor
file_gxe_results_snps <- "/home/pbrandes/20250701_GxE/gxe_results_snp_edi_snpxedi_results_only.txt"
file_gxe_results_all <- "/cluster/project2/DIVERGE/20250701_GxE/00_gwas_results.PHENO1.glm.logistic"
file_frq_data <- "/cluster/project2/DIVERGE/20250620_GWAS/QC/00_plink_files/08_hardy_final.frq"
file_pc_results <- "/cluster/project2/DIVERGE/20250701_GxE/covariates.txt"

# SNPs with significant interactions derived from other MDD GxE papers
snps_coleman <- "/home/pbrandes/20250701_GxE/snps_coleman_2020.txt"
snps_peterson <- "/home/pbrandes/20250701_GxE/snps_peterson_2018.txt"
snps_ye <- "/home/pbrandes/20250701_GxE/snps_ye_2021.txt"
snps_dunn_aa <-"snps_dunn_2016_AA.txt"
snps_dunn_hl <- "snps_dunn_2016_HL.txt"
snps_ArnauSoler <- "snps_ArnauSoler_2019.txt"



### Correctly format GxE results ####################################################################################################################################################################
gxe_results_snps <- read.table(file_gxe_results_snps, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(gxe_results_snps) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "OR", "LOGOR_SE", "Z_STAT", "P")

gxe_results_snps <- gxe_results_snps %>% 
	mutate(
		CHROM = case_when(CHROM == "Y" ~ "23",CHROM == "XY" ~ "24",CHROM == "MT" ~ "25",TRUE ~ as.character(CHROM)),
		CHROM = as.numeric(CHROM)
	) %>%
	filter(CHROM != 0 & CHROM <= 22 & !is.na(P))



### Get significiant GxE hits with MAF > 5% ####################################################################################################################################################################
frq_data <- read.table(file_frq_data, header = TRUE, stringsAsFactors = FALSE)

common_variants <- frq_data %>% 
	filter(MAF >= 0.05)



### Widen to one column per SNP ####################################################################################################################################################################
# Will have to adapt names here, when changing risk factors!
wider_df <- gxe_results_snps %>%
	filter(ID %in% common_variants$SNP) %>%
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
    		values_from = c(OR, LOGOR_SE, Z_STAT, P),
    		names_glue = "{.value}{suffix}"
	) %>%
  	mutate(
		log10_P_snp = log10(P_snp),
		log10_P_gxe = log10(P_gxe),
		LOG_OR_gxe = log(OR_gxe),
		significant_gxe = if_else(P_gxe <= 0.05, "1", "0"),
		significant_gxe = factor(significant_gxe, levels = c(0, 1)),
		highly_sign_gxe = if_else(P_gxe <= 1e-5, "1", "0"),
		highly_sign_gxe = factor(highly_sign_gxe, levels = c(0, 1)),
		significant_both = if_else(P_snp <= 1e-5 | P_gxe <= 1e-5, "1", "0"),
		significant_both = factor(significant_both, levels = c(0, 1)),
		log_or_gxe = log(OR_gxe),
		CI_lower_gxe = exp(log_or_gxe - 1.96 * LOGOR_SE_gxe),
		CI_upper_gxe = exp(log_or_gxe + 1.96 * LOGOR_SE_gxe)
	)



#### Get significant SNPs ##########################################################################################################################################
significant_snps <- wider_df %>%
	filter(P_snp < 0.00001) %>%
	arrange(P_snp)

significant_interaction <- wider_df %>%
	filter(P_gxe < 0.00001) %>%
	arrange(P_gxe)	



### plot beta of p-value of main effect vs beta of the interactions ###############################################################################################################
# For HPC env.: png("3007_beta_vs_p_2.png", width=1000, height=800)
ggplot(wider_df, aes(x = LOG_OR_gxe, y = log10_P_gxe)) +
	geom_point(aes(color = highly_sign_gxe, size = highly_sign_gxe), alpha = 0.7) + # Switch between significant_gxe and highly_sign_gxe
	geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
	geom_hline(yintercept = log10(1e-5), linetype = "dashed", color = "blue", linewidth = 0.5) +
	geom_hline(yintercept = log10(5e-8), linetype = "dashed", color = "red", linewidth = 0.5) +
	scale_y_reverse(
		limits = log10(c(1, 5e-9)),
		breaks = log10(c(1, 0.1, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)),
		labels = c("1", "1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"),
		expand = c(0, 0)
	) +
	scale_color_manual(
		values = c("1" = "#5768f1", "0" = "#969696"),
		name = "GxE Significance",
		labels = c("1" = "Significant Interaction (<= 0.05)", "0" = "No Significant Interaction")
	) +
	scale_size_manual(
		values = c("1" = 1, "0" = 1),  
		guide = "none"  # Hide size legend since redundant with colour
	) +
	theme_bw() +
	theme(
		text = element_text(family = "CMU Serif"),
		panel.grid.minor = element_blank(),
		plot.title = element_text(hjust = 0.5, size = 20),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.text = element_text(size = 16),
		legend.title = element_text(size = 18),
		legend.text = element_text(size = 16),
		legend.position = "none",
		panel.border = element_blank(),
		axis.line = element_line(color = "black")
	)
# For HPC env.: dev.off()



### Find SNPs which have been reported as having singificant interactions in other papers ###########################################################################
# All SNPs were manually extracted from papers, so processing is not standardised
df_coleman <- read.table(snps_coleman, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
df_coleman <- df_coleman %>%
	filter(V12 <= 1e-5) %>%
	select(V2)

df_ye <- read.table(snps_ye, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df_ye <- df_ye %>%
	select(SNP)

df_peterson <- read.table(snps_peterson, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

df_dunn_aa <- read.table(snps_dunn_aa, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

df_dunn_hl <- read.table(snps_dunn_hl, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

df_ArnauSoler <- read.table(snps_ArnauSoler, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

wider_df <- wider_df %>%
	mutate(
		reported_in = case_when(
			ID %in% df_coleman$V2 ~ "coleman",
			ID %in% df_peterson$V1 ~ "peterson",
			ID %in% df_dunn_aa$V1 ~ "dunn_aa",
			ID %in% df_dunn_hl$V1 ~ "dunn_hl",
			ID %in% df_ArnauSoler$V1 ~ "arnau_soler",
			ID %in% df_ye$SNP ~ "ye",
			TRUE ~ NA_character_  # for SNPs not found in any reference
		)
	)

wider_df %>% 
	filter(!is.na(reported_in)) %>% 
	select(ID, P_snp, P_gxe, reported_in)



### Plot main effect pval vs interaction pval #######################################################################################################################
# For HPC env.: png("1707_pgxe_vs_psnp_1.png", width=800, height=800)
ggplot(wider_df, aes(x = log10_P_gxe, y = log10_P_snp)) +
	geom_point(aes(color = significant_both), alpha = 0.7) + # Switch between significant_gxe and highly_sign_gxe
	geom_vline(xintercept = log10(1e-5), linetype = "dashed", color = "blue", linewidth = 0.5) +
	geom_vline(xintercept = log10(5e-8), linetype = "dashed", color = "red", linewidth = 0.5) +
	geom_hline(yintercept = log10(1e-5), linetype = "dashed", color = "blue", linewidth = 0.5) +
	geom_hline(yintercept = log10(5e-8), linetype = "dashed", color = "red", linewidth = 0.5) +
	scale_color_manual(
		values = c("1" = "#5768f1", "0" = "#969696")
	) +
	scale_y_reverse(
		limits = log10(c(1, 5e-9)),
		breaks = log10(c(1, 0.1, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9)),
		labels = c("1", "1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8", "1e-9"),
		expand = c(0, 0)
	) +
	scale_x_reverse(
		limits = log10(c(1, 5e-9)),
		breaks = log10(c(1, 0.1, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9)),
		labels = c("1", "1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8", "1e-9"),
		expand = c(0, 0)
	) +
	labs(
		x = "P-value of the GxE interaction",
		y = "P-value of the SNP main effect",
	) +
	theme_bw() +
	theme(
		text = element_text(family = "CMU Serif"),
		panel.grid.minor = element_blank(),
		plot.title = element_text(hjust = 0.5, size = 20),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.text = element_text(size = 16),
		legend.position = "none",
		panel.border = element_blank(),
		axis.line = element_line(color = "black")
	)
# For HPC env.: dev.off()

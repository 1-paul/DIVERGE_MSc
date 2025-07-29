library(dplyr)
library(ggplot2)
library(tidyr)
library(qqman)

gwas_results_snps <- read.table("/cluster/project2/DIVERGE/20250620_GWAS/GWAS/00_gwas_results_snp_results_only.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

colnames(gwas_results_snps) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "OR", "LOG(OR)_SE", "Z_STAT", "P")

gwas_results_snps <- gwas_results_snps %>% 
	mutate(CHROM = case_when(CHROM == "Y" ~ "23",CHROM == "XY" ~ "24",CHROM == "MT" ~ "25",TRUE ~ as.character(CHROM))) %>%
	filter(CHROM != 0) %>%
	mutate(CHROM = as.numeric(CHROM)) %>%
	filter(CHROM <= 22) %>% # Autosomes only
	filter(!is.na(P))


### Get minor allele frequencies of significant variants --------------------------------------------------
frq_data <- read.table("/cluster/project2/DIVERGE/20250620_GWAS/QC/00_plink_files/08_hardy_final.frq", header = TRUE, stringsAsFactors = FALSE)

significant_hits <- gwas_results_snps %>%
	filter(P < 0.00001) %>%
	select(ID)

significant_hits_maf <- frq_data %>% 
	filter(SNP %in% significant_hits$ID)

# Filter out rare variants 
common_variants <- frq_data %>% 
	filter(MAF >= 0.05)  # Filter out MAF < 5%
	
gwas_results_snps <- gwas_results_snps %>%
	filter(ID %in% common_variants$SNP)

# Get significant hits of common variants
significant_common_hits_maf <- gwas_results_snps %>%
	filter(P < 0.00001) 	



### Manhattan & QQ plot using qqman package ---------------------------------------------------------

# png("manhattan_plot.png", width=1200, height=600) 
manhattan(gwas_results_snps, chr="CHROM", bp="POS", snp="ID", p="P", col = c("#5768f1", "#9ea7f7"), ylim = c(0, 8))
# dev.off()


pvals <- gwas_results_snps$P
pvals <- sort(na.omit(pvals))
n <- length(pvals)

# Expected and observed -log10(p)
df <- data.frame(
  expected = -log10(ppoints(n)),
  observed = -log10(pvals),
  lower = -log10(qbeta(0.025, 1:n, n:1)),
  upper = -log10(qbeta(0.975, 1:n, n:1))
)

# Q-Q plot with shaded confidence interval
ggplot(df, aes(x = expected, y = observed)) +
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




### Manhattan plot using ggplot

don <- gwas_results_snps %>% 
	group_by(CHROM) %>% # Compute chromosome size
	summarise(chrlen=max(POS)) %>%
	mutate(tot=cumsum(chrlen) - chrlen) %>% # Calculate cumulative position of each chromosome
	dplyr::select(-chrlen) %>%
	left_join(gwas_results_snps, ., by=c("CHROM"="CHROM")) %>% # Add this info to the initial dataset
	arrange(CHROM, POS) %>% # Add a cumulative position of each SNP
	mutate(BPcum=POS+tot)

axisdf = don %>%
	group_by(CHROM) %>%
	summarize(center=(max(BPcum) + min(BPcum) ) / 2 )

ggplot(don, aes(x=BPcum, y=-log10(P))) +
    
    # Show all points
    geom_point( aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("#4dac26", "#d01c8b"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )



### Evaluate PCs ---------------------------------------------------------
## Summarise Mean, Median, Max, Min
gwas_results_all <- read.table("/cluster/project2/DIVERGE/20250620_GWAS/GWAS/00_gwas_results.PHENO1.glm.logistic")

colnames(gwas_results_all) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "OR", "LOG(OR)_SE", "Z_STAT", "P")

gwas_results_all$CHROM <- as.numeric(gwas_results_all$CHROM)

gwas_results_all <- gwas_results_all %>% 
	filter(!is.na(CHROM))

gwas_results_all <- gwas_results_all %>%
	filter(
		!is.na(P),          # Remove NA
	    P > 0,              # Remove negative/zero
	    !is.infinite(P)     # Remove infinite
	) %>%
	mutate(
		P = if_else(P > 1, 1, P)  # Cap at 1 if any P > 1 exists
	)

gwas_results_all %>% 
	group_by(TEST) %>% 
	summarize(
		mean = mean(P),
	    median = median(P),
		maximum = max(P),
		minimum = min(P)
	)



### Null model with PCs ---------------------------------------------------------
## Load PCs
pc_results_all <- read.table("/cluster/project2/DIVERGE/20250620_GWAS/GWAS/covariates.txt")
colnames(pc_results_all) <- c("FID", "IID", "SEX", "PHENO1", "FID2", "IID2", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")

## Match them with phenotypic data
combined_df <- merge(phenotype, pc_results_all, by.x = "subject_id", , by.y = "IID")

## Run the null model
model <- glm(subject_type_logical ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, data = combined_df, family = binomial)

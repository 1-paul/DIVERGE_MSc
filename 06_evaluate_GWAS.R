library(dplyr)
library(ggplot2)
library(tidyr)
library(qqman)

gwas_results_snps <- read.table("/cluster/project2/DIVERGE/20250605_GWAS/snp_results_only.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

colnames(gwas_results_snps) <- c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "OR", "LOG(OR)_SE", "Z_STAT", "P")

gwas_results_snps$CHROM <- as.numeric(gwas_results_snps$CHROM)

gwas_results_snps <- gwas_results_snps %>% 
	filter(!is.na(CHROM))

gwas_results_snps <- gwas_results_snps %>%
	filter(
		!is.na(P),          # Remove NA
	    P > 0,              # Remove negative/zero
	    !is.infinite(P)     # Remove infinite
	) %>%
	mutate(
		P = if_else(P > 1, 1, P)  # Cap at 1 if any P > 1 exists
	)



### Manhattan & QQ plot using qqman package

manhattan(gwas_results_snps, chr="CHROM", bp="POS", snp="ID", p="P")

qq(gwas_results_snps$P, main = "Q-Q plot of GWAS p-values")


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



### Evaluate PCs
## Summarise Mean, Median, Max, Min

gwas_results_all <- read.table("/cluster/project2/DIVERGE/20250605_GWAS/gwas_results.PHENO1.glm.logistic")

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



### Null model with PCs
## Load PCs
pc_results_all <- read.table("/cluster/project2/DIVERGE/20250605_GWAS/covariates.txt")
colnames(pc_results_all) <- c("FID", "IID", "SEX", "PHENO1", "FID2", "IID2", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")

## Match them with phenotypic data
combined_df <- merge(phenotype, pc_results_all, by.x = "subject_id", , by.y = "IID")

## Run the null model
model <- glm(subject_type_logical ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, data = combined_df, family = binomial)

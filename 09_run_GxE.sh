#$ -S /bin/bash
#$ -l h_rt=06:0:0
#$ -l tmem=30G
#$ -l h_vmem=30G
#$ -cwd
#$ -j y
#$ -N run_GxE

#/share/apps/plink-1.90-beta-6.10/plink --bfile DIVERGE_QCed --logistic --out gwas_results

/share/apps/genomics/plink-2.0/bin/plink2 \
  --pfile /cluster/project2/DIVERGE/20250620_GWAS/GWAS/DIVERGE_QCed \
  --pheno covariates2.txt \
  --pheno-name PHENO1 \
  --covar covariates2.txt \
  --covar-name SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,early_domestic_issues \
  --parameters 1, 2-11, 21 \
  --glm interaction\
  --out gxe_results

# 1, refers to variants additive effect, 2-11 to the covariates, and 21 to the interaction between variants and early_domestic_issues



### Filter for the main effects and interaction
awk '$7 == "ADD" || $7 == "early_domestic_issues" || $7 == "ADDxearly_domestic_issues"' gxe_results.PHENO1.glm.logistic > gxe_results_snp_edm_snpxedm_results_only.txt

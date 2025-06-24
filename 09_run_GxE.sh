#$ -S /bin/bash
#$ -l h_rt=06:0:0
#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -cwd
#$ -j y
#$ -N run_GWAS

#/share/apps/plink-1.90-beta-6.10/plink --bfile DIVERGE_QCed --logistic --out gwas_results

/share/apps/genomics/plink-2.0/bin/plink2 \
  --pfile /cluster/project2/DIVERGE/20250620_GWAS/GWAS/DIVERGE_QCed \
  --pheno covariates2.txt \
  --pheno-name PHENO1 \
  --covar covariates2.txt \
  --covar-name early_domestic_issues,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8 \
  --parameters 1, 3-12, 13 \
  --glm interaction\
  --out gxe_results

# 1, refers to Additive effect, 3-12 to the covariates, and 13 to the interaction between variants and early_domestic_issues

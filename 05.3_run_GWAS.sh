# Credit to Munim Husain

#$ -S /bin/bash
#$ -l h_rt=06:0:0
#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -cwd
#$ -j y
#$ -N run_GWAS

#/share/apps/plink-1.90-beta-6.10/plink --bfile DIVERGE_QCed --logistic --out gwas_results

/share/apps/genomics/plink-2.0/bin/plink2 \
  --pfile DIVERGE_QCed \
  --pheno covariates.txt \
  --pheno-name PHENO1 \
  --covar covariates.txt \
  --covar-name SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8 \
  --glm \
  --out gwas_results [pbrandes@pchuckle 20250605_GWAS]
  

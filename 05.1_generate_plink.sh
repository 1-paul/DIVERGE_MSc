# Credit to Munim Husain

#$ -S /bin/bash
#$ -l h_rt=01:0:0
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -cwd
#$ -j y
#$ -N generate_plink_files

inputfile="../munim_workspace/QC_pipeline_20250514/00_plink_files/05_sex_check_removed_xhcrom_abnormality"
hardyweinberg_snps="../munim_workspace/QC_pipeline_20250514/08_hardy_weinberg_20250602/meta_snps.txt"
removed_relatives="../munim_workspace/QC_pipeline_20250514/00_reports/06_relatedness_relatives_kingrobust.txt"

##the current input file has gone through QC up to the end of the sex checks, the following command removes snps 
##that have failed a meta-analysis of hardy weinberg on different population clusters
/share/apps/plink-1.90-beta-6.10/plink --bfile "$inputfile" --exclude "$hardyweinberg_snps" --make-bed --out ./DIVERGE_QCed_initial --allow-extra-chr

###removing close relatives here
/share/apps/plink-1.90-beta-6.10/plink --bfile DIVERGE_QCed_initial --remove "$removed_relatives" --make-bed --out DIVERGE_QCed --allow-extra-chr


###### populate the phenotype column
export PATH=/share/apps/python-3.9.5-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.9.5-shared/lib:$LD_LIBRARY_PATH
/share/apps/python-3.9.5-shared/bin/python3.9 02_populate_phenotype.py

rm DIVERGE_QCed.fam
mv DIVERGE_QCed_.fam DIVERGE_QCed.fam

##convert to plink2
/share/apps/genomics/plink-2.0/bin/plink2 --bfile DIVERGE_QCed --make-pgen --out DIVERGE_QCed

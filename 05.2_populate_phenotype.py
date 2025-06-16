# Credit to Munim Husain

import csv

# filenames
fam_file = 'DIVERGE_QCed.fam'
pheno_file = '../munim_workspace/QC_pipeline_20250514/07_PCA/DIVERGE-DIVERGEGenotypeQC_DATA_LABELS_2025-04-23_1334.csv'
output_fam = 'DIVERGE_QCed_.fam'

# Read phenotype file keyed by IID only
pheno_dict = {}
with open(pheno_file, 'r') as pf:
    reader = csv.reader(pf)
    next(reader)  # skip header row
    for row in reader:
        iid, label = row[0], row[2].strip().lower() 
        pheno = '2' if label == 'proband' else '1'
        pheno_dict[iid] = pheno

# Update .fam file using only IID
with open(fam_file, 'r') as fin, open(output_fam, 'w') as fout:
    for line in fin:
        parts = line.strip().split()
        if len(parts) == 6:
            iid = parts[1]
            if iid in pheno_dict:
                parts[5] = pheno_dict[iid]
            fout.write(' '.join(parts) + '\n')

print(f"Updated .fam file written to {output_fam}")

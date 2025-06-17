# Credit to Munim Husain

import csv

# filenames
fam_file = 'DIVERGE_QCed.fam'
pheno_file = '/cluster/project2/DIVERGE//munim_workspace/QC_pipeline_20250514/07_PCA/DIVERGE-DIVERGEGenotypeQC_DATA_LABELS_2025-04-23_1334.csv'
domestic_issues_file = '/SAN/ugi/ukhls/Paul_MS_proj/DIVERGE-PaulBrandesProject_DATA_2025-05-28_1315.csv'
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


### New Version with added early_domestic_issues ----------------------------------------------------------------------------------------------------------------------


import csv

# File paths
fam_file = 'DIVERGE_QCed.fam'
pheno_file = '/cluster/project2/DIVERGE/munim_workspace/QC_pipeline_20250514/07_PCA/DIVERGE-DIVERGEGenotypeQC_DATA_LABELS_2025-04-23_1334.csv'
domestic_issues_file = '/SAN/ugi/ukhls/Paul_MS_proj/DIVERGE-PaulBrandesProject_DATA_2025-05-28_1315.csv'
output_fam = 'DIVERGE_QCed_new.fam'

# Read phenotype file (case/control status)
pheno_dict = {}
with open(pheno_file, 'r') as pf:
    reader = csv.reader(pf)
    next(reader)  # Skip header
    for row in reader:
        iid = row[0]
        label = row[2].strip().lower()
        pheno = '2' if label == 'proband' else '1'  # 2=case, 1=control
        pheno_dict[iid] = pheno

# Read domestic issues file (early_domestic_issues)
domestic_issues_dict = {}
with open(domestic_issues_file, 'r') as dif:
    reader = csv.reader(dif)
    header = next(reader)  # Get header
    # Find column index for subject_id and early_domestic_issues
    try:
        id_col = header.index('subject_id')
        edi_col = header.index('early_domestic_issues')
    except ValueError as e:
        raise ValueError(f"Column not found in domestic issues file: {e}")
    
    for row in reader:
        subject_id = row[id_col]
        early_domestic = row[edi_col]
        domestic_issues_dict[subject_id] = early_domestic

# Update .fam file (adding phenotype and domestic issues) using only IID
with open(fam_file, 'r') as fin, open(output_fam, 'w') as fout:
    for line in fin:
        parts = line.strip().split()
        if len(parts) >= 6:  # Ensure standard .fam format
            iid = parts[1]
            
            # Update phenotype (PHENO1)
            if iid in pheno_dict:
                parts[5] = pheno_dict[iid]  # 6th column is phenotype
            
            # Write the line (you could add domestic issues as a new column if needed)
            fout.write(' '.join(parts) + '\n')

print(f"Updated .fam file written to {output_fam}")

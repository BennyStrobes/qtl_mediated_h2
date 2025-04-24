#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-32:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)







tissue_info_file="$1"
global_mesc_lasso_weight_stem="$2"
mesc_run_name="$3"
genotype_stem="$4"
gwas_genotype_stem="$5"
ldsc_baseline_ld_annotation_stem="$6"
tglr_expression_score_dir="$7"

if false; then
source ~/.bash_profile
fi

date

# Loop through tissues
sed 1d $tissue_info_file | while read Tissue_name sample_size sample_names_file; do
	echo ${Tissue_name}

	# loop through chromsomes
	mesc_lasso_weight_stem=${global_mesc_lasso_weight_stem}${Tissue_name}"_"
	for CHR in {1..22}
	do
		python3 create_TGLR_gene_level_ld_scores.py ${CHR} ${Tissue_name} ${mesc_lasso_weight_stem} ${mesc_run_name} ${genotype_stem} ${gwas_genotype_stem} ${ldsc_baseline_ld_annotation_stem} ${tglr_expression_score_dir}
	done
done

date
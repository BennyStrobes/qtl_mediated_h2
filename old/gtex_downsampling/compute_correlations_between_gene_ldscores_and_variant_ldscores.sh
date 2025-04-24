#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4G                         # Memory total in MiB (for all cores)

if false; then
source ~/.bash_profile
fi


mesc_run_name="$1"
expression_score_dir="$2"
ldsc_baseline_ld_annotation_stem="$3"
gene_ldscore_variant_ld_score_corr_dir="$4"



python3 compute_correlations_between_gene_ldscores_and_variant_ldscores.py ${mesc_run_name} ${expression_score_dir} ${ldsc_baseline_ld_annotation_stem} ${gene_ldscore_variant_ld_score_corr_dir}

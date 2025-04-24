#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-22:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=160G                         # Memory total in MiB (for all cores)






trans_protein_models_dir="$1"
gwas_genotype_stem="$2"
ldsc_baseline_ld_annotation_stem="$3"
tglr_trans_expression_score_dir="$4"
eqtl_class="$5"
cis_snp_dir="$6"


source ~/.bash_profile
module load python/3.7.4


date

echo $eqtl_class

method="sbayesrc"
python3 create_TGLR_trans_gene_level_ld_scores.py $eqtl_class ${gwas_genotype_stem} $trans_protein_models_dir ${ldsc_baseline_ld_annotation_stem} ${tglr_trans_expression_score_dir} $method $cis_snp_dir





date
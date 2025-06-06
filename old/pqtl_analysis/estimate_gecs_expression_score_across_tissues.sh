#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=80G                         # Memory total in MiB (for all cores)





protein_models_dir="$1"
gwas_genotype_stem="$2"
ldsc_baseline_ld_annotation_stem="$3"
gecs_expression_score_dir="$4"
CHR="$5"

if false; then
source ~/.bash_profile
fi

date

tissue_name="UKBB_pqtl"

module load python/3.7.4

method="susie_inf_pmces"

echo $CHR
python3 create_gecs_gene_level_ld_scores.py ${CHR} $tissue_name ${gwas_genotype_stem} $protein_models_dir ${ldsc_baseline_ld_annotation_stem} ${gecs_expression_score_dir} $method



if false; then
for CHR in {1..22}
do
echo $CHR
python3 create_TGLR_gene_level_ld_scores.py ${CHR} $tissue_name ${gwas_genotype_stem} $protein_models_dir ${ldsc_baseline_ld_annotation_stem} ${tglr_expression_score_dir} $method
done
fi



date
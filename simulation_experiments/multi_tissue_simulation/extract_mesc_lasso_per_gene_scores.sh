#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-14:50                         # Runtime in D-HH:MM format
#SBATCH -p medium                          # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)


simulation_number="$1"
chrom_string="$2"
simulation_name_string="$3"
simulation_genotype_dir="$4"
mesc_expression_score_dir="$5"


module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate

date
if false; then
python3 extract_mesc_lasso_per_gene_scores.py $simulation_number $chrom_string $simulation_name_string $simulation_genotype_dir $mesc_expression_score_dir
fi
date

python3 extract_mesc_lasso_standardized_per_gene_scores.py $simulation_number $chrom_string $simulation_name_string $simulation_genotype_dir $mesc_expression_score_dir
date
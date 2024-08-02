#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-0:45                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5GB                         # Memory total in MiB (for all cores)




simulation_number="$1"
simulation_name_string="$2"
simulated_trait_dir="$3"
simulated_gwas_dir="$4"
simulation_genotype_dir="$5"
simulated_learned_gene_models_dir="$6"
n_gwas_individuals="$7"
trait_med_h2_inference_dir="$8"

if false; then
module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate
fi


eqtl_ss="300"
python3 trait_mediated_h2_inference.py $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_ss $trait_med_h2_inference_dir

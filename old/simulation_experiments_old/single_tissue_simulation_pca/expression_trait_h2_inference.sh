#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:45                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5GB                         # Memory total in MiB (for all cores)






simulation_number="$1"
simulation_name_string="$2"
simulated_learned_gene_models_dir="$3"
simulation_genotype_dir="$4"
eqtl_sample_size="$5"
trait_h2_inference_dir="$6"

module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate


python3 expression_trait_h2_inference.py $simulation_number $simulation_name_string $simulated_learned_gene_models_dir $simulation_genotype_dir $eqtl_sample_size $trait_h2_inference_dir
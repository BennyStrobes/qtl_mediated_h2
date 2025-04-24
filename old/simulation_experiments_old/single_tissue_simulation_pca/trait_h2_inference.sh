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
n_gwas_individuals="$6"
trait_h2_inference_dir="$7"

module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate


python3 trait_h2_inference.py $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $n_gwas_individuals $trait_h2_inference_dir

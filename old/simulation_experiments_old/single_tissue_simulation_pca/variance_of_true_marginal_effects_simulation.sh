#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:30                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=25GB                         # Memory total in MiB (for all cores)







true_marginal_effect_simulation_dir="${1}"
chrom_num="${2}"
simulation_genotype_dir="${3}"


module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate




python3 variance_of_true_marginal_effects_simulation.py $true_marginal_effect_simulation_dir $chrom_num $simulation_genotype_dir
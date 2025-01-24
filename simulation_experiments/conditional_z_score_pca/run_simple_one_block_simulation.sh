#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-9:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=30GB                         # Memory total in MiB (for all cores)






simulation_genotype_dir="${1}"
output_root_dir="${2}"
simulation_number="${3}"




module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate


python3 run_simple_one_block_simulation.py $simulation_genotype_dir $output_root_dir $simulation_number
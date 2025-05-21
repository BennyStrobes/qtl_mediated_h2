#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-9:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)




ld_dir="$1"
output_dir="${2}"
gwas_ss="${3}"
eqtl_ss="${4}"
nm_h2="${5}"
med_h2="${6}"
eqtl_h2="${7}"
sim_sparsity="${8}"


if false; then
module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate
fi

date
python3 run_simulation.py $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity
date
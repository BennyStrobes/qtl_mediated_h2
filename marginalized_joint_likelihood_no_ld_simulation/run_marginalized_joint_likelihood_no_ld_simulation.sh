#!/bin/bash
#SBATCH -t 0-15:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)
#SBATCH -c 1                           # Partition to run in



simulation_number="$1"
eqtl_ss="$2"
gwas_ss="$3"
ge_h2="$4"
h2="$5"
frac_med_h2="$6"
simulation_input_dir="$7"
simulation_output_dir="$8"


if false; then
module load gcc/6.2.0
module load python/3.6.0
source /n/groups/price/ben/environments/tensor_flow_cpu/bin/activate
fi


echo "Simulation "${simulation_number}

python3 run_marginalized_joint_likelihood_no_ld_simulation.py $simulation_number $eqtl_ss $gwas_ss $ge_h2 $h2 $frac_med_h2 $simulation_input_dir $simulation_output_dir
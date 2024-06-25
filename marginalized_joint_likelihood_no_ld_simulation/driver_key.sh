#!/bin/bash
#SBATCH -t 0-15:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)
#SBATCH -c 1                           # Partition to run in









################################
# Data directories
################################
simulation_input_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/marginalized_joint_likelihood_no_ld_simulation/input_data/"

simulation_output_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/marginalized_joint_likelihood_no_ld_simulation/output_data/"







simulation_number="1"
eqtl_ss="500"
gwas_ss="20000"
ge_h2="0.075"
h2="0.5"
frac_med_h2="0.1"
for simulation_number in $(seq 1 10); do 
	sh run_marginalized_joint_likelihood_no_ld_simulation.sh $simulation_number $eqtl_ss $gwas_ss $ge_h2 $h2 $frac_med_h2 $simulation_input_dir $simulation_output_dir
done
#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-7:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=33GB                         # Memory total in MiB (for all cores)





################
# Command line args
#################
mediated_h2_results_dir="$1"
visualize_med_h2_results_dir="$2"




# Simulation parameters
n_genes="100"
fraction_genes_cis_h2="1.0"
ge_h2="0.05"
per_element_heritability="0.001"
total_heritability="0.3"
fraction_expression_mediated_heritability="0.1"


date
simulation_name_string="n_genes_"${n_genes}"_frac_genes_h2_"${fraction_genes_cis_h2}"_ele_h2_"${per_element_heritability}"_tot_h2_"${total_heritability}"_frac_mediated_"${fraction_expression_mediated_heritability}"_ge_h2_"${ge_h2}"_"


python3 organize_simulation_results.py $mediated_h2_results_dir $visualize_med_h2_results_dir $simulation_name_string


module load R/3.5.1
Rscript visualize_simulation.R $visualize_med_h2_results_dir

#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)




################
# Command line args
#################
simulation_number="$1"
processed_genotype_data_dir="$2"
simulated_eqtl_data_dir="$3"
simulated_expression_data_dir="$4"
simulated_gene_models_dir="$5"
simulated_gwas_data_dir="$6"
mediated_h2_results_dir="$7"


# Simulation parameters
n_genes="100"
fraction_genes_cis_h2="0.5"
ge_h2="0.05"
per_element_heritability="0.001"
total_heritability="0.3"
fraction_expression_mediated_heritability="0.1"



simulation_name_string="simulation_"${simulation_number}"_n_genes_"${n_genes}"_frac_genes_h2_"${fraction_genes_cis_h2}"_ele_h2_"${per_element_heritability}"_tot_h2_"${total_heritability}"_frac_mediated_"${fraction_expression_mediated_heritability}"_ge_h2_"${ge_h2}"_"



##########################
# Part 1: Simulate variant-gene links and causal eqtl effects
###########################
if false; then
python3 simulate_variant_gene_links_and_causal_eqtl_effects.py $simulation_name_string $processed_genotype_data_dir $simulated_eqtl_data_dir $n_genes $fraction_genes_cis_h2 $ge_h2
fi

##########################
# Part 2: Simulate gene expression
###########################
eqtl_sample_size_arr=( "100" "200" "300" "500" "1000" "5000")
if false; then
for eqtl_ss in "${eqtl_sample_size_arr[@]}"
do
	python3 simulate_gene_expression.py $simulation_name_string $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_expression_data_dir $eqtl_ss
done
fi





##########################
# Part 3: Learn gene models
###########################
# NOTE: currently ld pruning has a large effect on estimated h2 (particularly at low eqtl sample sizes..)
eqtl_sample_size_arr=( "100" "200" "300" "500" "1000" "5000")
if false; then
for eqtl_ss in "${eqtl_sample_size_arr[@]}"
do
	python3 learn_gene_models.py $simulation_name_string $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_expression_data_dir $simulated_gene_models_dir $eqtl_ss
done
fi




##########################
# Part 4: simulate gwas data
###########################
if false; then
python3 simulate_gwas_data.py $simulation_name_string $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_gwas_data_dir $total_heritability $fraction_expression_mediated_heritability $per_element_heritability
fi


##########################
# Part 5: Compute GWAS summary statistics
###########################
if false; then
python3 compute_gwas_summary_statistics.py $simulation_name_string $processed_genotype_data_dir $simulated_gwas_data_dir
fi



##########################
# Part 6: Perform heritability inference with various methods
###########################
eqtl_ss="500"
python3 perform_med_h2_inference_shell.py $simulation_name_string $simulated_gwas_data_dir $simulated_gene_models_dir $eqtl_ss $mediated_h2_results_dir $processed_genotype_data_dir





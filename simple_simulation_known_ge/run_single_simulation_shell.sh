#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-0:40                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=8GB                         # Memory total in MiB (for all cores)




################
# Command line args
#################
simulation_number="$1"
processed_genotype_data_dir="$2"
simulated_eqtl_data_dir="$3"
simulated_gwas_data_dir="$4"
mediated_h2_results_dir="$5"


# Simulation parameters
n_genes="400"
fraction_genes_cis_h2="1.0"
ge_h2="0.05"
per_element_heritability="0.0005"
total_heritability="0.3"
fraction_expression_mediated_heritability="0.3"


date


# for sub_simulation_number in $(seq 1 5); do 

#for sub_simulation_number in $(seq 1 10); do 
for sub_simulation_number in $(seq 1 10); do 

simulation_name_string="simulation_"${simulation_number}"_"${sub_simulation_number}"_n_genes_"${n_genes}"_frac_genes_h2_"${fraction_genes_cis_h2}"_ele_h2_"${per_element_heritability}"_tot_h2_"${total_heritability}"_frac_mediated_"${fraction_expression_mediated_heritability}"_ge_h2_"${ge_h2}"_"
echo $simulation_name_string
##########################
# Part 1: Simulate variant-gene links and causal eqtl effects
###########################
####**
if false; then
python3 simulate_variant_gene_links_and_causal_eqtl_effects.py $simulation_name_string $processed_genotype_data_dir $simulated_eqtl_data_dir $n_genes $fraction_genes_cis_h2 $ge_h2



##########################
# Part 2: simulate gwas data
###########################
####**
# NEGATIVE SELECTION
python3 simulate_gwas_data_negative_selection.py $simulation_name_string $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_gwas_data_dir $total_heritability $fraction_expression_mediated_heritability $per_element_heritability
# NEUTRAL SELECTION
python3 simulate_gwas_data_neutral_selection.py $simulation_name_string $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_gwas_data_dir $total_heritability $fraction_expression_mediated_heritability $per_element_heritability


##########################
# Part 5: Compute GWAS summary statistics
###########################
selection="negative"
python3 compute_gwas_summary_statistics.py $simulation_name_string $processed_genotype_data_dir $simulated_gwas_data_dir $selection
selection="neutral"
python3 compute_gwas_summary_statistics.py $simulation_name_string $processed_genotype_data_dir $simulated_gwas_data_dir $selection



##########################
# Part 8: Perform mediated heritability ldsc inference 
###########################
####**
selection="negative"
python3 perform_med_h2_ldsc_inference.py $simulation_name_string $simulated_gwas_data_dir $mediated_h2_results_dir $processed_genotype_data_dir $simulated_eqtl_data_dir $selection
selection="neutral"
python3 perform_med_h2_ldsc_inference.py $simulation_name_string $simulated_gwas_data_dir $mediated_h2_results_dir $processed_genotype_data_dir $simulated_eqtl_data_dir $selection
fi


selection="neutral"
python3 perform_med_h2_ldsc_inference_with_simulated_noise.py $simulation_name_string $simulated_gwas_data_dir $mediated_h2_results_dir $processed_genotype_data_dir $simulated_eqtl_data_dir $selection





date

done




















if false; then
eqtl_sample_size_arr=( "100" "300" "1000")
for eqtl_ss in "${eqtl_sample_size_arr[@]}"
do
	python3 perform_med_h2_ldsc_inference.py $simulation_name_string $simulated_gwas_data_dir $simulated_gene_models_dir $eqtl_ss $mediated_h2_results_dir $processed_genotype_data_dir $simulated_expression_data_dir
done
fi



if false; then
eqtl_ss="100"
python3 perform_med_h2_ldsc_inference.py $simulation_name_string $simulated_gwas_data_dir $simulated_gene_models_dir $eqtl_ss $mediated_h2_results_dir $processed_genotype_data_dir $simulated_expression_data_dir
fi












##########################
# Part 6: Perform heritability inference with various methods
###########################
if false; then
python3 perform_h2_inference.py $simulation_name_string $simulated_gwas_data_dir $simulated_gene_models_dir $mediated_h2_results_dir $processed_genotype_data_dir
fi



##########################
# Part 7: Perform mediated heritability inference with various methods
###########################
if false; then
eqtl_ss="500"
python3 perform_med_h2_inference.py $simulation_name_string $simulated_gwas_data_dir $simulated_gene_models_dir $eqtl_ss $mediated_h2_results_dir $processed_genotype_data_dir $simulated_expression_data_dir

eqtl_ss="1000"
python3 perform_med_h2_inference.py $simulation_name_string $simulated_gwas_data_dir $simulated_gene_models_dir $eqtl_ss $mediated_h2_results_dir $processed_genotype_data_dir $simulated_expression_data_dir
fi


eqtl_sample_size_arr=( "100" "300" "500" "1000" "5000")
if false; then
for eqtl_ss in "${eqtl_sample_size_arr[@]}"
do
	python3 perform_med_h2_inference_shell.py $simulation_name_string $simulated_gwas_data_dir $simulated_gene_models_dir $eqtl_ss $mediated_h2_results_dir $processed_genotype_data_dir $simulated_expression_data_dir
done
fi


##########################
# Part 8: Perform mediated heritability ldsc inference 
###########################
####**
eqtl_sample_size_arr=( "100" "300" "1000")
if false; then
for eqtl_ss in "${eqtl_sample_size_arr[@]}"
do
	python3 perform_med_h2_ldsc_inference.py $simulation_name_string $simulated_gwas_data_dir $simulated_gene_models_dir $eqtl_ss $mediated_h2_results_dir $processed_genotype_data_dir $simulated_expression_data_dir
done
fi












eqtl_sample_size_arr=( "100" "200" "300" "500" "1000" "5000")
if false; then
for eqtl_ss in "${eqtl_sample_size_arr[@]}"
do
	python3 perform_med_h2_inference_shell.py $simulation_name_string $simulated_gwas_data_dir $simulated_gene_models_dir $eqtl_ss $mediated_h2_results_dir $processed_genotype_data_dir $simulated_expression_data_dir
done
fi













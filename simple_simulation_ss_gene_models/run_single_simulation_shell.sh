#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-4:40                         # Runtime in D-HH:MM format
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
simulated_expression_data_dir="$6"
simulated_gene_models_dir="$7"


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
for sub_simulation_number in $(seq 1 4); do 

simulation_name_string="simulation_"${simulation_number}"_"${sub_simulation_number}"_n_genes_"${n_genes}"_frac_genes_h2_"${fraction_genes_cis_h2}"_ele_h2_"${per_element_heritability}"_tot_h2_"${total_heritability}"_frac_mediated_"${fraction_expression_mediated_heritability}"_ge_h2_"${ge_h2}"_"
echo $simulation_name_string
##########################
# Part 1: Simulate variant-gene links and causal eqtl effects
###########################
####**
python3 simulate_variant_gene_links_and_causal_eqtl_effects.py $simulation_name_string $processed_genotype_data_dir $simulated_eqtl_data_dir $n_genes $fraction_genes_cis_h2 $ge_h2



##########################
# Part 2: Simulate gene expression
###########################
####**
eqtl_sample_size_arr=( "100" "300" "1000")
for eqtl_ss in "${eqtl_sample_size_arr[@]}"
do
	echo $eqtl_ss
	python3 simulate_gene_expression.py $simulation_name_string $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_expression_data_dir $eqtl_ss
done



##########################
# Part 3: Learn gene models
###########################
# NOTE: currently ld pruning has a large effect on estimated h2 (particularly at low eqtl sample sizes..)
####**
eqtl_sample_size_arr=( "100" "300" "1000")
for eqtl_ss in "${eqtl_sample_size_arr[@]}"
do
	echo $eqtl_ss
	python3 learn_gene_models.py $simulation_name_string $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_expression_data_dir $simulated_gene_models_dir $eqtl_ss $ge_h2
done


##########################
# Part 4: simulate gwas data
###########################
####**
# NEUTRAL SELECTION
python3 simulate_gwas_data_neutral_selection.py $simulation_name_string $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_gwas_data_dir $total_heritability $fraction_expression_mediated_heritability $per_element_heritability


##########################
# Part 5: Compute GWAS summary statistics
###########################
selection="neutral"
python3 compute_gwas_summary_statistics.py $simulation_name_string $processed_genotype_data_dir $simulated_gwas_data_dir $selection

##########################
# Part 6: Perform mediated heritability ldsc inference using various approaches
###########################
selection="neutral"
eqtl_sample_size_arr=("100" "300" "1000")
for eqtl_ss in "${eqtl_sample_size_arr[@]}"
do
	python3 perform_med_h2_ldsc_inference.py $simulation_name_string $simulated_gwas_data_dir $mediated_h2_results_dir $processed_genotype_data_dir $simulated_eqtl_data_dir $selection $eqtl_ss $simulated_gene_models_dir
done


date

done



































##########################
# Part 6: Perform mediated heritability ldsc inference with known true noise
###########################
selection="neutral"
eqtl_sample_size_arr=("100" "300" "1000")
if false; then
for eqtl_ss in "${eqtl_sample_size_arr[@]}"
do
	python3 perform_med_h2_ldsc_inference_with_true_noise.py $simulation_name_string $simulated_gwas_data_dir $mediated_h2_results_dir $processed_genotype_data_dir $simulated_eqtl_data_dir $selection $eqtl_ss $simulated_gene_models_dir
done
fi

if false; then
##########################
# Part 7: Perform mediated heritability ldsc inference with estimated noise
###########################
selection="neutral"
eqtl_sample_size_arr=("100" "300" "1000")
for eqtl_ss in "${eqtl_sample_size_arr[@]}"
do
	python3 perform_med_h2_ldsc_inference_with_est_noise.py $simulation_name_string $simulated_gwas_data_dir $mediated_h2_results_dir $processed_genotype_data_dir $simulated_eqtl_data_dir $selection $eqtl_ss $simulated_gene_models_dir
done
fi
















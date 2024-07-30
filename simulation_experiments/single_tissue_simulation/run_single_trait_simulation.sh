#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=11GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
per_element_heritability="$8"
total_heritability="${9}"
fraction_expression_mediated_heritability="${10}"
ge_h2="${11}"
simulated_gene_expression_dir="${12}"
simulated_learned_gene_models_dir="${13}"
simulated_trait_dir="${14}"
simulated_gwas_dir="${15}"
eqtl_architecture="${16}"

source /home/bes710/.bash_profile
module load R/4.0.1
echo "Simulation"$simulation_number
date
echo $simulation_name_string
echo $eqtl_architecture



#######################################################
# Step 1: Simulate gene expression causal eqtl effects
#######################################################
echo "Simulation Step 1"
python3 simulate_gene_expression.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir $ge_h2 $eqtl_architecture $per_element_heritability $total_heritability $fraction_expression_mediated_heritability


#######################################################
# Step 2: Simulate trait values
#######################################################
echo "Simulation Step 2"
python3 simulate_trait_values.py $simulation_number $chrom_num $cis_window $simulated_gene_expression_dir $simulation_name_string $processed_genotype_data_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $simulated_trait_dir $n_gwas_individuals $eqtl_architecture $ge_h2


source /home/bes710/.bash_profile

#######################################################
# Step 3: Run GWAS on simulated trait on only snps in TGFM windows.
#######################################################
echo "Simulation Step 3"
python3 run_gwas_on_simulated_trait.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $simulated_gwas_dir


#######################################################
# Step 4: Get variant LD-scores
#######################################################
echo "Simulation Step 4"
python3 get_variant_ld_scores.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_gwas_dir

#######################################################
# Step 5: Simulate gene expression and compute eqtl summary statistics
#######################################################
echo "Simulation Step 5"
eqtl_sample_size_arr=( "100" "300" "500")
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 fit_gene_models.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir $ge_h2 $eqtl_architecture $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $eqtl_sample_size
done




date




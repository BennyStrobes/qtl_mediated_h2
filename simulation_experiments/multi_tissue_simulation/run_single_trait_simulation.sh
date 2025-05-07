#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-32:30                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_string="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
processed_genotype_data_dir="$6"
per_element_heritability="$7"
total_heritability="${8}"
fraction_expression_mediated_heritability="${9}"
ge_h2="${10}"
simulated_gene_expression_dir="${11}"
simulated_learned_gene_models_dir="${12}"
simulated_trait_dir="${13}"
simulated_gwas_dir="${14}"
eqtl_architecture="${15}"
n_tissues="${16}"
alt_simulated_learned_gene_models_dir="${17}"

echo "Simulation"$simulation_number
date
echo $simulation_name_string
echo $eqtl_architecture

module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate

date

#######################################################
# Step 1: Simulate gene expression causal eqtl effects
#######################################################
echo "Simulation Step 1"
python3 simulate_causal_eqtl_effect_sizes.py $simulation_number $chrom_string $cis_window $simulated_gene_expression_dir $simulation_name_string $processed_genotype_data_dir $ge_h2 $eqtl_architecture $n_tissues



gene_trait_architecture="linear"
gene_trait_architecture="stdExpr"
gene_trait_architecture_arr=( "linear" "stdExpr" ) 
for gene_trait_architecture in "${gene_trait_architecture_arr[@]}"
do
	echo $gene_trait_architecture
	gwas_simulation_name_string=${simulation_name_string}"_gt_arch_"${gene_trait_architecture}


	#######################################################
	# Step 2: Simulate trait values
	#######################################################
	echo "Simulation Step 2"
	python3 simulate_trait_values.py $simulation_number $chrom_string $cis_window $simulated_gene_expression_dir $simulation_name_string $gwas_simulation_name_string $processed_genotype_data_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $simulated_trait_dir $n_gwas_individuals $eqtl_architecture $ge_h2 $gene_trait_architecture

	#######################################################
	# Step 3: Run GWAS on simulated trait on only snps in TGFM windows.
	#######################################################
	echo "Simulation Step 3"
	python3 run_gwas_on_simulated_trait.py $simulation_number $chrom_string $gwas_simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $simulated_gwas_dir
done




#######################################################
# Step 4: Simulate gene expression and compute eqtl summary statistics
#######################################################
# TODO: Include up to 10K (will probably need to load in chunks of genotype at a time instead of whole thing at once)
# TODO: Add MESC CODE
echo "Simulation Step 4"
eqtl_sample_size_arr=( "100" "200" "300" "1000")
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	echo $eqtl_sample_size
	python3 simulate_gene_expression_and_compute_eqtl_ss.py $simulation_number $chrom_string $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir $eqtl_sample_size $n_tissues
done




















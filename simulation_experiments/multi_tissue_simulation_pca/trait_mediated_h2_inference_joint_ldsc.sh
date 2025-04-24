#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:10                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=8GB                         # Memory total in MiB (for all cores)




simulation_number="$1"
simulation_name_string="$2"
simulated_trait_dir="$3"
simulated_gwas_dir="$4"
simulation_genotype_dir="$5"
simulated_learned_gene_models_dir="$6"
n_gwas_individuals="$7"
trait_med_h2_inference_dir="$8"
simulated_gene_expression_dir="$9"






echo $simulation_number
echo $simulation_name_string



eqtl_sample_size_arr=( "100" "300" "1000" "10000")
date
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	echo $eqtl_sample_size

	date
	python3 trait_mediated_h2_inference_joint_ldsc.py $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $simulated_gene_expression_dir
	date
done
date




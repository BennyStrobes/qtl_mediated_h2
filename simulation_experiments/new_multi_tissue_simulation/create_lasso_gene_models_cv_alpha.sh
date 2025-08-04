#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-30:30                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)


simulation_number="${1}"
chrom_string="${2}"
cis_window="${3}"
simulation_name_string="${4}"
simulation_genotype_dir="${5}"
simulated_learned_gene_models_dir="${6}"
lasso_gene_models_dir="${7}"
calibrated_mesc_code_dir="${8}"
cis_window="${9}"
tissue_num="${10}"



echo "Simulation"$simulation_number
date
echo $simulation_name_string

chromosome_file="/n/groups/price/ben/qtl_mediated_h2/simulation_chromosomes.txt"


source ~/.bash_profile







eqtl_sample_size_arr=( "100" "200" "300" "1000")

sample_split_arr=( "full" "replicate1" "replicate2")
sample_split_arr=( "replicate1" "replicate2")


for eqtl_ss in "${eqtl_sample_size_arr[@]}"; do
for sample_split in "${sample_split_arr[@]}"; do


	# Print current iteration
	echo ${eqtl_ss}"_"${sample_split}"_"${tissue_num}"_CV"
	
	# Expression matrix for this iteration
	expr_file=${simulated_learned_gene_models_dir}${simulation_name_string}"_tissue"${tissue_num}"_"${eqtl_ss}"_"${sample_split}"_expression.txt"

	# Stem of bgen file
	bgen_file_stem=${simulation_genotype_dir}"simulated_eqtl_"${eqtl_ss}"_data_"

	output_file=${lasso_gene_models_dir}${simulation_name_string}"_tissue"${tissue_num}"_"${eqtl_ss}"_"${sample_split}"_lasso_CV_est_causal_effects.txt"

	date
	python3 ${calibrated_mesc_code_dir}fit_lasso_gene_models.py \
		--expr ${expr_file} \
		--bgen ${bgen_file_stem} \
		--chromosome-file ${chromosome_file} \
		--cis-window ${cis_window} \
		--cross-val-alpha \
		--output ${output_file}

done 
done



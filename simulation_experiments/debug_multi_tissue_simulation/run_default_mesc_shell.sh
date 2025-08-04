#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_string="$2"
simulation_name_string="$3"
simulation_genotype_dir="$4"
mesc_code_dir="$5"
mesc_expression_score_dir="$6"
mesc_processed_input_dir="$7"
mesc_results_dir="$8"
simulated_gwas_dir="$9"


module load python/2.7.12
conda activate mesc


sample_split="full"
eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")


##########################
# Linear version
##########################
gt_arch="linear"
orig_sumstats_file=${simulated_gwas_dir}${simulation_name_string}"_gt_arch_"${gt_arch}"_simualated_gwas_results.txt"
new_sumstats_file=${mesc_processed_input_dir}${simulation_name_string}"_gt_arch_"${gt_arch}"_simualated_gwas_results.txt"
python3 reformat_sumstats_for_mesc.py $orig_sumstats_file $new_sumstats_file


for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 preprocess_data_for_mesc_run.py $simulation_number $chrom_string $simulation_name_string $simulation_genotype_dir $mesc_expression_score_dir $eqtl_sample_size $sample_split $mesc_processed_input_dir $gt_arch


	expscore_prefix=${mesc_processed_input_dir}${simulation_name_string}"_"${gt_arch}"_"${eqtl_sample_size}"_"${sample_split}
	output_prefix=${mesc_results_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${sample_split}"_gt_arch_"${gt_arch}"_baselineLD"
	python ${mesc_code_dir}run_mesc.py --h2med $new_sumstats_file --exp-chr $expscore_prefix --out $output_prefix --chisq-max "100000000000000"
	output_prefix=${mesc_results_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${sample_split}"_gt_arch_"${gt_arch}"_genotypeIntercept"
	python ${mesc_code_dir}run_mesc.py --h2med $new_sumstats_file --exp-chr $expscore_prefix --out $output_prefix --chisq-max "100000000000000" --ref-ld-keep-annot "baseL2"
done



##########################
# stdExpr version
##########################
gt_arch="stdExpr"
orig_sumstats_file=${simulated_gwas_dir}${simulation_name_string}"_gt_arch_"${gt_arch}"_simualated_gwas_results.txt"
new_sumstats_file=${mesc_processed_input_dir}${simulation_name_string}"_gt_arch_"${gt_arch}"_simualated_gwas_results.txt"
python3 reformat_sumstats_for_mesc.py $orig_sumstats_file $new_sumstats_file


for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 preprocess_data_for_mesc_run.py $simulation_number $chrom_string $simulation_name_string $simulation_genotype_dir $mesc_expression_score_dir $eqtl_sample_size $sample_split $mesc_processed_input_dir $gt_arch

	expscore_prefix=${mesc_processed_input_dir}${simulation_name_string}"_"${gt_arch}"_"${eqtl_sample_size}"_"${sample_split}
	output_prefix=${mesc_results_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${sample_split}"_gt_arch_"${gt_arch}"_baselineLD"
	python ${mesc_code_dir}run_mesc.py --h2med $new_sumstats_file --exp-chr $expscore_prefix --out $output_prefix --chisq-max "100000000000000"
	output_prefix=${mesc_results_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${sample_split}"_gt_arch_"${gt_arch}"_genotypeIntercept"
	python ${mesc_code_dir}run_mesc.py --h2med $new_sumstats_file --exp-chr $expscore_prefix --out $output_prefix --chisq-max "100000000000000" --ref-ld-keep-annot "baseL2"

done
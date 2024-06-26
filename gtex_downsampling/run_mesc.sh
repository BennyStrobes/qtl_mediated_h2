#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)

source ~/.bash_profile





tissue_info_file="$1"
mesc_run_name="$2"
gwas_genotype_stem="$3"
expression_score_dir="$4"
sumstat_dir="$5"
non_redundent_summary_statistics_file="$6"
ldsc_genotype_intercept_annotation_stem="$7"
ldsc_baseline_ld_annotation_stem="$8"
sldsc_weights_stem="$9"
frq_file_stem="${10}"
mesc_results_dir="${11}"
mesc_code_dir="${12}"








if false; then
#############################
# Meta-analyze mesc gene scores across tissues
#############################
source ~/.bash_profile
python3 extract_input_list_of_tissues_to_meta_analyze_by_mesc.py $tissue_info_file $mesc_run_name $expression_score_dir


module load python/2.7.12
conda activate mesc
# Output prefix for meta-analyzed results
meta_analyzed_expr_score_prefix=${expression_score_dir}${mesc_run_name}"_meta_analyzed_scores"

# run meta-analysis in each chromosome
# Loop through chromosomes
for chr_num in {1..22}
do
	echo ${chr_num}
	list_to_meta_analyze=${expression_score_dir}${mesc_run_name}"_meta_analysis_list_"${chr_num}".txt"
	python ${mesc_code_dir}meta_analyze_weights.py --input-prefixes ${list_to_meta_analyze} --bfile ${gwas_genotype_stem}${chr_num} --chr ${chr_num} --out ${meta_analyzed_expr_score_prefix}
done



#############################
# Run MESC in each trait
#############################
sed 1d $non_redundent_summary_statistics_file | while read gwas_trait_name h2 h2_se h2_z intercept; do
	echo $gwas_trait_name
	# Summary stats
	gwas_file_name=${sumstat_dir}${gwas_trait_name}".sumstats"

	# MESC with genotype intercept
	mesc_output_stem=${mesc_results_dir}${gwas_trait_name}"_"${mesc_run_name}"_genotype_intercept"
	python ${mesc_code_dir}run_mesc.py --h2med ${gwas_file_name} --exp-chr ${meta_analyzed_expr_score_prefix} --ref-ld-chr ${ldsc_genotype_intercept_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --frqfile-chr ${frq_file_stem} --out ${mesc_output_stem}

	# MESC with baselineLD_no_qtl
	mesc_output_stem=${mesc_results_dir}${gwas_trait_name}"_"${mesc_run_name}"_baselineLD_no_qtl"
	python ${mesc_code_dir}run_mesc.py --h2med ${gwas_file_name} --exp-chr ${meta_analyzed_expr_score_prefix} --ref-ld-chr ${ldsc_baseline_ld_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --frqfile-chr ${frq_file_stem} --out ${mesc_output_stem}
done
fi


#############################
# Organize MESC results across traits
#############################
source ~/.bash_profile


mesc_run_identifier=${mesc_run_name}"_genotype_intercept"
xt_output_file=${mesc_results_dir}"cross_trait_med_h2_summary_"$mesc_run_identifier".txt"
python3 organize_mesc_results_across_traits.py $mesc_run_identifier $non_redundent_summary_statistics_file ${mesc_results_dir} ${xt_output_file}
echo $xt_output_file

mesc_run_identifier=${mesc_run_name}"_baselineLD_no_qtl"
xt_output_file=${mesc_results_dir}"cross_trait_med_h2_summary_"$mesc_run_identifier".txt"
python3 organize_mesc_results_across_traits.py $mesc_run_identifier $non_redundent_summary_statistics_file ${mesc_results_dir} ${xt_output_file}
echo $xt_output_file


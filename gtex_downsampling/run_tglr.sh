#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)


tissue_info_file="${1}"
mesc_run_name="${2}"
ldsc_genotype_intercept_annotation_stem="${3}"
ldsc_baseline_ld_annotation_stem="${4}"
tglr_expression_score_dir="${5}"
tglr_results_dir="${6}"
tglr_code_dir="$7"
sumstat_dir="$8"
non_redundent_summary_statistics_file="$9"
sldsc_weights_stem="${10}"


source ~/.bash_profile


########################################
# Organize TGLR expression scores across tissues
########################################
expression_score_annotation_stem=${tglr_expression_score_dir}"tglr_expression_scores_"${mesc_run_name}"."
python3 organize_tglr_expression_scores_across_tissues.py $tissue_info_file $mesc_run_name $ldsc_genotype_intercept_annotation_stem $tglr_expression_score_dir $expression_score_annotation_stem




########################################
# Run TGLR for each trait
########################################
sed 1d $non_redundent_summary_statistics_file | while read gwas_trait_name h2 h2_se h2_z intercept; do

	trait_file=${sumstat_dir}${gwas_trait_name}".sumstats"

	# TGLR with baseline ld
	tglr_output_stem=${tglr_results_dir}${gwas_trait_name}"_"${mesc_run_name}"_baselineLD_no_qtl_tglr"
	module load python/2.7.12
	python ${tglr_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_baseline_ld_annotation_stem}","${expression_score_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --print-delete-vals --print-coefficients --out ${tglr_output_stem}
	source ~/.bash_profile
	python3 organize_tglr_results.py $tglr_output_stem ${ldsc_baseline_ld_annotation_stem} ${expression_score_annotation_stem}

	# TGLR with genotype intercept
	tglr_output_stem=${tglr_results_dir}${gwas_trait_name}"_"${mesc_run_name}"_genotype_intercept_tglr"
	module load python/2.7.12
	python ${tglr_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_genotype_intercept_annotation_stem}","${expression_score_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --print-delete-vals --print-coefficients --out ${tglr_output_stem}
	source ~/.bash_profile
	python3 organize_tglr_results.py $tglr_output_stem ${ldsc_genotype_intercept_annotation_stem} ${expression_score_annotation_stem}
done



########################################
# Organize TGLR results across trait runs
########################################

mesc_run_identifier=${mesc_run_name}"_genotype_intercept"
xt_output_file=${tglr_results_dir}"cross_trait_med_h2_summary_"$mesc_run_identifier".txt"
python3 organize_tglr_results_results_across_traits.py $mesc_run_identifier $non_redundent_summary_statistics_file ${tglr_results_dir} ${xt_output_file}


mesc_run_identifier=${mesc_run_name}"_baselineLD_no_qtl"
xt_output_file=${tglr_results_dir}"cross_trait_med_h2_summary_"$mesc_run_identifier".txt"
python3 organize_tglr_results_results_across_traits.py $mesc_run_identifier $non_redundent_summary_statistics_file ${tglr_results_dir} ${xt_output_file}






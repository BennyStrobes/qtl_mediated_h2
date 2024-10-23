#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)


ldsc_genotype_intercept_annotation_stem="${1}"
ldsc_baseline_ld_annotation_stem="${2}"
tglr_expression_score_dir="${3}"
tglr_results_dir="${4}"
tglr_code_dir="$5"
sumstat_dir="$6"
non_redundent_summary_statistics_file="$7"
sldsc_weights_stem="${8}"

source ~/.bash_profile


########################################
# Organize TGLR expression scores across tissues
########################################

method="sbayesrc"
cis_trans_expression_score_annotation_stem=${tglr_expression_score_dir}"tglr_cis_trans_expression_scores_"${method}"."
cis_expression_score_annotation_stem=${tglr_expression_score_dir}"tglr_cis_expression_scores_"${method}"."
trans_expression_score_annotation_stem=${tglr_expression_score_dir}"tglr_trans_expression_scores_"${method}"."
if false; then
python3 organize_tglr_trans_expression_scores_across_tissues.py $ldsc_genotype_intercept_annotation_stem $tglr_expression_score_dir $method $cis_trans_expression_score_annotation_stem $cis_expression_score_annotation_stem $trans_expression_score_annotation_stem
fi

########################################
# Run TGLR for each trait
# Cis + trans
########################################
if false; then
sed 1d $non_redundent_summary_statistics_file | while read gwas_trait_name h2 h2_se h2_z intercept; do

	echo $gwas_trait_name
	trait_file=${sumstat_dir}${gwas_trait_name}".sumstats"

	# TGLR with baseline ld
	tglr_output_stem=${tglr_results_dir}${gwas_trait_name}"_"${method}"_cis_trans_expression_score_baselineLD_no_qtl_tglr"
	module load python/2.7.12
	python ${tglr_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_baseline_ld_annotation_stem}","${cis_trans_expression_score_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --print-delete-vals --print-coefficients --out ${tglr_output_stem}
	source ~/.bash_profile
	module load python/3.7.4
	python3 organize_tglr_results.py $tglr_output_stem ${ldsc_baseline_ld_annotation_stem} ${cis_trans_expression_score_annotation_stem}

	# TGLR with genotype intercept
	tglr_output_stem=${tglr_results_dir}${gwas_trait_name}"_"${method}"_cis_trans_expression_score_genotype_intercept_tglr"
	module load python/2.7.12
	python ${tglr_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_genotype_intercept_annotation_stem}","${cis_trans_expression_score_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --print-delete-vals --print-coefficients --out ${tglr_output_stem}
	source ~/.bash_profile
	module load python/3.7.4
	python3 organize_tglr_results.py $tglr_output_stem ${ldsc_genotype_intercept_annotation_stem} ${cis_trans_expression_score_annotation_stem}
done
fi

########################################
# Run TGLR for each trait
# Cis only
########################################
if false; then
sed 1d $non_redundent_summary_statistics_file | while read gwas_trait_name h2 h2_se h2_z intercept; do

	echo $gwas_trait_name
	trait_file=${sumstat_dir}${gwas_trait_name}".sumstats"

	# TGLR with baseline ld
	tglr_output_stem=${tglr_results_dir}${gwas_trait_name}"_"${method}"_cis_expression_score_baselineLD_no_qtl_tglr"
	module load python/2.7.12
	python ${tglr_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_baseline_ld_annotation_stem}","${cis_expression_score_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --print-delete-vals --print-coefficients --out ${tglr_output_stem}
	source ~/.bash_profile
	module load python/3.7.4
	python3 organize_tglr_results.py $tglr_output_stem ${ldsc_baseline_ld_annotation_stem} ${cis_expression_score_annotation_stem}

	# TGLR with genotype intercept
	tglr_output_stem=${tglr_results_dir}${gwas_trait_name}"_"${method}"_cis_expression_score_genotype_intercept_tglr"
	module load python/2.7.12
	python ${tglr_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_genotype_intercept_annotation_stem}","${cis_expression_score_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --print-delete-vals --print-coefficients --out ${tglr_output_stem}
	source ~/.bash_profile
	module load python/3.7.4
	python3 organize_tglr_results.py $tglr_output_stem ${ldsc_genotype_intercept_annotation_stem} ${cis_expression_score_annotation_stem}
done
fi

########################################
# Run TGLR for each trait
# Trans only
########################################
sed 1d $non_redundent_summary_statistics_file | while read gwas_trait_name h2 h2_se h2_z intercept; do

	echo $gwas_trait_name
	trait_file=${sumstat_dir}${gwas_trait_name}".sumstats"

	# TGLR with baseline ld
	tglr_output_stem=${tglr_results_dir}${gwas_trait_name}"_"${method}"_trans_expression_score_baselineLD_no_qtl_tglr"
	module load python/2.7.12
	python ${tglr_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_baseline_ld_annotation_stem}","${trans_expression_score_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --print-delete-vals --print-coefficients --out ${tglr_output_stem}
	source ~/.bash_profile
	module load python/3.7.4
	python3 organize_tglr_results.py $tglr_output_stem ${ldsc_baseline_ld_annotation_stem} ${trans_expression_score_annotation_stem}

	# TGLR with genotype intercept
	tglr_output_stem=${tglr_results_dir}${gwas_trait_name}"_"${method}"_trans_expression_score_genotype_intercept_tglr"
	module load python/2.7.12
	python ${tglr_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_genotype_intercept_annotation_stem}","${trans_expression_score_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --print-delete-vals --print-coefficients --out ${tglr_output_stem}
	source ~/.bash_profile
	module load python/3.7.4
	python3 organize_tglr_results.py $tglr_output_stem ${ldsc_genotype_intercept_annotation_stem} ${trans_expression_score_annotation_stem}
done


########################################
# Organize TGLR results across trait runs
########################################
if false; then
method="susie_inf_pmces"
mesc_run_identifier=${method}"_genotype_intercept"
xt_output_file=${tglr_results_dir}"cross_trait_med_h2_summary_"$mesc_run_identifier".txt"
python3 organize_tglr_results_results_across_traits.py $mesc_run_identifier $non_redundent_summary_statistics_file ${tglr_results_dir} ${xt_output_file}

mesc_run_identifier=${method}"_baselineLD_no_qtl"
xt_output_file=${tglr_results_dir}"cross_trait_med_h2_summary_"$mesc_run_identifier".txt"
python3 organize_tglr_results_results_across_traits.py $mesc_run_identifier $non_redundent_summary_statistics_file ${tglr_results_dir} ${xt_output_file}
fi







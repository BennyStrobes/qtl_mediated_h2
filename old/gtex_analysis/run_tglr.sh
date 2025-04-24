#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)




qtl_version="${1}"
ldsc_genotype_intercept_annotation_stem="${2}"
ldsc_baseline_ld_annotation_stem="${3}"
tglr_expression_score_dir="${4}"
tglr_results_dir="${5}"
sumstat_dir="${6}"
non_redundent_summary_statistics_file="${7}"
sldsc_weights_stem="${8}"
tglr_code_dir="${9}"


########################################
# Organize TGLR expression scores across tissues
########################################
expression_score_annotation_stem=${tglr_expression_score_dir}"tglr_gene_ld_scores_"${qtl_version}"."
if false; then
python3 organize_tglr_expression_scores.py $qtl_version $ldsc_genotype_intercept_annotation_stem $tglr_expression_score_dir $expression_score_annotation_stem
fi


########################################
# Run TGLR for each trait
########################################
gwas_trait_name="UKB_460K.lung_FEV1FVCzSMOKE"
gwas_trait_name="UKB_460K.blood_WHITE_COUNT"
gwas_trait_name="PASS_BrainstemVol_Satizabal2019"
gwas_trait_name="UKB_460K.biochemistry_TotalProtein"
gwas_trait_name="UKB_460K.biochemistry_IGF1"
gwas_trait_name="UKB_460K.bp_DIASTOLICadjMEDz"

	if false; then
	trait_file=${sumstat_dir}${gwas_trait_name}".sumstats"
	# TGLR with baseline ld
	tglr_output_stem=${tglr_results_dir}${gwas_trait_name}"_"${qtl_version}"_baselineLD_no_qtl_tglr"
	module load python/2.7.12
	python ${tglr_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_baseline_ld_annotation_stem}","${expression_score_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --print-delete-vals --print-coefficients --out ${tglr_output_stem}
	source ~/.bash_profile
	python3 organize_tglr_results.py $tglr_output_stem ${ldsc_baseline_ld_annotation_stem} ${expression_score_annotation_stem}
	fi

if false; then
sed 1d $non_redundent_summary_statistics_file | while read gwas_trait_name h2 h2_se h2_z intercept; do
	trait_file=${sumstat_dir}${gwas_trait_name}".sumstats"
	# TGLR with baseline ld
	tglr_output_stem=${tglr_results_dir}${gwas_trait_name}"_"${qtl_version}"_baselineLD_no_qtl_tglr"
	module load python/2.7.12
	python ${tglr_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_baseline_ld_annotation_stem}","${expression_score_annotation_stem} --w-ld-chr ${sldsc_weights_stem} --print-delete-vals --print-coefficients --out ${tglr_output_stem}
	source ~/.bash_profile
	python3 organize_tglr_results.py $tglr_output_stem ${ldsc_baseline_ld_annotation_stem} ${expression_score_annotation_stem}
done
fi



python3 meta_analyze_tglr_results_across_gwases.py $non_redundent_summary_statistics_file ${tglr_results_dir} $qtl_version





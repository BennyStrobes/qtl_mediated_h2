#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:20                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_string="$2"
simulation_name_string="$3"
simulation_genotype_dir="$4"
simulated_gene_expression_dir="$5"
calibrated_mesc_code_dir="$6"
estimated_cis_snp_h2_dir="$7"


echo $simulation_number

causal_eqtl_summary_file=${simulated_gene_expression_dir}${simulation_name_string}"_causal_eqtl_effect_summary.txt"

# Stem of bgen file
bgen_file_stem=${simulation_genotype_dir}"simulated_reference_genotype_data_"

study_names_file=${simulated_gene_expression_dir}${simulation_name_string}"_all_study_names.txt"

python3 extract_true_gene_causal_effects.py $causal_eqtl_summary_file $bgen_file_stem $simulated_gene_expression_dir $simulation_name_string $study_names_file


hm3_rsid_file_stem="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/genotype_processing/hm3_rsids_chr"

ldsc_annotation_file_stem="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/baselineLD."

gene_summary_file_stem=${simulation_genotype_dir}"gene_positions_chr"
chromosome_file="/n/groups/price/ben/qtl_mediated_h2/simulation_chromosomes.txt"

echo "START"
#STANDARDIZED FORM
python3 ${calibrated_mesc_code_dir}compute_predicted_marginal_summary_statistics.py \
	--study-names-file $study_names_file \
	--chromosome-file ${chromosome_file} \
	--hm3-file-stem $hm3_rsid_file_stem \
	--ldsc-annotation-file-stem ${ldsc_annotation_file_stem} \
	--bgen ${bgen_file_stem} \
	--gene-summary-file-stem ${gene_summary_file_stem} \
	--gene-summary-file-suffix ".txt" \
	--standardize \
	--dont-standardize-snp-effects \
	--output-dir $simulated_gene_expression_dir

echo "MID"

#UNSTANDARDIZED FORM
python3 ${calibrated_mesc_code_dir}compute_predicted_marginal_summary_statistics.py \
	--study-names-file $study_names_file \
	--chromosome-file ${chromosome_file} \
	--hm3-file-stem $hm3_rsid_file_stem \
	--ldsc-annotation-file-stem ${ldsc_annotation_file_stem} \
	--bgen ${bgen_file_stem} \
	--gene-summary-file-stem ${gene_summary_file_stem} \
	--gene-summary-file-suffix ".txt" \
	--dont-standardize-snp-effects \
	--output-dir $simulated_gene_expression_dir


if false; then

eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_100_true_est_cis_snp_h2_summary.txt"
# In sample STANDARDIZED FORM
python3 ${calibrated_mesc_code_dir}compute_predicted_marginal_summary_statistics.py \
	--study-names-file $study_names_file \
	--chromosome-file ${chromosome_file} \
	--hm3-file-stem $hm3_rsid_file_stem \
	--ldsc-annotation-file-stem ${ldsc_annotation_file_stem} \
	--bgen ${bgen_file_stem} \
	--gene-summary-file-stem ${gene_summary_file_stem} \
	--gene-summary-file-suffix ".txt" \
	--standardize \
	--gene-cis-h2-file $eqtl_cis_snp_h2_summary_file \
	--dont-standardize-snp-effects \
	--output-dir $simulated_gene_expression_dir
fi


echo "END"
python3 create_true_causal_eqtl_summary_files.py $simulation_number $simulation_name_string $simulated_gene_expression_dir









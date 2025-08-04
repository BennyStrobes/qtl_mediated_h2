#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-12:10                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_string="$2"
simulation_name_string="$3"
simulation_genotype_dir="$4"
simulated_gene_expression_dir="$5"
calibrated_mesc_code_dir="$6"

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




python3 create_true_causal_eqtl_summary_files.py $simulation_number $simulation_name_string $simulated_gene_expression_dir









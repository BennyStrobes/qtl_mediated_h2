#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:10                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)





simulation_number="${1}"
chrom_string="${2}"
simulation_name_string="${3}"
simulation_genotype_dir="${4}"
mesc_expression_score_dir="${5}"
lasso_gene_models_dir="${6}"
calibrated_mesc_code_dir="${7}"
alpha_0="${8}"



module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate




study_names_file=${lasso_gene_models_dir}${simulation_name_string}"_lasso_alpha_"${alpha_0}"_all_study_names.txt"
python3 extract_study_names_for_lasso_gene_ldscores.py $simulation_name_string $alpha_0 $lasso_gene_models_dir $study_names_file




# Names of chromosomes to run anlaysis on
chromosome_file="/n/groups/price/ben/qtl_mediated_h2/simulation_chromosomes.txt"


# Stem of bgen file
bgen_file_stem=${simulation_genotype_dir}"simulated_reference_genotype_data_"

hm3_rsid_file_stem="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/genotype_processing/hm3_rsids_chr"

ldsc_annotation_file_stem="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/baselineLD."

gene_summary_file_stem=${simulation_genotype_dir}"gene_level_ld_chr"


python3 ${calibrated_mesc_code_dir}compute_predicted_marginal_summary_statistics.py \
	--study-names-file $study_names_file \
	--chromosome-file ${chromosome_file} \
	--hm3-file-stem $hm3_rsid_file_stem \
	--ldsc-annotation-file-stem ${ldsc_annotation_file_stem} \
	--bgen ${bgen_file_stem} \
	--gene-summary-file-stem ${gene_summary_file_stem} \
	--gene-summary-file-suffix "_bins_10_summary_file.txt" \
	--standardize \
	--output-dir $lasso_gene_models_dir

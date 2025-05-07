#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-15:10                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
simulation_name_string="$2"
simulated_trait_dir="$3"
simulated_gwas_dir="$4"
simulation_genotype_dir="$5"
simulated_learned_gene_models_dir="$6"
n_gwas_individuals="$7"
trait_med_h2_inference_dir="$8"
simulated_gene_expression_dir="$9"
chrom_string="${10}"
joint_ldsc_code_dir="${11}"

echo "SIMULATION"$simulation_number


chrom_nums_file="/n/groups/price/ben/joint_ldsc_temp_data/simulation_chromosomes.txt"



gt_architecture="linear"
non_med_anno="genotype_intercept"


eqtl_sample_size_arr=( "100" "200" "300" "1000")
eqtl_snp_representation_arr=( "pca_90" "bins_10" "bins_20")
if false; then
for eqtl_snp_representation in "${eqtl_snp_representation_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do

	gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${gt_architecture}"_simualated_gwas_results.txt"
	variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	gene_ldscore_filestem=${simulation_genotype_dir}"gene_level_ld_chr"
	gene_ldscore_filesuffix="_"${eqtl_snp_representation}"_summary_file.txt"
	eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_eqtl_sumstats_xt_summary.txt"
	variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${gt_architecture}"_"${eqtl_sample_size}

	python3 ${joint_ldsc_code_dir}run_joint_ldsc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--gene-ldscore-filestem ${gene_ldscore_filestem} \
		--gene-ldscore-filesuffix ${gene_ldscore_filesuffix} \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-summary-file ${eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${gt_architecture} \
		--output-stem ${output_stem}

done
done
fi






eqtl_snp_representation="pca_90"
gt_architecture="linear"
non_med_anno="full_anno"
non_med_anno="genotype_intercept"
eqtl_sample_size="100"


	gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${gt_architecture}"_simualated_gwas_results.txt"
	variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	gene_ldscore_filestem=${simulation_genotype_dir}"gene_level_ld_chr"
	gene_ldscore_filesuffix="_"${eqtl_snp_representation}"_summary_file.txt"
	eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_eqtl_sumstats_xt_summary.txt"
	variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${gt_architecture}"_"${eqtl_sample_size}"_TMP"

	python3 ${joint_ldsc_code_dir}run_joint_ldsc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--gene-ldscore-filestem ${gene_ldscore_filestem} \
		--gene-ldscore-filesuffix ${gene_ldscore_filesuffix} \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-summary-file ${eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${gt_architecture} \
		--output-stem ${output_stem}




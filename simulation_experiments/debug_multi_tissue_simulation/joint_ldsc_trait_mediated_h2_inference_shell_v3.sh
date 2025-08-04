#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-13:10                         # Runtime in D-HH:MM format
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



if false; then

simulated_gt_architecture="linear"
inference_gt_architecture="linear"
inference_approach="mle"
inference_approach="gibbs"
echo $simulated_gt_architecture"_"$inference_gt_architecture


non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")

eqtl_snp_representation_arr=( "pca_90" "bins_10" "bins_20")
eqtl_snp_representation_arr=( "bins_10")

for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_snp_representation in "${eqtl_snp_representation_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do

	gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
	variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	gene_ldscore_filestem=${simulation_genotype_dir}"gene_level_ld_chr"
	gene_ldscore_filesuffix="_"${eqtl_snp_representation}"_summary_file.txt"
	eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_eqtl_sumstats_xt_summary.txt"
	variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${inference_approach}"_"${eqtl_sample_size}
	echo $non_med_anno"_"$eqtl_snp_representation"_"$eqtl_sample_size
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
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach $inference_approach \
		--output-stem ${output_stem}
done
done
done
fi



if false; then
simulated_gt_architecture="linear"
inference_gt_architecture="stdExpr"
inference_approach="mle"
inference_approach="gibbs"
echo $simulated_gt_architecture"_"$inference_gt_architecture


non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")

eqtl_snp_representation_arr=( "pca_90" "bins_10" "bins_20")
eqtl_snp_representation_arr=( "bines_10" "bins_20")

for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_snp_representation in "${eqtl_snp_representation_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do

	gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
	variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	gene_ldscore_filestem=${simulation_genotype_dir}"gene_level_ld_chr"
	gene_ldscore_filesuffix="_"${eqtl_snp_representation}"_summary_file.txt"
	eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_eqtl_sumstats_xt_summary.txt"
	variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${inference_approach}"_"${eqtl_sample_size}
	echo $non_med_anno"_"$eqtl_snp_representation"_"$eqtl_sample_size
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
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach $inference_approach \
		--output-stem ${output_stem}


done
done
done
fi



simulated_gt_architecture="linear"
inference_gt_architecture="random_bins"
inference_approach="gibbs"
inference_approach="mle"
echo $simulated_gt_architecture"_"$inference_gt_architecture


non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")

eqtl_snp_representation_arr=( "pca_90" "bins_10" "bins_20")
eqtl_snp_representation_arr=( "bins_10" "bins_20")

for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_snp_representation in "${eqtl_snp_representation_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do

	gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
	variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	gene_ldscore_filestem=${simulation_genotype_dir}"gene_level_ld_chr"
	gene_ldscore_filesuffix="_"${eqtl_snp_representation}"_summary_file.txt"
	eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_eqtl_sumstats_xt_summary.txt"
	variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_sample_size}
	echo $non_med_anno"_"$eqtl_snp_representation"_"$eqtl_sample_size
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
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach $inference_approach \
		--output-stem ${output_stem}


done
done
done





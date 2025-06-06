#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-4:50                         # Runtime in D-HH:MM format
#SBATCH -p short                          # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



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
calibrated_mesc_code_dir="${11}"
mesc_expression_score_dir="${12}"

echo "SIMULATION"$simulation_number


chrom_nums_file="/n/groups/price/ben/joint_ldsc_temp_data/simulation_chromosomes.txt"





simulated_gt_architecture="linear"
inference_gt_architecture="linear"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
gene_ld_score_type="squared_marginal_sumstats"
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")

eqtl_snp_representation_arr=( "pca_90" "pca_95" "bins_10" "bins_20")
eqtl_snp_representation_arr=( "bins_20")

eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "1000")

step1_gene_ldscores_arr=("mescLassoCorr" "mescLasso" "mescLassoPlusMarginalSS" "MarginalSS")

inference_approach_arr=("2SLS" "JIVE")


for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_snp_representation in "${eqtl_snp_representation_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do
for step1_gene_ldscores in "${step1_gene_ldscores_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"

gene_ldscore_filestem=${simulation_genotype_dir}"gene_level_ld_chr"
gene_ldscore_filesuffix="_"${eqtl_snp_representation}"_summary_file.txt"
eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${gene_ld_score_type}"_"${eqtl_beta_sq_filter}"_"${step1_gene_ldscores}"_"${inference_approach}"_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--gene-ldscore-filestem ${gene_ldscore_filestem} \
		--gene-ldscore-filesuffix ${gene_ldscore_filesuffix} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-summary-file ${eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--gene-ldscore-type ${gene_ld_score_type} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--step1-gene-ldscores ${step1_gene_ldscores} \
		--output-stem ${output_stem}
	date


done
done
done
done
done
done











if false; then

simulated_gt_architecture="stdExpr"
inference_gt_architecture="stdExpr"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
gene_ld_score_type="squared_marginal_sumstats"
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")

eqtl_snp_representation_arr=( "pca_90" "pca_95" "bins_10" "bins_20")
eqtl_snp_representation_arr=( "bins_20")

eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0")



for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_snp_representation in "${eqtl_snp_representation_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"


gene_ldscore_filestem=${simulation_genotype_dir}"gene_level_ld_chr"
gene_ldscore_filesuffix="_"${eqtl_snp_representation}"_summary_file.txt"
eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${gene_ld_score_type}"_"${eqtl_beta_sq_filter}"_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--gene-ldscore-filestem ${gene_ldscore_filestem} \
		--gene-ldscore-filesuffix ${gene_ldscore_filesuffix} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-summary-file ${eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--gene-ldscore-type ${gene_ld_score_type} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--output-stem ${output_stem}
	date



done
done
done
done




simulated_gt_architecture="linear"
inference_gt_architecture="random_bins"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
gene_ld_score_type="squared_marginal_sumstats"
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")

eqtl_snp_representation_arr=( "pca_90" "pca_95" "bins_10" "bins_20")
eqtl_snp_representation_arr=( "bins_20")

eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0")


for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_snp_representation in "${eqtl_snp_representation_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"


gene_ldscore_filestem=${simulation_genotype_dir}"gene_level_ld_chr"
gene_ldscore_filesuffix="_"${eqtl_snp_representation}"_summary_file.txt"
eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${gene_ld_score_type}"_"${eqtl_beta_sq_filter}"_"${eqtl_sample_size}"_weighted"

	step1_regression_method="non_zero_snps"
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${gene_ld_score_type}"_"${eqtl_beta_sq_filter}"_"${step1_regression_method}"_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--gene-ldscore-filestem ${gene_ldscore_filestem} \
		--gene-ldscore-filesuffix ${gene_ldscore_filesuffix} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-summary-file ${eqtl_summary_file} \
		--step1-regression-method $step1_regression_method \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--gene-ldscore-type ${gene_ld_score_type} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--output-stem ${output_stem}
	date

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${gene_ld_score_type}"_"${eqtl_beta_sq_filter}"_"${step1_regression_method}"_"${eqtl_sample_size}"_unweighted"
	echo $output_stem
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--gene-ldscore-filestem ${gene_ldscore_filestem} \
		--gene-ldscore-filesuffix ${gene_ldscore_filesuffix} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--step1-regression-method $step1_regression_method \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-summary-file ${eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--gene-ldscore-type ${gene_ld_score_type} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--output-stem ${output_stem}
	date

done
done
done
done



simulated_gt_architecture="linear"
inference_gt_architecture="stdExpr"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
gene_ld_score_type="squared_marginal_sumstats"
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")

eqtl_snp_representation_arr=( "pca_90" "pca_95" "bins_10" "bins_20")
eqtl_snp_representation_arr=( "bins_20")

eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0")


for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_snp_representation in "${eqtl_snp_representation_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"


gene_ldscore_filestem=${simulation_genotype_dir}"gene_level_ld_chr"
gene_ldscore_filesuffix="_"${eqtl_snp_representation}"_summary_file.txt"
eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${gene_ld_score_type}"_"${eqtl_beta_sq_filter}"_"${eqtl_sample_size}"_weighted"

	step1_regression_method="non_zero_snps"
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${gene_ld_score_type}"_"${eqtl_beta_sq_filter}"_"${step1_regression_method}"_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--gene-ldscore-filestem ${gene_ldscore_filestem} \
		--gene-ldscore-filesuffix ${gene_ldscore_filesuffix} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-summary-file ${eqtl_summary_file} \
		--step1-regression-method $step1_regression_method \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--gene-ldscore-type ${gene_ld_score_type} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--output-stem ${output_stem}
	date

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${gene_ld_score_type}"_"${eqtl_beta_sq_filter}"_"${step1_regression_method}"_"${eqtl_sample_size}"_unweighted"
	echo $output_stem
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--gene-ldscore-filestem ${gene_ldscore_filestem} \
		--gene-ldscore-filesuffix ${gene_ldscore_filesuffix} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--step1-regression-method $step1_regression_method \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-summary-file ${eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--gene-ldscore-type ${gene_ld_score_type} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--output-stem ${output_stem}
	date


done
done
done
done
fi















































if false; then




simulated_gt_architecture="linear"
inference_gt_architecture="linear"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
gene_ld_score_type="ashr_style_pred"
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")

eqtl_snp_representation_arr=( "pca_90" "pca_95" "bins_10" "bins_20")
eqtl_snp_representation_arr=( "bins_20")

eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" )

for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_snp_representation in "${eqtl_snp_representation_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"

gene_ldscore_filestem=${simulation_genotype_dir}"gene_level_ld_chr"
gene_ldscore_filesuffix="_"${eqtl_snp_representation}"_summary_file.txt"
eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${gene_ld_score_type}"_"${eqtl_beta_sq_filter}"_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--gene-ldscore-filestem ${gene_ldscore_filestem} \
		--gene-ldscore-filesuffix ${gene_ldscore_filesuffix} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-summary-file ${eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--gene-ldscore-type ${gene_ld_score_type} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--output-stem ${output_stem}
	date


done
done
done
done



fi









if false; then
#####################
gene_ld_score_type="ldsc_style_pred"

non_med_anno_arr=("genotype_intercept" "full_anno")
eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")
eqtl_snp_representation_arr=( "pca_90" "pca_95" "bins_10" "bins_20")


non_med_anno_arr=("genotype_intercept" "full_anno")
eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")
eqtl_snp_representation_arr=( "pca_90" "pca_95" "bins_10" "bins_20")

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
	
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${gene_ld_score_type}"_"${eqtl_sample_size}
	date
	echo $output_stem
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
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
		--gene-ldscore-type ${gene_ld_score_type} \
		--output-stem ${output_stem}
	date

done
done
done
fi



#####################
gene_ld_score_type="ashr_style_pred"


if false; then
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")
eqtl_sample_size_arr=( "100" "200")
eqtl_snp_representation_arr=("bins_20")

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
	
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${eqtl_snp_representation}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${gene_ld_score_type}"_"${eqtl_sample_size}
	date
	echo $output_stem
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
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
		--gene-ldscore-type ${gene_ld_score_type} \
		--output-stem ${output_stem}
	date

done
done
done
fi






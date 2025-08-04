#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:20                         # Runtime in D-HH:MM format
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
estimated_cis_snp_h2_dir="${13}"
lasso_gene_models_dir="${14}"

source ~/.bash_profile


echo "SIMULATION"$simulation_number


chrom_nums_file="/n/groups/price/ben/joint_ldsc_temp_data/simulation_chromosomes.txt"

####################
#**#**#*#
# Step 1: Create Cis snp heritability files
eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	echo $eqtl_sample_size
	eqtl_snp_representation="bins_10"
	gene_ldscore_filestem=${simulation_genotype_dir}"gene_level_ld_chr"
	gene_ldscore_filesuffix="_"${eqtl_snp_representation}"_summary_file.txt"
	eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
	variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"

	python3 create_cis_snp_heritability_files.py ${gene_ldscore_filestem} ${gene_ldscore_filesuffix} $eqtl_summary_file $chrom_nums_file $variant_stdev_filestem $eqtl_sample_size $estimated_cis_snp_h2_dir $simulation_name_string $mesc_expression_score_dir $simulated_trait_dir
done





simulated_gt_architecture="linear"
inference_gt_architecture="linear"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "greml" "ldsc" )
validation_data_gene_ldscores_type="MarginalSS"

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")

eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "1000")

training_data_gene_ldscores_type_arr=("MarginalSS" "mescLassoCorr" "mescLasso" "mescLassoPlusMarginalSS")
training_data_gene_ldscores_type_arr=("mescLasso" "MarginalSS" "mescLassoPlusMarginalSS")

inference_approach_arr=("2SLS" "JIVE")
inference_approach_arr=("2SLS")

if false; then
for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do
for training_data_gene_ldscores_type in "${training_data_gene_ldscores_type_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${cis_snp_h2_method}"_"${training_data_gene_ldscores_type}"_"${inference_approach}"_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--output-stem ${output_stem}
	date


done
done
done
done
done
done
fi



simulated_gt_architecture="linear"
inference_gt_architecture="linear"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "greml" )
validation_data_gene_ldscores_type="MarginalSS"

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")

eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0")

training_data_gene_ldscores_type_arr=("MarginalSS" "mescLassoCorr" "mescLasso" "mescLassoPlusMarginalSS")
training_data_gene_ldscores_type_arr=("mesc_lasso_corr")

inference_approach_arr=("2SLS" "JIVE")
inference_approach_arr=("2SLS")
if false; then
for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do
for training_data_gene_ldscores_type in "${training_data_gene_ldscores_type_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

marginal_eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$mesc_expression_score_dir${simulation_name_string}"_"${eqtl_sample_size}"_"${training_data_gene_ldscores_type}"_replicate_eqtl_sumstats_xt_summary.txt"
variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${cis_snp_h2_method}"_"${training_data_gene_ldscores_type}"_"${inference_approach}"_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${marginal_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--output-stem ${output_stem}
	date


done
done
done
done
done
done
fi


if false; then
simulated_gt_architecture="stdExpr"
inference_gt_architecture="stdExpr"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "greml" )
validation_data_gene_ldscores_type="MarginalSS"

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")

training_data_gene_ldscores_type_arr=("MarginalSS" "mescLassoCorr" "mescLasso" "mescLassoPlusMarginalSS")
training_data_gene_ldscores_type_arr=("mesc_lasso_corr_standardized")

inference_approach_arr=("2SLS" "JIVE")

for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do
for training_data_gene_ldscores_type in "${training_data_gene_ldscores_type_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

marginal_eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$mesc_expression_score_dir${simulation_name_string}"_"${eqtl_sample_size}"_"${training_data_gene_ldscores_type}"_replicate_eqtl_sumstats_xt_summary.txt"
variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${cis_snp_h2_method}"_"${training_data_gene_ldscores_type}"_"${inference_approach}"_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${marginal_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--standardize-gene-validation-scores \
		--output-stem ${output_stem}
	date


done
done
done
done
done
done
fi




if false; then
simulated_gt_architecture="stdExpr"
inference_gt_architecture="stdExpr"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "greml" )
validation_data_gene_ldscores_type="MarginalSS"

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")
eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")
eqtl_beta_sq_filter_arr=( "100" "1000.0")

training_data_gene_ldscores_type_arr=("MarginalSS" "mescLassoCorr" "mescLasso" "mescLassoPlusMarginalSS")
training_data_gene_ldscores_type_arr=("mesc_lasso_corr_standardized")
cis_h2_threshold_arr=(".005" ".01" ".025")

inference_approach_arr=("2SLS" "JIVE")

for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do
for training_data_gene_ldscores_type in "${training_data_gene_ldscores_type_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do
for cis_h2_threshold in "${cis_h2_threshold_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

marginal_eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$mesc_expression_score_dir${simulation_name_string}"_"${eqtl_sample_size}"_"${training_data_gene_ldscores_type}"_replicate_eqtl_sumstats_xt_summary.txt"
variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${cis_snp_h2_method}"_"${training_data_gene_ldscores_type}"_"${inference_approach}"_"${cis_h2_threshold}"_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${marginal_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--standardize-gene-validation-scores \
		--cis-h2-thresh ${cis_h2_threshold} \
		--output-stem ${output_stem}
	date


done
done
done
done
done
done
done
fi








simulated_gt_architecture="stdExpr"
inference_gt_architecture="stdExpr"

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "greml" )

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")
eqtl_beta_sq_filter_arr=( "1000.0")

training_data_gene_ldscores_type_arr=("MarginalSS" "mescLassoCorr" "mescLasso" "mescLassoPlusMarginalSS")
training_data_gene_ldscores_type_arr=("mesc_lasso_corr_standardized")
validation_data_gene_ldscores_type="mesc_lasso_corr_standardized"

inference_approach_arr=("2SLS" "JIVE")
inference_approach_arr=("2SLS" )

if false; then
for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do
for training_data_gene_ldscores_type in "${training_data_gene_ldscores_type_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do



gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

pred_eqtl_summary_file=$mesc_expression_score_dir${simulation_name_string}"_"${eqtl_sample_size}"_"${training_data_gene_ldscores_type}"_replicate_eqtl_sumstats_xt_summary.txt"
variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${cis_snp_h2_method}"_"${training_data_gene_ldscores_type}"_"${validation_data_gene_ldscores_type}"_"${inference_approach}"_cc3_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${pred_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--covariance-calibration \
		--output-stem ${output_stem}
	date


done
done
done
done
done
done
fi



if false; then
simulated_gt_architecture="stdExpr"
inference_gt_architecture="stdExpr"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "greml" )

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")
eqtl_beta_sq_filter_arr=( "1000.0")

training_data_gene_ldscores_type_arr=("MarginalSS" "mescLassoCorr" "mescLasso" "mescLassoPlusMarginalSS")
training_data_gene_ldscores_type_arr=("lasso_CV_corr_standardized")
validation_data_gene_ldscores_type="lasso_CV_corr_standardized"


inference_approach_arr=("2SLS" "JIVE")
inference_approach_arr=("2SLS" )

for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do
for training_data_gene_ldscores_type in "${training_data_gene_ldscores_type_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do



gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

pred_eqtl_summary_file=$lasso_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_"${training_data_gene_ldscores_type}"_replicate_eqtl_sumstats_xt_summary.txt"
variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${cis_snp_h2_method}"_"${training_data_gene_ldscores_type}"_"${validation_data_gene_ldscores_type}"_"${inference_approach}"_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${pred_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--output-stem ${output_stem}


	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${cis_snp_h2_method}"_"${training_data_gene_ldscores_type}"_"${validation_data_gene_ldscores_type}"_"${inference_approach}"_cc_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${pred_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--covariance-calibration \
		--output-stem ${output_stem}
	date


done
done
done
done
done
done
fi







##**#*#*#*#*#
if false; then
simulated_gt_architecture="stdExpr"
inference_gt_architecture="stdExpr"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "ldsc" "avgChisq" "greml" )
validation_data_gene_ldscores_type="MarginalSS"

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")
eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")
eqtl_beta_sq_filter_arr=( "100.0" "1000.0")

training_data_gene_ldscores_type_arr=("MarginalSS" "mescLassoCorr" "mescLasso" "mescLassoPlusMarginalSS")
training_data_gene_ldscores_type_arr=("lasso_CV_corr_standardized")
cis_h2_threshold_arr=(".01" ".025")

inference_approach_arr=("2SLS")

for non_med_anno in "${non_med_anno_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do
for training_data_gene_ldscores_type in "${training_data_gene_ldscores_type_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do
for cis_h2_threshold in "${cis_h2_threshold_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

marginal_eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$mesc_expression_score_dir${simulation_name_string}"_"${eqtl_sample_size}"_"${training_data_gene_ldscores_type}"_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$lasso_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_"${training_data_gene_ldscores_type}"_replicate_eqtl_sumstats_xt_summary.txt"

variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${cis_snp_h2_method}"_"${training_data_gene_ldscores_type}"_"${inference_approach}"_"${cis_h2_threshold}"_"${eqtl_sample_size}"_weighted"
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${marginal_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--standardize-gene-validation-scores \
		--cis-h2-thresh ${cis_h2_threshold} \
		--output-stem ${output_stem}
	date


done
done
done
done
done
done
done
fi





#######################
# Binned approach
# Weighted
if false; then
simulated_gt_architecture="stdExpr"
inference_gt_architecture="linear"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "true" "avgChisq" "ldsc" "greml" )

validation_data_gene_ldscores_type="MarginalSS"

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")
eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")
eqtl_beta_sq_filter_arr=( "1000.0")

training_data_gene_ldscores_type_arr=("MarginalSS" "mescLassoCorr" "mescLasso" "mescLassoPlusMarginalSS")
training_data_gene_ldscores_type_arr=("lasso_CV_corr" "lasso_CV_corr_standardized")

inference_approach_arr=("seperately_binned_2SLS" "binned_2SLS")



for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do
for non_med_anno in "${non_med_anno_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do
for training_data_gene_ldscores_type in "${training_data_gene_ldscores_type_arr[@]}"
do





gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

marginal_eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$mesc_expression_score_dir${simulation_name_string}"_"${eqtl_sample_size}"_"${training_data_gene_ldscores_type}"_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$lasso_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_"${training_data_gene_ldscores_type}"_replicate_eqtl_sumstats_xt_summary.txt"

variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${cis_snp_h2_method}"_"${training_data_gene_ldscores_type}"_"${validation_data_gene_ldscores_type}"_"${inference_approach}"_"${eqtl_sample_size}"_weighted"
	echo $output_stem
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_binned_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${marginal_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--output-stem ${output_stem}
	date


done
done
done
done
done
done
fi






#######################
# Binned approach
# unWeighted
if false; then
simulated_gt_architecture="stdExpr"
inference_gt_architecture="linear"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "true" "avgChisq" "ldsc" "greml" )
cis_snp_h2_method_arr=( "true" "avgChisq" "greml" )

validation_data_gene_ldscores_type="MarginalSS"

eqtl_sample_size_arr=( "100" "200" "300" "1000" "100-1000")
eqtl_sample_size_arr=( "1000" "100" "200" "300" "100-1000")


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")
eqtl_beta_sq_filter_arr=( "1000.0")

training_data_gene_ldscores_type_arr=("MarginalSS" "mescLassoCorr" "mescLasso" "mescLassoPlusMarginalSS")
training_data_gene_ldscores_type_arr=("lasso_CV_corr" "lasso_CV_corr_standardized")
training_data_gene_ldscores_type_arr=( "lasso_CV_corr_standardized")


inference_approach_arr=("seperately_binned_2SLS" "binned_2SLS")
n_bins_arr=("4" "6")



for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do
for non_med_anno in "${non_med_anno_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do
for training_data_gene_ldscores_type in "${training_data_gene_ldscores_type_arr[@]}"
do
for n_bins in "${n_bins_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

marginal_eqtl_summary_file=$simulated_learned_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$mesc_expression_score_dir${simulation_name_string}"_"${eqtl_sample_size}"_"${training_data_gene_ldscores_type}"_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$lasso_gene_models_dir${simulation_name_string}"_"${eqtl_sample_size}"_"${training_data_gene_ldscores_type}"_replicate_eqtl_sumstats_xt_summary.txt"

variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	

	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${cis_snp_h2_method}"_"${training_data_gene_ldscores_type}"_"${validation_data_gene_ldscores_type}"_"${inference_approach}"_"${n_bins}"_"${eqtl_sample_size}
	echo $output_stem
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_binned_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${marginal_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--n-expr-cis-h2-bins ${n_bins} \
		--output-stem ${output_stem}
	date
done
done
done
done
done
done
done
fi






# INFINITE SS
# Inverse variance arch

simulated_gt_architecture="stdExpr"
inference_gt_architecture="stdExpr"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "true" "avgChisq" "ldsc" "greml" )
cis_snp_h2_method_arr=( "true" )

validation_data_gene_ldscores_type="lasso_CV_corr_standardized"
training_data_gene_ldscores_type="lasso_CV_corr_standardized"


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")
eqtl_beta_sq_filter_arr=( "1000.0")


inference_approach_arr=( "2SLS")

eqtl_sample_size="100"
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do
for non_med_anno in "${non_med_anno_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

marginal_eqtl_summary_file=$simulated_gene_expression_dir${simulation_name_string}"_true_eqtls_corr_standardized_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$simulated_gene_expression_dir${simulation_name_string}"_true_eqtls_corr_standardized_replicate_eqtl_sumstats_xt_summary.txt"

echo $pred_eqtl_summary_file

variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${inference_approach}"_INF_SS_weighted"
	echo $output_stem
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--chromosome-file ${chrom_nums_file} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${marginal_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--output-stem ${output_stem}
	date
done
done
done
done


#######################
# Binned approach
# Infinite SAMPLE SIZE
simulated_gt_architecture="stdExpr"
inference_gt_architecture="linear"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "true" "avgChisq" "ldsc" "greml" )
cis_snp_h2_method_arr=( "true" )

validation_data_gene_ldscores_type="lasso_CV_corr"
training_data_gene_ldscores_type="lasso_CV_corr"


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")
eqtl_beta_sq_filter_arr=( "1000.0")


inference_approach_arr=( "binned_2SLS")
n_bins_arr=("1" "2" "3" "4" "5" "6" "7")


eqtl_sample_size="100"
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do
for non_med_anno in "${non_med_anno_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do
for n_bins in "${n_bins_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

marginal_eqtl_summary_file=$simulated_gene_expression_dir${simulation_name_string}"_true_eqtls_corr_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$simulated_gene_expression_dir${simulation_name_string}"_true_eqtls_corr_replicate_eqtl_sumstats_xt_summary.txt"

echo $pred_eqtl_summary_file

variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${cis_snp_h2_method}"_"${inference_approach}"_"${n_bins}"_INF_SS_weighted"
	echo $output_stem
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_binned_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--chromosome-file ${chrom_nums_file} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${marginal_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--n-expr-cis-h2-bins ${n_bins} \
		--output-stem ${output_stem}
	date
done
done
done
done
done



simulated_gt_architecture="stdExpr"
inference_gt_architecture="stdExpr"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "true" "avgChisq" "ldsc" "greml" )
cis_snp_h2_method_arr=( "true" )

validation_data_gene_ldscores_type="lasso_CV_corr_standardized"
training_data_gene_ldscores_type="lasso_CV_corr_standardized"


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")
eqtl_beta_sq_filter_arr=( "1000.0")


inference_approach_arr=( "2SLS")

if false; then
eqtl_sample_size="100"
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do
for non_med_anno in "${non_med_anno_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

marginal_eqtl_summary_file=$simulated_gene_expression_dir${simulation_name_string}"_true_eqtls_corr_standardized_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$simulated_gene_expression_dir${simulation_name_string}"_true_eqtls_corr_standardized_replicate_eqtl_sumstats_xt_summary.txt"

echo $pred_eqtl_summary_file

variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${inference_approach}"_INF_SS"
	echo $output_stem
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${marginal_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--output-stem ${output_stem}
	date
done
done
done
done
fi











simulated_gt_architecture="stdExpr"
inference_gt_architecture="linear"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "true" "avgChisq" "ldsc" "greml" )
cis_snp_h2_method_arr=( "true" )

validation_data_gene_ldscores_type="lasso_CV_corr"
training_data_gene_ldscores_type="lasso_CV_corr"


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")
eqtl_beta_sq_filter_arr=( "1000.0")


inference_approach_arr=( "2SLS")


eqtl_sample_size="100"
if false; then
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do
for non_med_anno in "${non_med_anno_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

marginal_eqtl_summary_file=$simulated_gene_expression_dir${simulation_name_string}"_true_eqtls_corr_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$simulated_gene_expression_dir${simulation_name_string}"_true_eqtls_corr_replicate_eqtl_sumstats_xt_summary.txt"

echo $pred_eqtl_summary_file

variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${inference_approach}"_INF_SS_weighted"
	echo $output_stem
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--chromosome-file ${chrom_nums_file} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${marginal_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--output-stem ${output_stem}
	date
done
done
done
done
fi



simulated_gt_architecture="linear"
inference_gt_architecture="linear"
echo $simulated_gt_architecture"_"$inference_gt_architecture

#####################
non_med_anno_arr=("genotype_intercept" "full_anno")
non_med_anno_arr=("genotype_intercept")
cis_snp_h2_method_arr=( "true" "avgChisq" "ldsc" "greml" )
cis_snp_h2_method_arr=( "true" )

validation_data_gene_ldscores_type="lasso_CV_corr"
training_data_gene_ldscores_type="lasso_CV_corr"


eqtl_beta_sq_filter_arr=( "0.5" "1.0" "2.0" "5.0" "10.0" "100.0" "10000.0")
eqtl_beta_sq_filter_arr=( "100.0" "500.0" "1000.0")
eqtl_beta_sq_filter_arr=( "1000.0")


inference_approach_arr=( "2SLS")


eqtl_sample_size="100"
if false; then
for cis_snp_h2_method in "${cis_snp_h2_method_arr[@]}"
do
for non_med_anno in "${non_med_anno_arr[@]}"
do
for inference_approach in "${inference_approach_arr[@]}"
do
for eqtl_beta_sq_filter in "${eqtl_beta_sq_filter_arr[@]}"
do


gwas_sumstat_file=$simulated_gwas_dir${simulation_name_string}"_gt_arch_"${simulated_gt_architecture}"_simualated_gwas_results.txt"
variant_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
variant_stdev_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_genotype_stdev_chrom"
regression_snp_ldscore_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_hm3_ldscores_chrom"
eqtl_cis_snp_h2_summary_file=${estimated_cis_snp_h2_dir}${simulation_name_string}"_"${eqtl_sample_size}"_"${cis_snp_h2_method}"_est_cis_snp_h2_summary.txt"

marginal_eqtl_summary_file=$simulated_gene_expression_dir${simulation_name_string}"_true_eqtls_corr_replicate_eqtl_sumstats_xt_summary.txt"
pred_eqtl_summary_file=$simulated_gene_expression_dir${simulation_name_string}"_true_eqtls_corr_replicate_eqtl_sumstats_xt_summary.txt"

echo $pred_eqtl_summary_file

variant_m_filestem=${simulation_genotype_dir}"variant_reference_genotype_data_ldscores_chrom"
	
	output_stem=${trait_med_h2_inference_dir}${simulation_name_string}"_"${non_med_anno}"_"${simulated_gt_architecture}"_"${inference_gt_architecture}"_"${eqtl_beta_sq_filter}"_"${inference_approach}"_INF_SS_weighted"
	echo $output_stem
	date
	python3 ${calibrated_mesc_code_dir}run_calibrated_mesc.py \
		--gwas-sumstat-file $gwas_sumstat_file \
		--gwas-N ${n_gwas_individuals} \
		--variant-ldscore-filestem ${variant_ldscore_filestem} \
		--variant-M-filestem ${variant_m_filestem} \
		--regression-snp-ldscore-filestem $regression_snp_ldscore_filestem \
		--variant-stdev-filestem ${variant_stdev_filestem} \
		--chromosome-file ${chrom_nums_file} \
		--eqtl-training-data-summary-file ${pred_eqtl_summary_file} \
		--eqtl-validation-data-summary-file ${marginal_eqtl_summary_file} \
		--non-mediated-annotation-version ${non_med_anno} \
		--gene-trait-architecture ${inference_gt_architecture} \
		--inference-approach ${inference_approach} \
		--squared-eqtl-effect-threshold ${eqtl_beta_sq_filter} \
		--mesc-expression-score-dir $mesc_expression_score_dir \
		--training-data-eqtl-ldscores-type ${training_data_gene_ldscores_type} \
		--validation-data-eqtl-ldscores-type ${validation_data_gene_ldscores_type} \
		--eqtl-cis-snp-h2-summary-file ${eqtl_cis_snp_h2_summary_file} \
		--output-stem ${output_stem}
	date
done
done
done
done
fi






























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






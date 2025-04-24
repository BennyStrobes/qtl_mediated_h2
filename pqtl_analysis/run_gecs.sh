#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=60G                         # Memory total in MiB (for all cores)





gecs_expression_score_dir="$1"
gecs_results_dir="$2"
sumstat_dir="$3"
non_redundent_summary_statistics_file="$4"
regresssion_gene_windows="$5"
tglr_expression_score_dir="$6"
tglr_results_dir="$7"

if false; then
source ~/.bash_profile
fi

########################################
# Run GECS for each trait
########################################
gwas_trait_name="PANUKBB_LegFatPercentageRight"

trait_file=${sumstat_dir}${gwas_trait_name}".sumstats"
if false; then
gecs_output_stem=${gecs_results_dir}${gwas_trait_name}"_"${regresssion_gene_windows}"_baselineLD_no_qtl_gecs"
trait_total_h2_est_stem=${tglr_results_dir}${gwas_trait_name}"_susie_inf_pmces_baselineLD_no_qtl_tglr_h2"
python3 run_gecs.py --gwas-data-file $trait_file --gecs-input-data-dir $gecs_expression_score_dir --regress-in-gene-windows --n-genes-dir $tglr_expression_score_dir --total-h2-input-stem ${trait_total_h2_est_stem} --out $gecs_output_stem
 
# GECS with genotype intercept
gecs_output_stem=${gecs_results_dir}${gwas_trait_name}"_"${regresssion_gene_windows}"_genotype_intercept_gecs"
trait_total_h2_est_stem=${tglr_results_dir}${gwas_trait_name}"_susie_inf_pmces_genotype_intercept_tglr_h2"
python3 run_gecs.py --gwas-data-file $trait_file --gecs-input-data-dir $gecs_expression_score_dir --genotype-intercept-only --regress-in-gene-windows --n-genes-dir $tglr_expression_score_dir --total-h2-input-stem ${trait_total_h2_est_stem} --out $gecs_output_stem
fi

if false; then
sed 1d $non_redundent_summary_statistics_file | while read gwas_trait_name h2 h2_se h2_z intercept; do


	trait_file=${sumstat_dir}${gwas_trait_name}".sumstats"

	echo $gwas_trait_name

	# GECS with genotype intercept
	gecs_output_stem=${gecs_results_dir}${gwas_trait_name}"_"${regresssion_gene_windows}"_genotype_intercept_gecs"
	trait_total_h2_est_stem=${tglr_results_dir}${gwas_trait_name}"_susie_inf_pmces_genotype_intercept_tglr_h2"
	python3 run_gecs.py --gwas-data-file $trait_file --gecs-input-data-dir $gecs_expression_score_dir --genotype-intercept-only --regress-in-gene-windows --n-genes-dir $tglr_expression_score_dir --total-h2-input-stem ${trait_total_h2_est_stem} --out $gecs_output_stem

	# GECS with baseline ld
	gecs_output_stem=${gecs_results_dir}${gwas_trait_name}"_"${regresssion_gene_windows}"_baselineLD_no_qtl_gecs"
	trait_total_h2_est_stem=${tglr_results_dir}${gwas_trait_name}"_susie_inf_pmces_baselineLD_no_qtl_tglr_h2"
	python3 run_gecs.py --gwas-data-file $trait_file --gecs-input-data-dir $gecs_expression_score_dir --regress-in-gene-windows --n-genes-dir $tglr_expression_score_dir --total-h2-input-stem ${trait_total_h2_est_stem} --out $gecs_output_stem

done
fi




python3 organize_gecs_results_across_traits.py $tglr_results_dir $gecs_results_dir $non_redundent_summary_statistics_file

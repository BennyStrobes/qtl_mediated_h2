#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-12:10                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)







simulation_number="$1"
chrom_string="$2"
simulation_name_string="$3"
simulation_genotype_dir="$4"
simulated_learned_gene_models_dir="$5"
mesc_code_dir="$6"
plink_executable="$7"
mesc_expression_score_dir="$8"


echo $simulation_genotype_dir

# RUNTIME ESTIMATE (50h?)
# Can 1000 ss run on 10GB




module load python/2.7.12
conda activate mesc

# Need to iterate across:
## 1. Chromosomes (1 and 2)
## 2. eQTL sample sizes (100, 200, 300, and 1000)
## 3. Sample splits (full data, replicate 1, and replicate 2)
## 4. Tissues (0, 1, 2, 3, 4)

################################
#NOTES!!!!!!!
# Consider removing --est-lasso-only flag
# Consider adding the --keep to specify which snps to make expression scores for
#################################




eqtl_sample_size_arr=( "100" "200" "300" "1000")
sample_split_arr=( "full" "replicate1" "replicate2")


for chrom_num in $(seq 1 2); do 
for eqtl_ss in "${eqtl_sample_size_arr[@]}"; do
for sample_split in "${sample_split_arr[@]}"; do
for tissue_num in $(seq 0 4); do


	# Print current iteration
	echo ${chrom_num}"_"${eqtl_ss}"_"${sample_split}"_"${tissue_num}

	# File containing HM3 ids
	hm3_rsid_file="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/genotype_processing/hm3_rsids_chr"${chrom_num}".txt"
	
	# Expression matrix for this iteration
	expr_file=${simulated_learned_gene_models_dir}${simulation_name_string}"_tissue"${tissue_num}"_"${eqtl_ss}"_"${sample_split}"_expression.txt"

	# Referene genotype data (used to make ldscors)
	ref_plink_stem=${simulation_genotype_dir}"simulated_reference_genotype_data_"${chrom_num}

	# Plink data corresponding to expression samples (though this data could contain too many samples)
	full_expr_plink_stem=${simulation_genotype_dir}"simulated_eqtl_"${eqtl_ss}"_data_"${chrom_num}

	# Where to write output to
	output_stem=${mesc_expression_score_dir}${simulation_name_string}"_tissue"${tissue_num}"_"${eqtl_ss}"_"${sample_split}"_"${chrom_num}

	# Where to write temporary output to
	tmp_output_dir=${mesc_expression_score_dir}"mesc_run_"$simulation_number"_"${chrom_num}"_"${eqtl_ss}"_"${sample_split}"_"${tissue_num}
	mkdir ${tmp_output_dir}


	#########################
	# Filter full_expr_plink_stem to only samples we have expression for
	#########################
	# First extract sample names we have expression for
	sample_names_expr_file=$tmp_output_dir"/"$simulation_number"_"${chrom_num}"_"${eqtl_ss}"_"${sample_split}"_"${tissue_num}"_gene_expression_sample_names.txt"
	python3 extract_sample_names_list_from_mesc_gene_expression_file.py $expr_file $sample_names_expr_file

	# Run plink
	expr_plink_stem=$tmp_output_dir"/"$simulation_number"_"${chrom_num}"_"${eqtl_ss}"_"${sample_split}"_"${tissue_num}"_geno_"${chrom_num}
	${plink_executable} --allow-no-sex --bfile ${full_expr_plink_stem} --keep ${sample_names_expr_file} --make-bed --out $expr_plink_stem --silent

	# Error check to make sure sample names match
	python3 error_check_that_sample_names_match.py $sample_names_expr_file $expr_plink_stem".fam"


	#########################
	# Run MESC
	#########################
	python ${mesc_code_dir}run_mesc.py --compute-expscore-indiv --plink-path $plink_executable --expression-matrix $expr_file --exp-bfile $expr_plink_stem --geno-bfile $ref_plink_stem --chr $chrom_num --tmp $tmp_output_dir --keep ${hm3_rsid_file} --out $output_stem


	#########################
	# Delete unnecessary files
	#########################
	rm $sample_names_expr_file
	rm ${expr_plink_stem}".bed"
	rm ${expr_plink_stem}".bim"
	rm ${expr_plink_stem}".fam"
	rm ${expr_plink_stem}".log"


done
done
done
done




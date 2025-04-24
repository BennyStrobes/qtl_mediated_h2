#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4G                         # Memory total in MiB (for all cores)

source ~/.bash_profile


downsampled_run_name="$1"
downsampled_tissue_info_file="$2"
processed_expression_dir="$3"
downsampled_expression_dir="$4"


echo $downsampled_run_name

sed 1d $downsampled_tissue_info_file | while read tissue_name sample_size tissue_sample_names_file; do

	echo $tissue_name
	# Non downsampled (full data) expression file
	full_data_expression_file=${processed_expression_dir}${tissue_name}"_normalized_expression.txt"

	# Downsampled expression (output) file
	downsampled_expression_file=${downsampled_expression_dir}${downsampled_run_name}"_"${tissue_name}"_normalized_expression.txt"
	downsampled_expression_pc_file=${downsampled_expression_dir}${downsampled_run_name}"_"${tissue_name}"_expression_pc.cov"

	# Downsample expression in this tissue
	python3 downsample_gene_expression.py $full_data_expression_file $tissue_sample_names_file $downsampled_expression_file $downsampled_expression_pc_file

done
















if false; then
tissue_name="Adipose_Subcutaneous"
tissue_sample_names_file="/n/scratch/users/b/bes710/qtl_mediated_h2/gtex_downsampling/downsampling_info/downsample_6536_Adipose_Subcutaneous_EUR_sample_names.txt"
fi
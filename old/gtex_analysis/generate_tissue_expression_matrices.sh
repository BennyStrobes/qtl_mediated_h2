#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4G                         # Memory total in MiB (for all cores)

source ~/.bash_profile



tissue_info_file="$1"
gtex_normalized_expression_dir="$2"
gene_annotation_file="$3"
pseudotissue_expression_dir="$4"


# Loop through tissues in tissue info file
sed 1d $tissue_info_file | while read tissue_name sample_size tissue_sample_names_file; do
	echo $tissue_name
	python3 generate_tissue_expression_matrix.py $tissue_name $gtex_normalized_expression_dir $gene_annotation_file $tissue_sample_names_file $pseudotissue_expression_dir
done




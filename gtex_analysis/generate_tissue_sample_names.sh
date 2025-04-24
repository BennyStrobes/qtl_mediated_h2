#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=2G                         # Memory total in MiB (for all cores)

gtex_v8_normalized_expression_matrices_dir="${1}"
gtex_v8_european_list_file="${2}"
tissue_sample_names_dir="${3}"
tissue_info_file="${4}"


python3 generate_tissue_sample_names.py $gtex_v8_normalized_expression_matrices_dir $gtex_v8_european_list_file $tissue_sample_names_dir $tissue_info_file
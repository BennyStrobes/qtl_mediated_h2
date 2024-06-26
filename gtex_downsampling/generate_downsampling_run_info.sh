#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=2G                         # Memory total in MiB (for all cores)

source ~/.bash_profile




tissue_info_file="${1}"
downsampling_info_dir="${2}"
downsample_run_summary_file="${3}"



python3 generate_downsampling_run_info.py $tissue_info_file $downsampling_info_dir $downsample_run_summary_file
#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)







downsample_run_summary_file="$1"
mesc_results_dir="$2"



variant_anno="genotype_intercept"
python3 meta_analyze_mesc_results_across_downsampling_runs.py $downsample_run_summary_file $mesc_results_dir $variant_anno



variant_anno="baselineLD_no_qtl"
python3 meta_analyze_mesc_results_across_downsampling_runs.py $downsample_run_summary_file $mesc_results_dir $variant_anno
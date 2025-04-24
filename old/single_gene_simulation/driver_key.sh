#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-22:40                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=12GB                         # Memory total in MiB (for all cores)





genotype_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simple_simulation_per_gene_calibration/processed_genotype/"


simulated_results_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/single_gene_simulation/simulated_results/"



python3 run_single_gene_simulation.py $genotype_dir $simulated_results_dir
#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:30                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=15GB                         # Memory total in MiB (for all cores)




raw_ld_file="$1"
mediated_h2_results="$2"
n_causal_eqtl_snps="$3"

echo $n_causal_eqtl_snps

python3 correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps_extra_nm_var_in_gene2.py $raw_ld_file $mediated_h2_results $n_causal_eqtl_snps

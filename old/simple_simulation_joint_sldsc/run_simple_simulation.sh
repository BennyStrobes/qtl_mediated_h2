#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-13:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



# COmmand line args
gwas_ss="$1"
n_snps="$2"
nm_h2="$3"
n_genes="$4"
med_h2="$5"
snps_per_gene="$6"
ge_h2="$7"
eqtl_ss="$8"
sim_iter="$9"
output_dir="${10}"

python3 run_simple_simulation.py $gwas_ss $n_snps $nm_h2 $n_genes $med_h2 $snps_per_gene $ge_h2 $eqtl_ss $sim_iter $output_dir
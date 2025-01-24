#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-1:35                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



sim_iter="$1"
n_sims="$2"
gwas_ss="$3"
n_snps="$4"
n_genes="$5"
snps_per_gene="$6"
med_h2="$7"
nm_h2="$8"
eqtl_ss="$9"
eqtl_arch="${10}"
mean_cis_h2="${11}"
frac_causal_genes="${12}"
output_stem="${13}"


module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate

echo $output_stem


python3 run_proportional_med_h2_no_ld_simulation.py $sim_iter $n_sims $gwas_ss $n_snps $n_genes $snps_per_gene $med_h2 $nm_h2 $eqtl_ss $eqtl_arch $mean_cis_h2 $frac_causal_genes $output_stem

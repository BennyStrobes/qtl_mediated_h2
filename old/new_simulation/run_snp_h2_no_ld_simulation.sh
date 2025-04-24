#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5GB                         # Memory total in MiB (for all cores)




n_snps="$1"
sim_h2="$2"
sample_size="$3"
trait_architecture="$4"
output_stem="$5"




module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate


echo ${trait_architecture}
echo ${sample_size}
echo ${n_snps}
echo ${sim_h2}
echo ${output_stem}


python3 run_snp_h2_no_ld_simulation.py $n_snps $sim_h2 $sample_size $trait_architecture $output_stem
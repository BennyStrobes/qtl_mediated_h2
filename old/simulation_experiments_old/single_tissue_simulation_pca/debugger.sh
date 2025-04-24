#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:45                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=3GB     


module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate


python3 debugger.py
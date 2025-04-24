#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-25:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=12G                         # Memory total in MiB (for all cores)




Tissue_name="$1"
sample_names_file="$2"
genotype_stem="$3"
expr_file="$4"
cov_file="$5"
gwas_genotype_stem="$6"
plink_executable="$7"
output_root="$8"


# Load in correct packages
module load python/3.7.4
module load R/4.0.1


tissue_output_dir=${output_root}"GTEx_v8_genotype_EUR_"${Tissue_name}"/"
mkdir $tissue_output_dir


# Loop through chromosomes
for CHR in {1..22}
do

	echo "CHROMOSOME "${CHR}

	# Filter GTEx data to samples from this tissue
    tmp_plink_stem=${tissue_output_dir}"GTEx_v8_genotype_EUR_tmp_"${Tissue_name}"."${CHR}  # Need to delete later
    plink2 \
        --bfile ${genotype_stem}${CHR} \
        --keep-allele-order \
        --chr $CHR \
        --keep ${sample_names_file} \
        --threads 1 \
        --make-bed \
        --out ${tmp_plink_stem}
    

    # Estimate causal eQTL effects
	eqtl_output_root=${tissue_output_dir}"GTEx_v8_eqtl_"${Tissue_name}"_"${CHR}
	python3 estimate_causal_eqtl_effects_and_eqtl_ss_for_single_tissue_chromosome_pairing.py $tmp_plink_stem $expr_file $cov_file $gwas_genotype_stem $CHR $eqtl_output_root


    # Delete unnessary, temporary plink file
    rm ${tmp_plink_stem}".bim"
    rm ${tmp_plink_stem}".fam"
    rm ${tmp_plink_stem}".bed"
    rm ${tmp_plink_stem}".log"
	
done
#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4G                         # Memory total in MiB (for all cores)

source ~/.bash_profile



gtex_v8_european_list_file="$1"
gtex_v8_genotype_dir="$2"
processed_genotype_dir="$3"
gwas_genotype_stem="$4"


# Loop through chromosomes
for CHR in {1..22}
do
    echo $CHR

	# Filter to only EUR individuals
    plink2 \
        --bfile ${gtex_v8_genotype_dir}"GTEx_v8_genotype."${CHR} \
        --keep-allele-order \
        --chr $CHR \
        --keep ${gtex_v8_european_list_file} \
        --threads 1 \
        --make-bed \
        --out ${processed_genotype_dir}"GTEx_v8_genotype_EUR_tmp."${CHR}

    # Impose MAF filter of .01
    plink2 \
        --bfile ${processed_genotype_dir}"GTEx_v8_genotype_EUR_tmp."${CHR} \
        --keep-allele-order \
        --chr $CHR \
        --maf 0.01 \
        --max-maf 0.99 \
        --threads 1 \
        --make-bed \
        --out ${processed_genotype_dir}"GTEx_v8_genotype_EUR_tmp2."${CHR}

    # Extract list of snp ids from GTEx data that have a corresponding rsid in 1KG data
    snp_id_file=${processed_genotype_dir}"GTEx_v8_genotype_EUR_tmp_rsid_"${CHR}".txt"
    python3 extract_list_of_gtex_snp_ids_in_gwas_rsid_bim.py ${processed_genotype_dir}"GTEx_v8_genotype_EUR_tmp2."${CHR}".bim" $gwas_genotype_stem${CHR}".bim" $snp_id_file

    # Filter to those snps
    plink2 \
        --bfile ${processed_genotype_dir}"GTEx_v8_genotype_EUR_tmp2."${CHR} \
        --keep-allele-order \
        --chr $CHR \
        --extract ${snp_id_file} \
        --threads 1 \
        --make-bed \
        --out ${processed_genotype_dir}"GTEx_v8_genotype_EUR."${CHR}

    # Hack to convert bim file from GTEx snpid to rsid
    cp ${processed_genotype_dir}"GTEx_v8_genotype_EUR."${CHR}".bim" ${processed_genotype_dir}"GTEx_v8_genotype_EUR_tmper."${CHR}".bim"
    python3 convert_gtex_snp_snp_id_bim_file_to_rsid_bim_file.py ${processed_genotype_dir}"GTEx_v8_genotype_EUR_tmper."${CHR}".bim" $gwas_genotype_stem${CHR}".bim" ${processed_genotype_dir}"GTEx_v8_genotype_EUR."${CHR}".bim"



    # Calculate allele frequencies
    plink2 \
        --bfile ${processed_genotype_dir}"GTEx_v8_genotype_EUR."${CHR} \
        --freq \
        --threads 1 \
        --out ${processed_genotype_dir}"GTEx_v8_genotype_EUR."${CHR}    

    # Remove unnessary file
    rm ${processed_genotype_dir}"GTEx_v8_genotype_EUR_tmp"*


done







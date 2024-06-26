#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-13:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)

source ~/.bash_profile




Tissue_name="$1"
sample_names_file="$2"
genotype_stem="$3"
expr_file="$4"
cov_file="$5"
gwas_genotype_stem="$6"
mesc_code_dir="$7"
plink_executable="$8"
mesc_run_name="${9}"
expression_score_dir="${10}"




module load python/2.7.12
conda activate mesc

echo $mesc_run_name
echo $Tissue_name
date


###############################
# Create MESC expression scores
################################
tmp_output_dir=${expression_score_dir}"GTEx_v8_genotype_EUR_tmp_"${mesc_run_name}"_"${Tissue_name}
mkdir ${tmp_output_dir}

# Loop through chromosomes
for CHR in {1..22}
do


    # Filter GTEx data to samples from this tissue
    tmp_plink_stem=${expression_score_dir}"GTEx_v8_genotype_EUR_tmp_"${mesc_run_name}"_"${Tissue_name}"."${CHR}  # Need to delete later
    plink2 \
        --bfile ${genotype_stem}${CHR} \
        --keep-allele-order \
        --chr $CHR \
        --keep ${sample_names_file} \
        --threads 1 \
        --make-bed \
        --out ${tmp_plink_stem}

    # Quick error checking
    python3 error_check_to_make_sure_expression_and_genotype_sample_names_match.py $expr_file $tmp_plink_stem".fam"


    # Run MESC analysis
	python ${mesc_code_dir}run_mesc.py --compute-expscore-indiv --plink-path $plink_executable --expression-matrix $expr_file --exp-bfile $tmp_plink_stem --geno-bfile $gwas_genotype_stem${CHR} --chr $CHR --covariates $cov_file --tmp $tmp_output_dir --out $expression_score_dir"overall_"${mesc_run_name}"_"${Tissue_name}"_"${CHR} --est-lasso-only
    
    # Delete unnessary, temporary plink file
    rm ${tmp_plink_stem}".bim"
    rm ${tmp_plink_stem}".fam"
    rm ${tmp_plink_stem}".bed"
    rm ${tmp_plink_stem}".log"
done


# Delete temporary output dir
rm -rf ${tmp_output_dir}

date


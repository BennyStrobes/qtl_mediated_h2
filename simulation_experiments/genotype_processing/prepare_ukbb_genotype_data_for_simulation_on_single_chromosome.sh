#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=15GB                         # Memory total in MiB (for all cores)

# First three parts ran at 200GB
# Needs to be run for longer (was just running small section)


ukbb_genotype_dir="$1"
processed_genotype_data_root_dir="$2"
chrom_num="$3"
n_gwas_individuals="$4"
ldsc_baseline_hg19_annotation_dir="$5"
kg_genotype_dir="$6"
hm3_snp_list_dir="$7"
quasi_independent_dir="$8"
gencode_gene_annotation_file="${9}"
calibrated_mesc_code_dir="${10}"
source ~/.bash_profile


###############################
# Extract list of hm3 variants and variants in ldsc baseline analysis
###############################
hm3_rs_id_file=${processed_genotype_data_root_dir}"hm3_rsids_chr"${chrom_num}".txt"
if false; then
python3 extract_list_of_hm3_rs_ids.py $hm3_snp_list_dir $chrom_num $kg_genotype_dir $hm3_rs_id_file
fi

baselineld_rs_id_file=${processed_genotype_data_root_dir}"baselineld_rsids_chr"${chrom_num}".txt"
if false; then
python3 extract_list_of_ldsc_annotation_rs_ids.py $ldsc_baseline_hg19_annotation_dir $chrom_num $kg_genotype_dir $baselineld_rs_id_file
fi

###############################
# Make genotype subdirectory for this gwas sample size
###############################
processed_genotype_data_dir=${processed_genotype_data_root_dir}"gwas_sample_size_"${n_gwas_individuals}"/"
mkdir $processed_genotype_data_dir
echo $processed_genotype_data_dir



###############################
# Filter UKBB genotype data to only include those variants in ldsc baseline analysis
###############################
if false; then
plink2 \
    --bgen /n/groups/price/UKBiobank/download_500K/ukb_imp_chr"${chrom_num}"_v3.bgen ref-unknown\
    --sample /n/groups/price/UKBiobank/download_500K/ukb14048_imp_chr1_v3_s487395.sample\
    --keep /n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/unrelated_337K.txt\
    --extract /n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/snp_info/snp_list_chr${chrom_num}.MAF_001_INFO_06.txt\
    --rm-dup force-first\
    --maj-ref\
    --geno 0.1\
    --maf 0.001\
    --hwe 1e-50\
    --make-pgen \
    --threads 1\
    --out ${processed_genotype_data_dir}"ukb_imp_chr"${chrom_num}"_tmper"


plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr"${chrom_num}"_tmper" --hwe .01 --maf .05 --extract ${baselineld_rs_id_file} --keep-allele-order --threads 1 --make-pgen --out ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num}
fi

###############################
# extract lists of Individuals for each data set
###############################
gwas_individual_file=${processed_genotype_data_dir}"gwas_individuals.txt"
eqtl_individual_stem=${processed_genotype_data_dir}"eqtl_individuals_"
ref_genotype_individual_file=${processed_genotype_data_dir}"ref_genotype_individuals.txt"
if false; then
python3 extract_lists_of_simulated_individuals.py ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} $n_gwas_individuals $gwas_individual_file $eqtl_individual_stem $ref_genotype_individual_file


###############################
# Filter UKBB genotype data to only include individuals in simulated gwas data set 
###############################
plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${gwas_individual_file} --make-pgen --keep-allele-order --threads 1 --out ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}


###############################
# Filter UKBB genotype data to only include individuals in eqtl gwas data set 
###############################
eqtl_sample_size="100"
plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-pgen --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="200"
plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-pgen --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="300"
plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-pgen --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="500"
plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-pgen --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="1000"
plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-pgen --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="10000"
plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-pgen --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

# Reference data set
plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${ref_genotype_individual_file} --make-pgen --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_reference_genotype_data_"${chrom_num}


###############################
# Convert to bgen
###############################
# GWAS 
plink2 --pfile ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num} --keep-allele-order --export bgen-1.2 --out ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}

# eQTL 
eqtl_sample_size="100"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --export bgen-1.2 --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}
eqtl_sample_size="200"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --export bgen-1.2 --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}
eqtl_sample_size="300"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --export bgen-1.2 --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}
eqtl_sample_size="500"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --export bgen-1.2 --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}
eqtl_sample_size="1000"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --export bgen-1.2 --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}
eqtl_sample_size="10000"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --export bgen-1.2 --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

# Reference genotype
plink2 --pfile ${processed_genotype_data_dir}"simulated_reference_genotype_data_"${chrom_num} --keep-allele-order --export bgen-1.2 --out ${processed_genotype_data_dir}"simulated_reference_genotype_data_"${chrom_num}
fi


###############################
# Convert to plink1 (needed for mesc)
###############################
# GWAS 
plink2 --pfile ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num} --keep-allele-order -make-bed --out ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}

# eQTL 
eqtl_sample_size="100"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --make-bed --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}
eqtl_sample_size="200"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --make-bed --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}
eqtl_sample_size="300"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --make-bed --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}
eqtl_sample_size="500"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --make-bed --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}
eqtl_sample_size="1000"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --make-bed --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}
eqtl_sample_size="10000"
plink2 --pfile ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num} --keep-allele-order --make-bed --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

# Reference genotype
plink2 --pfile ${processed_genotype_data_dir}"simulated_reference_genotype_data_"${chrom_num} --keep-allele-order -make-bed --out ${processed_genotype_data_dir}"simulated_reference_genotype_data_"${chrom_num}





if false; then
module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate
fi
genotype_version="reference_genotype_data"

# Extract variant level LD scores
variant_ld_score_file=${processed_genotype_data_dir}"variant_"${genotype_version}"_ldscores_chrom"${chrom_num}".txt"
variant_M_file=${processed_genotype_data_dir}"variant_"${genotype_version}"_ldscores_chrom"${chrom_num}"_M.txt"
if false; then
python3 ${calibrated_mesc_code_dir}extract_variant_ldscores.py --chrom $chrom_num --bgen-file ${processed_genotype_data_dir}"simulated_"${genotype_version}"_"${chrom_num}".bgen" --hm3-rsid-file $hm3_rs_id_file --sldsc-annotation-file ${ldsc_baseline_hg19_annotation_dir}"baselineLD."${chrom_num}".annot.gz" --variant-ld-score-file $variant_ld_score_file --variant-M-file $variant_M_file
fi

# Extract hm3 variant level LD scores
variant_ld_score_file=${processed_genotype_data_dir}"variant_"${genotype_version}"_hm3_ldscores_chrom"${chrom_num}".txt"
if false; then
python3 ${calibrated_mesc_code_dir}extract_hm3_variant_ldscores_for_weighting.py --chrom $chrom_num --bgen-file ${processed_genotype_data_dir}"simulated_"${genotype_version}"_"${chrom_num}".bgen" --cm-position-file ${ldsc_baseline_hg19_annotation_dir}"baselineLD."${chrom_num}".annot.gz" --hm3-rsid-file $hm3_rs_id_file --variant-ld-score-file $variant_ld_score_file
fi


genotype_version="reference_genotype_data"
# Extract variant level LD scores
variant_sdev_file=${processed_genotype_data_dir}"variant_"${genotype_version}"_genotype_stdev_chrom"${chrom_num}".txt"
if false; then
python3 ${calibrated_mesc_code_dir}extract_genotype_standard_deviations.py --chrom $chrom_num --bgen-file ${processed_genotype_data_dir}"simulated_"${genotype_version}"_"${chrom_num}".bgen" --output-file $variant_sdev_file
fi


# Extract gene list
# Genes are defined by actual tss
# Limit to protein coding genes
cis_window="500000"
simulated_gene_position_file=${processed_genotype_data_dir}"gene_positions_chr"${chrom_num}".txt"
if false; then
python3 prepare_simulated_gene_position_list.py $chrom_num $gencode_gene_annotation_file $simulated_gene_position_file $cis_window
fi






















###################
# OLD
###################


# Construct gene level LD matrices
if false; then
gene_snp_representation="bins_10"
gene_ld_output_root=${processed_genotype_data_dir}"gene_level_ld_chr"${chrom_num}"_"${gene_snp_representation}
python3 construct_gene_level_ld_matrices.py $chrom_num $simulated_gene_position_file ${processed_genotype_data_dir}"simulated_"${genotype_version}"_"${chrom_num} $variant_ld_score_file $kg_genotype_dir $gene_snp_representation $gene_ld_output_root
fi

# Construct gene level LD matrices
gene_snp_representation="pca"
gene_ld_output_root=${processed_genotype_data_dir}"gene_level_ld_chr"${chrom_num}"_"${gene_snp_representation}
if false; then
python3 construct_gene_level_ld_matrices.py $chrom_num $simulated_gene_position_file ${processed_genotype_data_dir}"simulated_"${genotype_version}"_"${chrom_num} $variant_ld_score_file $kg_genotype_dir $gene_snp_representation $gene_ld_output_root
fi


if false; then
gene_snp_representation="bins_20"
gene_ld_output_root=${processed_genotype_data_dir}"gene_level_ld_chr"${chrom_num}"_"${gene_snp_representation}
python3 construct_gene_level_ld_matrices.py $chrom_num $simulated_gene_position_file ${processed_genotype_data_dir}"simulated_"${genotype_version}"_"${chrom_num} $variant_ld_score_file $kg_genotype_dir $gene_snp_representation $gene_ld_output_root
fi



# Construct gene level LD matrices
if false; then
gene_snp_representation="pca2"
gene_ld_output_root=${processed_genotype_data_dir}"gene_level_ld_chr"${chrom_num}"_"${gene_snp_representation}
python3 construct_gene_level_ld_matrices.py $chrom_num $simulated_gene_position_file ${processed_genotype_data_dir}"simulated_"${genotype_version}"_"${chrom_num} $variant_ld_score_file $kg_genotype_dir $gene_snp_representation $gene_ld_output_root
fi














####################
# OLDer
####################



if false; then
#########################
# Get LD matrices
##########################
module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate

python3 get_quasi_independent_ld_matrices.py ${processed_genotype_data_dir} $quasi_independent_dir
fi




if false; then
echo "START"

python3 get_ld_matrices.py ${processed_genotype_data_dir}


echo "MID"


#########################
# Get variant LD scores
##########################
python3 get_variant_ld_scores.py ${processed_genotype_data_dir}

echo "DONE"
fi








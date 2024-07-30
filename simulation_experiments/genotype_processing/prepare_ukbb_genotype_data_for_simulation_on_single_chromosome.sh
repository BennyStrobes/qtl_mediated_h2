#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-25:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=200GB                         # Memory total in MiB (for all cores)

# First three parts ran at 200GB



ukbb_genotype_dir="$1"
processed_genotype_data_root_dir="$2"
chrom_num="$3"
n_gwas_individuals="$4"
ldsc_baseline_hg19_annotation_dir="$5"
kg_genotype_dir="$6"


source ~/.bash_profile

###############################
# Extract list of variants in ldsc baseline analysis
###############################
ldsc_annotation_rs_id_file=${processed_genotype_data_root_dir}"ldsc_annotation_rsids_chr"${chrom_num}".txt"
python3 extract_list_of_ldsc_annotation_rs_ids.py $ldsc_baseline_hg19_annotation_dir $chrom_num $kg_genotype_dir $ldsc_annotation_rs_id_file


###############################
# Make genotype subdirectory for this gwas sample size
###############################
processed_genotype_data_dir=${processed_genotype_data_root_dir}"gwas_sample_size_"${n_gwas_individuals}"/"
mkdir $processed_genotype_data_dir
echo $processed_genotype_data_dir

###############################
# Filter UKBB genotype data to only include those variants in ldsc baseline analysis
###############################
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

plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr"${chrom_num}"_tmper" --hwe .01 --extract ${ldsc_annotation_rs_id_file} --maf .05 --make-bed --keep-allele-order --threads 1 --out ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num}


###############################
# extract lists of Individuals for each data set
###############################
gwas_individual_file=${processed_genotype_data_dir}"gwas_individuals.txt"
eqtl_individual_stem=${processed_genotype_data_dir}"eqtl_individuals_"
ref_genotype_individual_file=${processed_genotype_data_dir}"ref_genotype_individuals.txt"
python3 extract_lists_of_simulated_individuals.py ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} $n_gwas_individuals $gwas_individual_file $eqtl_individual_stem $ref_genotype_individual_file


###############################
# Filter UKBB genotype data to only include individuals in simulated gwas data set 
###############################
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${gwas_individual_file} --make-bed --keep-allele-order --threads 1 --out ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}


###############################
# Filter UKBB genotype data to only include individuals in eqtl gwas data set 
###############################
eqtl_sample_size="100"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="200"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="300"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="500"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="1000"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}


###############################
# Filter UKBB genotype data to only include individuals in reference genotype data
###############################
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${ref_genotype_individual_file} --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_reference_genotype_data_"${chrom_num}


#########################
# Filter 1KG variants to only those in our analysis
##########################
plink2 --bfile ${kg_genotype_dir}"1000G.EUR.QC."${chrom_num} --extract ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}".bim" --threads 1 --make-bed --keep-allele-order --out ${processed_genotype_data_dir}"100G.EUR.QC.filtered."${chrom_num}


#########################
# Get variant LD scores
##########################
source ~/.bash_profile


python3 get_variant_ld_scores.py ${processed_genotype_data_dir}

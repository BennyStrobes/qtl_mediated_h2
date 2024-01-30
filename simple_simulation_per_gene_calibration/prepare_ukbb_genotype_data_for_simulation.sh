#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-13:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=200GB                         # Memory total in MiB (for all cores)

# First three parts ran at 160GB



ukbb_genotype_dir="$1"
processed_genotype_data_dir="$2"

echo $processed_genotype_data_dir

chrom_num="21"

source ~/.bash_profile


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
    --maf 0.05\
    --hwe 1e-50\
    --make-pgen \
    --threads 1\
    --out ${processed_genotype_data_dir}"ukb_imp_chr"${chrom_num}"_tmper"

plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr"${chrom_num}"_tmper" --hwe .01 --maf .05 --make-bed --keep-allele-order --threads 1 --out ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num}




###############################
# extract lists of Individuals for each data set
###############################
n_gwas_individuals="100000"
gwas_individual_file=${processed_genotype_data_dir}"gwas_individuals.txt"
eqtl_individual_stem=${processed_genotype_data_dir}"eqtl_individuals_"
python3 extract_lists_of_simulated_individuals.py ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} $n_gwas_individuals $gwas_individual_file $eqtl_individual_stem


###############################
# Filter UKBB genotype data to only include individuals in simulated gwas data set 
###############################
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${gwas_individual_file} --make-bed --keep-allele-order --threads 1 --out ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}


###############################
# Filter UKBB genotype data to only include individuals in eqtl gwas data set 
###############################
eqtl_sample_size="100"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}


eqtl_sample_size="300"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="500"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="1000"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="5000"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}



###############################
# Extract genotype data for a bunch of 4MB windows!
###############################
python3 extract_genotype_data_for_a_bunch_of_multi_mb_windows.py $chrom_num ${processed_genotype_data_dir}


python3 split_gwas_genotype_genotype_data_into_chunks.py ${processed_genotype_data_dir}










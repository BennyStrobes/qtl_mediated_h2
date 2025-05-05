#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-40:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=26GB                         # Memory total in MiB (for all cores)




chrom_num="${1}"
processed_genotype_data_root_dir="${2}"
ldsc_baseline_hg19_annotation_dir="${3}"
gene_snp_representation="${4}"
joint_ldsc_code_dir="${5}"
n_gwas_individuals="${6}"


echo ${chrom_num}"_"${gene_snp_representation}

module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate


processed_genotype_data_dir=${processed_genotype_data_root_dir}"gwas_sample_size_"${n_gwas_individuals}"/"


genotype_version="reference_genotype_data"
simulated_gene_position_file=${processed_genotype_data_dir}"gene_positions_chr"${chrom_num}".txt"
variant_ld_score_file=${processed_genotype_data_dir}"variant_"${genotype_version}"_ldscores_chrom"${chrom_num}".txt"
bgen_file=${processed_genotype_data_dir}"simulated_"${genotype_version}"_"${chrom_num}".bgen"
sldsc_anno_file=${ldsc_baseline_hg19_annotation_dir}"baselineLD."${chrom_num}".annot.gz"

gene_ld_output_root=${processed_genotype_data_dir}"gene_level_ld_chr"${chrom_num}"_"${gene_snp_representation}
python3 ${joint_ldsc_code_dir}construct_gene_level_ld_matrices.py --chrom ${chrom_num} --gene-position-file ${simulated_gene_position_file} --bgen-file ${bgen_file} --variant-ld-score-file ${variant_ld_score_file} --gene-snp-representation ${gene_snp_representation} --sldsc-annotation-file ${sldsc_anno_file} --gene-ld-output-root ${gene_ld_output_root}
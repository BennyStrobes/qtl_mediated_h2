#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-15:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)




chrom_num="$1"
causal_eqtl_summary_file="$2"
gwas_genotype_stem="$3"
ldsc_baseline_ld_annotation_stem="$4"
gene_annotation_file="$5"
tglr_expression_score_dir="$6"

date
version="tissue_stratification"
python3 estimate_gene_ld_scores_in_single_chromosome.py $chrom_num ${causal_eqtl_summary_file} ${gwas_genotype_stem} ${ldsc_baseline_ld_annotation_stem} ${gene_annotation_file} ${tglr_expression_score_dir} ${version}

date
version="tissue_aggregation"
python3 estimate_gene_ld_scores_in_single_chromosome.py $chrom_num ${causal_eqtl_summary_file} ${gwas_genotype_stem} ${ldsc_baseline_ld_annotation_stem} ${gene_annotation_file} ${tglr_expression_score_dir} ${version}

date

version="distance_stratification"
python3 estimate_gene_ld_scores_in_single_chromosome.py $chrom_num ${causal_eqtl_summary_file} ${gwas_genotype_stem} ${ldsc_baseline_ld_annotation_stem} ${gene_annotation_file} ${tglr_expression_score_dir} ${version}

date

version="qtl_rank_stratification"
python3 estimate_gene_ld_scores_in_single_chromosome.py $chrom_num ${causal_eqtl_summary_file} ${gwas_genotype_stem} ${ldsc_baseline_ld_annotation_stem} ${gene_annotation_file} ${tglr_expression_score_dir} ${version}
date
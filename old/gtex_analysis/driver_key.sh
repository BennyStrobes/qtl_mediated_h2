

##################
# Input data
##################

# GTEx v8 normalizd expression matrices directory
gtex_v8_normalized_expression_matrices_dir="/n/groups/price/GTEX/GTEX_V8/normalized_expression/"

# GTEx v8 genotype directory
gtex_v8_genotype_dir="/n/groups/price/GTEX/GTEX_V8/genotype/"

# File containing list of European ancestry GTEx individuals
gtex_v8_european_list_file="/n/groups/price/huwenbo/DATA/GTEx_v8/GTEx_v8_phenotypes/EUR_subject.txt"

# Code directory containing MESC
# ** mesc_code_dir="/n/groups/price/ben/tools/mesc/"
mesc_code_dir="/n/scratch/users/b/bes710/tmp_storage/mesc/mesc/"

# Code directory containing LDSC
ldsc_code_dir="/n/groups/price/ben/tools/ldsc/"
#ldsc_code_dir="/n/scratch/users/b/bes710/tmp_storage/ldsc_new_code/ldsc/"

# List of non-redundent summary statistics
non_redundent_summary_statistics_file="/n/groups/price/ldsc/sumstats_formatted_2024/non_redundent_traits_EUR_2024.txt"
# **non_redundent_summary_statistics_file="/n/scratch/users/b/bes710/tmp_storage/ldsc/sumstats_formatted_2024/non_redundent_traits_EUR_2024.txt"

# Directory containing summary statistics
sumstat_dir="/n/groups/price/ldsc/sumstats_formatted_2024/sumstats/"
# **sumstat_dir="/n/scratch/users/b/bes710/tmp_storage/ldsc/sumstats_formatted_2024/sumstats/"

# Directory containing plink files for 100G
kg_plink_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"
#  **kg_plink_dir="/n/scratch/users/b/bes710/tmp_storage/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"

# Directory containing baselineLD annotation
baselineLD_anno_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"
# baselineLD_anno_dir="/n/scratch/users/b/bes710/tmp_storage/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"

# Directory containing ldsc weights
ldsc_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/weights/"
#ldsc_weights_dir="/n/scratch/users/b/bes710/tmp_storage/ldsc/reference_files/1000G_EUR_Phase3_hg38/weights/"

# File containing hapmap3 rsid
hapmap3_rsid_file="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/w_hm3.noMHC.snplist"
#hapmap3_rsid_file="/n/scratch/users/b/bes710/tmp_storage/ldsc/reference_files/1000G_EUR_Phase3_hg38/w_hm3.noMHC.snplist"

# GTEx gencode gene annotation file
# Downloaded from https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf on Jan 19 2022
gene_annotation_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/input_data/gencode.v26.GRCh38.genes.gtf"
#gene_annotation_file="/n/scratch/users/b/bes710/tmp_storage/ge_misc/gencode.v26.GRCh38.genes.gtf"

# Plink executable
# https://www.cog-genomics.org/plink/ (Linux 64-bit)
plink_executable="/n/scratch/users/b/bes710/tmp_storage/plink/plink"

# Code directory for tglr
tglr_code_dir="/n/groups/price/ben/tools/tglr/tglr/"


##################
# Output data
##################
# Temporary output root
tmp_output_root="/n/scratch/users/b/bes710/qtl_mediated_h2/gtex_analysis/"
perm_output_root="/n/groups/price/ben/qtl_mediated_h2/gtex_analysis/"

# Directory containing tissue sample names
tissue_sample_names_dir=${tmp_output_root}"tissue_sample_names/"

# Directory containing processed gene expression
processed_expression_dir=${tmp_output_root}"processed_expression_dir/"

# Directory containing processed genotype data
processed_genotype_dir=${tmp_output_root}"processed_genotype/"

# Directory containing GTEx estimated causal eQTL effects
gtex_causal_eqtl_effects_dir=${perm_output_root}"gtex_causal_eqtl_effects/"

# Directory containing TGLR variant annotations
tglr_variant_annotations_dir=${tmp_output_root}"tglr_variant_annotations/"

# TGLR expression score directory
tglr_expression_score_dir=${tmp_output_root}"tglr_expression_scores/"

# TGLR results directory
tglr_results_dir=${tmp_output_root}"tglr_results/"

# Visualize tglr results
viz_results_dir=${perm_output_root}"visualize_tglr_results/"


##################
# Run analysis
##################


###################
# create ordered list of samples for each of the tissues
tissue_info_file=$tissue_sample_names_dir"tissue_info.txt"
if false; then
sh generate_tissue_sample_names.sh $gtex_v8_normalized_expression_matrices_dir $gtex_v8_european_list_file $tissue_sample_names_dir $tissue_info_file
fi

###################
# Generate expression matrices for each of tissues
if false; then
sh generate_tissue_expression_matrices.sh $tissue_info_file $gtex_v8_normalized_expression_matrices_dir $gene_annotation_file $processed_expression_dir
fi

###################
# Filter GTEx genotype data to only EUR donors
# ALso convert snp ids to rsids (as we need to filter on rsids in mesc)
gwas_genotype_stem=${kg_plink_dir}"1000G.EUR.hg38."
if false; then
sh filter_gtex_genotype_to_EUR_donors.sh $gtex_v8_european_list_file $gtex_v8_genotype_dir $processed_genotype_dir $gwas_genotype_stem
fi



###################
# Prepare variant annotations for TGLR
if false; then
sbatch preprocess_variant_annotation_for_tglr.sh $ldsc_code_dir $hapmap3_rsid_file $baselineLD_anno_dir $kg_plink_dir $tglr_variant_annotations_dir
fi



################## 
# Estimate causal eQTL effects
genotype_stem=$processed_genotype_dir"GTEx_v8_genotype_EUR."
gwas_genotype_stem=${kg_plink_dir}"1000G.EUR.hg38."

tissue_info_file3="/n/groups/price/ben/code_temp3/new_missing.txt"
if false; then
sed 1d $tissue_info_file3 | while read Tissue_name; do
	echo $Tissue_name
	sample_names_file="/n/scratch/users/b/bes710/qtl_mediated_h2/gtex_analysis/tissue_sample_names/"${Tissue_name}"_EUR_sample_names.txt"

	expr_file=${processed_expression_dir}${Tissue_name}"_normalized_expression.txt"
	cov_file=${processed_expression_dir}${Tissue_name}"_expression_pc.cov"
	sbatch estimate_causal_eqtl_effects_in_single_tissue.sh ${Tissue_name} ${sample_names_file} ${genotype_stem} ${expr_file} ${cov_file} ${gwas_genotype_stem} ${plink_executable} ${gtex_causal_eqtl_effects_dir}
done
fi

# Organize causal eqtl effects across tissues
tissue_info_file=$tissue_sample_names_dir"tissue_info.txt"
if false; then
python3 organize_causal_eqtl_effect_summaries_across_tissues.py ${tissue_info_file} ${gtex_causal_eqtl_effects_dir}
fi





################## 
# Estimate gene LD scores
ldsc_baseline_ld_annotation_stem=${tglr_variant_annotations_dir}"baselineLD_no_qtl."

causal_eqtl_summary_file=${gtex_causal_eqtl_effects_dir}"GTEx_v8_genotype_EUR_cross_tissue_gene_model_summary.txt"
if false; then
for chrom_num in $(seq 1 22); do 
	echo $chrom_num
	sbatch estimate_gene_ld_scores_in_single_chromosome.sh $chrom_num ${causal_eqtl_summary_file} ${gwas_genotype_stem} ${ldsc_baseline_ld_annotation_stem} ${gene_annotation_file} ${tglr_expression_score_dir}
done
fi




################## 
# Run TGLR
ldsc_baseline_ld_annotation_stem=${tglr_variant_annotations_dir}"baselineLD_no_qtl."
ldsc_genotype_intercept_annotation_stem=${tglr_variant_annotations_dir}"genotype_intercept."
sldsc_weights_stem=${tglr_variant_annotations_dir}"regression_weights."
if false; then
sqtl_version="tissue_aggregation"
sh run_tglr.sh ${qtl_version} ${ldsc_genotype_intercept_annotation_stem} ${ldsc_baseline_ld_annotation_stem} ${tglr_expression_score_dir} ${tglr_results_dir} $sumstat_dir ${non_redundent_summary_statistics_file} ${sldsc_weights_stem} ${tglr_code_dir}
fi

if false; then
qtl_version="tissue_stratification"
sh run_tglr.sh ${qtl_version} ${ldsc_genotype_intercept_annotation_stem} ${ldsc_baseline_ld_annotation_stem} ${tglr_expression_score_dir} ${tglr_results_dir} $sumstat_dir ${non_redundent_summary_statistics_file} ${sldsc_weights_stem} ${tglr_code_dir}

qtl_version="distance_stratification"
sh run_tglr.sh ${qtl_version} ${ldsc_genotype_intercept_annotation_stem} ${ldsc_baseline_ld_annotation_stem} ${tglr_expression_score_dir} ${tglr_results_dir} $sumstat_dir ${non_redundent_summary_statistics_file} ${sldsc_weights_stem} ${tglr_code_dir}

qtl_version="qtl_rank_stratification"
sh run_tglr.sh ${qtl_version} ${ldsc_genotype_intercept_annotation_stem} ${ldsc_baseline_ld_annotation_stem} ${tglr_expression_score_dir} ${tglr_results_dir} $sumstat_dir ${non_redundent_summary_statistics_file} ${sldsc_weights_stem} ${tglr_code_dir}
fi


Rscript visualize_tglr_results.R ${tglr_results_dir} $viz_results_dir















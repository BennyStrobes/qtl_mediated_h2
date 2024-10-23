
##################
# Input data
##################

# Code directory containing LDSC
ldsc_code_dir="/n/groups/price/ldsc/ldsc/"

# List of non-redundent summary statistics
non_redundent_summary_statistics_file="/n/groups/price/ldsc/sumstats_formatted_2024/non_redundent_traits_EUR_2024.txt"

# Directory containing summary statistics
sumstat_dir="/n/groups/price/ldsc/sumstats_formatted_2024/sumstats/"

# Directory containing plink files for 100G
kg_plink_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"

# Directory containing baselineLD annotation
baselineLD_anno_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/"

# Directory containing ldsc weights
ldsc_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/"

# File containing hapmap3 rsid
hapmap3_rsid_file="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/1000G_hm3_noMHC.snplist"

# GTEx gencode gene annotation file
# Downloaded from https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf on Jan 19 2022
gene_annotation_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/input_data/gencode.v26.GRCh38.genes.gtf"

# Sent from Ali 
# Generated using SbayesRC
trans_protein_models_dir="/n/groups/price/UKBiobank/proteomics_betas/beta/"

# Code directory for tglr
tglr_code_dir="/n/groups/price/ben/tools/tglr/tglr/"

# Directory containing models of protein dat
protein_models_dir="/n/groups/price/kangcheng/projects/ukb-ppp-data/protein-model/susie_inf_hm3/"

# Directory containing file for each gene of cis-snps relative to each gene
cis_snp_dir="/n/groups/price/ali/proteomic_proj/data/gene_annot/cis_snps_hapmap3/"



##################
# Output data
##################
# Temporary output root
tmp_output_root="/n/scratch/users/b/bes710/qtl_mediated_h2/protein_med_h2/"

# Directory containing TGLR variant annotations
tglr_variant_annotations_dir=${tmp_output_root}"tglr_variant_annotations/"

# TGLR expression score directory
tglr_expression_score_dir=${tmp_output_root}"tglr_expression_scores/"

# TGLR Trans expression score directory
tglr_trans_expression_score_dir=${tmp_output_root}"tglr_trans_expression_scores/"

# TGLR results directory
tglr_results_dir=${tmp_output_root}"tglr_results/"

# Visualize TGLR results directory
visualize_tglr_results_dir=${tmp_output_root}"visualize_tglr_results/"





###################
# Prepare variant annotations for TGLR
if false; then
sbatch preprocess_variant_annotation_for_tglr.sh $ldsc_code_dir $hapmap3_rsid_file $baselineLD_anno_dir $kg_plink_dir $tglr_variant_annotations_dir
fi




###################
# Estimate TGLR expression scores in each tissue, independently 
# Use MESC lasso weights as input
ldsc_baseline_ld_annotation_stem=${tglr_variant_annotations_dir}"baselineLD_no_qtl."
ldsc_genotype_intercept_annotation_stem=${tglr_variant_annotations_dir}"genotype_intercept."
sldsc_weights_stem=${tglr_variant_annotations_dir}"regression_weights."

gwas_genotype_stem=${kg_plink_dir}"1000G.EUR.QC."
if false; then
sh estimate_tglr_expression_score_across_tissues.sh $protein_models_dir ${gwas_genotype_stem} ${ldsc_baseline_ld_annotation_stem} ${tglr_expression_score_dir}
fi



###################
# Run TGLR
if false; then
sbatch run_tglr.sh ${ldsc_genotype_intercept_annotation_stem} ${ldsc_baseline_ld_annotation_stem} ${tglr_expression_score_dir} ${tglr_results_dir} $tglr_code_dir $sumstat_dir ${non_redundent_summary_statistics_file} ${sldsc_weights_stem}
fi



###################
# Estimate TGLR expression scores in each tissue, independently 
# Use cis+trans pqtl gene models from ali
ldsc_baseline_ld_annotation_stem=${tglr_variant_annotations_dir}"baselineLD_no_qtl."
ldsc_genotype_intercept_annotation_stem=${tglr_variant_annotations_dir}"genotype_intercept."
sldsc_weights_stem=${tglr_variant_annotations_dir}"regression_weights."
gwas_genotype_stem=${kg_plink_dir}"1000G.EUR.QC."

if false; then
eqtl_class="cis_trans_pqtl"
sbatch estimate_tglr_trans_expression_scores.sh $trans_protein_models_dir ${gwas_genotype_stem} ${ldsc_baseline_ld_annotation_stem} ${tglr_trans_expression_score_dir} $eqtl_class $cis_snp_dir
fi
if false; then
eqtl_class="cis_pqtl"
sbatch estimate_tglr_trans_expression_scores.sh $trans_protein_models_dir ${gwas_genotype_stem} ${ldsc_baseline_ld_annotation_stem} ${tglr_trans_expression_score_dir} $eqtl_class $cis_snp_dir
eqtl_class="trans_pqtl"
sbatch estimate_tglr_trans_expression_scores.sh $trans_protein_models_dir ${gwas_genotype_stem} ${ldsc_baseline_ld_annotation_stem} ${tglr_trans_expression_score_dir} $eqtl_class $cis_snp_dir
fi

###################
# Run TGLR
if false; then
sh run_tglr_w_trans_qtls.sh ${ldsc_genotype_intercept_annotation_stem} ${ldsc_baseline_ld_annotation_stem} ${tglr_trans_expression_score_dir} ${tglr_results_dir} $tglr_code_dir $sumstat_dir ${non_redundent_summary_statistics_file} ${sldsc_weights_stem}
fi



###################
# Visualize TGLR results
####################
if false; then
module load R/3.5.1
Rscript visualize_tglr_results.R ${tglr_results_dir} $visualize_tglr_results_dir
fi




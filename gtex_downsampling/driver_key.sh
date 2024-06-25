

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
# ** ldsc_code_dir="/n/groups/price/ldsc/ldsc/"
ldsc_code_dir="/n/scratch/users/b/bes710/tmp_storage/ldsc/ldsc/"

# List of non-redundent summary statistics
# ** non_redundent_summary_statistics_file="/n/groups/price/ldsc/sumstats_formatted_2024/non_redundent_traits_EUR_2024.txt"
non_redundent_summary_statistics_file="/n/scratch/users/b/bes710/tmp_storage/ldsc/sumstats_formatted_2024/non_redundent_traits_EUR_2024.txt"

# Directory containing summary statistics
# ** sumstat_dir="/n/groups/price/ldsc/sumstats_formatted_2024/sumstats/"
sumstat_dir="/n/scratch/users/b/bes710/tmp_storage/ldsc/sumstats_formatted_2024/sumstats/"

# Directory containing plink files for 100G
# ** kg_plink_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"
kg_plink_dir="/n/scratch/users/b/bes710/tmp_storage/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"

# Directory containing baselineLD annotation
# ** baselineLD_anno_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"
baselineLD_anno_dir="/n/scratch/users/b/bes710/tmp_storage/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"

# Directory containing ldsc weights
# ** ldsc_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/weights/"
ldsc_weights_dir="/n/scratch/users/b/bes710/tmp_storage/ldsc/reference_files/1000G_EUR_Phase3_hg38/weights/"

# GTEx gencode gene annotation file
# Downloaded from https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf on Jan 19 2022
# ** gene_annotation_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/input_data/gencode.v26.GRCh38.genes.gtf"
gene_annotation_file="/n/scratch/users/b/bes710/tmp_storage/ge_misc/gencode.v26.GRCh38.genes.gtf"



##################
# Output data
##################
# Temporary output root
tmp_output_root="/n/scratch/users/b/bes710/qtl_mediated_h2/gtex_downsampling/"

# Directory containing tissue sample names
tissue_sample_names_dir=${tmp_output_root}"tissue_sample_names/"

# Directory containing processed gene expression
processed_expression_dir=${tmp_output_root}"processed_expression_dir/"

# Directory containing Down-sampled gene expression
downsampled_expression_dir=${tmp_output_root}"downsampled_expression_dir/"

# Directory containing processed genotype data
processed_genotype_dir=${tmp_output_root}"processed_genotype/"

# Directory containing downsampled genotype data
downsampled_genotype_dir=${tmp_output_root}"downsampled_genotype/"

# Directory containing expression scores
expression_score_dir=${tmp_output_root}"expression_scores/"


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
# Estimate overall expression scores in each tissue, independently 
# See https://github.com/douglasyao/mesc/wiki/Estimating-overall-expression-scores
mesc_run_name="full_data"
genotype_stem=$processed_genotype_dir"GTEx_v8_genotype_EUR."
gwas_genotype_stem=${kg_plink_dir}"1000G.EUR.hg38."

# Loop
if false; then
sed 1d $tissue_info_file | while read Tissue_name sample_size sample_names_file; do
	expr_file=${processed_expression_dir}${Tissue_name}"_normalized_expression.txt"
	cov_file=${processed_expression_dir}${Tissue_name}"_expression_pc.cov"
	sbatch estimate_mesc_expression_score_in_single_tissue.sh ${Tissue_name} ${sample_names_file} ${genotype_stem} ${expr_file} ${cov_file} ${gwas_genotype_stem} ${mesc_code_dir} ${mesc_run_name} ${expression_score_dir}
done
fi









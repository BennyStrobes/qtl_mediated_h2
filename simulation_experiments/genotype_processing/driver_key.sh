############################
# Input data
############################

# Directory created by Martin containing UKBB genotype for 334K unrelated European individuals
ukbb_genotype_dir="/n/groups/price/UKBiobank/bgen_MAF001_500K_v3/"

# LDSC baseline annotations (hg19)
ldsc_baseline_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/"

# LDSC 1KG genotype files (hg19)
kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"

# Directory containing HM3 snp lists
hm3_snp_list_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/"

# Directory containing quasi indpendent ld blocks
quasi_independent_dir="/n/groups/price/ben/quasi_independent_ld_blocks/"

# Gencode hg19 gene annotation file
gencode_gene_annotation_file="/n/groups/price/ben/gene_annotation_files/gencode.v19.annotation.gtf.gz"

# Joint-LDSC code dir
joint_ldsc_code_dir="/n/groups/price/ben/joint_ldsc/"

# Calibrated mesc code dir
calibrated_mesc_code_dir="/n/groups/price/ben/calibrated_mesc_v3/"


############################
# Output data
############################
processed_genotype_data_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/genotype_processing/"



############################
# Simulation parameters
############################
# Chromosome to simulate on 



############################
# Prepare genotype data for analysis:
## 1. Filter number of individuals in original data
## 2. Filter sites to be those in LDSC annotation file
## 3. Convert to plink bed files
############################
# NOTE: THERE IS CURRENTLY A HACK IN HERE TO REMOVE 3 variants (out of 500000) on chrom 1 that have no variance across the 100-sample eqtl data set.
############################
# Needs to be 200GB and 25 h
if false; then
n_gwas_individuals="100000"
chrom_num="1"
sh prepare_ukbb_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir $hm3_snp_list_dir $quasi_independent_dir $gencode_gene_annotation_file $calibrated_mesc_code_dir

chrom_num="2"
sh prepare_ukbb_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir $hm3_snp_list_dir $quasi_independent_dir $gencode_gene_annotation_file $calibrated_mesc_code_dir
fi


# TO do
# 1. let above finish running
# 2. compare to original. make sure same number of rows and base ldscores are the same
# 3. Switch filenames so new is main and tmp is deleted
# 4. Run below on it.



gene_snp_representation_arr=( "bins_10" "bins_20" "pca_90" "pca_95" "pca_99")
if false; then
for gene_snp_representation in "${gene_snp_representation_arr[@]}"
do
	chrom_num="1"
	sbatch prepare_gene_ld_files.sh $chrom_num $processed_genotype_data_dir $ldsc_baseline_hg19_annotation_dir $gene_snp_representation $joint_ldsc_code_dir $n_gwas_individuals
	chrom_num="2"
	sbatch prepare_gene_ld_files.sh $chrom_num $processed_genotype_data_dir $ldsc_baseline_hg19_annotation_dir $gene_snp_representation $joint_ldsc_code_dir $n_gwas_individuals
done
fi


############################
# Input data
############################

# Directory created by Martin containing UKBB genotype for 334K unrelated European individuals
ukbb_genotype_dir="/n/groups/price/UKBiobank/bgen_MAF001_500K_v3/"

# LDSC baseline annotations (hg19)
ldsc_baseline_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline_v1.2/"

# LDSC 1KG genotype files (hg19)
kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"

# Directory containing HM3 snp lists
hm3_snp_list_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/"

# Directory containing quasi indpendent ld blocks
quasi_independent_dir="/n/groups/price/ben/quasi_independent_ld_blocks/"


############################
# Output data
############################
processed_genotype_data_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/genotype_processing/"



############################
# Simulation parameters
############################
# Chromosome to simulate on 
chrom_num="1"



############################
# Prepare genotype data for analysis:
## 1. Filter number of individuals in original data
## 2. Filter sites to be those in LDSC annotation file
## 3. Convert to plink bed files
############################
# NOTE: THERE IS CURRENTLY A HACK IN HERE TO REMOVE 3 variants (out of 500000) on chrom 1 that have no variance across the 100-sample eqtl data set.
############################
# Needs to be 200GB and 25 h
n_gwas_individuals="100000"
if false; then
sh prepare_ukbb_hm3_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir $hm3_snp_list_dir $quasi_independent_dir
fi







# No longer completely up to date
if false; then
sh prepare_ukbb_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir 
fi




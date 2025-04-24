#################################
# Simple (ie. small) simulation
#################################






############################
# Input data
############################

# Directory created by Martin containing UKBB genotype for 334K unrelated European individuals
ukbb_genotype_dir="/n/groups/price/UKBiobank/bgen_MAF001_500K_v3/"

# LDSC baseline annotations (hg19)
ldsc_baseline_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline_v1.2/"

# LDSC 1KG genotype files (hg19)
kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"


############################
# Output data
############################
# Root directory
output_root="/n/groups/price/ben/qtl_mediated_h2/simple_simulation_per_gene_calibration/"
temp_output_root="/n/scratch/users/b/bes710/qtl_mediated_h2/simple_simulation_per_gene_calibration/"

# Output directory for processed genotype
processed_genotype_data_dir=${temp_output_root}"processed_genotype/"




############################
# Run analysis
############################
# Prepare genotype data
if false; then
sbatch prepare_ukbb_genotype_data_for_simulation.sh $ukbb_genotype_dir $processed_genotype_data_dir
fi

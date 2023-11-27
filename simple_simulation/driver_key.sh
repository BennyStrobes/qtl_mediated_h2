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
output_root="/n/groups/price/ben/qtl_mediated_h2/simple_simulation/"
temp_output_root="/n/scratch3/users/b/bes710/qtl_mediated_h2/simple_simulation/"

# Output directory for processed genotype
processed_genotype_data_dir=${temp_output_root}"processed_genotype/"

# Output directory for simulated eqtl data
simulated_eqtl_data_dir=${temp_output_root}"simulated_eqtl_data/"

# Output directory for simulated expression data
simulated_expression_data_dir=${temp_output_root}"simulated_expression_data/"

# Output directory containing simulated learned gene models
simulated_gene_models_dir=${temp_output_root}"simulated_learned_gene_models/"

# Output directory for simulated gwas data
simulated_gwas_data_dir=${temp_output_root}"simulated_gwas_data/"

# Output directory for simulated gwas data
mediated_h2_results_dir=${temp_output_root}"mediated_h2_results/"





############################
# Run analysis
############################
# Prepare genotype data
if false; then
sh prepare_ukbb_genotype_data_for_simulation.sh $ukbb_genotype_dir $processed_genotype_data_dir
fi


# RUN SIMULATIONS
simulation_number="1"
sh run_single_simulation_shell.sh $simulation_number $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_expression_data_dir $simulated_gene_models_dir $simulated_gwas_data_dir $mediated_h2_results_dir



if false; then
for simulation_number in $(seq 2 20); do 
	sh run_single_simulation_shell.sh $simulation_number $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_expression_data_dir $simulated_gene_models_dir $simulated_gwas_data_dir $mediated_h2_results_dir
done
fi





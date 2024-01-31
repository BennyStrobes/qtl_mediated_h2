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
temp_output_root_old="/n/scratch/users/b/bes710/qtl_mediated_h2/simple_simulation_per_gene_calibration/"
temp_output_root="/n/scratch/users/b/bes710/qtl_mediated_h2/simple_simulation_multiple_eqtl_datasets/"


# Output directory for processed genotype
processed_genotype_data_dir=${temp_output_root_old}"processed_genotype/"

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

# Output directory for visualizing results
visualize_med_h2_results_dir=${temp_output_root}"visualize_mediated_h2_results/"





############################
# Run analysis
############################
# RUN SIMULATIONS
simulation_number="1"
if false; then
sbatch run_single_simulation_shell.sh $simulation_number $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_gwas_data_dir $mediated_h2_results_dir $simulated_expression_data_dir $simulated_gene_models_dir
fi

if false; then
for simulation_number in $(seq 2 100); do 
	sbatch run_single_simulation_shell.sh $simulation_number $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_gwas_data_dir $mediated_h2_results_dir $simulated_expression_data_dir $simulated_gene_models_dir
done
fi

if false; then
for simulation_number in $(seq 231 400); do 
	sbatch run_single_simulation_shell.sh $simulation_number $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_gwas_data_dir $mediated_h2_results_dir $simulated_expression_data_dir $simulated_gene_models_dir
done
fi





if false; then
sh organize_and_visualize_simulation_results.sh $mediated_h2_results_dir $visualize_med_h2_results_dir
fi






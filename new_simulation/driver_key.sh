#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=12GB                         # Memory total in MiB (for all cores)

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
output_root="/n/groups/price/ben/qtl_mediated_h2/new_simulation/"
temp_output_root="/n/scratch/users/b/bes710/qtl_mediated_h2/new_simulation/"

# Directory containing results of NO-LD simulation
no_ld_simulation_res_dir=${temp_output_root}"no_ld_simulation_res_dir/"

# Heritability analysis results
snp_h2_simulation_dir=${temp_output_root}"snp_h2_simulation/"

# mediated heritability analysis results
med_h2_simulation_dir=${temp_output_root}"med_h2_simulation/"




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

# Output directory for visualizing results
visualize_med_h2_results_dir=${temp_output_root}"visualize_mediated_h2_results/"





############################
# Run analysis
############################
if false; then
for simulation_number in $(seq 1 100); do 
	sbatch run_no_ld_simulations.sh $simulation_number $no_ld_simulation_res_dir
done
fi


if false; then
for simulation_number in $(seq 1 100); do 
	sbatch run_no_ld_simulations.sh $simulation_number $no_ld_simulation_res_dir
done
fi

if false; then
for simulation_number in $(seq 101 200); do 
	sbatch run_no_ld_simulations.sh $simulation_number $no_ld_simulation_res_dir
done
fi

if false; then
for simulation_number in $(seq 301 400); do 
	sbatch run_no_ld_simulations.sh $simulation_number $no_ld_simulation_res_dir
done
fi

if false; then
simulation_number="302"
sh run_no_ld_simulations.sh $simulation_number $no_ld_simulation_res_dir
fi



#################################
# Run snp h2 (no ld) simulation
#################################
n_snps="500"
sim_h2="0.1"
sample_size="100"
trait_architecture="sparse"
if false; then
sample_size_arr=("1000") 
trait_architecture_arr=("polygenic" "sparse")
sim_h2_arr=("0.005" "0.01" "0.05" "0.1" "0.2")
for sample_size in "${sample_size_arr[@]}"; do
for sim_h2 in "${sim_h2_arr[@]}"; do
for trait_architecture in "${trait_architecture_arr[@]}"; do
	output_stem=${snp_h2_simulation_dir}"snp_h2_simulation_n_snps_"${n_snps}"_sample_size_"${sample_size}"_h2_"${sim_h2}"_arch_"$trait_architecture
	sbatch run_snp_h2_no_ld_simulation.sh $n_snps $sim_h2 $sample_size $trait_architecture $output_stem
done
done
done
fi




#################################
# Run med h2 (no ld) simulation
# 1 causal tissue
#################################
n_sims="10"
gwas_ss="10000"
n_snps="40000"
n_genes="100"
snps_per_gene="200"
med_h2="0.1"
nm_h2="0.3"
eqtl_ss="500"
mean_cis_h2="0.05"
eqtl_arch="sparse"
frac_causal_genes="0.5"


frac_causal_genes_arr=("0.25" "0.5" "0.75" "1.0")
mean_cis_h2_arr=("0.01" "0.05" "0.1")


sim_iter="1"
if false; then
for sim_iter in $(seq 1 10); do 
for frac_causal_genes in "${frac_causal_genes_arr[@]}"; do
for mean_cis_h2 in "${mean_cis_h2_arr[@]}"; do
output_stem=${med_h2_simulation_dir}"proportional_med_h2_simulation_"${sim_iter}"_gwas_ss_"${gwas_ss}"_eqtl_ss_"${eqtl_ss}"_med_h2_"${med_h2}"_nm_h2_"${nm_h2}"_mean_cis_h2_"${mean_cis_h2}"_gene_frac_"${frac_causal_genes}"_arch_"$eqtl_arch
sbatch run_proportional_med_h2_no_ld_simulation.sh $sim_iter $n_sims $gwas_ss $n_snps $n_genes $snps_per_gene $med_h2 $nm_h2 $eqtl_ss $eqtl_arch $mean_cis_h2 $frac_causal_genes $output_stem
done
done
done
fi

sim_iter="1"
frac_causal_genes="0.5"
mean_cis_h2="0.05"
n_sims="1"

frac_causal_genes_arr=("0.5")
mean_cis_h2_arr=("0.01" "0.05" "0.1")

if false; then
for sim_iter in $(seq 11 100); do 
for frac_causal_genes in "${frac_causal_genes_arr[@]}"; do
for mean_cis_h2 in "${mean_cis_h2_arr[@]}"; do
output_stem=${med_h2_simulation_dir}"proportional_med_h2_simulation_"${sim_iter}"_gwas_ss_"${gwas_ss}"_eqtl_ss_"${eqtl_ss}"_med_h2_"${med_h2}"_nm_h2_"${nm_h2}"_mean_cis_h2_"${mean_cis_h2}"_gene_frac_"${frac_causal_genes}"_arch_"$eqtl_arch
sbatch run_proportional_med_h2_no_ld_simulation.sh $sim_iter $n_sims $gwas_ss $n_snps $n_genes $snps_per_gene $med_h2 $nm_h2 $eqtl_ss $eqtl_arch $mean_cis_h2 $frac_causal_genes $output_stem
done
done
done
fi


if false; then
python3 organize_proportional_med_h2_results.py $med_h2_simulation_dir



Rscript visualize_proportional_med_h2_results.R $med_h2_simulation_dir
fi


if false; then
Rscript visualize_no_ld_simulation_results.R $no_ld_simulation_res_dir
fi




if false; then
python3 temp_builder.py
fi

if false; then
module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/numpyro/bin/activate
fi

if false; then
eqtl_arch="sparse"
calibration_output_file=$no_ld_simulation_res_dir"calibration_analysis_"${eqtl_arch}".txt"
python3 calibration_analysis.py $eqtl_arch $calibration_output_file
fi
if false; then
eqtl_arch="infinitesimal"
calibration_output_file=$no_ld_simulation_res_dir"calibration_analysis_"${eqtl_arch}".txt"
python3 calibration_analysis.py $eqtl_arch $calibration_output_file
fi

if false; then
python3 organize_no_ld_simulation_results.py $no_ld_simulation_res_dir
fi

if false; then
Rscript visualize_no_ld_simulation_results.R $no_ld_simulation_res_dir
fi


# Prepare genotype data
if false; then
sbatch prepare_ukbb_genotype_data_for_simulation.sh $ukbb_genotype_dir $processed_genotype_data_dir
fi



simulation_number="1"
if false; then
sh run_single_simulation_shell.sh $simulation_number $processed_genotype_data_dir $simulated_eqtl_data_dir $simulated_expression_data_dir $simulated_gene_models_dir $simulated_gwas_data_dir $mediated_h2_results_dir
fi

if false; then
sh organize_and_visualize_simulation_results.sh $mediated_h2_results_dir $visualize_med_h2_results_dir
fi






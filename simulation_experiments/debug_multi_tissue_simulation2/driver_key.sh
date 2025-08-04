#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-30:45                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5GB   

############################
# Input data
############################
# Real genotype data directories
simulation_genotype_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/genotype_processing/gwas_sample_size_100000/"

# Joint-LDSC code dir
joint_ldsc_code_dir="/n/groups/price/ben/joint_ldsc/"
joint_ldsc_code_dir="/n/groups/price/ben/joint_ldsc_tmp/"
joint_ldsc_code_dir="/n/groups/price/ben/joint_ldsc_dbg/"

calibrated_mesc_code_dir="/n/groups/price/ben/calibrated_mesc/"
calibrated_mesc_code_dir="/n/groups/price/ben/calibrated_mesc_v2/"
calibrated_mesc_code_dir="/n/groups/price/ben/calibrated_mesc_v5/"


mesc_code_dir="/n/groups/price/ben/tools/mesc_2_chromosomes/"

plink_executable="/n/groups/price/ben/tools/plink_linux_x86_64_20210606/plink"


############################
# Output directories
############################
output_root_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/debug_multi_tissue_simulation2/"
perm_output_root_dir="/n/groups/price/ben/qtl_mediated_h2/simulation_experiments/debug_multi_tissue_simulation2/"


# Directory containing simulated gene expression (ie causal eQTL effect sizes in each tissue)
simulated_gene_expression_dir=$perm_output_root_dir"simulated_gene_expression/"

# Directory containing simulated gene expression (ie causal eQTL effect sizes in each tissue)
simulated_learned_gene_models_dir=$output_root_dir"simulated_learned_gene_models/"

alt_simulated_learned_gene_models_dir=${perm_output_root_dir}"simulated_learned_gene_models_alt/"

# Directory containing simulated trait values
simulated_trait_dir=$output_root_dir"simulated_trait/"

# Directory contaiing simulated gwas results
simulated_gwas_dir=$output_root_dir"simulated_gwas/"

# Directory containing results of trait h2 inference
trait_h2_inference_dir=$output_root_dir"trait_h2_inference/"

# Directory containing mesc expression scores
mesc_expression_score_dir=$output_root_dir"mesc_expression_scores/"

# Directory containing default mesc results
mesc_processed_input_dir=$output_root_dir"mesc_processed_input/"

# Directory containing default mesc results
mesc_results_dir=$output_root_dir"mesc_results/"

# Directory containing results of trait med h2 inference
trait_med_h2_inference_dir=$output_root_dir"trait_mediated_h2_inference/"

# Directory containing estimated cis-snp h2s
estimated_cis_snp_h2_dir=$output_root_dir"estimated_cis_snp_h2/"

# Directory containing lasso estimated gene models
lasso_gene_models_dir=$output_root_dir"lasso_gene_models_dir/"

# Directory containing organized results
organized_trait_med_h2_results_dir=${perm_output_root_dir}"organized_trait_med_h2_results/"

# trait h2 visualization dir
visualize_trait_med_h2_dir=$perm_output_root_dir"visualize_trait_med_h2/"




#####################
# Not sure these directories areused any more
expr_trait_h2_inference_dir=$output_root_dir"expression_trait_h2_inference/"
visualize_expression_trait_h2=$output_root_dir"visualize_expression_trait_h2/"




############################
# Main simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"

# Per genetic-element heritabilities
per_element_heritability="0.0005"

# Total heritability
total_heritability="0.3"

# Fraction of heritability mediated by gene expression
fraction_expression_mediated_heritability="0.1"

# Window size
cis_window="500000"

# Gene expression heritability
ge_h2="05"

# Number of tissues
n_tissues="5"

# eQTL architecture
eqtl_architecture="default"



#chrom_string
chrom_string="1_2"


############################
# Run main single simulation of simulating the trait data
############################
if false; then
for simulation_number in $(seq 1 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch run_single_trait_simulation.sh $simulation_number $chrom_string $cis_window $n_gwas_individuals $simulation_name_string $simulation_genotype_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $eqtl_architecture $n_tissues $alt_simulated_learned_gene_models_dir
done
fi

############################
# Convert true causal eqtl effects to gene ld scores
############################
if false; then
for simulation_number in $(seq 1 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch extract_true_gene_ld_scores.sh $simulation_number $chrom_string $simulation_name_string $simulation_genotype_dir $simulated_gene_expression_dir $calibrated_mesc_code_dir $estimated_cis_snp_h2_dir
done
fi



###############################
# Calibrated mesc approach
###############################
if false; then
for simulation_number in $(seq 1 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch calibrated_mesc_trait_mediated_h2_inference_shell.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $simulated_gene_expression_dir $chrom_string $calibrated_mesc_code_dir $mesc_expression_score_dir $estimated_cis_snp_h2_dir $lasso_gene_models_dir
done
fi


####################
# Organize results
####################
tmp_simulation_name_string="_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
if false; then
python3 organize_trait_med_h2_results.py $trait_med_h2_inference_dir $organized_trait_med_h2_results_dir $simulated_trait_dir $tmp_simulation_name_string
fi














############################
# Create lasso gene models with fixed alpha
############################
alpha_0="0.1"
if false; then
for simulation_number in $(seq 1 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch create_lasso_gene_models.sh $simulation_number $chrom_string $cis_window $simulation_name_string $simulation_genotype_dir $simulated_learned_gene_models_dir $lasso_gene_models_dir $calibrated_mesc_code_dir $cis_window $alpha_0
done
fi


############################
# Create lasso gene models with CV alpha
############################
if false; then
for simulation_number in $(seq 1 200); do 
for tissue_num in $(seq 0 4); do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch create_lasso_gene_models_cv_alpha.sh $simulation_number $chrom_string $cis_window $simulation_name_string $simulation_genotype_dir $simulated_learned_gene_models_dir $lasso_gene_models_dir $calibrated_mesc_code_dir $cis_window $tissue_num
done
done
fi

if false; then
filer="/n/groups/price/ben/code_temp7/1000_rep2_start.txt"
sed '1d' $filer \
| while IFS=$'\t' read -r simulation_number tissue_num; do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch create_lasso_gene_models_cv_alpha.sh $simulation_number $chrom_string $cis_window $simulation_name_string $simulation_genotype_dir $simulated_learned_gene_models_dir $lasso_gene_models_dir $calibrated_mesc_code_dir $cis_window $tissue_num
done
fi





############################
# Convert lasso runs to per gene effect estimates
############################
if false; then
alpha_0="CV"
for simulation_number in $(seq 1 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sh extract_lasso_gene_ld_scores.sh $simulation_number $chrom_string $simulation_name_string $simulation_genotype_dir $mesc_expression_score_dir $lasso_gene_models_dir $calibrated_mesc_code_dir $alpha_0
done
fi


############################
# Convert true causal eqtl effects to gene ld scores
############################
if false; then
for simulation_number in $(seq 1 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sh extract_true_gene_ld_scores.sh $simulation_number $chrom_string $simulation_name_string $simulation_genotype_dir $simulated_gene_expression_dir $calibrated_mesc_code_dir
done
fi




############################
# Use MESC to generate expression scores
############################
if false; then
for simulation_number in $(seq 1 200); do 
for tissue_num in $(seq 0 4); do
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch use_mesc_to_generate_expression_scores.sh $simulation_number $chrom_string $simulation_name_string $simulation_genotype_dir $simulated_learned_gene_models_dir $mesc_code_dir $plink_executable $mesc_expression_score_dir $tissue_num
done
done
fi


############################
# Run default mesc
############################
if false; then
for simulation_number in $(seq 1 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch run_default_mesc_shell.sh $simulation_number $chrom_string $simulation_name_string $simulation_genotype_dir $mesc_code_dir $mesc_expression_score_dir $mesc_processed_input_dir $mesc_results_dir $simulated_gwas_dir
done
fi

############################
# Convert mesc lasso runs to per gene effect estimates
############################
if false; then
for simulation_number in $(seq 1 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sh extract_mesc_lasso_per_gene_scores.sh $simulation_number $chrom_string $simulation_name_string $simulation_genotype_dir $mesc_expression_score_dir
done
fi



###############################
# Calibrated mesc approach
###############################
if false; then
for simulation_number in $(seq 2 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch calibrated_mesc_trait_mediated_h2_inference_shell.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $simulated_gene_expression_dir $chrom_string $calibrated_mesc_code_dir $mesc_expression_score_dir $estimated_cis_snp_h2_dir $lasso_gene_models_dir
done
fi

####################
# Organize results
####################
tmp_simulation_name_string="_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
if false; then
python3 organize_trait_med_h2_results.py $trait_med_h2_inference_dir $organized_trait_med_h2_results_dir $simulated_trait_dir $tmp_simulation_name_string
fi


tmp_simulation_name_string="_chrom"${chrom_string}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
if false; then
python3 organize_mesc_results.py $mesc_results_dir $organized_trait_med_h2_results_dir $simulated_trait_dir $tmp_simulation_name_string
fi


if false; then
module load R/3.5.1
Rscript visualize_trait_med_h2_results.R $trait_med_h2_inference_dir $organized_trait_med_h2_results_dir $visualize_trait_med_h2_dir
fi










































###############################
# Joint LDSC approach
###############################
if false; then
for simulation_number in $(seq 1 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch trait_mediated_h2_inference_no_pca_joint_ldsc.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $simulated_gene_expression_dir
done
fi


###############################
# Joint LDSC approach (relying on pca)
###############################
if false; then
for simulation_number in $(seq 4 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch trait_mediated_h2_inference_joint_ldsc.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $simulated_gene_expression_dir
done
fi




####################
# Organize results
####################
if false; then
python3 organize_trait_med_h2_results.py $trait_med_h2_inference_dir $organized_trait_med_h2_results_dir
fi


if false; then
module load R/3.5.1
fi
if false; then
Rscript visualize_trait_med_h2_results.R $trait_med_h2_inference_dir $organized_trait_med_h2_results_dir $visualize_trait_med_h2_dir
fi


















############
# OLDER
############










########################
# GECS trait mediated inference
###############################
if false; then
for simulation_number in $(seq 1 400); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch gecs_trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $alt_simulated_learned_gene_models_dir
done
fi





###############################
# LD4M trait mediated inference
###############################
simulation_number="1"

if false; then
for simulation_number in $(seq 2 200); do 
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
sbatch mixed_ld4m_trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $alt_simulated_learned_gene_models_dir
done
fi





########################
# Mediated trait h2 inference with PCA-marginal likelihood
###############################
cc_hyper_param_arr=( "0" "1e-6")
cc_hyper_param_arr=( "0")
if false; then
for cc_hyper_param in "${cc_hyper_param_arr[@]}"
do
for simulation_number in $(seq 1 50); do 
	eqtl_sample_size="100"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $cc_hyper_param

	eqtl_sample_size="300"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $cc_hyper_param

	eqtl_sample_size="1000"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $cc_hyper_param

	eqtl_sample_size="10000"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $cc_hyper_param
done
done
fi



########################
# Mediated trait h2 inference with PCA-marginal likelihood
###############################
if false; then
for simulation_number in $(seq 3 100); do 
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
sbatch trait_mediated_h2_inference_no_pca_joint_reml.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir
done
fi


########################
# TGLR trait mediated inference
###############################
if false; then
for simulation_number in $(seq 1 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch tglr_trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $alt_simulated_learned_gene_models_dir
done
fi




###############################
# Joint LDSC approach
###############################
# Using LDSC
if false; then
for simulation_number in $(seq 2 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch trait_mediated_h2_inference_no_pca_joint_ldsc.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $simulated_gene_expression_dir
done
fi


if false; then
for simulation_number in $(seq 201 500); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch trait_mediated_h2_inference_no_pca_joint_ldsc.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $simulated_gene_expression_dir
done
fi

if false; then
for simulation_number in $(seq 1 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch trait_mediated_h2_inference_no_pca_joint_ldsc.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $simulated_gene_expression_dir
done
fi

if false; then
for simulation_number in $(seq 201 500); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch trait_mediated_h2_inference_no_pca_joint_ldsc.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $simulated_gene_expression_dir
done
fi

###############################
# Dissattenuated LDSC approach
###############################
if false; then
for simulation_number in $(seq 4 200); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch trait_mediated_h2_inference_no_pca_dissaten_ldsc.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir $simulated_gene_expression_dir
done
fi


####################
# Organize results
####################
if false; then
python3 organize_trait_med_h2_results.py $trait_med_h2_inference_dir $visualize_trait_med_h2_dir
fi

if false; then
module load R/3.5.1
Rscript visualize_trait_med_h2_results.R $trait_med_h2_inference_dir $visualize_trait_med_h2_dir
fi






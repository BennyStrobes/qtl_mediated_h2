#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-30:45                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5GB   

############################
# Input data
############################
# Real genotype data directories
simulation_genotype_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/genotype_processing/hm3_gwas_sample_size_100000/"

# Gencode hg19 gene annotation file
gencode_gene_annotation_file="/n/groups/price/ben/gene_annotation_files/gencode.v19.annotation.gtf.gz"

# Directory containing quasi indpendent ld blocks
quasi_independent_dir="/n/groups/price/ben/quasi_independent_ld_blocks/"

############################
# Output directories
############################
output_root_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/multi_tissue_simulation_pca/"

# Directory containing simulated gene positions
simulated_gene_position_dir=$output_root_dir"simulated_gene_positions/"

# Directory containing simulated gene expression (ie causal eQTL effect sizes in each tissue)
simulated_gene_expression_dir=$output_root_dir"simulated_gene_expression/"

# Directory containing simulated gene expression (ie causal eQTL effect sizes in each tissue)
simulated_learned_gene_models_dir=$output_root_dir"simulated_learned_gene_models/"

alt_simulated_learned_gene_models_dir="/n/groups/price/ben/qtl_mediated_h2/simulation_experiments/multi_tissue_simulation_pca/simulated_learned_gene_models_alt/"

# Directory containing simulated trait values
simulated_trait_dir=$output_root_dir"simulated_trait/"

# Directory contaiing simulated gwas results
simulated_gwas_dir=$output_root_dir"simulated_gwas/"

# Directory containing results of trait h2 inference
trait_h2_inference_dir=$output_root_dir"trait_h2_inference/"

# Directory containing results of trait med h2 inference
trait_med_h2_inference_dir=$output_root_dir"trait_mediated_h2_inference/"

# Directory containing results of expression trait h2 inference
expr_trait_h2_inference_dir=$output_root_dir"expression_trait_h2_inference/"


# Expression trait h2 visualization dir
visualize_expression_trait_h2=$output_root_dir"visualize_expression_trait_h2/"

visualize_trait_med_h2_dir=$output_root_dir"visualize_trait_med_h2/"



############################
# Simulate data
############################

chrom_num="1"

############################
# Prepare gene file for simulation:
# Genes are defined by actual tss
# Limit to protein coding genes
# In simulation, I will assume gene is expressed in each tissue
############################
simulated_gene_position_file=${simulated_gene_position_dir}"gene_positions_chr"${chrom_num}".txt"
if false; then
sh prepare_simulated_gene_position_list.sh $chrom_num $gencode_gene_annotation_file $simulated_gene_position_file $quasi_independent_dir
fi



############################
# Main simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"

# cis window arround genes to define eQTLs
cis_window="100000"

# Per genetic-element heritabilities
per_element_heritability="0.0005"

# Total heritability
total_heritability="0.3"

# Fraction of heritability mediated by gene expression
fraction_expression_mediated_heritability="0.1"

# cis window arround genes to define eQTLs
cis_window="100000"

# Gene expression heritability
ge_h2="05"

# Number of tissues
n_tissues="5"

# eQTL architecture
eqtl_architecture="default"


############################
# Run main single simulation of simulating the trait data
############################
if false; then
for simulation_number in $(seq 2 201); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}"_n_tiss_"${n_tissues}
	sbatch run_single_trait_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $simulation_genotype_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $eqtl_architecture $n_tissues $alt_simulated_learned_gene_models_dir
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






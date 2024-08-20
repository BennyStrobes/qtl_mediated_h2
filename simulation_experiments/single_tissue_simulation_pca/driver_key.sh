

############################
# Input data
############################
# Real genotype data directories
simulation_genotype_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/genotype_processing/hm3_gwas_sample_size_100000/"

# Gencode hg19 gene annotation file
gencode_gene_annotation_file="/n/groups/price/ben/gene_annotation_files/gencode.v19.annotation.gtf.gz"

############################
# Output directories
############################
output_root_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/single_tissue_simulation_pca/"

# Directory containing simulated gene positions
simulated_gene_position_dir=$output_root_dir"simulated_gene_positions/"

# Directory containing simulated gene expression (ie causal eQTL effect sizes in each tissue)
simulated_gene_expression_dir=$output_root_dir"simulated_gene_expression/"

# Directory containing simulated gene expression (ie causal eQTL effect sizes in each tissue)
simulated_learned_gene_models_dir=$output_root_dir"simulated_learned_gene_models/"

# Directory containing simulated trait values
simulated_trait_dir=$output_root_dir"simulated_trait/"

# Directory contaiing simulated gwas results
simulated_gwas_dir=$output_root_dir"simulated_gwas/"

# Directory containing results of trait h2 inference
trait_h2_inference_dir=$output_root_dir"trait_h2_inference/"

# Directory containing results of trait med h2 inference
trait_med_h2_inference_dir=$output_root_dir"trait_mediated_h2_inference/"

# Directory containing quasi indpendent ld blocks
quasi_independent_dir="/n/groups/price/ben/quasi_independent_ld_blocks/"

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


# eQTL architecture
eqtl_architecture="default"


############################
# Run main single simulation of simulating the trait data
############################
if false; then
for simulation_number in $(seq 1 30); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_trait_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $simulation_genotype_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $eqtl_architecture
done
fi
if false; then
for simulation_number in $(seq 1 20); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
	sbatch run_single_trait_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $simulation_genotype_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $eqtl_architecture
done
fi

simulation_number="1"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
sh run_single_trait_simulation.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $simulation_genotype_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $eqtl_architecture



############################
# Run trait heritability inference
############################
if false; then
for simulation_number in $(seq 1 20); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
	sbatch trait_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $n_gwas_individuals $trait_h2_inference_dir
done
fi


############################
# Run expression trait heritability inference
############################
simulation_number="1"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
if false; then
eqtl_sample_size="100"
sbatch expression_trait_h2_inference.sh $simulation_number $simulation_name_string $simulated_learned_gene_models_dir $simulation_genotype_dir $eqtl_sample_size $trait_h2_inference_dir
eqtl_sample_size="300"
sbatch expression_trait_h2_inference.sh $simulation_number $simulation_name_string $simulated_learned_gene_models_dir $simulation_genotype_dir $eqtl_sample_size $trait_h2_inference_dir
eqtl_sample_size="500"
sbatch expression_trait_h2_inference.sh $simulation_number $simulation_name_string $simulated_learned_gene_models_dir $simulation_genotype_dir $eqtl_sample_size $trait_h2_inference_dir
fi



############################
# Run expression-mediated trait heritability inference
############################
simulation_number="1"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
eqtl_sample_size="300"
if false; then
for simulation_number in $(seq 1 30); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}


	cc_hyper_param="0.0"
	sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $cc_hyper_param

	cc_hyper_param="1e-16"
	sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $cc_hyper_param

	cc_hyper_param="1e-8"
	sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $cc_hyper_param

	cc_hyper_param="1e-5"
	sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $cc_hyper_param
done
fi



############################
# Run trait heritability inference with rss likelihood
############################
if false; then
for simulation_number in $(seq 1 30); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
	sbatch trait_h2_inference_w_rss_likelihood.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $n_gwas_individuals $trait_h2_inference_dir
done
fi





############################
# Run expression-mediated trait heritability inference with rss likelihood
############################
eqtl_sample_size="300"

simulation_number="1"
if false; then
for simulation_number in $(seq 1 10); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
	resid_var_bool="True"
	sbatch trait_mediated_h2_inference_w_rss_likelihood.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $resid_var_bool
	resid_var_bool="False"
	sbatch trait_mediated_h2_inference_w_rss_likelihood.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $resid_var_bool
done
fi

if false; then
for simulation_number in $(seq 11 20); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
	resid_var_bool="True"
	sbatch trait_mediated_h2_inference_w_rss_likelihood.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $resid_var_bool
	resid_var_bool="False"
	sbatch trait_mediated_h2_inference_w_rss_likelihood.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $resid_var_bool
done
fi

if false; then
simulation_number="21"
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
	resid_var_bool="True"
	sh trait_mediated_h2_inference_w_rss_likelihood.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $resid_var_bool
fi


if false; then
simulation_number="6"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
resid_var_bool="True"
sh trait_mediated_h2_inference_w_rss_likelihood_and_known_deltas.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $resid_var_bool $simulated_gene_expression_dir
fi

if false; then
simulation_number="6"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
resid_var_bool="True"
sh trait_mediated_h2_inference_w_rss_likelihood_in_sample_ld_for_gwas_eqtl.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir $resid_var_bool $simulated_gene_expression_dir
fi

if false; then
for simulation_number in $(seq 1 10); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
	sbatch trait_mediated_h2_inference_w_rss_likelihood_no_pca.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir

done
fi

if false; then
for simulation_number in $(seq 1 20); do 
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
sbatch trait_mediated_h2_inference_w_rss_likelihood_no_pca.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir
done
fi

simulation_number="1"
if false; then
sh trait_mediated_h2_inference_w_rss_likelihood_no_pca.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $eqtl_sample_size $trait_med_h2_inference_dir

fi













############################
# Run expression-mediated trait heritability inference
############################
if false; then
simulation_number="1"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir

simulation_number="2"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir

simulation_number="3"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir

simulation_number="4"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir

simulation_number="5"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
sbatch trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir
fi
if false; then
simulation_number="1"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
sh trait_mediated_h2_inference.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $simulated_learned_gene_models_dir $n_gwas_individuals $trait_med_h2_inference_dir
fi

############################
# Run trait heritability inference with rss likelihood
############################
if false; then
for simulation_number in $(seq 2 20); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
	sbatch trait_h2_inference_w_rss_likelihood.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $n_gwas_individuals $trait_h2_inference_dir
done
fi

if false; then
simulation_number="1"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
sh trait_h2_inference_w_rss_likelihood.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $n_gwas_individuals $trait_h2_inference_dir
fi

############################
# Run expression-mediated trait heritability inference with rss likelihood
############################
if false; then
simulation_number="2"
n_eqtl_individuals="300"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}"_qtl_arch_"${eqtl_architecture}
sh trait_mediated_h2_inference_w_rss_likelihood.sh $simulation_number $simulation_name_string $simulated_trait_dir $simulated_gwas_dir $simulation_genotype_dir $n_gwas_individuals $n_eqtl_individuals $simulated_learned_gene_models_dir $trait_med_h2_inference_dir
fi








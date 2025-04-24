import numpy as np
import os
import sys
import pdb


def simulate_nm_variant_to_trait_effect_sizes(n_snps, n_causal_snps, per_element_heritability,non_mediated_variant_causal_effect_sizes_file):
	# Initialize causal effects
	causal_effects = np.zeros(n_snps)

	# Randomly assign causal effects
	causal_effects[np.random.choice(np.arange(n_snps), size=n_causal_snps, replace=False)] = np.random.normal(loc=0.0, scale=np.sqrt(per_element_heritability),size=n_causal_snps)


	# Print causal effects
	t = open(non_mediated_variant_causal_effect_sizes_file,'w')
	t.write('variant_name\tcausal_effect\n')
	for variant_iter in range(n_snps):
		t.write('variant_' + str(variant_iter) + '\t' + str(causal_effects[variant_iter]) + '\n')
	t.close()

	return causal_effects 

def print_gene_causal_effect_sizes(gene_causal_effects, mediated_gene_causal_effect_sizes_file):
	t = open(mediated_gene_causal_effect_sizes_file,'w')
	t.write('gene_name\tcausal_effect\n')
	for g_iter in range(len(gene_causal_effects)):
		t.write('gene_' + str(g_iter) + '\t' + str(gene_causal_effects[g_iter]) + '\n')
	t.close()
	return

def get_n_genes(simulated_causal_eqtl_effect_summary_file):
	aa = np.loadtxt(simulated_causal_eqtl_effect_summary_file, dtype=str,delimiter='\t')
	return aa.shape[0] -1

def extract_genetic_component(processed_genotype_data_dir, variant_effect_sizes, n_bins=100):
	n_snps = len(variant_effect_sizes)

	tmp = np.load(processed_genotype_data_dir + 'gwas_genotype_1_0.npy')
	n_samp = tmp.shape[0]

	genetic_comp = np.zeros(n_samp)

	bin_start = 0
	for bin_iter in range(n_bins):
		bin_genotype = np.load(processed_genotype_data_dir + 'gwas_genotype_1_' + str(bin_iter) + '.npy')
		bin_end = bin_start + bin_genotype.shape[1]
		genetic_comp = genetic_comp + np.dot(bin_genotype, variant_effect_sizes[bin_start:bin_end])
		bin_start = int(np.copy(bin_end)*1.0)

	return genetic_comp




####################################
# Command line arguments
####################################
simulation_name_string = sys.argv[1]
processed_genotype_data_dir = sys.argv[2]
simulated_eqtl_data_dir = sys.argv[3]
simulated_gwas_data_dir = sys.argv[4]
total_heritability = float(sys.argv[5])
fraction_expression_mediated_heritability = float(sys.argv[6])
per_element_heritability = float(sys.argv[7])


# Load in gwas genotype data
#gwas_genotype = np.load(processed_genotype_data_dir + 'gwas_genotype_1.npy')
#n_snps = gwas_genotype.shape[1]

LD = np.load(processed_genotype_data_dir+ 'gwas_genotype_LD_1.npy')
n_snps = LD.shape[0]

# Extract amount of mediated and non-mediated h2
mediated_h2 = total_heritability*fraction_expression_mediated_heritability
nm_h2 = total_heritability*(1.0-fraction_expression_mediated_heritability)

# Get number of causal genetic elements
n_causal_genes = int(mediated_h2/per_element_heritability)
n_causal_snps = int(nm_h2/per_element_heritability)

simulated_causal_eqtl_effect_summary_file = simulated_eqtl_data_dir + simulation_name_string + 'causal_eqtl_effect_summary.txt'
n_genes = get_n_genes(simulated_causal_eqtl_effect_summary_file)


# Simulated causal non-mediated variant to trait effect sizes
non_mediated_variant_causal_effect_sizes_file = simulated_gwas_data_dir + simulation_name_string + 'negative_selection_nm_variant_to_trait_effect_sizes.txt'
nm_variant_effect_sizes = simulate_nm_variant_to_trait_effect_sizes(n_snps, n_causal_snps,per_element_heritability,non_mediated_variant_causal_effect_sizes_file)


# Now add mediated effects
# There are XX heritable genes
# Randomly select which ones have causal effects
gene_causal_effects = np.zeros(n_genes)

# simulated gene causal effect sizes of genes
gene_causal_effects[np.random.choice(np.arange(n_genes), size=n_causal_genes, replace=False)] = np.random.normal(loc=0.0, scale=np.sqrt(per_element_heritability),size=n_causal_genes)


# Print gene causal effects to output file
mediated_gene_causal_effect_sizes_file = simulated_gwas_data_dir + simulation_name_string + 'negative_selection_mediated_gene_to_trait_effect_sizes.txt'
print_gene_causal_effect_sizes(gene_causal_effects, mediated_gene_causal_effect_sizes_file)


# Now add gene causal effects to genetic component of trait
mediated_variant_effect_sizes = np.zeros(n_snps)
simulated_causal_eqtl_effect_summary_file = simulated_eqtl_data_dir + simulation_name_string + 'causal_eqtl_effect_summary.txt'
f = open(simulated_causal_eqtl_effect_summary_file)
#mediated_genetic_comp = np.copy(genetic_comp)*0.0
head_count = 0
gene_counter = 0
count = 0
aa = []
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gene_id = data[0]
	gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
	gene_snp_causal_effects = np.asarray(data[2].split(',')).astype(float)
	gene_trait_causal_effect = gene_causal_effects[gene_counter]

	if gene_trait_causal_effect == 0.0:
		gene_counter = gene_counter + 1
		continue
	'''
	# Predicted gene expression
	pred_ge = np.dot(gwas_genotype[:, gene_snp_indices], gene_snp_causal_effects)

	# Standardize genetically predicted gene expression
	std_pred_ge = (pred_ge - np.mean(pred_ge))/np.std(pred_ge)

	# Add to genetic component of trait
	gene_effect = std_pred_ge*gene_trait_causal_effect
	genetic_comp = genetic_comp + gene_effect
	mediated_genetic_comp = mediated_genetic_comp + gene_effect
	'''

	gene_variance = np.dot(np.dot(gene_snp_causal_effects,LD[:, gene_snp_indices][gene_snp_indices,:]), gene_snp_causal_effects)
	mediated_variant_effect_sizes[gene_snp_indices] = mediated_variant_effect_sizes[gene_snp_indices] + gene_trait_causal_effect*gene_snp_causal_effects*(1.0/np.sqrt(gene_variance))

	gene_counter = gene_counter + 1


f.close()


nm_genetic_comp = extract_genetic_component(processed_genotype_data_dir, nm_variant_effect_sizes)
mediated_genetic_comp = extract_genetic_component(processed_genotype_data_dir, mediated_variant_effect_sizes)
genetic_comp = nm_genetic_comp + mediated_genetic_comp

total_nm_genetic_var = np.var(nm_genetic_comp)
total_mediated_genetic_var = np.var(mediated_genetic_comp)
total_genetic_var = np.var(genetic_comp)

total_genetic_var_output = simulated_gwas_data_dir + simulation_name_string + 'negative_selection_simulated_genetic_var.npy'
np.save(total_genetic_var_output, total_genetic_var)

total_nm_genetic_var_output = simulated_gwas_data_dir + simulation_name_string + 'negative_selection_simulated_nm_genetic_var.npy'
np.save(total_nm_genetic_var_output, total_nm_genetic_var)

total_mediated_genetic_var_output = simulated_gwas_data_dir + simulation_name_string + 'negative_selection_simulated_mediated_genetic_var.npy'
np.save(total_mediated_genetic_var_output, total_mediated_genetic_var)

# Residual variance
residual_var = 1.0 - total_genetic_var

# Draw trait values
trait_values = np.random.normal(loc=(genetic_comp), scale=np.sqrt(residual_var))


# Standardize trait values
standardized_trait_values = (trait_values - np.mean(trait_values))/np.std(trait_values)

# Save to output
trait_values_output = simulated_gwas_data_dir + simulation_name_string + 'negative_selection_simulated_trait.txt'
np.savetxt(trait_values_output, standardized_trait_values, fmt="%s", delimiter='\n')


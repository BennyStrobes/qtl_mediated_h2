import numpy as np 
import os
import sys
import pdb




def assign_snps_to_each_gene(n_snps, n_genes, n_snps_per_gene, gene_to_snp_mapping_file):
	# Open output file handle
	t = open(gene_to_snp_mapping_file,'w')
	# Print header of output file handle
	t.write('gene_name\tassigned_snp_indices\n')

	# Initialize dictionary to keep track of gene to snp mapping
	gene_to_snp_mapping = {}

	# Loop through genes
	for gene_iter in range(n_genes):
		# Randomly select gene position 
		gene_center = np.random.choice(np.arange(n_snps)[42:-43]) # Subsetting done to avoid weird edge cases
		# Assign snps to gene
		gene_snps = np.arange(gene_center-(n_snps_per_gene/2), gene_center+(n_snps_per_gene/2)).astype(int)

		# name of the gene
		gene_name = 'gene_' + str(gene_iter)

		# Print to output file
		t.write(gene_name + '\t' + ','.join(gene_snps.astype(str)) + '\n')

		# Add gene to dictionary of gene to snp mapping
		gene_to_snp_mapping[gene_name] = gene_snps

	# Close output file handle
	t.close()

	return gene_to_snp_mapping


def simulate_causal_eqtl_effect_sizes(gene_to_snp_mapping, simulation_name_string, simulated_causal_eqtl_effect_summary_file, fraction_genes_cis_h2, ge_h2, n_causal_snps_per_gene):
	# Number of genes
	n_genes = len(gene_to_snp_mapping)

	# Open output file handle
	t = open(simulated_causal_eqtl_effect_summary_file,'w')
	t.write('gene_id\tcis_snp_indices\tcis_snp_causal_effects\n')

	# Loop through genes
	counter = 0
	for gene_iter in range(n_genes):
		# Get name of gene
		gene_name = 'gene_' + str(gene_iter)

		# Indices of snps corresponding to the gene
		snp_indices = gene_to_snp_mapping[gene_name]

		# Simulate causal eqtl effects
		causal_eqtl_effects = np.zeros(len(snp_indices))
		if gene_iter < n_genes*fraction_genes_cis_h2:
			counter = counter + 1
			causal_eqtl_effects[np.random.choice(np.arange(len(snp_indices)), size=5, replace=False)] = np.random.normal(loc=0.0, scale=np.sqrt(ge_h2/n_causal_snps_per_gene),size=n_causal_snps_per_gene)
		else:
			print('hi')

		'''
		gene_cis_h2_boolean = np.random.binomial(n=1, p=fraction_genes_cis_h2, size=1)[0]
		if gene_cis_h2_boolean == 0:
			causal_eqtl_effects = causal_eqtl_effects*0.0
		'''

		t.write(gene_name + '\t' + str(','.join(snp_indices.astype(str))) + '\t' + str(','.join(causal_eqtl_effects.astype(str))) + '\n')
	t.close()
	return

####################################
# Command line arguments
####################################
simulation_name_string = sys.argv[1]
processed_genotype_data_dir = sys.argv[2]
simulated_eqtl_data_dir = sys.argv[3]
n_genes = int(sys.argv[4])
fraction_genes_cis_h2 = float(sys.argv[5])
ge_h2 = float(sys.argv[6])

# Number of snps assigned to each gene
n_snps_per_gene=80
n_causal_snps_per_gene=5


# Load in genotype data in order to get number of snps
gwas_genotype = np.load(processed_genotype_data_dir + 'gwas_genotype_1.npy')
n_snps = gwas_genotype.shape[1]


# Create file containing assignment of snps to each gene
gene_to_snp_mapping_file = simulated_eqtl_data_dir + simulation_name_string + 'gene_to_snp_mapping_file.txt'
gene_to_snp_mapping = assign_snps_to_each_gene(n_snps, n_genes, n_snps_per_gene, gene_to_snp_mapping_file)


############################
# Simulate causal eQTL effect sizes
############################
# Create file to keep track of causal eqtl effect sizes across genes
simulated_causal_eqtl_effect_summary_file = simulated_eqtl_data_dir + simulation_name_string + 'causal_eqtl_effect_summary.txt'
simulate_causal_eqtl_effect_sizes(gene_to_snp_mapping, simulation_name_string, simulated_causal_eqtl_effect_summary_file, fraction_genes_cis_h2, ge_h2, n_causal_snps_per_gene)













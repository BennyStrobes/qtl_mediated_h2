import numpy as np 
import os
import sys
import pdb


def simulate_gene_expression(genotype, gene_snp_indices, gene_snp_causal_effects):
	genetic_gene_expression = np.dot(genotype[:, gene_snp_indices], gene_snp_causal_effects)
	gene_heritability = np.var(genetic_gene_expression,ddof=1)

	if gene_heritability > 1:
		print('assumption eroror!')
		pdb.set_trace()

	ge = np.random.normal(loc=genetic_gene_expression, scale=np.sqrt(1.0-gene_heritability))


	return ge, genetic_gene_expression





######################
# Command line args
######################
simulation_name_string = sys.argv[1]
processed_genotype_data_dir = sys.argv[2]
simulated_eqtl_data_dir = sys.argv[3]
simulated_expression_data_dir = sys.argv[4]
eqtl_ss = sys.argv[5]


# Previously generated file
simulated_causal_eqtl_effect_file = simulated_eqtl_data_dir + simulation_name_string + 'causal_eqtl_effect_summary.txt'
genotype_file = processed_genotype_data_dir + 'eqtl_' + str(eqtl_ss) + '_genotype_1.npy'

# Load in genotype data
genotype = np.load(genotype_file)

# Loop through genes
simulated_gene_expression = []
simulated_genetic_gene_expression = []

f = open(simulated_causal_eqtl_effect_file)
head_count = 0
aa = []
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract relevent fields
	gene_name = data[0]
	gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
	gene_snp_causal_effects = np.asarray(data[2].split(',')).astype(float)

	# Simulate gene expression
	gene_expression, genetic_gene_expression = simulate_gene_expression(genotype, gene_snp_indices, gene_snp_causal_effects)

	# Standardized gene expression
	std_gene_expression = (gene_expression - np.mean(gene_expression))/np.std(gene_expression)

	# Append to global array
	simulated_gene_expression.append(std_gene_expression)
	simulated_genetic_gene_expression.append(genetic_gene_expression)
f.close()

# Convert to numpy format
simulated_gene_expression = np.asarray(simulated_gene_expression)
simulated_genetic_gene_expression = np.asarray(simulated_genetic_gene_expression)


# SAVE TO OUTPUT
output_expression_file = simulated_expression_data_dir + simulation_name_string + 'eqtl_ss_' + str(eqtl_ss) + '.npy'
np.save(output_expression_file, simulated_gene_expression)

output_genetic_expression_file = simulated_expression_data_dir + simulation_name_string + 'eqtl_ss_' + str(eqtl_ss) + '_genetic_ge.npy'
np.save(output_genetic_expression_file, simulated_genetic_gene_expression)


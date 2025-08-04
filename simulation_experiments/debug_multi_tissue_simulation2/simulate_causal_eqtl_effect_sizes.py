import sys
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from bgen import BgenReader

def generate_tissue_covariance_structure_across_causal_effects(ge_h2, n_tiss):
	ge_per_snp_h2 = ge_h2/5
	cov_mat = np.zeros((n_tiss, n_tiss)) + 0.789*ge_per_snp_h2
	np.fill_diagonal(cov_mat, ge_per_snp_h2)
	return cov_mat


def simulate_causal_eqtl_effect_sizes_across_tissues(n_cis_snps, ge_h2, n_tiss=10):
	ge_per_snp_h2 = ge_h2/5
	# Initialize matrix of causal eqtl effect sizes across tissues
	causal_effect_sizes = np.zeros((n_cis_snps, n_tiss))


	# Generate cross tissue covariance structure matrix
	tissue_covariance_mat = generate_tissue_covariance_structure_across_causal_effects(ge_h2, n_tiss)

	# First simulate 3 variants with shared effects
	shared_variant_indices = np.random.choice(np.arange(n_cis_snps), size=3, replace=False, p=None)
	causal_effect_sizes[shared_variant_indices, :] =  np.random.multivariate_normal(np.zeros(n_tiss), tissue_covariance_mat,size=3)
	
	# Now simulate 2 variants with tissue-specific effects
	remaining_variant_indices = np.delete(np.arange(n_cis_snps), shared_variant_indices)
	for tiss_iter in range(n_tiss):
		tissue_specific_indices = np.random.choice(remaining_variant_indices, size=2, replace=False, p=None)
		causal_effect_sizes[tissue_specific_indices, tiss_iter] = np.random.normal(loc=0.0, scale=np.sqrt(ge_per_snp_h2),size=2)

	return causal_effect_sizes



def simulate_causal_eqtl_effect_sizes_across_tissues_shell(n_cis_snps, fraction_genes_cis_h2, ge_h2, eqtl_architecture, n_tissues):
	min_h2 = 0
	max_h2 = .9
	while min_h2 < .01 or max_h2 > .6:
		if eqtl_architecture == 'default':
			causal_eqtl_effects = simulate_causal_eqtl_effect_sizes_across_tissues(n_cis_snps, ge_h2, n_tiss=n_tissues)
		else:
			print('assumption error: Eqtl architecture: ' + eqtl_architecture + ' not currently implemented')
			pdb.set_trace()
		min_h2 = np.min(np.sum(np.square(causal_eqtl_effects),axis=0))
		max_h2 = np.max(np.sum(np.square(causal_eqtl_effects),axis=0))
	# Make (1.0 - fraction_genes_cis_h2)% of genes not cis heritable
	gene_cis_h2_boolean = np.random.binomial(n=1, p=fraction_genes_cis_h2, size=causal_eqtl_effects.shape[1])
	for tiss_iter, boolean_value in enumerate(gene_cis_h2_boolean):
		if boolean_value == 0:
			causal_eqtl_effects[:, tiss_iter] = (causal_eqtl_effects[:, tiss_iter])*0.0
	return causal_eqtl_effects



def simulate_causal_eqtl_effect_sizes(cis_window, simulated_gene_expression_dir,simulation_name_string, ref_eqtl_sample_size, processed_genotype_data_dir, chrom_string, simulated_causal_eqtl_effect_summary_file, fraction_genes_cis_h2, ge_h2, eqtl_architecture, n_tissues):

	# Get chromosomes corresponding to chrom_string
	if chrom_string == '1_2':
		chrom_arr = np.asarray([1,2])

	t = open(simulated_causal_eqtl_effect_summary_file,'w')
	t.write('gene_id\tchr\ttss\tcausal_eqtl_effect_file\tcis_snp_id_file\tcis_snp_index_file\ttotal_n_cis_snps\tsim_ge_h2\n')

	for chrom_num in chrom_arr:

		# Load in genotype data across chromosome for single eQTL data set (note: it doesn't matter what eqtl data set we are using because we are just using snp positions here and all eqtl data sets have the same snps)
		genotype_stem = processed_genotype_data_dir + 'simulated_eqtl_' + str(ref_eqtl_sample_size) + '_data_' + str(chrom_num)
		# Load in genotype object
		bfile = BgenReader(genotype_stem + '.bgen')
		G_obj_pos = np.asarray(bfile.positions())
		G_obj_rsids = np.asarray(bfile.rsids())


		gene_summary_file = processed_genotype_data_dir + 'gene_level_ld_chr' + str(chrom_num) + '_bins_10_summary_file.txt'
		f = open(gene_summary_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			ensamble_id = data[1]
			line_chrom_num = data[0]
			if line_chrom_num != str(chrom_num):
				print('asssumptioneroror')
				pdb.set_trace()
			tss = int(data[2])

			# Get snp names in cis window fo this gene
			cis_window_start = tss - cis_window
			cis_window_end = tss + cis_window
			cis_snp_indices = (G_obj_pos >= cis_window_start) & (G_obj_pos <= cis_window_end)
			# Quick error check
			if np.sum(np.load(data[8])) != np.sum(cis_snp_indices):
				print('assumption error')
				pdb.set_trace()
			cis_rsids = G_obj_rsids[cis_snp_indices]
			n_cis_snps = len(cis_rsids)

			# Randomly simulate gene expression h2
			gene_ge_h2 = np.random.choice([ge_h2/2.0, ge_h2, ge_h2*2.0])

			# Simulate causal eqtl effects across tissues
			causal_eqtl_effects = simulate_causal_eqtl_effect_sizes_across_tissues_shell(n_cis_snps, fraction_genes_cis_h2, gene_ge_h2, eqtl_architecture, n_tissues)

			# Save results to output
			# Causal eqtl effects
			gene_causal_effect_file = simulated_gene_expression_dir + simulation_name_string + '_' + ensamble_id + '_causal_eqtl_effects.npy'
			np.save(gene_causal_effect_file, causal_eqtl_effects)
			# SNP ids
			gene_cis_snpid_file = simulated_gene_expression_dir + simulation_name_string + '_' + ensamble_id + '_cis_snpids.npy'
			np.save(gene_cis_snpid_file,cis_rsids)
			# SNP indices
			gene_cis_snp_indices_file = simulated_gene_expression_dir + simulation_name_string + '_' + ensamble_id + '_cis_snp_indices.npy'
			np.save(gene_cis_snp_indices_file, np.where(cis_snp_indices==True)[0])

			# Write to output file
			t.write(ensamble_id + '\t' + str(chrom_num) + '\t' + str(tss) + '\t' + gene_causal_effect_file + '\t' + gene_cis_snpid_file + '\t' + gene_cis_snp_indices_file + '\t' + str(n_cis_snps) + '\t' + str(gene_ge_h2) + '\n')

		f.close()
	t.close()
	print(simulated_causal_eqtl_effect_summary_file)

	return










############################
# COMMAND LINE ARGUMENTS
############################
simulation_number = int(sys.argv[1])
chrom_string = sys.argv[2]
cis_window = int(sys.argv[3])
simulated_gene_expression_dir = sys.argv[4]
simulation_name_string = sys.argv[5]
processed_genotype_data_dir = sys.argv[6]
ge_h2_str = sys.argv[7]
eqtl_architecture = sys.argv[8]
n_tissues = int(sys.argv[9])


# Get true gene expression h2
ge_h2 = float('.' + ge_h2_str)

# Set seed
np.random.seed(simulation_number + 3000)

# Define vector of eQTL sample sizes
eqtl_sample_sizes = np.asarray([100, 300, 500, 1000])
eqtl_sample_sizes = np.asarray([100])


# Fraction of genes cis heritable for a given tissue
fraction_genes_cis_h2 = 0.8


############################
# Simulate causal eQTL effect sizes across tissues
############################
# Default and random_n
# Create file to keep track of causal eqtl effect sizes across genes
simulated_causal_eqtl_effect_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'
simulate_causal_eqtl_effect_sizes(cis_window, simulated_gene_expression_dir, simulation_name_string, eqtl_sample_sizes[0], processed_genotype_data_dir, chrom_string, simulated_causal_eqtl_effect_summary_file, fraction_genes_cis_h2, ge_h2, eqtl_architecture, n_tissues)


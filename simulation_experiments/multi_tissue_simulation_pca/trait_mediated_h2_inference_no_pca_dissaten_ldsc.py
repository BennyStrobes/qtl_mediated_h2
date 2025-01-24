import sys
import numpy as np 
import pandas as pd
import os
import pdb
import tensorflow as tf
import gzip
import time
import statsmodels.api as sm
import tensorflow_probability as tfp
import time
import pickle



def load_in_gwas_data(gwas_summary_file):
	rsids = []
	betas = []
	beta_ses = []

	head_count = 0

	f = open(gwas_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count = head_count + 1
			continue

		rsids.append(data[0])
		betas.append(float(data[1]))
		beta_ses.append(float(data[2]))

	f.close()

	return np.asarray(rsids), np.asarray(betas), np.asarray(beta_ses)


def get_gwas_variant_ld_scores(gwas_beta, gwas_rsids, quasi_ld_window_summary_file):
	var_ld_score = []
	f = open(quasi_ld_window_summary_file)
	window_to_ld_files = {}
	head_count = 0
	tmp_rsids = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		window_snp_indices = np.load(data[4])
		window_rsids = np.load(data[5])
		# QUick error checking
		if np.array_equal(window_rsids, gwas_rsids[window_snp_indices]) == False:
			print('assumption eroror')
			pdb.set_trace()
		ld_file = data[3]
		ld_mat = np.load(ld_file)
		tmp_rsids.append(window_rsids)

		squared_ld_mat = np.square(ld_mat)
		squared_adj_ld_mat = squared_ld_mat - ((1.0-squared_ld_mat)/(100000-2.0))

		window_ld_scores = np.sum(squared_adj_ld_mat,axis=0)

		var_ld_score.append(window_ld_scores)

		if window_name in window_to_ld_files:
			print('assumption eororor')
			pdb.set_trace()
		window_to_ld_files[window_name] = (ld_file)


	if np.array_equal(np.hstack(tmp_rsids), gwas_rsids) == False:
		print('assumption erorror')
		pdb.set_trace()



	return np.hstack(var_ld_score), window_to_ld_files

def create_mapping_from_rsid_to_position(rsids):
	rsid_to_position = {}

	for ii, rsid in enumerate(rsids):
		if rsid in rsid_to_position:
			print('assumption eroror')
			pdb.set_trace()
		rsid_to_position[rsid] = ii
	return rsid_to_position

def extract_eqtl_sumstats_for_specific_gene(sumstats_file, gene_name):
	f = open(sumstats_file)
	sumstats = []
	sumstat_ses = []
	cis_snps = []
	window_names = []
	rsids = []
	class_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != gene_name:
			continue
		sumstats.append(data[5])
		cis_snps.append(data[4])
		sumstat_ses.append(data[6])
		window_name = data[3]
		window_names.append(window_name)
		rsids.append(data[2])
		class_names.append(data[1])
	f.close()

	unique_window_names = np.unique(window_names)
	if len(unique_window_names) != 1:
		print('assumption eroror')
		pdb.set_trace()
	class_names = np.unique(class_names)
	if len(class_names) != 1:
		print('assumption eroror')
		pdb.set_trace()

	gene_window_name = unique_window_names[0]

	return np.asarray(sumstats).astype(float), np.asarray(sumstat_ses).astype(float), np.asarray(cis_snps).astype(float), gene_window_name, np.asarray(rsids), class_names[0]


def load_in_eqtl_data_with_bins_per_gene(eqtl_sumstat_file, snp_name_to_position, window_name_to_ld_files, eqtl_sample_size, eqtl_ld='out_of_sample', n_bins=5):
	# First get list of gene names
	gene_names = []
	f = open(eqtl_sumstat_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_names.append(data[0])
	f.close()
	gene_names = np.unique(gene_names)

	print(len(gene_names))

	gene_arr = []
	beta_arr = []
	gene_class_arr = []

	beta_se_arr = []
	ldscore_arr = []
	index_arr = []
	cis_snp_index_arr = []
	n_cis_snp_arr = []
	gene_indexes = []
	# Now loop through genes
	for gg_iter, gene_name in enumerate(gene_names):
		# Extract summary stats for specific gene
		eqtl_gene_beta, eqtl_gene_beta_se, gene_cis_snp_indices, gene_window_name, gene_rsids, gene_class_name = extract_eqtl_sumstats_for_specific_gene(eqtl_sumstat_file, gene_name)

		# Load in q_mat and w_mat for this window
		ld_file = window_name_to_ld_files[gene_window_name]
		if eqtl_ld == 'out_of_sample':
			ld_mat = np.load(ld_file)
			squared_ld_mat = np.square(ld_mat)
			squared_adj_ld_mat = squared_ld_mat - ((1.0 - squared_ld_mat)/(100000.0-2.0))
			eqtl_ld_scores = np.sum(squared_adj_ld_mat[gene_cis_snp_indices==1,:], axis=0)
		elif eqtl_ld == 'in_sample_adjusted':
			new_ld_file = ld_file.split('ref_geno_gwas_')[0] + 'ref_geno_eqtl_' + str(eqtl_sample_size) + '_' + ld_file.split('ref_geno_gwas_')[1]
			ld_mat = np.load(new_ld_file)
			squared_ld_mat = np.square(ld_mat)
			squared_adj_ld_mat = squared_ld_mat - ((1.0 - squared_ld_mat)/(eqtl_sample_size-2.0))
			eqtl_ld_scores = np.sum(squared_adj_ld_mat[gene_cis_snp_indices==1,:], axis=0)


		# Get names of snps in gene
		gene_snp_positions = []
		for gene_rsid in gene_rsids:
			gene_snp_positions.append(snp_name_to_position[gene_rsid])
		gene_snp_positions = np.asarray(gene_snp_positions)

		cis_snps_boolean = gene_cis_snp_indices == 1
		temp_vec = np.arange(len(gene_cis_snp_indices))
		snp_chunks = np.array_split(temp_vec[cis_snps_boolean],n_bins)

		for bin_iter in range(n_bins):
			gene_bin_cis_snp_indices = np.zeros(len(gene_cis_snp_indices))
			gene_bin_cis_snp_indices[snp_chunks[bin_iter]] = 1

			# compute eqtl ld scores for just this bin
			eqtl_ld_scores = np.sum(squared_adj_ld_mat[gene_bin_cis_snp_indices==1,:], axis=0)

			# Add to global array
			gene_class_arr.append(gene_class_name)
			gene_arr.append(gene_name + ';' + 'bin' + str(bin_iter))
			beta_arr.append(eqtl_gene_beta)
			beta_se_arr.append(eqtl_gene_beta_se)
			ldscore_arr.append(eqtl_ld_scores)
			index_arr.append(gene_snp_positions)
			cis_snp_index_arr.append(gene_snp_positions[gene_bin_cis_snp_indices==1])
			n_cis_snp_arr.append(np.sum(gene_bin_cis_snp_indices))
			gene_indexes.append(gg_iter)

	return np.asarray(gene_arr), beta_arr, beta_se_arr, ldscore_arr, index_arr, cis_snp_index_arr, np.asarray(n_cis_snp_arr), np.asarray(gene_class_arr), np.asarray(gene_indexes)


def load_in_eqtl_data_ld_scores(eqtl_sumstat_file, snp_name_to_position, window_name_to_ld_files, eqtl_sample_size, n_total_snps, n_tissues, n_subsamples=10):
	# Initialize eQTL LD scores
	eqtl_ld_scores = np.zeros((n_total_snps, n_tissues))
	subsampled_eqtl_ld_scores = []
	for subsample_iter in range(n_subsamples):
		subsampled_eqtl_ld_scores.append((np.zeros((n_total_snps, n_tissues)), np.zeros((n_total_snps, n_tissues))))

	# Keep track of number of genes per tissue
	per_tissue_genes = []
	for tissue_iter in range(n_tissues):
		per_tissue_genes.append({})

	# Loop through eQTL Ld score file
	f = open(eqtl_sumstat_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[0]
		tissue_index = int(data[1].split('issue')[1])
		rsid = data[2]
		snp_position = snp_name_to_position[rsid]
		sq_sumstats = np.square(np.asarray(data[5:]).astype(float))

		# Add (non-subsampled eqtl ld scores)
		eqtl_ld_scores[snp_position, tissue_index] = eqtl_ld_scores[snp_position, tissue_index] + sq_sumstats[0] - sq_sumstats[1]

		per_tissue_genes[tissue_index][gene_name] = 0

		# Add subsampled eqtl ld scores
		for subsample_iter in range(n_subsamples):
			base_index = 2 + 4*subsample_iter
			if np.isnan(sq_sumstats[base_index]) or np.isnan(sq_sumstats[base_index+2]):
				if subsample_iter == 0:
					base_index = 2 + 4*1
				else:
					base_index = 2 + 4*(subsample_iter-1)
			if np.isnan(sq_sumstats[base_index]) or np.isnan(sq_sumstats[base_index+2]):
				if subsample_iter == 0:
					base_index = 2 + 4*2
				elif subsample_iter == 1:
					base_index = 2 + 4*(subsample_iter+1)
				else:
					base_index = 2 + 4*(subsample_iter-2)
			# Term a
			subsampled_eqtl_ld_scores[subsample_iter][0][snp_position, tissue_index] = subsampled_eqtl_ld_scores[subsample_iter][0][snp_position, tissue_index] + sq_sumstats[base_index] - sq_sumstats[base_index+1]
			# Term b
			subsampled_eqtl_ld_scores[subsample_iter][1][snp_position, tissue_index] = subsampled_eqtl_ld_scores[subsample_iter][1][snp_position, tissue_index] + sq_sumstats[base_index+2] - sq_sumstats[base_index+3]
	f.close()

	# Keep track of n_genes per tissues
	n_genes_per_tissue = np.zeros(n_tissues)
	for tissue_iter in range(n_tissues):
		n_genes_per_tissue[tissue_iter] = len(per_tissue_genes[tissue_iter])

	return eqtl_ld_scores, subsampled_eqtl_ld_scores, n_genes_per_tissue




def compute_correlations_in_real_vs_permuted_eqtls(eqtl_sumstat_file, perm_eqtl_sumstat_file, snp_name_to_position, window_name_to_ld_files, eqtl_sample_size,simulation_name_string, perm_simulation_name_string, simulated_gene_expression_dir, eqtl_ld='out_of_sample', n_bins=5):
	# First get list of gene names
	gene_names = []
	f = open(eqtl_sumstat_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_names.append(data[0])
	f.close()
	gene_names = np.unique(gene_names)

	gene_arr = []
	beta_arr = []
	gene_class_arr = []

	beta_se_arr = []
	ldscore_arr = []
	index_arr = []
	cis_snp_index_arr = []
	n_cis_snp_arr = []
	gene_indexes = []
	# Now loop through genes
	ldscore_beta_sq_arr = []
	R_beta_sq_arr = []
	genetic_value_corr_squared_arr = []
	genetic_corr_squared_arr = []
	R_beta_beta_t_R_corr_arr = []
	beta_sq_corr_arr = []
	for gg_iter, gene_name in enumerate(gene_names):

		if gene_name.endswith('tissue0') == False:
			continue

		# Extract summary stats for specific gene
		eqtl_gene_beta, eqtl_gene_beta_se, gene_cis_snp_indices, gene_window_name, gene_rsids, gene_class_name = extract_eqtl_sumstats_for_specific_gene(eqtl_sumstat_file, gene_name)

		# Load in q_mat and w_mat for this window
		ld_file = window_name_to_ld_files[gene_window_name]
		if eqtl_ld == 'out_of_sample':
			ld_mat = np.load(ld_file)
			squared_ld_mat = np.square(ld_mat)
			squared_adj_ld_mat = squared_ld_mat - ((1.0 - squared_ld_mat)/(100000.0-2.0))
			eqtl_ld_scores = np.sum(squared_adj_ld_mat[gene_cis_snp_indices==1,:], axis=0)
		elif eqtl_ld == 'in_sample_adjusted':
			new_ld_file = ld_file.split('ref_geno_gwas_')[0] + 'ref_geno_eqtl_' + str(eqtl_sample_size) + '_' + ld_file.split('ref_geno_gwas_')[1]
			ld_mat = np.load(new_ld_file)
			squared_ld_mat = np.square(ld_mat)
			squared_adj_ld_mat = squared_ld_mat - ((1.0 - squared_ld_mat)/(eqtl_sample_size-2.0))
			eqtl_ld_scores = np.sum(squared_adj_ld_mat[gene_cis_snp_indices==1,:], axis=0)

		# Get names of snps in gene
		gene_snp_positions = []
		for gene_rsid in gene_rsids:
			gene_snp_positions.append(snp_name_to_position[gene_rsid])
		gene_snp_positions = np.asarray(gene_snp_positions)


		ensamble_id = gene_name.split(':')[0]

		causal_eqtl_file = simulated_gene_expression_dir + simulation_name_string + '_' + ensamble_id + '_causal_eqtl_effects.npy'
		causal_eqtl = np.load(causal_eqtl_file)
		perm_causal_eqtl_file = simulated_gene_expression_dir + perm_simulation_name_string + '_' + ensamble_id + '_causal_eqtl_effects.npy'
		perm_causal_eqtl = np.load(perm_causal_eqtl_file)

		aa = np.dot(np.square(ld_mat[:, gene_cis_snp_indices==1]), np.square(causal_eqtl[:,0]))
		bb = np.dot(np.square(ld_mat[:, gene_cis_snp_indices==1]), np.square(perm_causal_eqtl[:,0]))

		cc = np.square(np.dot(ld_mat[:, gene_cis_snp_indices==1],causal_eqtl[:,0]))
		dd = np.square(np.dot(ld_mat[:, gene_cis_snp_indices==1],perm_causal_eqtl[:,0]))

		small_ld = ld_mat[gene_cis_snp_indices==1,:][:,gene_cis_snp_indices==1]

		numerator = np.dot(np.dot(causal_eqtl[:,0], small_ld), perm_causal_eqtl[:,0])
		denom = np.sqrt(np.dot(np.dot(causal_eqtl[:,0], small_ld), causal_eqtl[:,0]))*np.sqrt(np.dot(np.dot(perm_causal_eqtl[:,0], small_ld), perm_causal_eqtl[:,0]))

		tmp_mat = np.dot(np.dot(small_ld, np.diag(np.square(causal_eqtl[:,0]))), small_ld)
		tmp_mat_perm = np.dot(np.dot(small_ld, np.diag(np.square(perm_causal_eqtl[:,0]))), small_ld)

		tmp_mat2 = np.dot(np.dot(small_ld, np.dot(causal_eqtl[:,0].reshape(-1,1), causal_eqtl[:,0].reshape(1,-1))), small_ld)
		tmp_mat_perm2 = np.dot(np.dot(small_ld, np.dot(perm_causal_eqtl[:,0].reshape(-1,1), perm_causal_eqtl[:,0].reshape(1,-1))), small_ld)


		ldscore_beta_sq_arr.append(np.corrcoef(aa,bb)[0,1])
		R_beta_sq_arr.append(np.corrcoef(cc,dd)[0,1])
		genetic_value_corr_squared_arr.append(np.square(numerator/denom))
		genetic_corr_squared_arr.append(np.square(np.corrcoef(causal_eqtl[:,0], perm_causal_eqtl[:,0]))[0,1])
		R_beta_beta_t_R_corr_arr.append(np.corrcoef(tmp_mat2[np.triu_indices(tmp_mat2.shape[0])], tmp_mat_perm2[np.triu_indices(tmp_mat_perm2.shape[0])])[0,1])
		beta_sq_corr_arr.append(np.corrcoef(np.square(causal_eqtl[:,0]), np.square(perm_causal_eqtl[:,0]))[0,1])


	return np.asarray(ldscore_beta_sq_arr), np.asarray(R_beta_sq_arr), np.asarray(genetic_value_corr_squared_arr), np.asarray(genetic_corr_squared_arr), np.asarray(R_beta_beta_t_R_corr_arr), np.asarray(beta_sq_corr_arr)







#####################
# Command line args
#####################
simulation_number= int(sys.argv[1])
simulation_name_string= sys.argv[2]
simulated_trait_dir=sys.argv[3]
simulated_gwas_dir=sys.argv[4]
simulation_genotype_dir=sys.argv[5]
simulated_learned_gene_models_dir=sys.argv[6]
n_gwas_individuals=int(sys.argv[7])
eqtl_sample_size=int(sys.argv[8])
trait_med_h2_inference_dir=sys.argv[9]
simulated_gene_expression_dir=sys.argv[10]

N_gwas = n_gwas_individuals

####################
# Load in data
####################
# Load in true simulated data parameters
genetic_trait_expr_med_file = simulated_trait_dir + simulation_name_string +'_expression_mediated_trait_values.txt'
genetic_trait_nm_file = simulated_trait_dir +simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'
sim_med_h2 = np.var(np.loadtxt(genetic_trait_expr_med_file))
sim_nm_h2 = np.var(np.loadtxt(genetic_trait_nm_file))
sim_h2 = np.var(np.loadtxt(genetic_trait_nm_file) + np.loadtxt(genetic_trait_expr_med_file))

print(sim_h2)
print(sim_med_h2)
print(sim_nm_h2)
# Load in GWAS summary statistics
gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results.txt'
print(gwas_summary_file)
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(gwas_summary_file)

# Get variant ld scores
#quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_full_space_ld_ld_summary.txt'
quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_summary.txt'
#quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_summary.txt'
gwas_variant_ld_scores, window_to_ld_files = get_gwas_variant_ld_scores(gwas_beta, gwas_rsids, quasi_ld_window_summary_file)

# Create mapping from rsid to position
snp_name_to_position = create_mapping_from_rsid_to_position(gwas_rsids)



# load in eqtl data
# Out of sample eqtl ld
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats_with_subsampling.txt'
#genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, eqtl_cis_snp_position, eqtl_n_cis_snps, eqtl_classes, gene_indexer = load_in_eqtl_data_with_bins_per_gene(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, eqtl_ld='out_of_sample', n_bins=5)
eqtl_ld_scores, subsampled_eqtl_ld_scores, n_genes_per_tissue = load_in_eqtl_data_ld_scores(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files, eqtl_sample_size, len(gwas_beta),5, n_subsamples=10)


#######################################
# Run standard S-LDSC
#######################################
X = sm.add_constant(gwas_variant_ld_scores)
Y = np.square(gwas_beta/gwas_beta_se)
model = sm.OLS(Y,X).fit()
model_constrained_intercept = sm.OLS(Y-1, gwas_variant_ld_scores).fit()
ldsc_snp_h2 = len(gwas_beta)*(model.params[1]/N_gwas)
ldsc_constrained_intercept_snp_h2 = len(gwas_beta)*(model_constrained_intercept.params[0]/N_gwas)

#######################################
# Med-h2 Inference (no dissattenuation)
########################################
Y = np.square(gwas_beta) - np.square(gwas_beta_se)
X = np.hstack((gwas_variant_ld_scores.reshape(-1,1), eqtl_ld_scores))
# Manually
#model = sm.OLS(Y,X).fit()
# By hand
xtx = np.dot(np.transpose(X), X)
med_h2_no_atten_model_raw_params = np.dot(np.dot(np.transpose(X), Y), np.linalg.inv(xtx))
med_h2_no_atten_model = {}
med_h2_no_atten_model['est_nm'] = med_h2_no_atten_model_raw_params[0]*len(gwas_variant_ld_scores)
med_h2_no_atten_model['est_per_tissue_med'] = med_h2_no_atten_model_raw_params[1:]*n_genes_per_tissue*.06
med_h2_no_atten_model['est_total_med'] = np.sum(med_h2_no_atten_model['est_per_tissue_med'])
med_h2_no_atten_model['est_total_h2'] = med_h2_no_atten_model['est_total_med'] + med_h2_no_atten_model['est_nm']


#######################################
# Med-h2 Inference (dissattenuation with 10 subsamples)
########################################
n_subsamples = 10
Y = np.square(gwas_beta) - np.square(gwas_beta_se)
X = np.hstack((gwas_variant_ld_scores.reshape(-1,1), eqtl_ld_scores))
# Initialize xtx
xtx = np.zeros((X.shape[1], X.shape[1]))
# Fill in xtx
for subsample_iter in range(n_subsamples):
	X1 = np.hstack((gwas_variant_ld_scores.reshape(-1,1), subsampled_eqtl_ld_scores[subsample_iter][0]))
	X2 = np.hstack((gwas_variant_ld_scores.reshape(-1,1), subsampled_eqtl_ld_scores[subsample_iter][1]))
	xtx = xtx + np.dot(np.transpose(X1), X2)
xtx = xtx/n_subsamples
# Run regression
med_h2_atten_subsample10_model_raw_params = np.dot(np.dot(np.transpose(X), Y), np.linalg.inv(xtx))
med_h2_atten_subsample10_model = {}
med_h2_atten_subsample10_model['est_nm'] = med_h2_atten_subsample10_model_raw_params[0]*len(gwas_variant_ld_scores)
med_h2_atten_subsample10_model['est_per_tissue_med'] = med_h2_atten_subsample10_model_raw_params[1:]*n_genes_per_tissue*.06
med_h2_atten_subsample10_model['est_total_med'] = np.sum(med_h2_atten_subsample10_model['est_per_tissue_med'])
med_h2_atten_subsample10_model['est_total_h2'] = med_h2_atten_subsample10_model['est_total_med'] + med_h2_atten_subsample10_model['est_nm']


#######################################
# Med-h2 Inference (dissattenuation with 1 subsamples)
########################################
n_subsamples = 1
Y = np.square(gwas_beta) - np.square(gwas_beta_se)
X = np.hstack((gwas_variant_ld_scores.reshape(-1,1), eqtl_ld_scores))
# Initialize xtx
xtx = np.zeros((X.shape[1], X.shape[1]))
# Fill in xtx
for subsample_iter in range(n_subsamples):
	X1 = np.hstack((gwas_variant_ld_scores.reshape(-1,1), subsampled_eqtl_ld_scores[subsample_iter][0]))
	X2 = np.hstack((gwas_variant_ld_scores.reshape(-1,1), subsampled_eqtl_ld_scores[subsample_iter][1]))
	xtx = xtx + np.dot(np.transpose(X1), X2)
xtx = xtx/n_subsamples
# Run regression
med_h2_atten_subsample1_model_raw_params = np.dot(np.dot(np.transpose(X), Y), np.linalg.inv(xtx))
med_h2_atten_subsample1_model = {}
med_h2_atten_subsample1_model['est_nm'] = med_h2_atten_subsample1_model_raw_params[0]*len(gwas_variant_ld_scores)
med_h2_atten_subsample1_model['est_per_tissue_med'] = med_h2_atten_subsample1_model_raw_params[1:]*n_genes_per_tissue*.06
med_h2_atten_subsample1_model['est_total_med'] = np.sum(med_h2_atten_subsample1_model['est_per_tissue_med'])
med_h2_atten_subsample1_model['est_total_h2'] = med_h2_atten_subsample1_model['est_total_med'] + med_h2_atten_subsample1_model['est_nm']

output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + str(eqtl_sample_size) + '_dissaten_ldsc_multimethod1.txt'
t = open(output_file,'w')
t.write('method\teQTL_SS\tsim_h2\tsim_med_h2\tsim_nm_h2\test_med_h2\test_med_h2_per_tissue\test_nm_h2\test_h2_ldsc\test_h2_ldsc_constrained_intercept\n')
t.write('med_no_dissaten_ldsc\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(med_h2_no_atten_model['est_total_med']) + '\t' + ','.join(med_h2_no_atten_model['est_per_tissue_med'].astype(str)) + '\t' + str(med_h2_no_atten_model['est_nm']) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
t.write('med_dissaten_ldsc_1_subsample\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(med_h2_atten_subsample1_model['est_total_med']) + '\t' + ','.join(med_h2_atten_subsample1_model['est_per_tissue_med'].astype(str)) + '\t' + str(med_h2_atten_subsample1_model['est_nm']) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
t.write('med_dissaten_ldsc_10_subsample\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(med_h2_atten_subsample10_model['est_total_med']) + '\t' + ','.join(med_h2_atten_subsample10_model['est_per_tissue_med'].astype(str)) + '\t' + str(med_h2_atten_subsample10_model['est_nm']) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')




#######################################
# Med-h2 Inference (no dissattenuation) (1 causal tissue only)
########################################
Y = np.square(gwas_beta) - np.square(gwas_beta_se)
X = np.hstack((gwas_variant_ld_scores.reshape(-1,1), eqtl_ld_scores[:,0:1]))
# Manually
#model = sm.OLS(Y,X).fit()
# By hand
xtx = np.dot(np.transpose(X), X)
med_h2_no_atten_model_raw_params = np.dot(np.dot(np.transpose(X), Y), np.linalg.inv(xtx))
med_h2_no_atten_model = {}
med_h2_no_atten_model['est_nm'] = med_h2_no_atten_model_raw_params[0]*len(gwas_variant_ld_scores)
med_h2_no_atten_model['est_per_tissue_med'] = med_h2_no_atten_model_raw_params[1:]*n_genes_per_tissue[0]*.06
med_h2_no_atten_model['est_total_med'] = np.sum(med_h2_no_atten_model['est_per_tissue_med'])
med_h2_no_atten_model['est_total_h2'] = med_h2_no_atten_model['est_total_med'] + med_h2_no_atten_model['est_nm']



#######################################
# Med-h2 Inference (dissattenuation with 10 subsamples) (1 causal tissu eonly)
########################################
n_subsamples = 10
Y = np.square(gwas_beta) - np.square(gwas_beta_se)
X = np.hstack((gwas_variant_ld_scores.reshape(-1,1), eqtl_ld_scores[:,0:1]))
# Initialize xtx
xtx = np.zeros((X.shape[1], X.shape[1]))
# Fill in xtx
for subsample_iter in range(n_subsamples):
	X1 = np.hstack((gwas_variant_ld_scores.reshape(-1,1), subsampled_eqtl_ld_scores[subsample_iter][0][:,0:1]))
	X2 = np.hstack((gwas_variant_ld_scores.reshape(-1,1), subsampled_eqtl_ld_scores[subsample_iter][1][:,0:1]))
	xtx = xtx + np.dot(np.transpose(X1), X2)
xtx = xtx/n_subsamples
# Run regression
med_h2_atten_subsample10_model_raw_params = np.dot(np.dot(np.transpose(X), Y), np.linalg.inv(xtx))
med_h2_atten_subsample10_model = {}
med_h2_atten_subsample10_model['est_nm'] = med_h2_atten_subsample10_model_raw_params[0]*len(gwas_variant_ld_scores)
med_h2_atten_subsample10_model['est_per_tissue_med'] = med_h2_atten_subsample10_model_raw_params[1:]*n_genes_per_tissue[0]*.06
med_h2_atten_subsample10_model['est_total_med'] = np.sum(med_h2_atten_subsample10_model['est_per_tissue_med'])
med_h2_atten_subsample10_model['est_total_h2'] = med_h2_atten_subsample10_model['est_total_med'] + med_h2_atten_subsample10_model['est_nm']

#######################################
# Med-h2 Inference (dissattenuation with 1 subsamples) ( 1 causal tissue only)
########################################
n_subsamples = 1
Y = np.square(gwas_beta) - np.square(gwas_beta_se)
X = np.hstack((gwas_variant_ld_scores.reshape(-1,1), eqtl_ld_scores[:,0:1]))
# Initialize xtx
xtx = np.zeros((X.shape[1], X.shape[1]))
# Fill in xtx
for subsample_iter in range(n_subsamples):
	X1 = np.hstack((gwas_variant_ld_scores.reshape(-1,1), subsampled_eqtl_ld_scores[subsample_iter][0][:,0:1]))
	X2 = np.hstack((gwas_variant_ld_scores.reshape(-1,1), subsampled_eqtl_ld_scores[subsample_iter][1][:,0:1]))
	xtx = xtx + np.dot(np.transpose(X1), X2)
xtx = xtx/n_subsamples
# Run regression
med_h2_atten_subsample1_model_raw_params = np.dot(np.dot(np.transpose(X), Y), np.linalg.inv(xtx))
med_h2_atten_subsample1_model = {}
med_h2_atten_subsample1_model['est_nm'] = med_h2_atten_subsample1_model_raw_params[0]*len(gwas_variant_ld_scores)
med_h2_atten_subsample1_model['est_per_tissue_med'] = med_h2_atten_subsample1_model_raw_params[1:]*n_genes_per_tissue[0]*.06
med_h2_atten_subsample1_model['est_total_med'] = np.sum(med_h2_atten_subsample1_model['est_per_tissue_med'])
med_h2_atten_subsample1_model['est_total_h2'] = med_h2_atten_subsample1_model['est_total_med'] + med_h2_atten_subsample1_model['est_nm']


t.write('med_no_dissaten_ldsc_1_tissue\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(med_h2_no_atten_model['est_total_med']) + '\t' + ','.join(med_h2_no_atten_model['est_per_tissue_med'].astype(str)) + '\t' + str(med_h2_no_atten_model['est_nm']) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
t.write('med_dissaten_ldsc_1_subsample_1_tissue\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(med_h2_atten_subsample1_model['est_total_med']) + '\t' + ','.join(med_h2_atten_subsample1_model['est_per_tissue_med'].astype(str)) + '\t' + str(med_h2_atten_subsample1_model['est_nm']) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
t.write('med_dissaten_ldsc_10_subsample_1_tissue\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(med_h2_atten_subsample10_model['est_total_med']) + '\t' + ','.join(med_h2_atten_subsample10_model['est_per_tissue_med'].astype(str)) + '\t' + str(med_h2_atten_subsample10_model['est_nm']) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

t.close()

print(output_file)

'''
##########################
# Temp saving
np.save('gwas_beta.npy', gwas_beta)
np.save('gwas_beta_se.npy', gwas_beta_se)
np.save('gwas_variant_ld_scores.npy', gwas_variant_ld_scores)
np.save('eqtl_beta_mat.npy', eqtl_beta_mat)
np.save('eqtl_beta_se_mat.npy', eqtl_beta_se_mat)
np.save('eqtl_ld_score_mat.npy', eqtl_ld_score_mat)
np.save('eqtl_n_cis_snps.npy', eqtl_n_cis_snps)
np.save('eqtl_classes.npy', eqtl_classes)
np.save('gene_indexer.npy', gene_indexer)
f = open('eqtl_cis_snp_position.pkl', 'wb')
pickle.dump(eqtl_cis_snp_position, f, pickle.HIGHEST_PROTOCOL)
f.close()
##########################
# temp loading
gwas_beta = np.load('gwas_beta.npy')
gwas_beta_se = np.load('gwas_beta_se.npy')
gwas_variant_ld_scores = np.load('gwas_variant_ld_scores.npy')
eqtl_beta_mat = np.load('eqtl_beta_mat.npy')
eqtl_beta_se_mat = np.load('eqtl_beta_se_mat.npy')
eqtl_ld_score_mat = np.load('eqtl_ld_score_mat.npy')
eqtl_n_cis_snps = np.load('eqtl_n_cis_snps.npy')
eqtl_classes = np.load('eqtl_classes.npy')
gene_indexer = np.load('gene_indexer.npy')
f = open('eqtl_cis_snp_position.pkl', 'rb')
eqtl_cis_snp_position = pickle.load(f)
f.close()
'''

'''
# Run standard S-LDSC
X = sm.add_constant(gwas_variant_ld_scores)
Y = np.square(gwas_beta/gwas_beta_se)
model = sm.OLS(Y,X).fit()
model_constrained_intercept = sm.OLS(Y-1, gwas_variant_ld_scores).fit()
ldsc_snp_h2 = len(gwas_beta)*(model.params[1]/N_gwas)
ldsc_constrained_intercept_snp_h2 = len(gwas_beta)*(model_constrained_intercept.params[0]/N_gwas)



output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod12.txt'
t = open(output_file,'w')
t.write('method\teQTL_SS\tsim_h2\tsim_med_h2\tsim_nm_h2\test_med_h2_joint_reml\test_med_h2_per_tissue_joint_reml\test_nm_h2_joint_reml\test_mean_eqtl_h2_joint_reml\test_h2_ldsc\test_h2_ldsc_constrained_intercept\n')


est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X_init = med_h2_with_sumstat_ldsc_two_step_multivariate_gene_bins(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, burn_in_iters=10, total_iters=20, gibbs=False)
t.write('eqtl_5_binned_no_intercept_two_step\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_data_set_variance(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=995, total_iters=1000, gibbs=False, eqtl_only=False)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
output_file2 = trait_med_h2_inference_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod12_X.npy'
np.save(output_file2, X)


# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_data_set_variance(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=995, total_iters=1000, gibbs=False, eqtl_only=True)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance_cis_eqtl_only\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


t.close()
print(output_file)
'''



#####################
# Run analysis permuting the eqtls we use
#####################
'''
# load in eqtl data
# Out of sample eqtl ld
perm_simulation_name_string = 'simulation_' + str(simulation_number+1) + '_chrom' + simulation_name_string.split('_chrom')[1]
eqtl_sumstat_file = simulated_learned_gene_models_dir + perm_simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, eqtl_cis_snp_position, eqtl_n_cis_snps, eqtl_classes, gene_indexer = load_in_eqtl_data_with_bins_per_gene(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, eqtl_ld='out_of_sample', n_bins=5)
# Generate matrix form of eqtl data
eqtl_beta_mat = get_matrix_form_of_eqtl_data(eqtl_beta, eqtl_position, 0.0, len(gwas_variant_ld_scores))
eqtl_beta_se_mat = get_matrix_form_of_eqtl_data(eqtl_beta_se, eqtl_position, 1.0/np.sqrt(eqtl_sample_size), len(gwas_variant_ld_scores))
eqtl_ld_score_mat = get_matrix_form_of_eqtl_data(eqtl_ldscore, eqtl_position, 0.0, len(gwas_variant_ld_scores))



# Run standard S-LDSC
X = sm.add_constant(gwas_variant_ld_scores)
Y = np.square(gwas_beta/gwas_beta_se)
model = sm.OLS(Y,X).fit()
model_constrained_intercept = sm.OLS(Y-1, gwas_variant_ld_scores).fit()
ldsc_snp_h2 = len(gwas_beta)*(model.params[1]/N_gwas)
ldsc_constrained_intercept_snp_h2 = len(gwas_beta)*(model_constrained_intercept.params[0]/N_gwas)



output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod12_permuted_eqtls.txt'
t = open(output_file,'w')
t.write('method\teQTL_SS\tsim_h2\tsim_med_h2\tsim_nm_h2\test_med_h2_joint_reml\test_med_h2_per_tissue_joint_reml\test_nm_h2_joint_reml\test_mean_eqtl_h2_joint_reml\test_h2_ldsc\test_h2_ldsc_constrained_intercept\n')


est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X_init = med_h2_with_sumstat_ldsc_two_step_multivariate_gene_bins(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, burn_in_iters=10, total_iters=20, gibbs=False)
t.write('eqtl_5_binned_no_intercept_two_step\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_data_set_variance(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=995, total_iters=1000, gibbs=False, eqtl_only=False)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
output_file2 = trait_med_h2_inference_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod12_permuted_eqtls_X.npy'
np.save(output_file2, X)


# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_data_set_variance(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=995, total_iters=1000, gibbs=False, eqtl_only=True)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance_cis_eqtl_only\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


t.close()
print(output_file)
'''












########################################
#*************************************
# Only single causal tissue included
########################################
'''
# load in eqtl data
# Out of sample eqtl ld
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, eqtl_cis_snp_position, eqtl_n_cis_snps, eqtl_classes, gene_indexer = load_in_eqtl_data_with_bins_per_gene(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, eqtl_ld='out_of_sample', n_bins=5)
# Generate matrix form of eqtl data
eqtl_beta_mat = get_matrix_form_of_eqtl_data(eqtl_beta, eqtl_position, 0.0, len(gwas_variant_ld_scores))
eqtl_beta_se_mat = get_matrix_form_of_eqtl_data(eqtl_beta_se, eqtl_position, 1.0/np.sqrt(eqtl_sample_size), len(gwas_variant_ld_scores))
eqtl_ld_score_mat = get_matrix_form_of_eqtl_data(eqtl_ldscore, eqtl_position, 0.0, len(gwas_variant_ld_scores))

###########################################
# Remove all tissues except causal tissue
###########################################
indices = eqtl_classes == 'tissue0'
causal_eqtl_beta_mat = eqtl_beta_mat[indices,:]
causal_eqtl_beta_se_mat = eqtl_beta_se_mat[indices,:]
causal_eqtl_ld_score_mat = eqtl_ld_score_mat[indices,:]
causal_eqtl_n_cis_snps = eqtl_n_cis_snps[indices]
causal_eqtl_classes = eqtl_classes[indices]
causal_gene_indexer = gene_indexer[indices]


# Run standard S-LDSC
X = sm.add_constant(gwas_variant_ld_scores)
Y = np.square(gwas_beta/gwas_beta_se)
model = sm.OLS(Y,X).fit()
model_constrained_intercept = sm.OLS(Y-1, gwas_variant_ld_scores).fit()
ldsc_snp_h2 = len(gwas_beta)*(model.params[1]/N_gwas)
ldsc_constrained_intercept_snp_h2 = len(gwas_beta)*(model_constrained_intercept.params[0]/N_gwas)



output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + str(eqtl_sample_size) + '_joint_ldsc_single_tissue_multimethod12.txt'
t = open(output_file,'w')
t.write('method\teQTL_SS\tsim_h2\tsim_med_h2\tsim_nm_h2\test_med_h2_joint_reml\test_med_h2_per_tissue_joint_reml\test_nm_h2_joint_reml\test_mean_eqtl_h2_joint_reml\test_h2_ldsc\test_h2_ldsc_constrained_intercept\n')


est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X_init = med_h2_with_sumstat_ldsc_two_step_multivariate_gene_bins(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, causal_eqtl_beta_mat, causal_eqtl_beta_se_mat, causal_eqtl_ld_score_mat, causal_eqtl_n_cis_snps, causal_eqtl_classes, causal_gene_indexer, eqtl_sample_size, burn_in_iters=10, total_iters=20, gibbs=False)
t.write('eqtl_5_binned_no_intercept_two_step\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_data_set_variance(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, causal_eqtl_beta_mat, causal_eqtl_beta_se_mat, causal_eqtl_ld_score_mat, causal_eqtl_n_cis_snps, causal_eqtl_classes, causal_gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=995, total_iters=1000, gibbs=False, eqtl_only=False)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

t.close()
print(output_file)



#####################
# Run analysis permuting the eqtls we use
#####################
# load in eqtl data
# Out of sample eqtl ld
perm_simulation_name_string = 'simulation_' + str(simulation_number+1) + '_chrom' + simulation_name_string.split('_chrom')[1]
eqtl_sumstat_file = simulated_learned_gene_models_dir + perm_simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, eqtl_cis_snp_position, eqtl_n_cis_snps, eqtl_classes, gene_indexer = load_in_eqtl_data_with_bins_per_gene(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, eqtl_ld='out_of_sample', n_bins=5)
# Generate matrix form of eqtl data
eqtl_beta_mat = get_matrix_form_of_eqtl_data(eqtl_beta, eqtl_position, 0.0, len(gwas_variant_ld_scores))
eqtl_beta_se_mat = get_matrix_form_of_eqtl_data(eqtl_beta_se, eqtl_position, 1.0/np.sqrt(eqtl_sample_size), len(gwas_variant_ld_scores))
eqtl_ld_score_mat = get_matrix_form_of_eqtl_data(eqtl_ldscore, eqtl_position, 0.0, len(gwas_variant_ld_scores))


###########################################
# Remove all tissues except causal tissue
###########################################
indices = eqtl_classes == 'tissue0'
causal_eqtl_beta_mat = eqtl_beta_mat[indices,:]
causal_eqtl_beta_se_mat = eqtl_beta_se_mat[indices,:]
causal_eqtl_ld_score_mat = eqtl_ld_score_mat[indices,:]
causal_eqtl_n_cis_snps = eqtl_n_cis_snps[indices]
causal_eqtl_classes = eqtl_classes[indices]
causal_gene_indexer = gene_indexer[indices]

# Run standard S-LDSC
X = sm.add_constant(gwas_variant_ld_scores)
Y = np.square(gwas_beta/gwas_beta_se)
model = sm.OLS(Y,X).fit()
model_constrained_intercept = sm.OLS(Y-1, gwas_variant_ld_scores).fit()
ldsc_snp_h2 = len(gwas_beta)*(model.params[1]/N_gwas)
ldsc_constrained_intercept_snp_h2 = len(gwas_beta)*(model_constrained_intercept.params[0]/N_gwas)


output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + str(eqtl_sample_size) + '_joint_ldsc_single_tissue_multimethod12_permuted_eqtls.txt'
t = open(output_file,'w')
t.write('method\teQTL_SS\tsim_h2\tsim_med_h2\tsim_nm_h2\test_med_h2_joint_reml\test_med_h2_per_tissue_joint_reml\test_nm_h2_joint_reml\test_mean_eqtl_h2_joint_reml\test_h2_ldsc\test_h2_ldsc_constrained_intercept\n')


est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X_init = med_h2_with_sumstat_ldsc_two_step_multivariate_gene_bins(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, causal_eqtl_beta_mat, causal_eqtl_beta_se_mat, causal_eqtl_ld_score_mat, causal_eqtl_n_cis_snps, causal_eqtl_classes, causal_gene_indexer, eqtl_sample_size, burn_in_iters=10, total_iters=20, gibbs=False)
t.write('eqtl_5_binned_no_intercept_two_step\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_data_set_variance(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, causal_eqtl_beta_mat, causal_eqtl_beta_se_mat, causal_eqtl_ld_score_mat, causal_eqtl_n_cis_snps, causal_eqtl_classes, causal_gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=995, total_iters=1000, gibbs=False, eqtl_only=False)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

t.close()
print(output_file)





#####################
# Emperically evaluate correlation between genes across real and permuted runs
#####################
perm_simulation_name_string = 'simulation_' + str(simulation_number+1) + '_chrom' + simulation_name_string.split('_chrom')[1]
perm_eqtl_sumstat_file = simulated_learned_gene_models_dir + perm_simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
ldscore_beta_sq_arr, R_beta_sq_arr, genetic_value_corr_squared_arr,genetic_corr_squared_arr, R_beta_beta_t_R_corr_arr, beta_sq_corr_arr = compute_correlations_in_real_vs_permuted_eqtls(eqtl_sumstat_file, perm_eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, simulation_name_string, perm_simulation_name_string, simulated_gene_expression_dir, eqtl_ld='out_of_sample', n_bins=5)



output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + str(eqtl_sample_size) + '_joint_ldsc_single_tissue_permuted_eqtls_correlation_values.txt'
t = open(output_file,'w')
t.write('method\tcorrelation\n')

t.write('ldscore_beta_sq\t' + str(np.nanmean(ldscore_beta_sq_arr)) + '\n')
ldscore_beta_sq_arr[np.isnan(ldscore_beta_sq_arr)] = 0.0
t.write('ldscore_beta_sq_nan_zero\t' + str(np.mean(ldscore_beta_sq_arr)) + '\n')

t.write('beta_sq_corr\t' + str(np.nanmean(beta_sq_corr_arr)) + '\n')
beta_sq_corr_arr[np.isnan(beta_sq_corr_arr)] = 0.0
t.write('beta_sq_corr_nan_zero\t' + str(np.mean(beta_sq_corr_arr)) + '\n')

t.write('R_beta_sq\t' + str(np.nanmean(R_beta_sq_arr)) + '\n')
R_beta_sq_arr[np.isnan(R_beta_sq_arr)] = 0.0
t.write('R_beta_sq_nan_zero\t' + str(np.mean(R_beta_sq_arr)) + '\n')

t.write('genetic_value_corr_sq\t' + str(np.nanmean(genetic_value_corr_squared_arr)) + '\n')
genetic_value_corr_squared_arr[np.isnan(genetic_value_corr_squared_arr)] = 0.0
t.write('genetic_value_corr_sq_nan_zero\t' + str(np.mean(genetic_value_corr_squared_arr)) + '\n')


t.write('genetic_corr_sq\t' + str(np.nanmean(genetic_corr_squared_arr)) + '\n')
genetic_corr_squared_arr[np.isnan(genetic_corr_squared_arr)] = 0.0
t.write('genetic_corr_sq_nan_zero\t' + str(np.mean(genetic_corr_squared_arr)) + '\n')

t.write('R_beta_beta_t_R_corr\t' + str(np.nanmean(R_beta_beta_t_R_corr_arr)) + '\n')
R_beta_beta_t_R_corr_arr[np.isnan(R_beta_beta_t_R_corr_arr)] = 0.0
t.write('R_beta_beta_t_R_corr_nan_zero\t' + str(np.mean(R_beta_beta_t_R_corr_arr)) + '\n')

t.close()
print(output_file)

'''




























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


def get_gwas_gene_window_variant_ld_scores(gwas_beta, gwas_rsids, quasi_ld_window_summary_file, single_bin_genes, eqtl_cis_snp_position_single_bin):
	anno_vec = np.zeros(len(gwas_beta))
	used_genes = {}
	for gene_iter, full_gene_name in enumerate(single_bin_genes):
		gene_name = full_gene_name.split(':')[0]
		if gene_name in used_genes:
			continue
		used_genes[gene_name] = 1
		anno_vec[eqtl_cis_snp_position_single_bin[gene_iter]] = anno_vec[eqtl_cis_snp_position_single_bin[gene_iter]] + 1
	
	gene_window_ld_scores = []
	f = open(quasi_ld_window_summary_file)
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
		squared_ld_mat = np.square(ld_mat)
		squared_adj_ld_mat = squared_ld_mat - ((1.0-squared_ld_mat)/(100000-2.0))

		window_ld_scores = np.sum(squared_adj_ld_mat*anno_vec[window_snp_indices],axis=1)

		gene_window_ld_scores.append(window_ld_scores)


	return np.hstack(gene_window_ld_scores), anno_vec

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



def load_in_eqtl_data(eqtl_sumstat_file, snp_name_to_position, window_name_to_ld_files, eqtl_sample_size, eqtl_ld='out_of_sample'):
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
	# Now loop through genes
	for gene_name in gene_names:
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
		elif eqtl_ld == 'in_sample_unadjusted':
			new_ld_file = ld_file.split('ref_geno_gwas_')[0] + 'ref_geno_eqtl_' + str(eqtl_sample_size) + '_' + ld_file.split('ref_geno_gwas_')[1]
			ld_mat = np.load(new_ld_file)
			squared_ld_mat = np.square(ld_mat)
			eqtl_ld_scores = np.sum(squared_ld_mat[gene_cis_snp_indices==1,:], axis=0)

		# Get names of snps in gene
		gene_snp_positions = []
		for gene_rsid in gene_rsids:
			gene_snp_positions.append(snp_name_to_position[gene_rsid])
		gene_snp_positions = np.asarray(gene_snp_positions)

		# Add to global array
		gene_class_arr.append(gene_class_name)
		gene_arr.append(gene_name)
		beta_arr.append(eqtl_gene_beta)
		beta_se_arr.append(eqtl_gene_beta_se)
		ldscore_arr.append(eqtl_ld_scores)
		index_arr.append(gene_snp_positions)
		cis_snp_index_arr.append(gene_snp_positions[gene_cis_snp_indices==1])
		n_cis_snp_arr.append(np.sum(gene_cis_snp_indices))

	return np.asarray(gene_arr), beta_arr, beta_se_arr, ldscore_arr, index_arr, cis_snp_index_arr, np.asarray(n_cis_snp_arr), np.asarray(gene_class_arr)


def get_matrix_form_of_eqtl_data(input_data_arr, input_data_positions, null_value, n_snps):
	n_genes = len(input_data_arr)
	output_mat = np.zeros((n_genes, n_snps)) + null_value


	for gene_iter, gene_arr in enumerate(input_data_arr):
		gene_input_position = input_data_positions[gene_iter]

		output_mat[gene_iter, gene_input_position] = gene_arr


	return output_mat


def sumstat_proportional_joint_reml_loss_softplus_link(gwas_beta_sq, gwas_beta_var, gwas_ld_scores, eqtl_beta_sq, eqtl_beta_var, eqtl_ld_scores,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, gene_to_class_index_matrix, gwas_noise, eqtl_noise):


	eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable

	per_gene_alpha_sq = tf.linalg.matmul(gene_to_class_index_matrix, tf.reshape(alpha_sq_variable, [-1,1]))[:,0]

	big_eqtl_beta_sq_variable1 = eqtl_ld_scores*(tf.tile(tf.reshape(per_gene_alpha_sq*eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps]))
	big_eqtl_beta_sq_variable2 = eqtl_ld_scores*(tf.tile(tf.reshape(eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps]))

	pred_per_snp_gwas_h2 = (gwas_ld_scores*beta_sq_variable) + (tf.math.reduce_sum((big_eqtl_beta_sq_variable1)*eqtl_mask,axis=0))

	gwas_log_like = compute_univariate_gaussian_log_like(gwas_beta_sq, gwas_beta_var + pred_per_snp_gwas_h2,  tf.math.softplus(gwas_noise))


	eqtl_indices = eqtl_mask==1.0
	eqtl_log_like = compute_univariate_gaussian_log_like(eqtl_beta_sq[eqtl_indices], eqtl_beta_var[eqtl_indices] + big_eqtl_beta_sq_variable2[eqtl_indices], tf.math.softplus(eqtl_noise))

	loss = -tf.reduce_sum(gwas_log_like) - tf.reduce_sum(eqtl_log_like)

	return loss


def sumstat_proportional_joint_reml_loss_w_intercept(gwas_beta_sq, gwas_beta_var, gwas_ld_scores, eqtl_beta_sq, eqtl_beta_var, eqtl_ld_scores,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, gene_to_class_index_matrix, gwas_noise, eqtl_noise, gwas_intercept, eqtl_intercept):


	eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable

	per_gene_alpha_sq = tf.linalg.matmul(gene_to_class_index_matrix, tf.reshape(alpha_sq_variable, [-1,1]))[:,0]

	big_eqtl_beta_sq_variable1 = eqtl_ld_scores*(tf.tile(tf.reshape(per_gene_alpha_sq*eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps]))
	big_eqtl_beta_sq_variable2 = eqtl_ld_scores*(tf.tile(tf.reshape(eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps]))

	pred_per_snp_gwas_h2 = (gwas_ld_scores*beta_sq_variable) + (tf.math.reduce_sum((big_eqtl_beta_sq_variable1)*eqtl_mask,axis=0))

	gwas_log_like = compute_univariate_gaussian_log_like(gwas_beta_sq, gwas_intercept + gwas_beta_var + pred_per_snp_gwas_h2,  gwas_noise)


	eqtl_indices = eqtl_mask==1.0
	eqtl_log_like = compute_univariate_gaussian_log_like(eqtl_beta_sq[eqtl_indices], eqtl_intercept + eqtl_beta_var[eqtl_indices] + big_eqtl_beta_sq_variable2[eqtl_indices], eqtl_noise)

	loss = -tf.reduce_sum(gwas_log_like) - tf.reduce_sum(eqtl_log_like)

	return loss



def sumstat_proportional_joint_reml_loss(gwas_beta_sq, gwas_beta_var, gwas_ld_scores, eqtl_beta_sq, eqtl_beta_var, eqtl_ld_scores,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, gene_to_class_index_matrix, gwas_noise, eqtl_noise):


	eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable

	per_gene_alpha_sq = tf.linalg.matmul(gene_to_class_index_matrix, tf.reshape(alpha_sq_variable, [-1,1]))[:,0]

	big_eqtl_beta_sq_variable1 = eqtl_ld_scores*(tf.tile(tf.reshape(per_gene_alpha_sq*eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps]))
	big_eqtl_beta_sq_variable2 = eqtl_ld_scores*(tf.tile(tf.reshape(eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps]))

	pred_per_snp_gwas_h2 = (gwas_ld_scores*beta_sq_variable) + (tf.math.reduce_sum((big_eqtl_beta_sq_variable1)*eqtl_mask,axis=0))

	gwas_log_like = compute_univariate_gaussian_log_like(gwas_beta_sq, gwas_beta_var + pred_per_snp_gwas_h2,  gwas_noise)


	eqtl_indices = eqtl_mask==1.0
	eqtl_log_like = compute_univariate_gaussian_log_like(eqtl_beta_sq[eqtl_indices], eqtl_beta_var[eqtl_indices] + big_eqtl_beta_sq_variable2[eqtl_indices], eqtl_noise)

	loss = -tf.reduce_sum(gwas_log_like) - tf.reduce_sum(eqtl_log_like)

	return loss

def compute_univariate_gaussian_log_like(xx, mean_value, variance_vec):
	log_like = -(tf.math.log(variance_vec)/2) - tf.math.divide(tf.square(xx-mean_value), (2.0*variance_vec))
	
	return log_like


def med_h2_with_sumstat_ldsc_w_intercept(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, max_epochs=40000, conv_thresh=1e-12, intercept_variables=False, learning_rate=5e-6):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)


	optimizer = tf.keras.optimizers.Adam(learning_rate=1e-5)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)


	# Convert variabels to tf tensors
	gwas_beta_sq = tf.convert_to_tensor(np.square(gwas_beta), dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	gwas_ld_scores = tf.convert_to_tensor(gwas_ld_scores, dtype=tf.float32)
	eqtl_beta_sq = tf.convert_to_tensor(np.square(eqtl_beta), dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)
	eqtl_ld_scores = tf.convert_to_tensor(eqtl_ld_scores, dtype=tf.float32)
	gene_to_class_index = tf.convert_to_tensor(gene_to_class_index.astype(int))
	gene_to_class_index_matrix = tf.convert_to_tensor(gene_to_class_index_matrix, dtype=tf.float32)




	# Initialize variables to optimize over
	beta_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='beta_sq', dtype=tf.float32)
	alpha_sq_variable = tf.Variable(initial_value=np.ones(n_eqtl_classes)*0.0000001,trainable=True, name='alpha_sq',dtype=tf.float32)
	eqtl_beta_sq_variable = tf.Variable(initial_value=np.ones(n_genes)*0.000000001,trainable=True, name='eqtl_beta_sq', dtype=tf.float32)

	gwas_noise = tf.Variable(initial_value=1e-2,trainable=True, name='gwas_noise')
	eqtl_noise = tf.Variable(initial_value=1e-2,trainable=True, name='eqtl_noise')

	gwas_intercept = tf.Variable(initial_value=0.0,trainable=True, name='gwas_intercept')
	eqtl_intercept = tf.Variable(initial_value=0.0,trainable=True, name='eqtl_intercept')


	converged = False
	prev_est_alpha_sq=10000
	best_loss = 1e10
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = sumstat_proportional_joint_reml_loss_w_intercept(gwas_beta_sq, gwas_beta_var,gwas_ld_scores, eqtl_beta_sq, eqtl_beta_var,eqtl_ld_scores,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, gene_to_class_index_matrix, gwas_noise, eqtl_noise, gwas_intercept, eqtl_intercept)

		trainable_variables = []
		trainable_variables.append(alpha_sq_variable)
		trainable_variables.append(eqtl_beta_sq_variable)
		trainable_variables.append(beta_sq_variable)
		trainable_variables.append(gwas_intercept)
		trainable_variables.append(eqtl_intercept)
		#if intercept_variables:
			#trainable_variables.append(gwas_noise)
			#trainable_variables.append(eqtl_noise)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est = np.asarray(alpha_sq_variable)*1.0

		diff = np.abs(prev_est_alpha_sq -cur_est)
		'''
		if diff < conv_thresh:
			converged = True
			break
		'''

		prev_est_alpha_sq = cur_est


		cur_loss = np.asmatrix(loss_value)[0,0]
		if cur_loss < best_loss:
			best_loss = cur_loss
			eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable
			best_gene_cis_h2 = np.asarray(eqtl_beta_sq_variable_scaled)*eqtl_n_cis_snps
			per_gene_alpha_sq_variable = np.asarray(alpha_sq_variable)[gene_to_class_index]
			per_gene_med_h2 = per_gene_alpha_sq_variable*best_gene_cis_h2
			best_per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
			best_total_med_h2 = np.sum(per_gene_med_h2)
			best_nm_h2 = np.sum(beta_sq_variable)*n_snps


		if np.mod(epoch_iter, 100) == 0.0:
			eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable 
			gene_cis_h2 = np.asarray(eqtl_beta_sq_variable_scaled)*eqtl_n_cis_snps
			per_gene_alpha_sq_variable = np.asarray(alpha_sq_variable)[gene_to_class_index]
			per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
			per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
			total_med_h2 = np.sum(per_gene_med_h2)
			nm_h2 = np.sum(beta_sq_variable)*n_snps
			print(epoch_iter)
			print('loss: ' + str(loss_value))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('gwas_noise: ' + str(np.asmatrix(gwas_noise)[0,0]))
			print('eqtl noise: ' + str(np.asmatrix(eqtl_noise)[0,0]))
			print('gwas_intercept: ' + str(np.asmatrix(gwas_intercept)[0,0]))
			print('eqtl intercept: ' + str(np.asmatrix(eqtl_intercept)[0,0]))

	print(best_loss)


	return best_total_med_h2, best_per_class_med_h2, best_nm_h2, best_gene_cis_h2



def med_h2_with_sumstat_ldsc_bayesian_gibbs_per_data_set_noise(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, eqtl_cis_snp_position, burn_in_iters=4600, total_iters=5000):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)
	class_indices = np.sort(np.unique(gene_to_class_index))

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1

	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2

	gwas_resid_var = 1e-4
	eqtl_resid_var = np.ones(n_eqtl_classes)*1e-4

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		params = np.random.multivariate_normal(mean=mean, cov=S)
		sig_sq = params[0]
		alpha_sq = params[1:]
		# Run OLS
		#model = sm.OLS(Y,X).fit()
		#sig_sq = model.params[0]
		#alpha_sq = model.params[1:]

		###############################################################
		# Second update psi_g
		###############################################################
		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - np.dot(X, params)
		for gene_iter in np.random.permutation(np.arange(n_genes)):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid + ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
			gwas_beta_sq_X = eqtl_ld_scores[gene_iter,:]*gene_alpha_sq
			gwas_beta_sq_resid_var = np.ones(len(gwas_beta_sq_X))*gwas_resid_var

			gene_indices = eqtl_beta[gene_iter,:] != 0.0
			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]
			eqtl_beta_sq_resid_var = np.ones(len(eqtl_beta_sq_X))*eqtl_resid_var[gene_to_class_index[gene_iter]]

			gene_Y = np.hstack((gwas_beta_sq_resid, eqtl_beta_sq_resid))
			gene_X = np.hstack((gwas_beta_sq_X, eqtl_beta_sq_X))
			gene_resid_var = np.hstack((gwas_beta_sq_resid_var, eqtl_beta_sq_resid_var))

			# Get sampling distribution
			S = 1.0/np.sum(np.square(gene_X)/gene_resid_var)
			mean = S*np.sum(gene_X*gene_Y/gene_resid_var)

			# Sample
			psi_g[gene_iter] = np.random.normal(loc=mean, scale=np.sqrt(S))

			# Remove current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)


		###############################################################
		# Update residual variances
		###############################################################
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)
		eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))
		masker = eqtl_mask==1

		for class_index in class_indices:
			gene_indices = gene_to_class_index==class_index
			class_eqtl_beta_sq_resid = (eqtl_beta_sq_resid[gene_indices, :])[masker[gene_indices,:]]
			eqtl_resid_var[class_index] = np.sum(np.square(class_eqtl_beta_sq_resid))/len(class_eqtl_beta_sq_resid)



		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 20) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + ','.join(eqtl_resid_var.astype(str)))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)

	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2



def med_h2_with_sumstat_ldsc_bayesian_gibbs_with_intercepts(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, eqtl_cis_snp_position, burn_in_iters=4600, total_iters=5000):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1

	# Initialize intercept variables
	gwas_intercept_variable = 0.0
	eqtl_intercept_variables = np.zeros(n_genes)

	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2

	gwas_resid_var = 1e-4
	eqtl_resid_var = 1e-4

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(np.ones(len(gwas_ld_scores)).reshape(1,-1)))
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		params = np.random.multivariate_normal(mean=mean, cov=S)
		gwas_intercept_variable = params[0]
		sig_sq = params[1]
		alpha_sq = params[2:]


		###############################################################
		# Second update psi_g
		###############################################################
		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - np.dot(X, params)
		for gene_iter in np.random.permutation(np.arange(n_genes)):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid + ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
			gwas_beta_sq_X = eqtl_ld_scores[gene_iter,:]*gene_alpha_sq
			gwas_beta_sq_resid_var = np.ones(len(gwas_beta_sq_X))*gwas_resid_var

			gene_indices = eqtl_beta[gene_iter,:] != 0.0
			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices]) - eqtl_intercept_variables[gene_iter]
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]
			eqtl_beta_sq_resid_var = np.ones(len(eqtl_beta_sq_X))*eqtl_resid_var

			gene_Y = np.hstack((gwas_beta_sq_resid, eqtl_beta_sq_resid))
			gene_X = np.hstack((gwas_beta_sq_X, eqtl_beta_sq_X))
			gene_resid_var = np.hstack((gwas_beta_sq_resid_var, eqtl_beta_sq_resid_var))

			# Get sampling distribution
			S = 1.0/np.sum(np.square(gene_X)/gene_resid_var)
			mean = S*np.sum(gene_X*gene_Y/gene_resid_var)

			# Sample
			psi_g[gene_iter] = np.random.normal(loc=mean, scale=np.sqrt(S))

			# Remove current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)

		###############################################################
		# Update eqtl intercepts
		###############################################################
		for gene_iter in np.random.permutation(np.arange(n_genes)):
			
			gene_indices = eqtl_beta[gene_iter,:] != 0.0
			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices]) - (psi_g[gene_iter]*eqtl_ld_scores[gene_iter, gene_indices])
			eqtl_beta_sq_X = np.ones(len(eqtl_beta_sq_resid))
			eqtl_beta_sq_resid_var = np.ones(len(eqtl_beta_sq_X))*eqtl_resid_var

			# Get sampling distribution
			S = 1.0/np.sum(np.square(eqtl_beta_sq_X)/eqtl_beta_sq_resid_var)
			mean = S*np.sum(eqtl_beta_sq_X*eqtl_beta_sq_resid/eqtl_beta_sq_resid_var)

			eqtl_intercept_variables[gene_iter] = np.random.normal(loc=mean, scale=np.sqrt(S))


		###############################################################
		# Update residual variances
		###############################################################
		tmp_beta_sq = np.transpose(np.transpose(np.square(eqtl_beta)) - eqtl_intercept_variables)
		eqtl_beta_sq_resid = (tmp_beta_sq - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
		eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)


		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 20) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(eqtl_resid_var))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)

	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2


def med_h2_with_sumstat_ldsc_bayesian_fixed_sample_size_resid_var(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, eqtl_sample_size, burn_in_iters=4600, total_iters=5000, gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	'''
	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1
	'''

	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2


	gwas_resid_var = np.square(gwas_beta)
	eqtl_resid_vars = np.square(eqtl_beta)
	gwas_resid_var = gwas_resid_var*0.0 + np.square(1.0/100000)
	eqtl_resid_vars = eqtl_resid_vars*0.0 + np.square(1.0/eqtl_sample_size)


	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X)/gwas_resid_var, X))
		mean = np.dot(S, np.dot(Y/gwas_resid_var,X))
		# Draw from distribution
		#params = np.copy(mean)
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]


		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = np.dot(alpha_weighted_eqtl_ld_scores/gwas_resid_var, np.transpose(alpha_weighted_eqtl_ld_scores))
		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])

		for gene_iter in np.arange(n_genes):

			gene_indices = eqtl_beta[gene_iter,:] != 0.0
			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]

			xt_x_term[gene_iter, gene_iter] = xt_x_term[gene_iter, gene_iter] + np.dot(eqtl_beta_sq_X/eqtl_resid_vars[gene_iter, gene_indices], eqtl_beta_sq_X)

			xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid/eqtl_resid_vars[gene_iter, gene_indices], eqtl_beta_sq_X)) + np.dot(gwas_beta_sq_resid/gwas_resid_var, alpha_weighted_eqtl_ld_scores[gene_iter,:])

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)


		'''
		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)

		eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))
		for gene_iter in np.arange(n_genes):
			gene_eqtl_beta_sq_resid = eqtl_beta_sq_resid[gene_iter,eqtl_mask[gene_iter,:] == 1]
			eqtl_resid_vars[gene_iter]= np.sum(np.square(gene_eqtl_beta_sq_resid))/len(gene_eqtl_beta_sq_resid)

		#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
		#eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)
		'''

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 20) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(np.mean(eqtl_resid_vars)))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2


def med_h2_with_sumstat_ldsc_bayesian_fixed_sample_size_plus_h2_resid_var_only_cis_eqtl_non_central_chi_sq(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=4600, total_iters=5000, gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1

	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2


	gwas_resid_var = np.square(gwas_beta)
	eqtl_resid_vars = np.square(eqtl_beta)
	# variance is (2 + 4*mean^2/variance)*variance^2
	gwas_mean_sq = (gwas_ld_scores*.3/n_snps)
	gwas_variance = np.square(gwas_beta_se)
	gwas_resid_var = gwas_resid_var*0.0 + (2.0 + (4.0*gwas_mean_sq/gwas_variance))*np.square(gwas_variance)

	eqtl_mean_sq = np.transpose(np.transpose(eqtl_ld_scores)*.05/eqtl_n_cis_snps)
	eqtl_variance = np.square(eqtl_beta_se)
	eqtl_resid_vars = eqtl_resid_vars*0.0 + (2.0 + (4.0*eqtl_mean_sq/eqtl_variance))*np.square(eqtl_variance)


	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X)/gwas_resid_var, X))
		mean = np.dot(S, np.dot(Y/gwas_resid_var,X))
		# Draw from distribution
		#params = np.copy(mean)
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]


		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = np.dot(alpha_weighted_eqtl_ld_scores/gwas_resid_var, np.transpose(alpha_weighted_eqtl_ld_scores))
		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])

		for gene_iter in np.arange(n_genes):

			#gene_indices = eqtl_beta[gene_iter,:] != 0.0
			gene_indices = cis_eqtl_mask[gene_iter,:] == 1.0

			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]

			xt_x_term[gene_iter, gene_iter] = xt_x_term[gene_iter, gene_iter] + np.dot(eqtl_beta_sq_X/eqtl_resid_vars[gene_iter, gene_indices], eqtl_beta_sq_X)

			xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid/eqtl_resid_vars[gene_iter, gene_indices], eqtl_beta_sq_X)) + np.dot(gwas_beta_sq_resid/gwas_resid_var, alpha_weighted_eqtl_ld_scores[gene_iter,:])

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)


		'''
		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)

		eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))
		for gene_iter in np.arange(n_genes):
			gene_eqtl_beta_sq_resid = eqtl_beta_sq_resid[gene_iter,eqtl_mask[gene_iter,:] == 1]
			eqtl_resid_vars[gene_iter]= np.sum(np.square(gene_eqtl_beta_sq_resid))/len(gene_eqtl_beta_sq_resid)

		#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
		#eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)
		'''

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 20) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(np.mean(eqtl_resid_vars)))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2


def med_h2_with_sumstat_ldsc_bayesian_fixed_sample_size_plus_h2_resid_var_only_cis_eqtl(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=4600, total_iters=5000, gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1

	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2


	gwas_resid_var = np.square(gwas_beta)
	eqtl_resid_vars = np.square(eqtl_beta)
	gwas_resid_var = gwas_resid_var*0.0 + np.square((1.0/100000) +  (gwas_ld_scores*.3/n_snps))
	eqtl_resid_vars = eqtl_resid_vars*0.0 + np.square((1.0/eqtl_sample_size) + np.transpose(np.transpose(eqtl_ld_scores)*.05/eqtl_n_cis_snps))


	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X)/gwas_resid_var, X))
		mean = np.dot(S, np.dot(Y/gwas_resid_var,X))
		# Draw from distribution
		#params = np.copy(mean)
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]


		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = np.dot(alpha_weighted_eqtl_ld_scores/gwas_resid_var, np.transpose(alpha_weighted_eqtl_ld_scores))
		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])

		for gene_iter in np.arange(n_genes):

			#gene_indices = eqtl_beta[gene_iter,:] != 0.0
			gene_indices = cis_eqtl_mask[gene_iter,:] == 1.0

			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]

			xt_x_term[gene_iter, gene_iter] = xt_x_term[gene_iter, gene_iter] + np.dot(eqtl_beta_sq_X/eqtl_resid_vars[gene_iter, gene_indices], eqtl_beta_sq_X)

			xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid/eqtl_resid_vars[gene_iter, gene_indices], eqtl_beta_sq_X)) + np.dot(gwas_beta_sq_resid/gwas_resid_var, alpha_weighted_eqtl_ld_scores[gene_iter,:])

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)


		'''
		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)

		eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))
		for gene_iter in np.arange(n_genes):
			gene_eqtl_beta_sq_resid = eqtl_beta_sq_resid[gene_iter,eqtl_mask[gene_iter,:] == 1]
			eqtl_resid_vars[gene_iter]= np.sum(np.square(gene_eqtl_beta_sq_resid))/len(gene_eqtl_beta_sq_resid)

		#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
		#eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)
		'''

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 20) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(np.mean(eqtl_resid_vars)))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2

def med_h2_with_sumstat_ldsc_bayesian_fixed_sample_size_plus_h2_resid_var_non_central_chi_sq(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, eqtl_sample_size, burn_in_iters=4600, total_iters=5000, gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	'''
	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1
	'''

	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2


	gwas_resid_var = np.square(gwas_beta)
	eqtl_resid_vars = np.square(eqtl_beta)


	gwas_mean_sq = (gwas_ld_scores*.3/n_snps)
	gwas_variance = np.square(gwas_beta_se)
	#gwas_variance = 1.0/100000.0
	gwas_resid_var = gwas_resid_var*0.0 + (2.0 + (4.0*gwas_mean_sq/gwas_variance))*np.square(gwas_variance)

	eqtl_mean_sq = np.transpose(np.transpose(eqtl_ld_scores)*.05/eqtl_n_cis_snps)
	eqtl_variance = 1.0/eqtl_sample_size
	eqtl_variance = np.square(eqtl_beta_se)
	eqtl_resid_vars = eqtl_resid_vars*0.0 + (2.0 + (4.0*eqtl_mean_sq/eqtl_variance))*np.square(eqtl_variance)



	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X)/gwas_resid_var, X))
		mean = np.dot(S, np.dot(Y/gwas_resid_var,X))
		# Draw from distribution
		#params = np.copy(mean)
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]


		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = np.dot(alpha_weighted_eqtl_ld_scores/gwas_resid_var, np.transpose(alpha_weighted_eqtl_ld_scores))
		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])

		for gene_iter in np.arange(n_genes):

			gene_indices = eqtl_beta[gene_iter,:] != 0.0
			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]

			xt_x_term[gene_iter, gene_iter] = xt_x_term[gene_iter, gene_iter] + np.dot(eqtl_beta_sq_X/eqtl_resid_vars[gene_iter, gene_indices], eqtl_beta_sq_X)

			xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid/eqtl_resid_vars[gene_iter, gene_indices], eqtl_beta_sq_X)) + np.dot(gwas_beta_sq_resid/gwas_resid_var, alpha_weighted_eqtl_ld_scores[gene_iter,:])

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)


		'''
		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)

		eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))
		for gene_iter in np.arange(n_genes):
			gene_eqtl_beta_sq_resid = eqtl_beta_sq_resid[gene_iter,eqtl_mask[gene_iter,:] == 1]
			eqtl_resid_vars[gene_iter]= np.sum(np.square(gene_eqtl_beta_sq_resid))/len(gene_eqtl_beta_sq_resid)

		#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
		#eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)
		'''

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 20) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(np.mean(eqtl_resid_vars)))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2

def med_h2_with_sumstat_ldsc_bayesian_fixed_sample_size_plus_h2_resid_var(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, eqtl_sample_size, burn_in_iters=4600, total_iters=5000, gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	'''
	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1
	'''

	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2


	gwas_resid_var = np.square(gwas_beta)
	eqtl_resid_vars = np.square(eqtl_beta)
	gwas_resid_var = gwas_resid_var*0.0 + np.square((1.0/100000) +  (gwas_ld_scores*.3/n_snps))
	eqtl_resid_vars = eqtl_resid_vars*0.0 + np.square((1.0/eqtl_sample_size) + np.transpose(np.transpose(eqtl_ld_scores)*.05/eqtl_n_cis_snps))


	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X)/gwas_resid_var, X))
		mean = np.dot(S, np.dot(Y/gwas_resid_var,X))
		# Draw from distribution
		#params = np.copy(mean)
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]


		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = np.dot(alpha_weighted_eqtl_ld_scores/gwas_resid_var, np.transpose(alpha_weighted_eqtl_ld_scores))
		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])

		for gene_iter in np.arange(n_genes):

			gene_indices = eqtl_beta[gene_iter,:] != 0.0
			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]

			xt_x_term[gene_iter, gene_iter] = xt_x_term[gene_iter, gene_iter] + np.dot(eqtl_beta_sq_X/eqtl_resid_vars[gene_iter, gene_indices], eqtl_beta_sq_X)

			xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid/eqtl_resid_vars[gene_iter, gene_indices], eqtl_beta_sq_X)) + np.dot(gwas_beta_sq_resid/gwas_resid_var, alpha_weighted_eqtl_ld_scores[gene_iter,:])

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)


		'''
		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)

		eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))
		for gene_iter in np.arange(n_genes):
			gene_eqtl_beta_sq_resid = eqtl_beta_sq_resid[gene_iter,eqtl_mask[gene_iter,:] == 1]
			eqtl_resid_vars[gene_iter]= np.sum(np.square(gene_eqtl_beta_sq_resid))/len(gene_eqtl_beta_sq_resid)

		#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
		#eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)
		'''

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 20) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(np.mean(eqtl_resid_vars)))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2


def med_h2_with_sumstat_ldsc_bayesian_per_gene_resid_var_only_cis_eqtl_snps(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, eqtl_cis_snp_position, burn_in_iters=4600, total_iters=5000, gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1

	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2


	gwas_resid_var = 1e-4
	eqtl_resid_vars = np.ones(n_genes)*1e-4

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		#params = np.copy(mean)
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]


		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = np.dot(alpha_weighted_eqtl_ld_scores, np.transpose(alpha_weighted_eqtl_ld_scores))/gwas_resid_var
		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])

		for gene_iter in np.arange(n_genes):

			#gene_indices = eqtl_beta[gene_iter,:] != 0.0
			gene_indices = cis_eqtl_mask[gene_iter,:] == 1
			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]

			xt_x_term[gene_iter, gene_iter] = xt_x_term[gene_iter, gene_iter] + np.dot(eqtl_beta_sq_X, eqtl_beta_sq_X)/eqtl_resid_vars[gene_iter]

			xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid, eqtl_beta_sq_X)/eqtl_resid_vars[gene_iter]) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_iter,:])/gwas_resid_var

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)


		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)

		eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))
		for gene_iter in np.arange(n_genes):
			gene_eqtl_beta_sq_resid = eqtl_beta_sq_resid[gene_iter,cis_eqtl_mask[gene_iter,:] == 1]
			eqtl_resid_vars[gene_iter]= np.sum(np.square(gene_eqtl_beta_sq_resid))/len(gene_eqtl_beta_sq_resid)

		#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
		#eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)


		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 20) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(np.mean(eqtl_resid_vars)))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2

def med_h2_with_sumstat_ldsc_bayesian_per_gene_resid_var(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, burn_in_iters=4600, total_iters=5000, gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	'''
	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1
	'''

	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2


	gwas_resid_var = 1e-4
	eqtl_resid_vars = np.ones(n_genes)*1e-4

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		#params = np.copy(mean)
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]


		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = np.dot(alpha_weighted_eqtl_ld_scores, np.transpose(alpha_weighted_eqtl_ld_scores))/gwas_resid_var
		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])

		for gene_iter in np.arange(n_genes):

			gene_indices = eqtl_beta[gene_iter,:] != 0.0
			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]

			xt_x_term[gene_iter, gene_iter] = xt_x_term[gene_iter, gene_iter] + np.dot(eqtl_beta_sq_X, eqtl_beta_sq_X)/eqtl_resid_vars[gene_iter]

			xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid, eqtl_beta_sq_X)/eqtl_resid_vars[gene_iter]) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_iter,:])/gwas_resid_var

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)


		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)

		eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))
		for gene_iter in np.arange(n_genes):
			gene_eqtl_beta_sq_resid = eqtl_beta_sq_resid[gene_iter,eqtl_mask[gene_iter,:] == 1]
			eqtl_resid_vars[gene_iter]= np.sum(np.square(gene_eqtl_beta_sq_resid))/len(gene_eqtl_beta_sq_resid)

		#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
		#eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)


		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 20) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(np.mean(eqtl_resid_vars)))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2

def med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, burn_in_iters=4600, total_iters=5000,gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	'''
	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1
	'''
	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2

	'''
	######***********#########
	# Organize data
	eqtl_block_beta_sq = []
	eqtl_block_intercept = []
	eqtl_block_eqtl_ld_scores = []
	for gene_iter in range(n_genes):
		gene_indices = eqtl_beta[gene_iter,:] != 0.0
		eqtl_block_beta_sq.append(np.square(eqtl_beta[gene_iter, gene_indices]))
		eqtl_block_intercept.append(np.square(eqtl_beta_se[gene_iter, gene_indices]))
		eqtl_block_ld_score = np.zeros((np.sum(gene_indices), n_genes))
		eqtl_block_ld_score[:, gene_iter] = eqtl_ld_scores[gene_iter, gene_indices]
		eqtl_block_eqtl_ld_scores.append(eqtl_block_ld_score)
	eqtl_block_beta_sq = np.hstack(eqtl_block_beta_sq)
	eqtl_block_intercept = np.hstack(eqtl_block_intercept)
	eqtl_block_eqtl_ld_scores = np.vstack(eqtl_block_eqtl_ld_scores)
	#####************#########
	'''

	precomputed_xt_x = np.dot(eqtl_ld_scores, np.transpose(eqtl_ld_scores))

	gwas_resid_var = 1e-4
	eqtl_resid_var = 1e-4

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]

		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = (precomputed_xt_x*np.dot(alpha_sq[gene_to_class_index].reshape(-1,1), alpha_sq[gene_to_class_index].reshape(1,-1)))/gwas_resid_var

		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])

		for gene_iter in np.arange(n_genes):

			gene_indices = eqtl_beta[gene_iter,:] != 0.0
			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]

			xt_x_term[gene_iter, gene_iter] = xt_x_term[gene_iter, gene_iter] + np.dot(eqtl_beta_sq_X, eqtl_beta_sq_X)/eqtl_resid_var

			xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid, eqtl_beta_sq_X)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_iter,:])/gwas_resid_var
		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		'''
		# Error checking
		if itera > 3:
			# V2
			Y_top = gwas_block_beta_sq - gwas_block_intercept - (gwas_ld_scores*sig_sq)
			Y_bottom = eqtl_block_beta_sq - eqtl_block_intercept
			X_top = np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index]
			# Note that X_bottom is simply eqtl_block_eqtl_ld_scores
			Y = np.hstack((Y_top, Y_bottom))
			X = np.vstack((X_top, eqtl_block_eqtl_ld_scores))

			weights = np.hstack((np.ones(len(Y_top))*gwas_resid_var, np.ones(len(Y_bottom))*eqtl_resid_var))

			mod_wls = sm.WLS(Y, X, weights=1.0 / weights).fit()
			pdb.set_trace()
		'''
		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)
		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)


		eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
		eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 1) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(eqtl_resid_var))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2



def med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_data_set_variance(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, nm_snp_count, burn_in_iters=4600, total_iters=5000,gibbs=True, eqtl_only=False):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1

	unique_gene_indices = np.sort(np.unique(gene_indexer))
	gene_index_to_gene_bin_pair_indices = []
	for unique_gene_index in unique_gene_indices:
		gene_index_to_gene_bin_pair_indices.append(np.where(gene_indexer==unique_gene_index)[0])


	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)

	for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
		representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
		snp_indices = eqtl_ld_scores[representative_gene_bin_pair_index,:] != 0.0 
		Y = np.square(eqtl_beta[representative_gene_bin_pair_index, snp_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, snp_indices])
		X = np.transpose(eqtl_ld_scores[:, snp_indices][gene_bin_pair_indices,:])

		model = sm.OLS(Y,X).fit()
		#cis_snp_h2 = model.params[0]
		psi_g[gene_bin_pair_indices] = model.params


	precomputed_xt_x = np.dot(eqtl_ld_scores, np.transpose(eqtl_ld_scores))

	gwas_resid_var = np.square(1/100000.0)
	eqtl_resid_var = np.ones(n_eqtl_classes)*np.square(1.0/eqtl_sample_size)

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		print(itera)
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(gwas_ld_scores)
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)

		n_nm_vars = gwas_ld_scores.shape[1]
		sig_sq = params[0:n_nm_vars]
		alpha_sq = params[n_nm_vars:]
		#sig_sq = params[0]
		#alpha_sq = params[1:]




		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = (precomputed_xt_x*np.dot(alpha_sq[gene_to_class_index].reshape(-1,1), alpha_sq[gene_to_class_index].reshape(1,-1)))/gwas_resid_var

		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - np.dot(gwas_ld_scores, sig_sq)


		for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
			representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same

			variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0

			tmp_eqtl_class = gene_to_class_index[representative_gene_bin_pair_index]
			
			eqtl_beta_sq_resid = np.square(eqtl_beta[representative_gene_bin_pair_index, variant_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, variant_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_bin_pair_indices, :][:, variant_indices]

			tmp_gene_mat = np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_var[tmp_eqtl_class]

			
			for ii, gene_iter in enumerate(gene_bin_pair_indices):
				xt_x_term[gene_iter, gene_bin_pair_indices] = xt_x_term[gene_iter, gene_bin_pair_indices] + tmp_gene_mat[ii,:]
			
			#xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] = xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] + np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_var

			#xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid, eqtl_beta_sq_X)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_iter,:])/gwas_resid_var
			#xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:])/gwas_resid_var
			xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_var[tmp_eqtl_class]) + np.dot(alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:], gwas_beta_sq_resid)/gwas_resid_var

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)
		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)

		for eqtl_class_index in range(n_eqtl_classes):
			eqtl_beta_sq_resid = []
			for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
				representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
				tmp_eqtl_class = gene_to_class_index[representative_gene_bin_pair_index]
				if eqtl_class_index != tmp_eqtl_class:
					continue
				if eqtl_only:
					variant_indices = np.sum(cis_eqtl_mask[gene_bin_pair_indices,:],axis=0) == 1.0
				else:
					variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
				tmp_eqtl_resid = (np.square(eqtl_beta[representative_gene_bin_pair_index, :]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, :])) - np.dot(np.transpose(eqtl_ld_scores[gene_bin_pair_indices,:]), psi_g[gene_bin_pair_indices])
				eqtl_beta_sq_resid.append(tmp_eqtl_resid[variant_indices])
			eqtl_beta_sq_resid = np.hstack(eqtl_beta_sq_resid)
			eqtl_resid_var[eqtl_class_index] = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = np.dot(nm_snp_count, sig_sq)
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 1) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(np.mean(eqtl_resid_var)))
			print(eqtl_resid_var)
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)


	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2, X


def med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, burn_in_iters=4600, total_iters=5000,gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	'''
	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1
	'''

	unique_gene_indices = np.sort(np.unique(gene_indexer))
	gene_index_to_gene_bin_pair_indices = []
	for unique_gene_index in unique_gene_indices:
		gene_index_to_gene_bin_pair_indices.append(np.where(gene_indexer==unique_gene_index)[0])


	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)

	for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
		representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
		snp_indices = eqtl_ld_scores[representative_gene_bin_pair_index,:] != 0.0 
		Y = np.square(eqtl_beta[representative_gene_bin_pair_index, snp_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, snp_indices])
		X = np.transpose(eqtl_ld_scores[:, snp_indices][gene_bin_pair_indices,:])

		model = sm.OLS(Y,X).fit()
		#cis_snp_h2 = model.params[0]
		psi_g[gene_bin_pair_indices] = model.params


	precomputed_xt_x = np.dot(eqtl_ld_scores, np.transpose(eqtl_ld_scores))

	gwas_resid_var = np.square(1/100000.0)
	eqtl_resid_var = np.square(1.0/eqtl_sample_size)

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		print(itera)
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]




		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = (precomputed_xt_x*np.dot(alpha_sq[gene_to_class_index].reshape(-1,1), alpha_sq[gene_to_class_index].reshape(1,-1)))/gwas_resid_var

		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])


		for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
			representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same

			variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
			
			eqtl_beta_sq_resid = np.square(eqtl_beta[representative_gene_bin_pair_index, variant_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, variant_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_bin_pair_indices, :][:, variant_indices]

			tmp_gene_mat = np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_var

			
			for ii, gene_iter in enumerate(gene_bin_pair_indices):
				xt_x_term[gene_iter, gene_bin_pair_indices] = xt_x_term[gene_iter, gene_bin_pair_indices] + tmp_gene_mat[ii,:]
			
			#xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] = xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] + np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_var

			#xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid, eqtl_beta_sq_X)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_iter,:])/gwas_resid_var
			#xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:])/gwas_resid_var
			xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_var) + np.dot(alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:], gwas_beta_sq_resid)/gwas_resid_var

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)
		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)

		eqtl_beta_sq_resid = []
		for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
			representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
			variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
			tmp_eqtl_resid = (np.square(eqtl_beta[representative_gene_bin_pair_index, :]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, :])) - np.dot(np.transpose(eqtl_ld_scores[gene_bin_pair_indices,:]), psi_g[gene_bin_pair_indices])
			eqtl_beta_sq_resid.append(tmp_eqtl_resid[variant_indices])

		#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
		eqtl_beta_sq_resid = np.hstack(eqtl_beta_sq_resid)
		eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 1) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(eqtl_resid_var))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)


	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2, X


def med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_single_init_variance(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size,eqtl_cis_snp_position, burn_in_iters=4600, total_iters=5000,gibbs=True, eqtl_only=False):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1

	unique_gene_indices = np.sort(np.unique(gene_indexer))
	gene_index_to_gene_bin_pair_indices = []
	for unique_gene_index in unique_gene_indices:
		gene_index_to_gene_bin_pair_indices.append(np.where(gene_indexer==unique_gene_index)[0])


	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)

	for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
		representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
		snp_indices = eqtl_ld_scores[representative_gene_bin_pair_index,:] != 0.0 
		Y = np.square(eqtl_beta[representative_gene_bin_pair_index, snp_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, snp_indices])
		X = np.transpose(eqtl_ld_scores[:, snp_indices][gene_bin_pair_indices,:])

		model = sm.OLS(Y,X).fit()
		#cis_snp_h2 = model.params[0]
		psi_g[gene_bin_pair_indices] = model.params


	precomputed_xt_x = np.dot(eqtl_ld_scores, np.transpose(eqtl_ld_scores))

	gwas_resid_var = np.square(1/100000.0)

	'''
	gwas_resid_var = np.square(1/eqtl_sample_size)

	eqtl_resid_var = np.square(1.0/eqtl_sample_size)
	'''

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	###############################################################
	# First update sig_sq and alpha_sq (based two step approach)
	###############################################################
	# Prepare data
	Y = gwas_block_beta_sq - gwas_block_intercept
	X = []
	X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
	X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
	X = np.hstack(X)
	S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
	mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
	# Draw from distribution
	if gibbs:
		params = np.random.multivariate_normal(mean=mean, cov=S)
	else:
		params = np.copy(mean)
	sig_sq = params[0]
	alpha_sq = params[1:]



	###############################################################
	# Set residual variances
	###############################################################
	gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])

	for gene_iter in np.arange(n_genes):
		gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
		# Re-include current effect
		gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
	gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)

	eqtl_beta_sq_resid = []
	for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
		representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
		if eqtl_only:
			variant_indices = np.sum(cis_eqtl_mask[gene_bin_pair_indices,:],axis=0) == 1.0
		else:
			variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
		tmp_eqtl_resid = (np.square(eqtl_beta[representative_gene_bin_pair_index, :]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, :])) - np.dot(np.transpose(eqtl_ld_scores[gene_bin_pair_indices,:]), psi_g[gene_bin_pair_indices])
		eqtl_beta_sq_resid.append(tmp_eqtl_resid[variant_indices])

	#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
	eqtl_beta_sq_resid = np.hstack(eqtl_beta_sq_resid)
	eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)



	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		print(itera)
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]




		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = (precomputed_xt_x*np.dot(alpha_sq[gene_to_class_index].reshape(-1,1), alpha_sq[gene_to_class_index].reshape(1,-1)))/gwas_resid_var

		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])


		for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
			representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same

			variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
			
			eqtl_beta_sq_resid = np.square(eqtl_beta[representative_gene_bin_pair_index, variant_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, variant_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_bin_pair_indices, :][:, variant_indices]

			tmp_gene_mat = np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_var

			
			for ii, gene_iter_alt in enumerate(gene_bin_pair_indices):
				xt_x_term[gene_iter_alt, gene_bin_pair_indices] = xt_x_term[gene_iter_alt, gene_bin_pair_indices] + tmp_gene_mat[ii,:]
			
			#xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] = xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] + np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_var

			#xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid, eqtl_beta_sq_X)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_iter,:])/gwas_resid_var
			#xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:])/gwas_resid_var
			xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_var) + np.dot(alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:], gwas_beta_sq_resid)/gwas_resid_var

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)
		###############################################################
		# Update residual variances
		###############################################################
		'''
		if itera == 0:
			for gene_iter in np.arange(n_genes):
				gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
				# Re-include current effect
				gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
			gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)

			eqtl_beta_sq_resid = []
			for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
				representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
				variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
				tmp_eqtl_resid = (np.square(eqtl_beta[representative_gene_bin_pair_index, :]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, :])) - np.dot(np.transpose(eqtl_ld_scores[gene_bin_pair_indices,:]), psi_g[gene_bin_pair_indices])
				eqtl_beta_sq_resid.append(tmp_eqtl_resid[variant_indices])

			#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
			eqtl_beta_sq_resid = np.hstack(eqtl_beta_sq_resid)
			eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		'''

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 1) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(eqtl_resid_var))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)


	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2, X




def med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_single_init_variance_per_gene(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size,eqtl_cis_snp_position, burn_in_iters=4600, total_iters=5000,gibbs=True, eqtl_only=False):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1

	unique_gene_indices = np.sort(np.unique(gene_indexer))
	gene_index_to_gene_bin_pair_indices = []
	for unique_gene_index in unique_gene_indices:
		gene_index_to_gene_bin_pair_indices.append(np.where(gene_indexer==unique_gene_index)[0])


	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)

	for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
		representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
		snp_indices = eqtl_ld_scores[representative_gene_bin_pair_index,:] != 0.0 
		Y = np.square(eqtl_beta[representative_gene_bin_pair_index, snp_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, snp_indices])
		X = np.transpose(eqtl_ld_scores[:, snp_indices][gene_bin_pair_indices,:])

		model = sm.OLS(Y,X).fit()
		#cis_snp_h2 = model.params[0]
		psi_g[gene_bin_pair_indices] = model.params


	precomputed_xt_x = np.dot(eqtl_ld_scores, np.transpose(eqtl_ld_scores))

	gwas_resid_var = np.square(1/100000.0)

	'''
	gwas_resid_var = np.square(1/eqtl_sample_size)

	eqtl_resid_var = np.square(1.0/eqtl_sample_size)
	'''

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	###############################################################
	# First update sig_sq and alpha_sq (based two step approach)
	###############################################################
	# Prepare data
	Y = gwas_block_beta_sq - gwas_block_intercept
	X = []
	X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
	X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
	X = np.hstack(X)
	S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
	mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
	# Draw from distribution
	if gibbs:
		params = np.random.multivariate_normal(mean=mean, cov=S)
	else:
		params = np.copy(mean)
	sig_sq = params[0]
	alpha_sq = params[1:]



	###############################################################
	# Set residual variances
	###############################################################
	gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])

	for gene_iter in np.arange(n_genes):
		gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
		# Re-include current effect
		gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
	gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)

	#eqtl_beta_sq_resid = []
	eqtl_resid_vars = np.zeros(len(gene_index_to_gene_bin_pair_indices))
	for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
		representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
		if eqtl_only:
			variant_indices = np.sum(cis_eqtl_mask[gene_bin_pair_indices,:],axis=0) == 1.0
		else:
			variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
		tmp_eqtl_resid = (np.square(eqtl_beta[representative_gene_bin_pair_index, :]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, :])) - np.dot(np.transpose(eqtl_ld_scores[gene_bin_pair_indices,:]), psi_g[gene_bin_pair_indices])
		#eqtl_beta_sq_resid.append(tmp_eqtl_resid[variant_indices])
		eqtl_beta_sq_resid = tmp_eqtl_resid[variant_indices]
		eqtl_resid_vars[gene_iter] = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)

	#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
	#eqtl_beta_sq_resid = np.hstack(eqtl_beta_sq_resid)
	#eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)



	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		print(itera)
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]




		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = (precomputed_xt_x*np.dot(alpha_sq[gene_to_class_index].reshape(-1,1), alpha_sq[gene_to_class_index].reshape(1,-1)))/gwas_resid_var

		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])


		for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
			representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same

			variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
			
			eqtl_beta_sq_resid = np.square(eqtl_beta[representative_gene_bin_pair_index, variant_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, variant_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_bin_pair_indices, :][:, variant_indices]

			tmp_gene_mat = np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_vars[gene_iter]

			
			for ii, gene_iter_alt in enumerate(gene_bin_pair_indices):
				xt_x_term[gene_iter_alt, gene_bin_pair_indices] = xt_x_term[gene_iter_alt, gene_bin_pair_indices] + tmp_gene_mat[ii,:]
			
			#xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] = xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] + np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_var

			#xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid, eqtl_beta_sq_X)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_iter,:])/gwas_resid_var
			#xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:])/gwas_resid_var
			xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_vars[gene_iter]) + np.dot(alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:], gwas_beta_sq_resid)/gwas_resid_var

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)
		###############################################################
		# Update residual variances
		###############################################################
		'''
		if itera == 0:
			for gene_iter in np.arange(n_genes):
				gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
				# Re-include current effect
				gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
			gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)

			eqtl_beta_sq_resid = []
			for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
				representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
				variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
				tmp_eqtl_resid = (np.square(eqtl_beta[representative_gene_bin_pair_index, :]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, :])) - np.dot(np.transpose(eqtl_ld_scores[gene_bin_pair_indices,:]), psi_g[gene_bin_pair_indices])
				eqtl_beta_sq_resid.append(tmp_eqtl_resid[variant_indices])

			#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
			eqtl_beta_sq_resid = np.hstack(eqtl_beta_sq_resid)
			eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		'''

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 1) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(np.mean(eqtl_resid_vars)))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)


	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2, X




def med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_fixed_resid_vars(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, burn_in_iters=4600, total_iters=5000,gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	'''
	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1
	'''

	unique_gene_indices = np.sort(np.unique(gene_indexer))
	gene_index_to_gene_bin_pair_indices = []
	for unique_gene_index in unique_gene_indices:
		gene_index_to_gene_bin_pair_indices.append(np.where(gene_indexer==unique_gene_index)[0])


	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)

	for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
		representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
		snp_indices = eqtl_ld_scores[representative_gene_bin_pair_index,:] != 0.0 
		Y = np.square(eqtl_beta[representative_gene_bin_pair_index, snp_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, snp_indices])
		X = np.transpose(eqtl_ld_scores[:, snp_indices][gene_bin_pair_indices,:])

		model = sm.OLS(Y,X).fit()
		#cis_snp_h2 = model.params[0]
		psi_g[gene_bin_pair_indices] = model.params


	precomputed_xt_x = np.dot(eqtl_ld_scores, np.transpose(eqtl_ld_scores))

	gwas_resid_var = np.square(1/100000.0)
	eqtl_resid_var = np.square(1.0/eqtl_sample_size)

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		print(itera)
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]




		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = (precomputed_xt_x*np.dot(alpha_sq[gene_to_class_index].reshape(-1,1), alpha_sq[gene_to_class_index].reshape(1,-1)))/gwas_resid_var

		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])


		for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
			representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same

			variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
			
			eqtl_beta_sq_resid = np.square(eqtl_beta[representative_gene_bin_pair_index, variant_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, variant_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_bin_pair_indices, :][:, variant_indices]

			tmp_gene_mat = np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_var

			
			for ii, gene_iter in enumerate(gene_bin_pair_indices):
				xt_x_term[gene_iter, gene_bin_pair_indices] = xt_x_term[gene_iter, gene_bin_pair_indices] + tmp_gene_mat[ii,:]
			
			#xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] = xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] + np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_var

			#xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid, eqtl_beta_sq_X)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_iter,:])/gwas_resid_var
			#xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:])/gwas_resid_var
			xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_var) + np.dot(alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:], gwas_beta_sq_resid)/gwas_resid_var

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)
		###############################################################
		# Update residual variances
		###############################################################
		if itera == 0:
			for gene_iter in np.arange(n_genes):
				gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
				# Re-include current effect
				gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
			gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)

			eqtl_beta_sq_resid = []
			for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
				representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
				variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
				tmp_eqtl_resid = (np.square(eqtl_beta[representative_gene_bin_pair_index, :]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, :])) - np.dot(np.transpose(eqtl_ld_scores[gene_bin_pair_indices,:]), psi_g[gene_bin_pair_indices])
				eqtl_beta_sq_resid.append(tmp_eqtl_resid[variant_indices])

			#eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
			eqtl_beta_sq_resid = np.hstack(eqtl_beta_sq_resid)
			eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 1) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(eqtl_resid_var))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)


	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2, X


def med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_gene_resid_var(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, burn_in_iters=4600, total_iters=5000,gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	'''
	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1
	'''

	unique_gene_indices = np.sort(np.unique(gene_indexer))
	gene_index_to_gene_bin_pair_indices = []
	for unique_gene_index in unique_gene_indices:
		gene_index_to_gene_bin_pair_indices.append(np.where(gene_indexer==unique_gene_index)[0])


	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)

	for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
		representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
		snp_indices = eqtl_ld_scores[representative_gene_bin_pair_index,:] != 0.0 
		Y = np.square(eqtl_beta[representative_gene_bin_pair_index, snp_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, snp_indices])
		X = np.transpose(eqtl_ld_scores[:, snp_indices][gene_bin_pair_indices,:])

		model = sm.OLS(Y,X).fit()
		#cis_snp_h2 = model.params[0]
		psi_g[gene_bin_pair_indices] = model.params

	#psi_g = np.random.normal(loc=np.mean(psi_g), scale=np.std(psi_g), size=len(psi_g))


	precomputed_xt_x = np.dot(eqtl_ld_scores, np.transpose(eqtl_ld_scores))

	gwas_resid_var = np.square(1/100000.0)
	#eqtl_resid_var = np.square(1.0/eqtl_sample_size)
	eqtl_resid_vars = np.square(1.0/eqtl_sample_size)*np.ones(len(gene_index_to_gene_bin_pair_indices))


	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		print(itera)
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]




		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = (precomputed_xt_x*np.dot(alpha_sq[gene_to_class_index].reshape(-1,1), alpha_sq[gene_to_class_index].reshape(1,-1)))/gwas_resid_var

		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])


		for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
			representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same

			variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
			
			eqtl_beta_sq_resid = np.square(eqtl_beta[representative_gene_bin_pair_index, variant_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, variant_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_bin_pair_indices, :][:, variant_indices]

			tmp_gene_mat = np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_vars[gene_iter]

			
			for ii, local_gene_iter in enumerate(gene_bin_pair_indices):
				xt_x_term[local_gene_iter, gene_bin_pair_indices] = xt_x_term[local_gene_iter, gene_bin_pair_indices] + tmp_gene_mat[ii,:]
			
			#xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] = xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] + np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_var

			#xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid, eqtl_beta_sq_X)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_iter,:])/gwas_resid_var
			#xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:])/gwas_resid_var
			xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_vars[gene_iter]) + np.dot(alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:], gwas_beta_sq_resid)/gwas_resid_var

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)
		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)

		#eqtl_beta_sq_resid = []
		for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
			representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
			variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
			tmp_eqtl_resid = (np.square(eqtl_beta[representative_gene_bin_pair_index, :]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, :])) - np.dot(np.transpose(eqtl_ld_scores[gene_bin_pair_indices,:]), psi_g[gene_bin_pair_indices])
			#eqtl_beta_sq_resid.append(tmp_eqtl_resid[variant_indices])
			eqtl_resid_vars[gene_iter] = np.sum(np.square(tmp_eqtl_resid[variant_indices]))/len(tmp_eqtl_resid[variant_indices])
		#eqtl_beta_sq_resid = np.hstack(eqtl_beta_sq_resid)
		#eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		#eqtl_resid_vars = eqtl_resid_vars*0.0 + eqtl_resid_var

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 1) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(np.mean(eqtl_resid_vars)))
			#print(np.sort(eqtl_resid_vars))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)


	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2, X


def med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_gene_resid_var_eqtl_snp_only(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=4600, total_iters=5000,gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1

	unique_gene_indices = np.sort(np.unique(gene_indexer))
	gene_index_to_gene_bin_pair_indices = []
	for unique_gene_index in unique_gene_indices:
		gene_index_to_gene_bin_pair_indices.append(np.where(gene_indexer==unique_gene_index)[0])


	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)

	for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
		representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
		snp_indices = eqtl_ld_scores[representative_gene_bin_pair_index,:] != 0.0 
		Y = np.square(eqtl_beta[representative_gene_bin_pair_index, snp_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, snp_indices])
		X = np.transpose(eqtl_ld_scores[:, snp_indices][gene_bin_pair_indices,:])

		model = sm.OLS(Y,X).fit()
		#cis_snp_h2 = model.params[0]
		psi_g[gene_bin_pair_indices] = model.params

	#psi_g = np.random.normal(loc=np.mean(psi_g), scale=np.std(psi_g), size=len(psi_g))


	precomputed_xt_x = np.dot(eqtl_ld_scores, np.transpose(eqtl_ld_scores))

	gwas_resid_var = np.square(1/100000.0)
	#eqtl_resid_var = np.square(1.0/eqtl_sample_size)
	eqtl_resid_vars = np.square(1.0/eqtl_sample_size)*np.ones(len(gene_index_to_gene_bin_pair_indices))


	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		print(itera)
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		if gibbs:
			params = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			params = np.copy(mean)
		sig_sq = params[0]
		alpha_sq = params[1:]




		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = (precomputed_xt_x*np.dot(alpha_sq[gene_to_class_index].reshape(-1,1), alpha_sq[gene_to_class_index].reshape(1,-1)))/gwas_resid_var

		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])


		for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
			representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same

			variant_indices = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
			
			eqtl_beta_sq_resid = np.square(eqtl_beta[representative_gene_bin_pair_index, variant_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, variant_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_bin_pair_indices, :][:, variant_indices]

			tmp_gene_mat = np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_vars[gene_iter]

			
			for ii, local_gene_iter in enumerate(gene_bin_pair_indices):
				xt_x_term[local_gene_iter, gene_bin_pair_indices] = xt_x_term[local_gene_iter, gene_bin_pair_indices] + tmp_gene_mat[ii,:]
			
			#xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] = xt_x_term[gene_bin_pair_indices, :][:, gene_bin_pair_indices] + np.dot(eqtl_beta_sq_X, np.transpose(eqtl_beta_sq_X))/eqtl_resid_var

			#xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid, eqtl_beta_sq_X)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_iter,:])/gwas_resid_var
			#xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:])/gwas_resid_var
			xt_y_term[gene_bin_pair_indices] = (np.dot(eqtl_beta_sq_X, eqtl_beta_sq_resid)/eqtl_resid_vars[gene_iter]) + np.dot(alpha_weighted_eqtl_ld_scores[gene_bin_pair_indices,:], gwas_beta_sq_resid)/gwas_resid_var

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		if gibbs:
			psi_g = np.random.multivariate_normal(mean=mean, cov=S)
		else:
			psi_g = np.copy(mean)
		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)

		#eqtl_beta_sq_resid = []
		for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
			representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
			#variant_indices_old = eqtl_beta[representative_gene_bin_pair_index,:] != 0.0
			variant_indices = np.sum(cis_eqtl_mask[gene_bin_pair_indices,:],axis=0) == 1.0
			tmp_eqtl_resid = (np.square(eqtl_beta[representative_gene_bin_pair_index, :]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, :])) - np.dot(np.transpose(eqtl_ld_scores[gene_bin_pair_indices,:]), psi_g[gene_bin_pair_indices])
			#eqtl_beta_sq_resid.append(tmp_eqtl_resid[variant_indices])
			eqtl_resid_vars[gene_iter] = np.sum(np.square(tmp_eqtl_resid[variant_indices]))/len(tmp_eqtl_resid[variant_indices])
		#eqtl_beta_sq_resid = np.hstack(eqtl_beta_sq_resid)
		#eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		#eqtl_resid_vars = eqtl_resid_vars*0.0 + eqtl_resid_var

		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 1) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(np.mean(eqtl_resid_vars)))
			#print(np.sort(eqtl_resid_vars))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)


	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2, X





def med_h2_with_sumstat_ldsc_two_step_multivariate_gene_bins(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, nm_var_snp_size_vec, burn_in_iters=4600, total_iters=5000,gibbs=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	'''
	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1
	'''

	unique_gene_indices = np.sort(np.unique(gene_indexer))
	gene_index_to_gene_bin_pair_indices = []
	for unique_gene_index in unique_gene_indices:
		gene_index_to_gene_bin_pair_indices.append(np.where(gene_indexer==unique_gene_index)[0])


	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)

	for gene_iter, gene_bin_pair_indices in enumerate(gene_index_to_gene_bin_pair_indices):
		representative_gene_bin_pair_index = gene_bin_pair_indices[0] # Doesn't have to be zero index all the same
		snp_indices = eqtl_ld_scores[representative_gene_bin_pair_index,:] != 0.0 
		Y = np.square(eqtl_beta[representative_gene_bin_pair_index, snp_indices]) - np.square(eqtl_beta_se[representative_gene_bin_pair_index, snp_indices])
		X = np.transpose(eqtl_ld_scores[:, snp_indices][gene_bin_pair_indices,:])

		model = sm.OLS(Y,X).fit()
		#cis_snp_h2 = model.params[0]
		psi_g[gene_bin_pair_indices] = model.params


	precomputed_xt_x = np.dot(eqtl_ld_scores, np.transpose(eqtl_ld_scores))

	gwas_resid_var = np.square(1/100000.0)
	eqtl_resid_var = np.square(1.0/eqtl_sample_size)

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []


	###############################################################
	# First update sig_sq and alpha_sq (based on only gwas block)
	###############################################################
	# Prepare data
	Y = gwas_block_beta_sq - gwas_block_intercept
	X = []
	X.append(gwas_ld_scores)
	X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
	X = np.hstack(X)
	S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
	mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
	# Draw from distribution
	if gibbs:
		params = np.random.multivariate_normal(mean=mean, cov=S)
	else:
		params = np.copy(mean)

	n_nm_vars = gwas_ld_scores.shape[1]
	sig_sq = params[0:n_nm_vars]
	alpha_sq = params[n_nm_vars:]



	gene_cis_h2 = psi_g*eqtl_n_cis_snps
	per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
	per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
	per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
	nm_h2 = np.dot(sig_sq, nm_var_snp_size_vec)
	total_med_h2 = np.sum(per_gene_med_h2)



	return total_med_h2, per_class_med_h2, nm_h2, gene_cis_h2, X



def med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_only_cis_eqtl_snp_pos(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, eqtl_cis_snp_position, burn_in_iters=4600, total_iters=5000):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1

	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2



	gwas_resid_var = 1e-4
	eqtl_resid_var = 1e-4

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		params = np.random.multivariate_normal(mean=mean, cov=S)
		sig_sq = params[0]
		alpha_sq = params[1:]


		###############################################################
		# Second update psi_g (multivariate update)
		###############################################################
		alpha_weighted_eqtl_ld_scores = np.transpose(np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index])
		# Create base X^TX (a matrix that will be used multiple times throughout analysis)
		xt_x_term = np.dot(alpha_weighted_eqtl_ld_scores, np.transpose(alpha_weighted_eqtl_ld_scores))/gwas_resid_var
		xt_y_term = np.zeros(n_genes)

		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - (sig_sq*X[:,0])

		for gene_iter in np.arange(n_genes):

			#gene_indices = eqtl_beta[gene_iter,:] != 0.0
			gene_indices = cis_eqtl_mask[gene_iter,:] == 1

			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]

			xt_x_term[gene_iter, gene_iter] = xt_x_term[gene_iter, gene_iter] + np.dot(eqtl_beta_sq_X, eqtl_beta_sq_X)/eqtl_resid_var

			xt_y_term[gene_iter] = (np.dot(eqtl_beta_sq_resid, eqtl_beta_sq_X)/eqtl_resid_var) + np.dot(gwas_beta_sq_resid, alpha_weighted_eqtl_ld_scores[gene_iter,:])/gwas_resid_var

		S = np.linalg.inv(xt_x_term)
		mean = np.dot(S, xt_y_term)

		psi_g = np.random.multivariate_normal(mean=mean, cov=S)

		###############################################################
		# Update residual variances
		###############################################################
		for gene_iter in np.arange(n_genes):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)


		eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[cis_eqtl_mask==1]
		eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)


		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 20) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(eqtl_resid_var))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2



def med_h2_with_sumstat_ldsc_bayesian_gibbs(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, burn_in_iters=4600, total_iters=5000):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	'''
	cis_eqtl_mask = np.copy(eqtl_beta)*0.0
	for gene_iter in range(n_genes):
		cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] = cis_eqtl_mask[gene_iter,eqtl_cis_snp_position[gene_iter]] + 1
	'''
	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2

	gwas_resid_var = 1e-4
	eqtl_resid_var = 1e-4

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)

	sampled_total_med = []
	sampled_nm = []
	sampled_per_tissue_med = []
	sampled_eqtl_h2 = []

	# Now loop through iterations
	for itera in range(total_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		S = np.linalg.inv(np.dot(np.transpose(X), X)/gwas_resid_var)
		mean = (1.0/gwas_resid_var)*np.dot(S, np.dot(Y,X))
		# Draw from distribution
		params = np.random.multivariate_normal(mean=mean, cov=S)
		sig_sq = params[0]
		alpha_sq = params[1:]


		###############################################################
		# Second update psi_g
		###############################################################
		gwas_beta_sq_resid = gwas_block_beta_sq - gwas_block_intercept - np.dot(X, params)
		for gene_iter in np.random.permutation(np.arange(n_genes)):
			gene_alpha_sq = alpha_sq[gene_to_class_index[gene_iter]]
			# Re-include current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid + ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)
			gwas_beta_sq_X = eqtl_ld_scores[gene_iter,:]*gene_alpha_sq
			gwas_beta_sq_resid_var = np.ones(len(gwas_beta_sq_X))*gwas_resid_var

			gene_indices = eqtl_beta[gene_iter,:] != 0.0
			eqtl_beta_sq_resid = np.square(eqtl_beta[gene_iter, gene_indices]) - np.square(eqtl_beta_se[gene_iter, gene_indices])
			eqtl_beta_sq_X = eqtl_ld_scores[gene_iter, gene_indices]
			eqtl_beta_sq_resid_var = np.ones(len(eqtl_beta_sq_X))*eqtl_resid_var

			gene_Y = np.hstack((gwas_beta_sq_resid, eqtl_beta_sq_resid))
			gene_X = np.hstack((gwas_beta_sq_X, eqtl_beta_sq_X))
			gene_resid_var = np.hstack((gwas_beta_sq_resid_var, eqtl_beta_sq_resid_var))

			# Get sampling distribution
			S = 1.0/np.sum(np.square(gene_X)/gene_resid_var)
			mean = S*np.sum(gene_X*gene_Y/gene_resid_var)

			# Sample
			psi_g[gene_iter] = np.random.normal(loc=mean, scale=np.sqrt(S))

			# Remove current effect
			gwas_beta_sq_resid = gwas_beta_sq_resid - ((eqtl_ld_scores[gene_iter,:])*(psi_g[gene_iter])*gene_alpha_sq)


		###############################################################
		# Update residual variances
		###############################################################
		eqtl_beta_sq_resid = (np.square(eqtl_beta) - np.square(eqtl_beta_se) - np.transpose(np.transpose(eqtl_ld_scores)*psi_g))[eqtl_mask==1]
		eqtl_resid_var = np.sum(np.square(eqtl_beta_sq_resid))/len(eqtl_beta_sq_resid)
		gwas_resid_var = np.sum(np.square(gwas_beta_sq_resid))/len(gwas_beta_sq_resid)


		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)


		if np.mod(itera, 20) ==0:
			print('Iteration: ' + str(itera))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('GWAS variance: ' + str(gwas_resid_var))
			print('eQTL variance: ' + str(eqtl_resid_var))
			print('')

		if itera > burn_in_iters:
			sampled_total_med.append(total_med_h2)
			sampled_nm.append(nm_h2)
			sampled_per_tissue_med.append(per_class_med_h2)
			sampled_eqtl_h2.append(gene_cis_h2)

	avg_total_med_h2 = np.mean(sampled_total_med)
	avg_per_class_med_h2 = np.mean(np.asarray(sampled_per_tissue_med),axis=0)
	avg_nm_h2 = np.mean(sampled_nm)
	avg_gene_cis_h2 = np.mean(sampled_eqtl_h2, axis=0)


	return avg_total_med_h2, avg_per_class_med_h2, avg_nm_h2, avg_gene_cis_h2


def med_h2_with_sumstat_ldsc_two_step(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, max_iters=1):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	# Initialize per gene snp h2
	psi_g = np.zeros(n_genes)
	for gene_iter in range(n_genes):
		indices = eqtl_ld_scores[gene_iter,:] != 0.0
		Y = np.square(eqtl_beta[gene_iter, indices]) - np.square(eqtl_beta_se[gene_iter, indices])
		X = eqtl_ld_scores[gene_iter, indices]
		model = sm.OLS(Y,X).fit()
		cis_snp_h2 = model.params[0]
		psi_g[gene_iter] = cis_snp_h2

	# Initialize non-mediated per snp h2 (doesn't matter what they are. all that matters is psi_g)
	sig_sq = 0.0
	alpha_sq = np.zeros(n_eqtl_classes)

	# Organize data
	gwas_block_beta_sq = np.square(gwas_beta)
	gwas_block_intercept = np.square(gwas_beta_se)
	eqtl_block_beta_sq = []
	eqtl_block_intercept = []
	eqtl_block_eqtl_ld_scores = []
	for gene_iter in range(n_genes):
		gene_indices = eqtl_beta[gene_iter,:] != 0.0
		eqtl_block_beta_sq.append(np.square(eqtl_beta[gene_iter, gene_indices]))
		eqtl_block_intercept.append(np.square(eqtl_beta_se[gene_iter, gene_indices]))
		eqtl_block_ld_score = np.zeros((np.sum(gene_indices), n_genes))
		eqtl_block_ld_score[:, gene_iter] = eqtl_ld_scores[gene_iter, gene_indices]
		eqtl_block_eqtl_ld_scores.append(eqtl_block_ld_score)
	eqtl_block_beta_sq = np.hstack(eqtl_block_beta_sq)
	eqtl_block_intercept = np.hstack(eqtl_block_intercept)
	eqtl_block_eqtl_ld_scores = np.vstack(eqtl_block_eqtl_ld_scores)

	# Now loop through iterations
	for itera in range(max_iters):
		###############################################################
		# First update sig_sq and alpha_sq (based on only gwas block)
		###############################################################
		# Prepare data
		Y = gwas_block_beta_sq - gwas_block_intercept
		X = []
		X.append(np.transpose(gwas_ld_scores.reshape(1,-1)))
		X.append(np.dot(np.transpose(eqtl_ld_scores)*psi_g, gene_to_class_index_matrix))
		X = np.hstack(X)
		# Run OLS
		model = sm.OLS(Y,X).fit()
		sig_sq = model.params[0]
		alpha_sq = model.params[1:]

		###############################################################
		# Second update psi_g
		###############################################################
		'''
		# Prepare data
		Y_top = gwas_block_beta_sq - gwas_block_intercept - (gwas_ld_scores*sig_sq)
		Y_bottom = eqtl_block_beta_sq - eqtl_block_intercept
		X_top = np.transpose(eqtl_ld_scores)*alpha_sq[gene_to_class_index]
		# Note that X_bottom is simply eqtl_block_eqtl_ld_scores
		Y = np.hstack((Y_top, Y_bottom))
		X = np.vstack((X_top, eqtl_block_eqtl_ld_scores))
		# Run OLS
		model = sm.OLS(Y,X).fit()
		psi_g = model.params
		'''
		gene_cis_h2 = psi_g*eqtl_n_cis_snps
		per_gene_alpha_sq_variable = np.asarray(alpha_sq)[gene_to_class_index]
		per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
		per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
		nm_h2 = sig_sq*n_snps
		total_med_h2 = np.sum(per_gene_med_h2)

		print('med: ' + str(total_med_h2))
		print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
		print('nm: ' + str(nm_h2))
		print('eqtl: ' + str(np.mean(gene_cis_h2)))


	return total_med_h2, per_class_med_h2, nm_h2, gene_cis_h2


def med_h2_with_sumstat_ldsc(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, max_epochs=40000, conv_thresh=1e-12, intercept_variables=False, learning_rate=5e-6):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)


	optimizer = tf.keras.optimizers.Adam(learning_rate=1e-5)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)


	# Convert variabels to tf tensors
	gwas_beta_sq = tf.convert_to_tensor(np.square(gwas_beta), dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	gwas_ld_scores = tf.convert_to_tensor(gwas_ld_scores, dtype=tf.float32)
	eqtl_beta_sq = tf.convert_to_tensor(np.square(eqtl_beta), dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)
	eqtl_ld_scores = tf.convert_to_tensor(eqtl_ld_scores, dtype=tf.float32)
	gene_to_class_index = tf.convert_to_tensor(gene_to_class_index.astype(int))
	gene_to_class_index_matrix = tf.convert_to_tensor(gene_to_class_index_matrix, dtype=tf.float32)




	# Initialize variables to optimize over
	beta_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='beta_sq', dtype=tf.float32)
	alpha_sq_variable = tf.Variable(initial_value=np.ones(n_eqtl_classes)*0.0000001,trainable=True, name='alpha_sq',dtype=tf.float32)
	eqtl_beta_sq_variable = tf.Variable(initial_value=np.ones(n_genes)*0.000000001,trainable=True, name='eqtl_beta_sq', dtype=tf.float32)

	gwas_noise = tf.Variable(initial_value=1e-2,trainable=True, name='gwas_noise')
	eqtl_noise = tf.Variable(initial_value=1e-2,trainable=True, name='eqtl_noise')

	converged = False
	prev_est_alpha_sq=10000
	best_loss = 1e10
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = sumstat_proportional_joint_reml_loss(gwas_beta_sq, gwas_beta_var,gwas_ld_scores, eqtl_beta_sq, eqtl_beta_var,eqtl_ld_scores,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, gene_to_class_index_matrix, gwas_noise, eqtl_noise)

		trainable_variables = []
		trainable_variables.append(alpha_sq_variable)
		trainable_variables.append(eqtl_beta_sq_variable)
		trainable_variables.append(beta_sq_variable)
		#trainable_variables.append(gwas_noise)
		#trainable_variables.append(eqtl_noise)
		#if intercept_variables:
			#trainable_variables.append(gwas_noise)
			#trainable_variables.append(eqtl_noise)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est = np.asarray(alpha_sq_variable)*1.0

		diff = np.abs(prev_est_alpha_sq -cur_est)
		'''
		if diff < conv_thresh:
			converged = True
			break
		'''

		prev_est_alpha_sq = cur_est


		cur_loss = np.asmatrix(loss_value)[0,0]
		if cur_loss <= best_loss:
			best_loss = cur_loss
			eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable
			best_gene_cis_h2 = np.asarray(eqtl_beta_sq_variable_scaled)*eqtl_n_cis_snps
			per_gene_alpha_sq_variable = np.asarray(alpha_sq_variable)[gene_to_class_index]
			per_gene_med_h2 = per_gene_alpha_sq_variable*best_gene_cis_h2
			best_per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
			best_total_med_h2 = np.sum(per_gene_med_h2)
			best_nm_h2 = np.sum(beta_sq_variable)*n_snps


		if np.mod(epoch_iter, 100) == 0.0:
			eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable 
			gene_cis_h2 = np.asarray(eqtl_beta_sq_variable_scaled)*eqtl_n_cis_snps
			per_gene_alpha_sq_variable = np.asarray(alpha_sq_variable)[gene_to_class_index]
			per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
			per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
			total_med_h2 = np.sum(per_gene_med_h2)
			nm_h2 = np.sum(beta_sq_variable)*n_snps
			print(epoch_iter)
			print('loss: ' + str(loss_value))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('gwas_noise: ' + str(np.asmatrix(gwas_noise)[0,0]))
			print('eqtl noise: ' + str(np.asmatrix(eqtl_noise)[0,0]))


	print(best_loss)


	return best_total_med_h2, best_per_class_med_h2, best_nm_h2, best_gene_cis_h2



def med_h2_with_sumstat_ldsc_learn_resid_vars(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, max_epochs=40000, conv_thresh=1e-12, intercept_variables=False, learning_rate=5e-6):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)


	optimizer = tf.keras.optimizers.Adam(learning_rate=1e-3)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)


	# Convert variabels to tf tensors
	gwas_beta_sq = tf.convert_to_tensor(np.square(gwas_beta), dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	gwas_ld_scores = tf.convert_to_tensor(gwas_ld_scores, dtype=tf.float32)
	eqtl_beta_sq = tf.convert_to_tensor(np.square(eqtl_beta), dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)
	eqtl_ld_scores = tf.convert_to_tensor(eqtl_ld_scores, dtype=tf.float32)
	gene_to_class_index = tf.convert_to_tensor(gene_to_class_index.astype(int))
	gene_to_class_index_matrix = tf.convert_to_tensor(gene_to_class_index_matrix, dtype=tf.float32)




	# Initialize variables to optimize over
	beta_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='beta_sq', dtype=tf.float32)
	alpha_sq_variable = tf.Variable(initial_value=np.ones(n_eqtl_classes)*0.0000001,trainable=True, name='alpha_sq',dtype=tf.float32)
	eqtl_beta_sq_variable = tf.Variable(initial_value=np.ones(n_genes)*0.000000001,trainable=True, name='eqtl_beta_sq', dtype=tf.float32)

	gwas_noise = tf.Variable(initial_value=-12.0,trainable=True, name='gwas_noise')
	eqtl_noise = tf.Variable(initial_value=-12.0,trainable=True, name='eqtl_noise')

	converged = False
	prev_est_alpha_sq=10000
	best_loss = 1e10
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = sumstat_proportional_joint_reml_loss_softplus_link(gwas_beta_sq, gwas_beta_var,gwas_ld_scores, eqtl_beta_sq, eqtl_beta_var,eqtl_ld_scores,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, gene_to_class_index_matrix, gwas_noise, eqtl_noise)

		trainable_variables = []
		trainable_variables.append(alpha_sq_variable)
		trainable_variables.append(eqtl_beta_sq_variable)
		trainable_variables.append(beta_sq_variable)
		trainable_variables.append(gwas_noise)
		trainable_variables.append(eqtl_noise)
		#if intercept_variables:
			#trainable_variables.append(gwas_noise)
			#trainable_variables.append(eqtl_noise)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est = np.asarray(alpha_sq_variable)*1.0

		diff = np.abs(prev_est_alpha_sq -cur_est)
		'''
		if diff < conv_thresh:
			converged = True
			break
		'''

		prev_est_alpha_sq = cur_est


		cur_loss = np.asmatrix(loss_value)[0,0]
		if cur_loss < best_loss:
			best_loss = cur_loss
			eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable
			best_gene_cis_h2 = np.asarray(eqtl_beta_sq_variable_scaled)*eqtl_n_cis_snps
			per_gene_alpha_sq_variable = np.asarray(alpha_sq_variable)[gene_to_class_index]
			per_gene_med_h2 = per_gene_alpha_sq_variable*best_gene_cis_h2
			best_per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
			best_total_med_h2 = np.sum(per_gene_med_h2)
			best_nm_h2 = np.sum(beta_sq_variable)*n_snps


		if np.mod(epoch_iter, 100) == 0.0:
			eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable 
			gene_cis_h2 = np.asarray(eqtl_beta_sq_variable_scaled)*eqtl_n_cis_snps
			per_gene_alpha_sq_variable = np.asarray(alpha_sq_variable)[gene_to_class_index]
			per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
			per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
			total_med_h2 = np.sum(per_gene_med_h2)
			nm_h2 = np.sum(beta_sq_variable)*n_snps
			print(epoch_iter)
			print('loss: ' + str(loss_value))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('gwas_noise: ' + str(np.asmatrix(tf.math.softplus(gwas_noise))[0,0]))
			print('eqtl noise: ' + str(np.asmatrix(tf.math.softplus(eqtl_noise))[0,0]))


	print(best_loss)


	return best_total_med_h2, best_per_class_med_h2, best_nm_h2, best_gene_cis_h2



def med_h2_with_sumstat_ldsc_w_per_study_residual_variances(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, eqtl_classes, max_epochs=40000, conv_thresh=1e-12, intercept_variables=False, learning_rate=5e-6):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Ordered eqtl classes
	ordered_eqtl_classes = np.sort(np.unique(eqtl_classes))
	n_eqtl_classes = len(ordered_eqtl_classes)
	class_mapping = {}
	for ii, eqtl_class in enumerate(ordered_eqtl_classes):
		class_mapping[eqtl_class] = ii
	gene_to_class_index = []
	gene_to_class_index_matrix = np.zeros((n_genes, n_eqtl_classes))
	for gene_iter,eqtl_class in enumerate(eqtl_classes):
		gene_to_class_index.append(class_mapping[eqtl_class])
		gene_to_class_index_matrix[gene_iter, class_mapping[eqtl_class]] = 1
	gene_to_class_index = np.asarray(gene_to_class_index)


	optimizer = tf.keras.optimizers.Adam(learning_rate=1e-5)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)


	# Convert variabels to tf tensors
	gwas_beta_sq = tf.convert_to_tensor(np.square(gwas_beta), dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	gwas_ld_scores = tf.convert_to_tensor(gwas_ld_scores, dtype=tf.float32)
	eqtl_beta_sq = tf.convert_to_tensor(np.square(eqtl_beta), dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)
	eqtl_ld_scores = tf.convert_to_tensor(eqtl_ld_scores, dtype=tf.float32)
	gene_to_class_index = tf.convert_to_tensor(gene_to_class_index.astype(int))
	gene_to_class_index_matrix = tf.convert_to_tensor(gene_to_class_index_matrix, dtype=tf.float32)




	# Initialize variables to optimize over
	beta_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='beta_sq', dtype=tf.float32)
	alpha_sq_variable = tf.Variable(initial_value=np.ones(n_eqtl_classes)*0.0000001,trainable=True, name='alpha_sq',dtype=tf.float32)
	eqtl_beta_sq_variable = tf.Variable(initial_value=np.ones(n_genes)*0.000000001,trainable=True, name='eqtl_beta_sq', dtype=tf.float32)

	#gwas_noise = tf.Variable(initial_value=1.0,trainable=True, name='gwas_noise')
	#eqtl_noise = tf.Variable(initial_value=1.0,trainable=True, name='eqtl_noise')
	gwas_noise = tf.convert_to_tensor(np.square(np.mean(np.square(gwas_beta))), dtype=tf.float32)
	eqtl_noise = tf.convert_to_tensor(np.square(np.mean(np.square(eqtl_beta[eqtl_beta != 0.0]))), dtype=tf.float32)

	converged = False
	prev_est_alpha_sq=10000
	best_loss = 1e10
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = sumstat_proportional_joint_reml_loss(gwas_beta_sq, gwas_beta_var,gwas_ld_scores, eqtl_beta_sq, eqtl_beta_var,eqtl_ld_scores,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, gene_to_class_index_matrix, gwas_noise, eqtl_noise)

		trainable_variables = []
		trainable_variables.append(alpha_sq_variable)
		trainable_variables.append(eqtl_beta_sq_variable)
		trainable_variables.append(beta_sq_variable)


		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est = np.asarray(alpha_sq_variable)*1.0

		diff = np.abs(prev_est_alpha_sq -cur_est)
		'''
		if diff < conv_thresh:
			converged = True
			break
		'''

		prev_est_alpha_sq = cur_est


		cur_loss = np.asmatrix(loss_value)[0,0]
		if cur_loss < best_loss:
			best_loss = cur_loss
			eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable
			best_gene_cis_h2 = np.asarray(eqtl_beta_sq_variable_scaled)*eqtl_n_cis_snps
			per_gene_alpha_sq_variable = np.asarray(alpha_sq_variable)[gene_to_class_index]
			per_gene_med_h2 = per_gene_alpha_sq_variable*best_gene_cis_h2
			best_per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
			best_total_med_h2 = np.sum(per_gene_med_h2)
			best_nm_h2 = np.sum(beta_sq_variable)*n_snps


		if np.mod(epoch_iter, 100) == 0.0:
			eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable 
			gene_cis_h2 = np.asarray(eqtl_beta_sq_variable_scaled)*eqtl_n_cis_snps
			per_gene_alpha_sq_variable = np.asarray(alpha_sq_variable)[gene_to_class_index]
			per_gene_med_h2 = per_gene_alpha_sq_variable*gene_cis_h2
			per_class_med_h2 = np.dot(per_gene_med_h2, np.asarray(gene_to_class_index_matrix))
			total_med_h2 = np.sum(per_gene_med_h2)
			nm_h2 = np.sum(beta_sq_variable)*n_snps
			print(epoch_iter)
			print('loss: ' + str(loss_value))
			print('med: ' + str(total_med_h2))
			print('per-tissue med: ' + ','.join(per_class_med_h2.astype(str)))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('gwas_noise: ' + str(np.asmatrix(gwas_noise)[0,0]))
			print('eqtl noise: ' + str(np.asmatrix(eqtl_noise)[0,0]))


	print(best_loss)


	return best_total_med_h2, best_per_class_med_h2, best_nm_h2, best_gene_cis_h2


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

# Load in gene window variant ld scores file
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
single_bin_genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, eqtl_cis_snp_position_single_bin, eqtl_n_cis_snps, eqtl_classes, gene_indexer = load_in_eqtl_data_with_bins_per_gene(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, eqtl_ld='out_of_sample', n_bins=1)
gwas_gene_window_variant_ld_scores, gwas_gene_window_variant_anno_vec = get_gwas_gene_window_variant_ld_scores(gwas_beta, gwas_rsids, quasi_ld_window_summary_file, single_bin_genes, eqtl_cis_snp_position_single_bin)

# load in eqtl data
# Out of sample eqtl ld
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, eqtl_cis_snp_position, eqtl_n_cis_snps, eqtl_classes, gene_indexer = load_in_eqtl_data_with_bins_per_gene(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, eqtl_ld='out_of_sample', n_bins=5)
# Generate matrix form of eqtl data
eqtl_beta_mat = get_matrix_form_of_eqtl_data(eqtl_beta, eqtl_position, 0.0, len(gwas_variant_ld_scores))
eqtl_beta_se_mat = get_matrix_form_of_eqtl_data(eqtl_beta_se, eqtl_position, 1.0/np.sqrt(eqtl_sample_size), len(gwas_variant_ld_scores))
eqtl_ld_score_mat = get_matrix_form_of_eqtl_data(eqtl_ldscore, eqtl_position, 0.0, len(gwas_variant_ld_scores))



'''
##########################
# Temp saving
np.save('gwas_beta.npy', gwas_beta)
np.save('gwas_beta_se.npy', gwas_beta_se)
np.save('gwas_variant_ld_scores.npy', gwas_variant_ld_scores)
np.save('gwas_gene_window_variant_ld_scores.npy', gwas_gene_window_variant_ld_scores)
np.save('eqtl_beta_mat.npy', eqtl_beta_mat)
np.save('eqtl_beta_se_mat.npy', eqtl_beta_se_mat)
np.save('eqtl_ld_score_mat.npy', eqtl_ld_score_mat)
np.save('eqtl_n_cis_snps.npy', eqtl_n_cis_snps)
np.save('eqtl_classes.npy', eqtl_classes)
np.save('gene_indexer.npy', gene_indexer)
f = open('eqtl_cis_snp_position.pkl', 'wb')
pickle.dump(eqtl_cis_snp_position, f, pickle.HIGHEST_PROTOCOL)
f.close()
'''
'''
##########################
# temp loading
gwas_beta = np.load('gwas_beta.npy')
gwas_beta_se = np.load('gwas_beta_se.npy')
gwas_variant_ld_scores = np.load('gwas_variant_ld_scores.npy')
gwas_gene_window_variant_ld_scores = np.load('gwas_gene_window_variant_ld_scores.npy')
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


# Run standard S-LDSC
X = sm.add_constant(gwas_variant_ld_scores)
Y = np.square(gwas_beta/gwas_beta_se)
model = sm.OLS(Y,X).fit()
model_constrained_intercept = sm.OLS(Y-1, gwas_variant_ld_scores).fit()
ldsc_snp_h2 = len(gwas_beta)*(model.params[1]/N_gwas)
ldsc_constrained_intercept_snp_h2 = len(gwas_beta)*(model_constrained_intercept.params[0]/N_gwas)


output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod13.txt'
t = open(output_file,'w')
t.write('method\teQTL_SS\tsim_h2\tsim_med_h2\tsim_nm_h2\test_med_h2_joint_reml\test_med_h2_per_tissue_joint_reml\test_nm_h2_joint_reml\test_mean_eqtl_h2_joint_reml\test_h2_ldsc\test_h2_ldsc_constrained_intercept\n')


est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X_init = med_h2_with_sumstat_ldsc_two_step_multivariate_gene_bins(gwas_beta, gwas_beta_se, gwas_variant_ld_scores.reshape(-1,1), eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, np.asarray([len(gwas_beta)]), burn_in_iters=10, total_iters=20, gibbs=False)
t.write('eqtl_5_binned_no_intercept_two_step\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_data_set_variance(gwas_beta, gwas_beta_se, gwas_variant_ld_scores.reshape(-1,1), eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, np.asarray([len(gwas_beta)]), burn_in_iters=995, total_iters=1000, gibbs=False, eqtl_only=False)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# With nm-variant gene window
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X_init = med_h2_with_sumstat_ldsc_two_step_multivariate_gene_bins(gwas_beta, gwas_beta_se, np.transpose(np.vstack((gwas_variant_ld_scores,gwas_gene_window_variant_ld_scores))), eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, np.asarray([len(gwas_beta), np.sum(gwas_gene_window_variant_anno_vec)]), burn_in_iters=10, total_iters=20, gibbs=False)
t.write('eqtl_5_binned_no_intercept_w_nm_gene_window_two_step\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_data_set_variance(gwas_beta, gwas_beta_se, np.transpose(np.vstack((gwas_variant_ld_scores,gwas_gene_window_variant_ld_scores))), eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, np.asarray([len(gwas_beta), np.sum(gwas_gene_window_variant_anno_vec)]), burn_in_iters=995, total_iters=1000, gibbs=False, eqtl_only=False)
t.write('eqtl_5_binned_no_intercept_w_nm_gene_window_bayesian_gibbs_resid_var_multivariate_per_data_set_variance\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

t.close()
print(output_file)


#####################
# Run analysis permuting the eqtls we use
#####################


# Load in gene window variant ld scores file
perm_simulation_name_string = 'simulation_' + str(simulation_number+1) + '_chrom' + simulation_name_string.split('_chrom')[1]
eqtl_sumstat_file = simulated_learned_gene_models_dir + perm_simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
single_bin_genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, eqtl_cis_snp_position_single_bin, eqtl_n_cis_snps, eqtl_classes, gene_indexer = load_in_eqtl_data_with_bins_per_gene(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, eqtl_ld='out_of_sample', n_bins=1)
gwas_gene_window_variant_ld_scores, gwas_gene_window_variant_anno_vec = get_gwas_gene_window_variant_ld_scores(gwas_beta, gwas_rsids, quasi_ld_window_summary_file, single_bin_genes, eqtl_cis_snp_position_single_bin)


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



output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod13_permuted_eqtls.txt'
t = open(output_file,'w')
t.write('method\teQTL_SS\tsim_h2\tsim_med_h2\tsim_nm_h2\test_med_h2_joint_reml\test_med_h2_per_tissue_joint_reml\test_nm_h2_joint_reml\test_mean_eqtl_h2_joint_reml\test_h2_ldsc\test_h2_ldsc_constrained_intercept\n')


est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X_init = med_h2_with_sumstat_ldsc_two_step_multivariate_gene_bins(gwas_beta, gwas_beta_se, gwas_variant_ld_scores.reshape(-1,1), eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, np.asarray([len(gwas_beta)]), burn_in_iters=10, total_iters=20, gibbs=False)
t.write('eqtl_5_binned_no_intercept_two_step\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_data_set_variance(gwas_beta, gwas_beta_se, gwas_variant_ld_scores.reshape(-1,1), eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, np.asarray([len(gwas_beta)]), burn_in_iters=995, total_iters=1000, gibbs=False, eqtl_only=False)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


# With nm-variant gene window
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X_init = med_h2_with_sumstat_ldsc_two_step_multivariate_gene_bins(gwas_beta, gwas_beta_se, np.transpose(np.vstack((gwas_variant_ld_scores,gwas_gene_window_variant_ld_scores))), eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, np.asarray([len(gwas_beta), np.sum(gwas_gene_window_variant_anno_vec)]), burn_in_iters=10, total_iters=20, gibbs=False)
t.write('eqtl_5_binned_no_intercept_w_nm_gene_window_two_step\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_data_set_variance(gwas_beta, gwas_beta_se, np.transpose(np.vstack((gwas_variant_ld_scores,gwas_gene_window_variant_ld_scores))), eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, np.asarray([len(gwas_beta), np.sum(gwas_gene_window_variant_anno_vec)]), burn_in_iters=995, total_iters=1000, gibbs=False, eqtl_only=False)
t.write('eqtl_5_binned_no_intercept_w_nm_gene_window_bayesian_gibbs_resid_var_multivariate_per_data_set_variance\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


t.close()
print(output_file)












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
'''




#####################
# Emperically evaluate correlation between genes across real and permuted runs
#####################
'''
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





























#############
# OLD
#############




'''
output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod11.txt'
t = open(output_file,'w')
t.write('method\teQTL_SS\tsim_h2\tsim_med_h2\tsim_nm_h2\test_med_h2_joint_reml\test_med_h2_per_tissue_joint_reml\test_nm_h2_joint_reml\test_mean_eqtl_h2_joint_reml\test_h2_ldsc\test_h2_ldsc_constrained_intercept\n')


est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X_init = med_h2_with_sumstat_ldsc_two_step_multivariate_gene_bins(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, burn_in_iters=10, total_iters=20, gibbs=False)
t.write('eqtl_5_binned_no_intercept_two_step\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
#output_file = trait_med_h2_inference_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod11_X_init.npy'
#np.save(output_file, X_init)



# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_single_init_variance(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=246, total_iters=250, gibbs=False, eqtl_only=False)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_single_init_variance\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_single_init_variance(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=246, total_iters=250, gibbs=False, eqtl_only=True)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_single_init_variance_cis_eqtl_only\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_single_init_variance_per_gene(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=246, total_iters=250, gibbs=False, eqtl_only=False)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_single_init_per_gene_variance\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_single_init_variance_per_gene(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=246, total_iters=250, gibbs=False, eqtl_only=True)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_single_init_per_gene_variance_cis_eqtl_only\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

'''


'''
# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_gene_resid_var(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, burn_in_iters=995, total_iters=1000, gibbs=False)
#est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, burn_in_iters=1500, total_iters=1505, gibbs=False)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_gene_resid_var\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
output_file = trait_med_h2_inference_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod11_X.npy'
np.save(output_file, X)


# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2, X = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_gene_bins_per_gene_resid_var_eqtl_snp_only(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, gene_indexer, eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=995, total_iters=1000, gibbs=False)
t.write('eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_gene_resid_var_eqtl_snp_only\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
output_file = trait_med_h2_inference_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod11_X2.npy'
np.save(output_file, X)


'''

'''


# load in eqtl data
# Out of sample eqtl ld
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, eqtl_cis_snp_position, eqtl_n_cis_snps, eqtl_classes = load_in_eqtl_data(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, eqtl_ld='out_of_sample')
# Generate matrix form of eqtl data
eqtl_beta_mat = get_matrix_form_of_eqtl_data(eqtl_beta, eqtl_position, 0.0, len(gwas_variant_ld_scores))
eqtl_beta_se_mat = get_matrix_form_of_eqtl_data(eqtl_beta_se, eqtl_position, 1.0/np.sqrt(eqtl_sample_size), len(gwas_variant_ld_scores))
eqtl_ld_score_mat = get_matrix_form_of_eqtl_data(eqtl_ldscore, eqtl_position, 0.0, len(gwas_variant_ld_scores))


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
f = open('eqtl_cis_snp_position.pkl', 'wb')
pickle.dump(eqtl_cis_snp_position, f, pickle.HIGHEST_PROTOCOL)
f.close()



'''
'''
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

print(ldsc_constrained_intercept_snp_h2)
'''
# Iterative equal weight
#est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, max_epochs=9000)
#t.write('no_intercept_equal_weight\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Two step (implicitely equal weight)
#est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_two_step(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes)
#t.write('no_intercept_equal_weight_two_step\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Iterative (learn unequal weights bayesian)
#est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_gibbs(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, eqtl_cis_snp_position, burn_in_iters=400, total_iters=600)
#t.write('no_intercept_bayesian_gibbs_resid_var\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Iterative (learn unequal weights bayesian)
#est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, eqtl_cis_snp_position, burn_in_iters=400, total_iters=600)
#t.write('no_intercept_bayesian_gibbs_resid_var_multivariate\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')



# Iterative (learn unequal weights bayesian)
#est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate_only_cis_eqtl_snp_pos(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, eqtl_cis_snp_position, burn_in_iters=400, total_iters=600)
#t.write('no_intercept_bayesian_gibbs_resid_var_multivariate_only_cis_eqtl_snps\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

# Iterative (learn unequal weights bayesian)
#est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_per_gene_resid_var_only_cis_eqtl_snps(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, eqtl_cis_snp_position, burn_in_iters=400, total_iters=600, gibbs=False)
#t.write('no_intercept_bayesian_gibbs_per_gene_resid_var_multivariate_only_cis_eqtl_snps\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

'''
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_fixed_sample_size_plus_h2_resid_var_non_central_chi_sq(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes,eqtl_sample_size, burn_in_iters=400, total_iters=600, gibbs=True)
t.write('no_intercept_bayesian_gibbs_sample_size_plus_h2_based_resid_var_multivariate_non_central_chi_sq\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_fixed_sample_size_plus_h2_resid_var_only_cis_eqtl_non_central_chi_sq(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes,eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=400, total_iters=600, gibbs=True)
t.write('no_intercept_bayesian_gibbs_sample_size_plus_h2_based_only_cis_eqtl_resid_var_multivariate_non_central_chi_sq\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
'''
'''
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_fixed_sample_size_plus_h2_resid_var_only_cis_eqtl(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes,eqtl_sample_size, eqtl_cis_snp_position, burn_in_iters=400, total_iters=600, gibbs=True)
t.write('no_intercept_bayesian_gibbs_sample_size_plus_h2_based_only_cis_eqtl_resid_var_multivariate\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_fixed_sample_size_plus_h2_resid_var(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes,eqtl_sample_size, burn_in_iters=400, total_iters=600, gibbs=True)
t.write('no_intercept_bayesian_gibbs_sample_size_plus_h2_based_resid_var_multivariate\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_fixed_sample_size_resid_var(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes,eqtl_sample_size, burn_in_iters=400, total_iters=600, gibbs=True)
t.write('no_intercept_bayesian_gibbs_sample_size_based_resid_var_multivariate\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
'''



'''
# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_per_gene_resid_var(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, burn_in_iters=400, total_iters=600, gibbs=True)
t.write('no_intercept_bayesian_gibbs_per_gene_resid_var_multivariate\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_per_gene_resid_var(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, burn_in_iters=400, total_iters=600, gibbs=False)
t.write('no_intercept_bayesian_als_per_gene_resid_var_multivariate\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
'''
# Iterative (learn unequal weights bayesian)
#est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_gibbs_with_intercepts(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, eqtl_cis_snp_position, burn_in_iters=9500, total_iters=10000)
#t.write('no_intercept_bayesian_gibbs_resid_var_w_intercepts\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


# Iterative (learn unequal weights bayesian)
#est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_gibbs_per_data_set_noise(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, eqtl_cis_snp_position, burn_in_iters=900, total_iters=1000)
#t.write('no_intercept_bayesian_gibbs_resid_var_per_data_set\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

###########################################
# Remove all tissues except causal tissue
###########################################
'''
indices = eqtl_classes == 'tissue0'
causal_eqtl_beta_mat = eqtl_beta_mat[indices,:]
causal_eqtl_beta_se_mat = eqtl_beta_se_mat[indices,:]
causal_eqtl_ld_score_mat = eqtl_ld_score_mat[indices,:]
causal_eqtl_n_cis_snps = eqtl_n_cis_snps[indices]
causal_eqtl_classes = eqtl_classes[indices]
'''

'''
# Two step (implicitely equal weight)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_two_step(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, causal_eqtl_beta_mat, causal_eqtl_beta_se_mat, causal_eqtl_ld_score_mat, causal_eqtl_n_cis_snps, causal_eqtl_classes)
t.write('causal_tissue_only_no_intercept_equal_weight_two_step\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_gibbs(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, causal_eqtl_beta_mat, causal_eqtl_beta_se_mat, causal_eqtl_ld_score_mat, causal_eqtl_n_cis_snps, causal_eqtl_classes, burn_in_iters=400, total_iters=600)
t.write('causal_tissue_only_no_intercept_bayesian_gibbs_resid_var\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


# Iterative (learn unequal weights bayesian)
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_gibbs_multivariate(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, causal_eqtl_beta_mat, causal_eqtl_beta_se_mat, causal_eqtl_ld_score_mat, causal_eqtl_n_cis_snps, causal_eqtl_classes, burn_in_iters=400, total_iters=600)
t.write('causal_tissue_only_no_intercept_bayesian_gibbs_resid_var_multivariate\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')


est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_per_gene_resid_var(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, causal_eqtl_beta_mat, causal_eqtl_beta_se_mat, causal_eqtl_ld_score_mat, causal_eqtl_n_cis_snps, causal_eqtl_classes, burn_in_iters=400, total_iters=600, gibbs=True)
t.write('causal_tissue_only_no_intercept_bayesian_gibbs_per_gene_resid_var_multivariate\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_bayesian_per_gene_resid_var(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, causal_eqtl_beta_mat, causal_eqtl_beta_se_mat, causal_eqtl_ld_score_mat, causal_eqtl_n_cis_snps, causal_eqtl_classes, burn_in_iters=400, total_iters=600, gibbs=False)
t.write('causal_tissue_only_no_intercept_bayesian_als_per_gene_resid_var_multivariate\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
'''






# load in eqtl data
# In-sample (adjusted) eqtl ld 
'''
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, eqtl_cis_snp_position, eqtl_n_cis_snps, eqtl_classes = load_in_eqtl_data(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, eqtl_ld='in_sample_adjusted')
# Generate matrix form of eqtl data
eqtl_beta_mat = get_matrix_form_of_eqtl_data(eqtl_beta, eqtl_position, 0.0, len(gwas_variant_ld_scores))
eqtl_beta_se_mat = get_matrix_form_of_eqtl_data(eqtl_beta_se, eqtl_position, 1.0/np.sqrt(eqtl_sample_size), len(gwas_variant_ld_scores))
eqtl_ld_score_mat = get_matrix_form_of_eqtl_data(eqtl_ldscore, eqtl_position, 0.0, len(gwas_variant_ld_scores))

est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, max_epochs=9000)
t.write('no_intercept_equal_weight_in_sample_adjusted_eqtl_ldscores\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
'''

'''
# load in eqtl data
# In-sample (unadjusted) eqtl ld 
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, eqtl_cis_snp_position, eqtl_n_cis_snps, eqtl_classes = load_in_eqtl_data(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, eqtl_ld='in_sample_unadjusted')
# Generate matrix form of eqtl data
eqtl_beta_mat = get_matrix_form_of_eqtl_data(eqtl_beta, eqtl_position, 0.0, len(gwas_variant_ld_scores))
eqtl_beta_se_mat = get_matrix_form_of_eqtl_data(eqtl_beta_se, eqtl_position, 1.0/np.sqrt(eqtl_sample_size), len(gwas_variant_ld_scores))
eqtl_ld_score_mat = get_matrix_form_of_eqtl_data(eqtl_ldscore, eqtl_position, 0.0, len(gwas_variant_ld_scores))

est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, max_epochs=9000)
t.write('no_intercept_equal_weight_in_sample_unadjusted_eqtl_ldscores\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')
'''









# NO longer used
'''
est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_w_intercept(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, max_epochs=9000)
t.write('intercept_equal_weight\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

est_total_med_h2, est_per_tissue_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_ldsc_w_per_study_residual_variances(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, eqtl_classes, max_epochs=9000)
t.write('no_intercept_predefined_weight\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_total_med_h2) + '\t' + ','.join(est_per_tissue_med_h2.astype(str)) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

'''


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
'''








# DEBUGGING A BIT: run ld score regression on eqtl data
'''
n_genes = eqtl_beta_mat.shape[0]
est_gene_h2s = []
est_gene_h2s2 = []
sim_gene_h2s = []

for gene_iter in range(n_genes):
	indices = eqtl_ld_score_mat[gene_iter,:] != 0.0

	gene_tissue_name = genes[gene_iter]
	gene_name = gene_tissue_name.split(':')[0]
	tissue_name = int(gene_tissue_name.split(':tissue')[1])
	causal_eqtl_effect_file = '/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/multi_tissue_simulation_pca/simulated_gene_expression/' + simulation_name_string + '_' + gene_name + '_causal_eqtl_effects.npy'
	causal_eqtls = np.load(causal_eqtl_effect_file)
	sim_gene_h2 =  np.sum(np.square(causal_eqtls),axis=0)[tissue_name]
	sim_gene_h2s.append(sim_gene_h2)

	Y = np.square(eqtl_beta_mat[gene_iter, indices]) - np.square(eqtl_beta_se_mat[gene_iter, indices])
	X = eqtl_ld_score_mat[gene_iter, indices]
	model = sm.OLS(Y,X).fit()
	cis_h2 = eqtl_n_cis_snps[gene_iter]*model.params[0]
	est_gene_h2s.append(cis_h2)

	Y = np.square(eqtl_beta_mat[gene_iter, indices]/eqtl_beta_se_mat[gene_iter, indices]) - 1.0
	X = eqtl_ld_score_mat[gene_iter, indices]
	model = sm.OLS(Y,X).fit()
	cis_h2 = eqtl_n_cis_snps[gene_iter]*model.params[0]/eqtl_sample_size
	est_gene_h2s2.append(cis_h2)

pdb.set_trace()

output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + str(eqtl_sample_size) + '_eqtl_h2_debug.txt'
t = open(output_file,'w')
t.write(str(simulation_number) + '\t' + str(eqtl_sample_size) + '\t' + str(np.mean(est_gene_h2s)) + '\t' + str(np.mean(est_gene_h2s2)) + '\t' + str(np.mean(sim_gene_h2s)) + '\n')
t.close()

print(output_file)
'''

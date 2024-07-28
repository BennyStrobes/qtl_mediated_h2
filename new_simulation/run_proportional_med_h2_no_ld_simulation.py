import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm

'''
import tensorflow as tf
import gzip
import time
import statsmodels.api as sm
import tensorflow_probability as tfp
'''



def simulate_data(gwas_ss, n_snps, eqtl_ss_1, med_h2_1, ge_h2_options, n_genes, snps_per_gene, nm_h2, eqtl_architecture, frac_causal_genes):
	gene_midpoint_lb = int(snps_per_gene/2) +1
	gene_midpoint_ub = n_snps - (int(snps_per_gene/2) +1)

	gene_causal_eqtl_effects = []
	gene_snp_indices_arr = []
	sim_ge_h2s = []
	gene_counts = np.zeros(n_snps)
	for gene_iter in range(n_genes):

		gene_midpoint = np.random.choice(np.arange(gene_midpoint_lb, gene_midpoint_ub))
		gene_start = gene_midpoint - int(snps_per_gene/2)
		gene_end = gene_midpoint + int(snps_per_gene/2)
		gene_snp_indices = np.arange(gene_start, gene_end)
		gene_counts[gene_snp_indices] = gene_counts[gene_snp_indices] + 1.0/len(gene_snp_indices)

		# Choose gene expression h2
		ge_h2 = np.random.choice(ge_h2_options)
		sim_ge_h2s.append(ge_h2)


		if eqtl_architecture == 'polygenic':
			gene_snp_causal_eqtl_effects = np.zeros(n_snps)
			gene_snp_causal_eqtl_effects[gene_snp_indices] = np.random.normal(loc=0, scale=np.sqrt(ge_h2/snps_per_gene), size=snps_per_gene)
		elif eqtl_architecture == 'sparse':
			gene_snp_causal_eqtl_effects = np.zeros(n_snps)
			causal_indices = np.random.choice(gene_snp_indices, size=10, replace=False)
			gene_snp_causal_eqtl_effects[causal_indices] = np.random.normal(loc=0, scale=np.sqrt(ge_h2/10), size=10)
		else:
			print('assumption error: eqtl architecture ' + eqtl_architecture + ' currently not implemented')
			pdb.set_trace()

		gene_causal_eqtl_effects.append(gene_snp_causal_eqtl_effects)
		gene_snp_indices_arr.append(gene_snp_indices)
	gene_causal_eqtl_effects = np.asarray(gene_causal_eqtl_effects)

	sim_alpha_1_sq = med_h2_1/np.mean(sim_ge_h2s)

	# Simulate eQTL genotype data
	eqtl_geno_1 = np.random.normal(loc=0,scale=1.0,size=(eqtl_ss_1,n_snps))
	for snp_iter in range(n_snps):
		eqtl_geno_1[:,snp_iter] = (eqtl_geno_1[:,snp_iter] - np.mean(eqtl_geno_1[:,snp_iter]))/np.std(eqtl_geno_1[:,snp_iter])

	# Simulate GWAS genotype data
	gwas_geno = np.random.normal(loc=0,scale=1.0,size=(gwas_ss,n_snps))
	for snp_iter in range(n_snps):
		gwas_geno[:,snp_iter] = (gwas_geno[:,snp_iter] - np.mean(gwas_geno[:,snp_iter]))/np.std(gwas_geno[:,snp_iter])

	# Simulate Gene trait value
	Es = []
	for gene_iter in range(n_genes):
		genetic_gene_trait_1 = np.dot(eqtl_geno_1, gene_causal_eqtl_effects[gene_iter, :])
		E_1 = np.random.normal(loc=genetic_gene_trait_1, scale=np.sqrt(1.0-np.var(genetic_gene_trait_1)))
		E_1 = (E_1 - np.mean(E_1))/np.std(E_1)
		Es.append(E_1)

	# Simulate gwas trait value
	n_causal_genes = int(np.floor(n_genes*frac_causal_genes))
	sim_alphas = np.zeros(n_genes)
	causal_gene_indices = np.random.choice(np.arange(n_genes), size=n_causal_genes, replace=False)
	sim_alphas[causal_gene_indices] = np.random.normal(loc=0, scale=np.sqrt(sim_alpha_1_sq/n_causal_genes), size=n_causal_genes)
	#sim_alphas = np.random.normal(loc=0, scale=np.sqrt(sim_alpha_1_sq/n_genes), size=n_genes)
	genetic_trait = np.zeros(gwas_ss)
	for gene_iter in range(n_genes):
		genetic_gene = np.dot(gwas_geno, gene_causal_eqtl_effects[gene_iter, :])
		genetic_trait =genetic_trait + genetic_gene*sim_alphas[gene_iter]
	print(np.var(genetic_trait))
	# Causal nm var-trait effect
	n_causal_snps = 3000
	sim_beta = np.zeros(n_snps)
	causal_indices = np.random.choice(np.arange(n_snps), size=n_causal_snps, replace=False)
	sim_beta[causal_indices] = np.random.normal(loc=0, scale=np.sqrt(nm_h2/n_causal_snps), size=n_causal_snps)
	genetic_trait = genetic_trait + np.dot(gwas_geno, sim_beta)
	print(np.var(genetic_trait))

	Y = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0-np.var(genetic_trait)))
	Y = (Y - np.mean(Y))/np.std(Y)

	return Y, Es, gwas_geno, eqtl_geno_1, gene_causal_eqtl_effects, sim_beta, sim_alphas, gene_snp_indices_arr, np.asarray(sim_ge_h2s), gene_counts





def get_marginal_summary_statistics(Y, G):
	n_snps = G.shape[1]
	beta = []
	beta_se = []
	for snp_iter in range(n_snps):
		olser = sm.OLS(Y, sm.add_constant(G[:,snp_iter])).fit()
		beta.append(olser.params[1])
		beta_se.append(olser.bse[1])
	return np.asarray(beta), np.asarray(beta_se)

def get_marginal_summary_statistics_across_genes(Es, G, gene_snp_indices_arr):
	n_genes = len(Es)
	n_snps = E_geno.shape[1]

	beta_arr = []
	beta_se_arr = []
	# Loop through genes
	for gene_iter in range(n_genes):
		beta = np.zeros(n_snps)
		beta_se = np.ones(n_snps)*1e-7
		gene_snp_indices = gene_snp_indices_arr[gene_iter]

		for gene_snp_index in gene_snp_indices:
			olser = sm.OLS(Es[gene_iter], sm.add_constant(G[:,gene_snp_index])).fit()
			beta[gene_snp_index] = olser.params[1]
			beta_se[gene_snp_index] = olser.bse[1]

		# Add to global array
		beta_arr.append(beta)
		beta_se_arr.append(beta_se)

	beta_arr = np.asarray(beta_arr)
	beta_se_arr = np.asarray(beta_se_arr)

	return beta_arr, beta_se_arr

def compute_univariate_gaussian_log_like(xx, mean_value, variance_vec):
	log_like = -(tf.math.log(variance_vec)/2) - tf.math.divide(tf.square(xx-mean_value), (2.0*variance_vec))
	
	return log_like


def sumstat_proportional_joint_reml_chunked_no_ld_loss(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, per_snp_eqtl_h2_est_lb, per_snp_gwas_h2_est_lb, chunk_to_snp_mappings):

	eqtl_beta_sq_variable_scaled = tf.nn.relu(eqtl_beta_sq_variable) + (per_snp_eqtl_h2_est_lb)


	#big_eqtl_beta_sq_variable = tf.convert_to_tensor(np.zeros(eqtl_mask.shape), dtype=tf.float32)
	#big_eqtl_beta_sq_variable = tf.tile(tf.reshape(eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps])

	pred_per_snp_gwas_h2 = beta_sq_variable

	eqtl_log_like = 0.0
	for gene_iter, chunk_to_snp_mapping in enumerate(chunk_to_snp_mappings):
		gene_eqtl_beta_sq = tf.matmul(eqtl_beta_sq_variable_scaled[gene_iter:(gene_iter)+1,:], chunk_to_snp_mapping)[0,:]
		pred_per_snp_gwas_h2 = pred_per_snp_gwas_h2 + (alpha_sq_variable*gene_eqtl_beta_sq*eqtl_mask[gene_iter,:])

		indices = eqtl_mask[gene_iter,:] != 0.0
		temp_eqtl_log_like = compute_univariate_gaussian_log_like(eqtl_beta[gene_iter, :][indices], 0.0, eqtl_beta_var[gene_iter,:][indices] + gene_eqtl_beta_sq[indices])

		eqtl_log_like = eqtl_log_like + tf.math.reduce_sum(temp_eqtl_log_like)


	#pred_per_snp_gwas_h2 = beta_sq_variable + alpha_sq_variable*(tf.math.reduce_sum((big_eqtl_beta_sq_variable)*eqtl_mask,axis=0))


	gwas_log_like = compute_univariate_gaussian_log_like(gwas_beta, 0.0, gwas_beta_var + pred_per_snp_gwas_h2)


	#eqtl_indices = eqtl_mask==1.0
	#eqtl_log_like = compute_univariate_gaussian_log_like(eqtl_beta[eqtl_indices], 0.0, eqtl_beta_var[eqtl_indices] + big_eqtl_beta_sq_variable[eqtl_indices])

	loss = -tf.reduce_sum(gwas_log_like) - eqtl_log_like

	return loss


def sumstat_proportional_joint_reml_no_ld_loss(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, per_snp_eqtl_h2_est_lb, per_snp_gwas_h2_est_lb):

	eqtl_beta_sq_variable_scaled = tf.nn.relu(eqtl_beta_sq_variable) + (per_snp_eqtl_h2_est_lb)

	big_eqtl_beta_sq_variable = tf.tile(tf.reshape(eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps])



	pred_per_snp_gwas_h2 = beta_sq_variable + alpha_sq_variable*(tf.math.reduce_sum((big_eqtl_beta_sq_variable)*eqtl_mask,axis=0))

	#thresholded_pred_per_snp_gwas_h2 = tf.maximum(pred_per_snp_gwas_h2, per_snp_gwas_h2_est_lb)
	#pred_per_snp_gwas_h2[pred_per_snp_gwas_h2 < per_snp_gwas_h2_est_lb] = per_snp_gwas_h2_est_lb

	gwas_log_like = compute_univariate_gaussian_log_like(gwas_beta, 0.0, gwas_beta_var + pred_per_snp_gwas_h2)


	eqtl_indices = eqtl_mask==1.0
	eqtl_log_like = compute_univariate_gaussian_log_like(eqtl_beta[eqtl_indices], 0.0, eqtl_beta_var[eqtl_indices] + big_eqtl_beta_sq_variable[eqtl_indices])


	loss = -tf.reduce_sum(gwas_log_like) - tf.reduce_sum(eqtl_log_like)

	return loss

def sumstat_proportional_joint_reml_per_snp_no_ld_loss(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, per_snp_eqtl_h2_est_lb, per_snp_gwas_h2_est_lb, eqtl_indices):

	big_eqtl_beta_sq_variable = tf.nn.relu(eqtl_beta_sq_variable) + (per_snp_eqtl_h2_est_lb)

	#big_eqtl_beta_sq_variable = tf.tile(tf.reshape(eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps])

	pred_per_snp_gwas_h2 = beta_sq_variable + alpha_sq_variable*(tf.math.reduce_sum((big_eqtl_beta_sq_variable)*eqtl_mask,axis=0))

	#thresholded_pred_per_snp_gwas_h2 = tf.maximum(pred_per_snp_gwas_h2, per_snp_gwas_h2_est_lb)
	#pred_per_snp_gwas_h2[pred_per_snp_gwas_h2 < per_snp_gwas_h2_est_lb] = per_snp_gwas_h2_est_lb

	gwas_log_like = compute_univariate_gaussian_log_like(gwas_beta, 0.0, tf.nn.relu(gwas_beta_var + pred_per_snp_gwas_h2)+1e-12)


	eqtl_log_like = compute_univariate_gaussian_log_like(eqtl_beta[eqtl_indices], 0.0, eqtl_beta_var[eqtl_indices] + big_eqtl_beta_sq_variable[eqtl_indices])


	loss = -tf.reduce_sum(gwas_log_like) - tf.reduce_sum(eqtl_log_like)

	return loss


def med_h2_with_sumstat_reml_no_ld_two_step(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, max_epochs=60000, conv_thresh=1e-12):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)


	optimizer = tf.keras.optimizers.Adam(learning_rate=1e-5)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)

	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)

	# Convert variabels to tf tensors
	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)

	per_snp_eqtl_h2_est_lb = -np.min(eqtl_beta_var[eqtl_mask==1])*.99
	per_snp_gwas_h2_est_lb = -np.min(gwas_beta_var)*.99



	# Step 1: estimate beta-sq
	eqtl_beta_sq_init = []
	for gene_iter in range(n_genes):
		gene_eqtl_beta = eqtl_beta[gene_iter,:]
		gene_eqtl_beta_var = eqtl_beta_var[gene_iter, :]

		indices = gene_eqtl_beta != 0.0
		gene_eqtl_beta2 = gene_eqtl_beta[indices]
		gene_eqtl_beta_var2 = gene_eqtl_beta_var[indices]
		eqtl_beta_sq_init.append(np.mean(np.square(gene_eqtl_beta2/np.sqrt(gene_eqtl_beta_var2)) - 1.0)*np.mean(gene_eqtl_beta_var2))
	eqtl_beta_sq_init = np.asarray(eqtl_beta_sq_init)


	# Initialize variables to optimize over
	beta_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='beta_sq')
	alpha_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='alpha_sq')
	eqtl_beta_sq_variable = tf.Variable(initial_value=eqtl_beta_sq_init- per_snp_eqtl_h2_est_lb,trainable=False, name='eqtl_beta_sq', dtype=tf.float32)


	converged = False
	prev_est_alpha_sq=10000
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = sumstat_proportional_joint_reml_no_ld_loss(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, per_snp_eqtl_h2_est_lb, per_snp_gwas_h2_est_lb)

		trainable_variables = []
		trainable_variables.append(alpha_sq_variable)
		trainable_variables.append(beta_sq_variable)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est = np.asarray(alpha_sq_variable)*1.0

		diff = np.abs(prev_est_alpha_sq -cur_est)
		if diff < conv_thresh:
			converged = True
			break

		prev_est_alpha_sq = cur_est

		print(cur_est)
		print(loss_value)

	if converged == False:
		print('did not converge')
		print(diff)

	#print(np.mean(np.sum(eqtl_beta_sq_var,axis=1)))
	#print(np.mean(np.sum(eqtl_beta_sq_var,axis=1))*cur_est)

	#print(loss_value)
	eqtl_beta_sq_variable_scaled = tf.nn.relu(eqtl_beta_sq_variable) + (per_snp_eqtl_h2_est_lb)
	big_eqtl_beta_sq_variable = tf.tile(tf.reshape(eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps])

	med_h2 = np.sum(alpha_sq_variable*(tf.math.reduce_sum((big_eqtl_beta_sq_variable)*eqtl_mask,axis=0)))
	nm_h2 = np.sum(beta_sq_variable)*n_snps

	gene_cis_h2 = eqtl_beta_sq_variable_scaled*snps_per_gene_arr

	return med_h2, nm_h2, gene_cis_h2



def med_h2_with_sumstat_reml_chunked_no_ld(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, chunks_per_gene=10, max_epochs=80000, conv_thresh=1e-12):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)


	optimizer = tf.keras.optimizers.Adam(learning_rate=1e-5)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)




	# Convert variabels to tf tensors
	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)

	per_snp_eqtl_h2_est_lb = -np.min(eqtl_beta_var[eqtl_mask==1])*.99
	per_snp_gwas_h2_est_lb = -np.min(gwas_beta_var)*.99



	# Initialize variables to optimize over
	beta_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='beta_sq')
	alpha_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='alpha_sq')
	eqtl_beta_sq_variable = tf.Variable(initial_value=np.ones((n_genes,chunks_per_gene))*0.0 - per_snp_eqtl_h2_est_lb,trainable=True, name='eqtl_beta_sq', dtype=tf.float32)


	big_eqtl_beta_sq_variable = tf.convert_to_tensor(np.zeros(eqtl_mask.shape), dtype=tf.float32)

	chunk_to_snp_mappings = []
	for gene_iter in range(n_genes):
		chunk_to_snp_mapping = np.zeros((chunks_per_gene, n_snps))
		gene_snps = eqtl_mask[gene_iter,:]
		
		chunk_index_array = np.array_split(np.where(gene_snps==1)[0], chunks_per_gene)


		if len(chunk_index_array) != chunks_per_gene:
			print('assumption eroror')
			pdb.set_trace()

		for chunk_iter, chunk_indices in enumerate(chunk_index_array):
			chunk_to_snp_mapping[chunk_iter, chunk_indices] = 1.0
			#big_eqtl_beta_sq_variable[gene_iter, :] = tf.gather(big_eqtl_beta_sq_variable[gene_iter, :], tf.expand_dims(chunk_indices,axis=1), eqtl_beta_sq_variable[gene_iter,chunk_iter])[:,0]

		chunk_to_snp_mapping = tf.convert_to_tensor(chunk_to_snp_mapping, dtype=tf.float32)
		chunk_to_snp_mappings.append(chunk_to_snp_mapping)




	converged = False
	prev_est_alpha_sq=10000
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = sumstat_proportional_joint_reml_chunked_no_ld_loss(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, per_snp_eqtl_h2_est_lb, per_snp_gwas_h2_est_lb, chunk_to_snp_mappings)

		trainable_variables = []
		trainable_variables.append(alpha_sq_variable)
		trainable_variables.append(eqtl_beta_sq_variable)
		trainable_variables.append(beta_sq_variable)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est = np.asarray(alpha_sq_variable)*1.0

		diff = np.abs(prev_est_alpha_sq -cur_est)
		if diff < conv_thresh:
			converged = True
			break

		prev_est_alpha_sq = cur_est

		print(cur_est)
		print(loss_value)

	if converged == False:
		print('did not converge')
		print(diff)

	eqtl_beta_sq_variable_scaled = tf.nn.relu(eqtl_beta_sq_variable) + (per_snp_eqtl_h2_est_lb)
	big_eqtl_beta_sq_variable = tf.tile(tf.reshape(eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps])

	med_h2 = np.sum(alpha_sq_variable*(tf.math.reduce_sum((big_eqtl_beta_sq_variable)*eqtl_mask,axis=0)))
	nm_h2 = np.sum(beta_sq_variable)*n_snps

	gene_cis_h2 = eqtl_beta_sq_variable_scaled*snps_per_gene_arr



	return med_h2, nm_h2, gene_cis_h2


def med_h2_with_sumstat_reml_no_ld(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, max_epochs=80000, conv_thresh=1e-12):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)


	optimizer = tf.keras.optimizers.Adam(learning_rate=1e-5)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)


	# Convert variabels to tf tensors
	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)

	per_snp_eqtl_h2_est_lb = -np.min(eqtl_beta_var[eqtl_mask==1])*.99
	per_snp_gwas_h2_est_lb = -np.min(gwas_beta_var)*.99



	# Initialize variables to optimize over
	beta_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='beta_sq')
	alpha_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='alpha_sq')
	eqtl_beta_sq_variable = tf.Variable(initial_value=np.ones(n_genes)*0.0 - per_snp_eqtl_h2_est_lb,trainable=True, name='eqtl_beta_sq', dtype=tf.float32)


	converged = False
	prev_est_alpha_sq=10000
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = sumstat_proportional_joint_reml_no_ld_loss(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, per_snp_eqtl_h2_est_lb, per_snp_gwas_h2_est_lb)

		trainable_variables = []
		trainable_variables.append(alpha_sq_variable)
		trainable_variables.append(eqtl_beta_sq_variable)
		trainable_variables.append(beta_sq_variable)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est = np.asarray(alpha_sq_variable)*1.0

		diff = np.abs(prev_est_alpha_sq -cur_est)
		if diff < conv_thresh:
			converged = True
			break

		prev_est_alpha_sq = cur_est

		print(cur_est)
		print(loss_value)

	if converged == False:
		print('did not converge')
		print(diff)

	eqtl_beta_sq_variable_scaled = tf.nn.relu(eqtl_beta_sq_variable) + (per_snp_eqtl_h2_est_lb)
	big_eqtl_beta_sq_variable = tf.tile(tf.reshape(eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps])

	med_h2 = np.sum(alpha_sq_variable*(tf.math.reduce_sum((big_eqtl_beta_sq_variable)*eqtl_mask,axis=0)))
	nm_h2 = np.sum(beta_sq_variable)*n_snps

	gene_cis_h2 = eqtl_beta_sq_variable_scaled*snps_per_gene_arr



	return med_h2, nm_h2, gene_cis_h2





def med_h2_with_sumstat_reml_per_snp_no_ld(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, max_epochs=80000, conv_thresh=1e-12):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)


	optimizer = tf.keras.optimizers.Adam(learning_rate=1e-5)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)
	eqtl_indices = eqtl_mask==1.0



	# Convert variabels to tf tensors
	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)

	per_snp_eqtl_h2_est_lb = -np.min(eqtl_beta_var[eqtl_mask==1])*.99
	per_snp_gwas_h2_est_lb = -np.min(gwas_beta_var)*.99



	# Initialize variables to optimize over
	beta_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='beta_sq')
	alpha_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='alpha_sq')
	eqtl_beta_sq_variable = tf.Variable(initial_value=np.ones(eqtl_mask.shape)*0.0 - per_snp_eqtl_h2_est_lb,trainable=True, name='eqtl_beta_sq', dtype=tf.float32)

	converged = False
	prev_est_alpha_sq=10000
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = sumstat_proportional_joint_reml_per_snp_no_ld_loss(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, per_snp_eqtl_h2_est_lb, per_snp_gwas_h2_est_lb, eqtl_indices)

		trainable_variables = []
		trainable_variables.append(alpha_sq_variable)
		trainable_variables.append(eqtl_beta_sq_variable)
		trainable_variables.append(beta_sq_variable)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est = np.asarray(alpha_sq_variable)*1.0

		diff = np.abs(prev_est_alpha_sq -cur_est)
		if diff < conv_thresh:
			converged = True
			#break

		prev_est_alpha_sq = cur_est

		if np.mod(epoch_iter,100) == 0.0:

			big_eqtl_beta_sq_variable = tf.nn.relu(eqtl_beta_sq_variable) + (per_snp_eqtl_h2_est_lb)

			med_h2 = np.sum(alpha_sq_variable*(tf.math.reduce_sum((big_eqtl_beta_sq_variable)*eqtl_mask,axis=0)))
			nm_h2 = np.sum(beta_sq_variable)*n_snps
			print(med_h2)
			print(nm_h2)

			print(loss_value)

	if converged == False:
		print('did not converge')
		print(diff)

	big_eqtl_beta_sq_variable = tf.nn.relu(eqtl_beta_sq_variable) + (per_snp_eqtl_h2_est_lb)

	med_h2 = np.sum(alpha_sq_variable*(tf.math.reduce_sum((big_eqtl_beta_sq_variable)*eqtl_mask,axis=0)))
	nm_h2 = np.sum(beta_sq_variable)*n_snps

	gene_cis_h2 = np.sum(big_eqtl_beta_sq_variable,axis=1)
	pdb.set_trace()



	return med_h2, nm_h2, gene_cis_h2


def update_alpha(gwas_beta_resid, gwas_beta_var, alpha, alpha_var, deltas, eqtl_indices):
	n_genes = len(alpha)
	# Loop through genes
	for gg in np.random.permutation(range(n_genes)):
		# Re include effects of current gene
		cis_gwas_beta = gwas_beta_resid[eqtl_indices[gg]] + deltas[gg]*alpha[gg]
		cis_gwas_beta_var = gwas_beta_var[eqtl_indices[gg]]

		# Compute posterior distribution
		posterior_var = 1.0/(np.sum(np.square(deltas[gg])/cis_gwas_beta_var) + (1.0/alpha_var))
		posterior_mean = np.sum(cis_gwas_beta*deltas[gg]/cis_gwas_beta_var)*posterior_var

		# Sample
		alpha[gg] = np.random.normal(loc=posterior_mean, scale=np.sqrt(posterior_var))


		# Remove updated effects of this gene
		gwas_beta_resid[eqtl_indices[gg]] = cis_gwas_beta - deltas[gg]*alpha[gg]
		
	return alpha, gwas_beta_resid


def update_deltas(gwas_beta_resid, gwas_beta_var, eqtl_betas, eqtl_beta_vars, alpha, deltas, delta_vars, eqtl_indices):
	n_genes = len(alpha)
	# Loop through genes
	for gg in np.random.permutation(range(n_genes)):

		# Re include effects of current gene
		cis_gwas_beta = gwas_beta_resid[eqtl_indices[gg]] + deltas[gg]*alpha[gg]
		cis_gwas_beta_var = gwas_beta_var[eqtl_indices[gg]]

		# Compute posterior distribution
		posterior_vars = 1.0/((np.square(alpha[gg])/cis_gwas_beta_var) + (1.0/eqtl_beta_vars[gg]) + (1.0/delta_vars[gg]))
		posterior_means = ((cis_gwas_beta*alpha[gg]/cis_gwas_beta_var) + (eqtl_betas[gg]/eqtl_beta_vars[gg]))*posterior_vars

		# Sample from posterior
		deltas[gg] = np.random.normal(loc=posterior_means, scale=np.sqrt(posterior_vars))

		# Remove updated effects of this gene
		gwas_beta_resid[eqtl_indices[gg]] = cis_gwas_beta - deltas[gg]*alpha[gg]

	return deltas, gwas_beta_resid

def update_gamma(gwas_beta_resid, gwas_beta_var, gamma, gamma_var):
	# Re-include current effects
	gwas_beta_resid = gwas_beta_resid + gamma

	# Compute posterior distribution
	posterior_vars = 1.0/((1.0/gwas_beta_var) + (1.0/gamma_var))
	posterior_means = (gwas_beta_resid/gwas_beta_var)*posterior_vars


	# Sample from posterior distribution
	gamma = np.random.normal(loc=posterior_means, scale=np.sqrt(posterior_vars))


	gwas_beta_resid = gwas_beta_resid - gamma

	return gamma, gwas_beta_resid

def update_gamma_var(gamma):
	vv = len(gamma)
	tau_sq = np.sum(np.square(gamma))/len(gamma)


	# Initialize inverse gamma distribution
	invgamma_dist = invgamma(vv/2 + 1e-5, scale=vv*tau_sq/2 + 1e-5)
	# Sample from it
	gamma_var = invgamma_dist.rvs(size=1)[0]

	return gamma_var



def med_h2_with_sumstat_bayesian_no_ld(gwas_beta, gwas_beta_se, eqtl_beta_raw, eqtl_beta_se_raw, max_iters=18000, burn_in_iter=15000):
	# Initialize variables
	KK = len(gwas_beta)
	GG = eqtl_beta_raw.shape[0]
	# Reorganize gwas data
	gwas_beta_var = np.square(gwas_beta_se)
	# Reorganize eqtl data
	eqtl_betas = []
	eqtl_beta_vars = []
	eqtl_indices = []
	for gene_iter in range(GG):
		gene_indices = eqtl_beta_raw[gene_iter,:] != 0.0
		eqtl_betas.append(eqtl_beta_raw[gene_iter, gene_indices])
		eqtl_beta_vars.append(np.square(eqtl_beta_se_raw[gene_iter, gene_indices]))
		eqtl_indices.append(gene_indices)

	# Initialize causal effcts
	gamma = np.zeros(KK)
	alpha = np.zeros(GG)
	deltas = []
	for gene_iter in range(GG):
		deltas.append(eqtl_betas[gene_iter])
	# Initialize variance parameters
	gamma_var = 1e-6
	alpha_var = 1e-6
	delta_vars = np.ones(GG)*1e-6

	# Remove causal effects from gwas beta
	# Only works because initialize all causal effects to zero
	gwas_beta_resid = np.copy(gwas_beta)

	sampled_gamma_vars = []
	sampled_alpha_vars = []
	sampled_delta_vars = []
	sampled_nm_h2s = []
	sampled_med_h2s = []

	for itera in range(max_iters):
		# Update gamma
		gamma, gwas_beta_resid = update_gamma(gwas_beta_resid, gwas_beta_var, gamma, gamma_var)

		# Update alphas
		alpha, gwas_beta_resid = update_alpha(gwas_beta_resid, gwas_beta_var, alpha, alpha_var, deltas, eqtl_indices)

		# Update deltas
		deltas, gwas_beta_resid = update_deltas(gwas_beta_resid, gwas_beta_var, eqtl_betas, eqtl_beta_vars, alpha, deltas, delta_vars, eqtl_indices)

		# Update gamma_var
		gamma_var = update_gamma_var(gamma)

		# Update alpha var
		alpha_var = update_gamma_var(alpha)

		# Update delta var
		for gg in range(GG):
			delta_vars[gg] = update_gamma_var(deltas[gg])

		if np.mod(itera, 1000) == 0.0:
			print(itera)

		# Update delta var
		if itera > burn_in_iter:
			sampled_gamma_vars.append(gamma_var)
			sampled_alpha_vars.append(alpha_var)
			sampled_delta_vars.append(delta_vars)

			med_h2 = 0.0
			for gg in range(GG):
				med_h2 = med_h2 + np.sum(np.square(deltas[gg]*alpha[gg]))
			nm_h2 = np.sum(np.square(gamma))

			sampled_nm_h2s.append(nm_h2)
			sampled_med_h2s.append(med_h2)


	sampled_gamma_vars = np.asarray(sampled_gamma_vars)
	sampled_alpha_vars = np.asarray(sampled_alpha_vars)
	sampled_med_h2s = np.asarray(sampled_med_h2s)
	sampled_nm_h2s = np.asarray(sampled_nm_h2s)
	sampled_delta_vars = np.asarray(sampled_delta_vars)


	return np.mean(sampled_med_h2s), np.mean(sampled_nm_h2s), np.mean(sampled_delta_vars,axis=0)*200




########################
# Command line args
########################
global_sim_iter = int(sys.argv[1])
n_sims = int(sys.argv[2])
gwas_ss = int(sys.argv[3])
n_snps = int(sys.argv[4])
n_genes = int(sys.argv[5])
snps_per_gene = int(sys.argv[6])
med_h2 = float(sys.argv[7])
nm_h2 = float(sys.argv[8])
eqtl_ss = int(sys.argv[9])
eqtl_architecture = sys.argv[10]
mean_cis_h2 = float(sys.argv[11])
frac_causal_genes = float(sys.argv[12])
output_stem = sys.argv[13]


# Set seed
np.random.seed(global_sim_iter)

# Options for ge_h2
ge_h2_options = [mean_cis_h2/2.0, mean_cis_h2, mean_cis_h2*2.0]


t = open(output_stem + '_est_summary.txt','w')
t.write('sim_iter\tmethod\tsim_nm_h2\tsim_med_h2\tsim_mean_gene_cis_h2\test_nm_h2\test_med_h2\test_mean_gene_cis_h2\n')


# Loop through independent simulations
for sim_iter in range(n_sims):
	# Simulate data
	Y, Es, Y_geno, E_geno, gene_causal_eqtl_effects, nm_var_causal_effects, gene_trait_effects, gene_snp_indices_arr, sim_ge_h2s, gene_window_based_ld_score = simulate_data(gwas_ss, n_snps, eqtl_ss, med_h2, ge_h2_options, n_genes, snps_per_gene, nm_h2, eqtl_architecture, frac_causal_genes)


	# Next get summary statistics
	gwas_beta, gwas_beta_se = get_marginal_summary_statistics(Y, Y_geno)
	eqtl_beta, eqtl_beta_se = get_marginal_summary_statistics_across_genes(Es, E_geno, gene_snp_indices_arr)

	# Temp loading of data
	'''
	np.save('gwas_beta.npy',gwas_beta)
	np.save('gwas_beta_se.npy',gwas_beta_se)
	np.save('eqtl_beta.npy',eqtl_beta)
	np.save('eqtl_beta_se.npy',eqtl_beta_se)

	gwas_beta = np.load('gwas_beta.npy')
	gwas_beta_se = np.load('gwas_beta_se.npy')
	eqtl_beta = np.load('eqtl_beta.npy')
	eqtl_beta_se = np.load('eqtl_beta_se.npy')
	'''

	sim_med_h2 = np.sum(np.square(np.dot(np.transpose(gene_causal_eqtl_effects), gene_trait_effects)))
	sim_nm_h2 = np.sum(np.square(nm_var_causal_effects))
	sim_eqtl_h2 = np.sum(np.square(gene_causal_eqtl_effects),axis=1)


	# Goes to nans
	#med_h2_sumstat_reml_per_snps, nm_h2_sumstat_reml_per_snp, eqtl_h2_sumstat_reml_per_snp = med_h2_with_sumstat_reml_per_snp_no_ld(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se)
	#med_h2_sumstat_reml_chunked, nm_h2_sumstat_reml_chunked, eqtl_h2_sumstat_reml_chunked = med_h2_with_sumstat_reml_chunked_no_ld(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se)

	# Bayesian approach
	med_h2_sumstat_bayesian, nm_h2_sumstat_bayesian, eqtl_h2_sumstat_bayesian = med_h2_with_sumstat_bayesian_no_ld(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se)


	t.write(str(sim_iter) + '\t' + 'sumstat_bayesian\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(np.mean(sim_eqtl_h2)) + '\t' + str(nm_h2_sumstat_bayesian) + '\t' + str(med_h2_sumstat_bayesian) + '\t' +  str(np.mean(eqtl_h2_sumstat_bayesian)) + '\n')

	'''

	med_h2_sumstat_reml_two_step, nm_h2_sumstat_reml_two_step, eqtl_h2_sumstat_reml_two_step = med_h2_with_sumstat_reml_no_ld_two_step(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se)
	med_h2_sumstat_reml, nm_h2_sumstat_reml, eqtl_h2_sumstat_reml = med_h2_with_sumstat_reml_no_ld(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se)

	t.write(str(sim_iter) + '\t' + 'sumstat_reml\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(np.mean(sim_eqtl_h2)) + '\t' + str(nm_h2_sumstat_reml) + '\t' + str(med_h2_sumstat_reml) + '\t' +  str(np.mean(eqtl_h2_sumstat_reml)) + '\n')
	t.write(str(sim_iter) + '\t' + 'sumstat_reml_two_step\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(np.mean(sim_eqtl_h2)) + '\t' + str(nm_h2_sumstat_reml_two_step) + '\t' + str(med_h2_sumstat_reml_two_step) + '\t' +  str(np.mean(eqtl_h2_sumstat_reml_two_step)) + '\n')

	t.flush()
	'''

t.close()



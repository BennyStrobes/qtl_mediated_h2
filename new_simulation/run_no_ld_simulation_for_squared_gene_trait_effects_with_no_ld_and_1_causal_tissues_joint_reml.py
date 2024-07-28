import sys
import numpy as np 
import pandas as pd
import os
import pdb
import tensorflow as tf
import gzip
import time
#from scipy.stats import gamma
#os.environ['OPENBLAS_NUM_THREADS'] = '20'
#os.environ['MKL_NUM_THREADS'] = '20'
#from joblib.externals.loky import set_loky_pickler
#from joblib import wrap_non_picklable_objects
#from tensorflow.keras.layers import Conv1D, Dense, BatchNormalization, Activation, GlobalAveragePooling1D, Add, Input, Flatten
#from tensorflow.keras import Model
#from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
#from scipy import odr
import tensorflow_probability as tfp
#from tensorflow_probability.python.math import bessel
#from tensorflow_probability.python.distributions import noncentral_chi2
#from scipy import stats
#tfd = tfp.distributions
#from tensorflow_probability.python.math import generic as tfp_math


def simulate_data(gwas_ss, n_snps, eqtl_ss_1, med_h2_1, ge_h2_options, n_genes, snps_per_gene, nm_h2, eqtl_architecture):
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
	sim_alphas = np.random.normal(loc=0, scale=np.sqrt(med_h2_1/n_genes), size=n_genes)
	genetic_trait = np.zeros(gwas_ss)
	for gene_iter in range(n_genes):
		genetic_gene = np.dot(gwas_geno, gene_causal_eqtl_effects[gene_iter, :])
		standardized_genetic_gene = genetic_gene/np.std(genetic_gene)
		genetic_trait =genetic_trait + standardized_genetic_gene*sim_alphas[gene_iter]
	# Causal nm var-trait effect
	n_causal_snps = 3000
	sim_beta = np.zeros(n_snps)
	causal_indices = np.random.choice(np.arange(n_snps), size=n_causal_snps, replace=False)
	sim_beta[causal_indices] = np.random.normal(loc=0, scale=np.sqrt(nm_h2/n_causal_snps), size=n_causal_snps)
	genetic_trait = genetic_trait + np.dot(gwas_geno, sim_beta)

	Y = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0-np.var(genetic_trait)))
	Y = (Y - np.mean(Y))/np.std(Y)

	return Y, Es, gwas_geno, eqtl_geno_1, gene_causal_eqtl_effects, sim_beta, sim_alphas, gene_snp_indices_arr, np.asarray(sim_ge_h2s), gene_counts


def get_marginal_summary_statistics_no_intercept(Y, G):
	n_snps = G.shape[1]
	beta = []
	beta_se = []
	for snp_iter in range(n_snps):
		olser = sm.OLS(Y, G[:,snp_iter]).fit()
		beta.append(olser.params[0])
		beta_se.append(olser.bse[0])
	return np.asarray(beta), np.asarray(beta_se)

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


def posterior_distribution_on_causal_eqtl_effects_closed_form(expression_vec, genotype_mat, ge_h2, resid_var):
	n_snps = genotype_mat.shape[1]
	#resid_var = 1.0 - ge_h2
	per_snp_h2 = ge_h2/n_snps 

	tmp_prec = (np.eye(n_snps)/per_snp_h2) + (np.dot(np.transpose(genotype_mat), genotype_mat)/(resid_var))
	posterior_var = np.linalg.inv(tmp_prec)

	posterior_mean = (1.0/resid_var)*np.dot(np.dot(posterior_var, np.transpose(genotype_mat)), expression_vec)

	return posterior_mean, posterior_var

def compute_univariate_gaussian_log_like(xx, mean_value, variance_vec):
	log_like = -(tf.math.log(variance_vec)/2) - tf.math.divide(tf.square(xx-mean_value), (2.0*variance_vec))
	
	return log_like


def sumstat_reml_no_ld_loss(gwas_beta, gwas_beta_var, h2_variable, num_snps):
	log_like = compute_univariate_gaussian_log_like(gwas_beta, 0.0, gwas_beta_var + (h2_variable/num_snps))

	loss = -tf.reduce_sum(log_like)
	
	return loss


def snp_h2_with_sumstat_reml_no_ld(gwas_beta, gwas_beta_se, max_epochs=10000, conv_thresh=1e-15):
	num_snps = len(gwas_beta)


	# Initialize variables to optimize over
	h2_variable = tf.Variable(initial_value=.00001,trainable=True, name='beta_sq')

	# Set up optimization function
	optimizer = tf.keras.optimizers.Adam()

	# Convert summary stats to tf tensors
	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)



	# Lopp through windows
	prev_est_h2 = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = sumstat_reml_no_ld_loss(gwas_beta, gwas_beta_var, h2_variable, num_snps)

		trainable_variables = []
		trainable_variables.append(h2_variable)


		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est_h2 = np.asarray(h2_variable)*1.0

		diff = np.abs(prev_est_h2 -cur_est_h2)
		if diff < conv_thresh:
			converged = True
			break

		prev_est_h2 = cur_est_h2

	if converged == False:
		print('did not converge')
		print(diff)

	return cur_est_h2









	snp_h2 = 0.0
	return snp_h2






#########################
# Command line args
#########################
global_simulation_number = sys.argv[1]
n_sims = int(sys.argv[2])
gwas_ss = int(sys.argv[3])
n_snps = int(sys.argv[4])
eqtl_ss_1 = int(sys.argv[5])
med_h2_1 = float(sys.argv[6])
output_root = sys.argv[7]
global_simulation_number = int(sys.argv[8])
n_genes = int(sys.argv[9])
snps_per_gene = int(sys.argv[10])
nm_h2 = float(sys.argv[11])
eqtl_architecture = sys.argv[12] # currently implemented for sparse or polygenic


# Set seed
np.random.seed(global_simulation_number)


# Open and print header to output file
output_file = output_root + '_effect_est_res_summary.txt'

# Options for ge_h2
ge_h2_options = [.05,.1, .2]

# Loop through sims
for sim_iter in range(n_sims):
	print(sim_iter)

	# First simulate data
	Y, Es, Y_geno, E_geno, gene_causal_eqtl_effects, nm_var_causal_effects, gene_trait_effects, gene_snp_indices_arr, sim_ge_h2s, gene_window_based_ld_score = simulate_data(gwas_ss, n_snps, eqtl_ss_1, med_h2_1, ge_h2_options, n_genes, snps_per_gene, nm_h2, eqtl_architecture)

	# Sim h2s
	sim_nm_var_h2 = np.sum(np.square(nm_var_causal_effects))
	sim_med_h2 = np.sum(np.square(gene_trait_effects))

	# Next get summary statistics
	gwas_beta, gwas_beta_se = get_marginal_summary_statistics(Y, Y_geno)
	#eqtl_beta, eqtl_beta_se = get_marginal_summary_statistics_across_genes(Es, E_geno, gene_snp_indices_arr)
	# Compute gene level z-scores
	#gwas_z = gwas_beta/gwas_beta_se

	# Temp loading of data
	#np.save('gwas_beta.npy',gwas_beta)
	#np.save('gwas_beta_se.npy',gwas_beta_se)



	# Part 1: evaluate gwas snp h2 using sumstat reml
	gwas_snp_h2_sumstat_reml = snp_h2_with_sumstat_reml_no_ld(gwas_beta, gwas_beta_se)
	print(gwas_snp_h2_sumstat_reml)



t.close()
t_gene.close()


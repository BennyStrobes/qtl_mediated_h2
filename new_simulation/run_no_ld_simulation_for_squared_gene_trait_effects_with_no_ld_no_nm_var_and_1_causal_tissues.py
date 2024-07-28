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
from tensorflow.keras.layers import Conv1D, Dense, BatchNormalization, Activation, GlobalAveragePooling1D, Add, Input, Flatten
from tensorflow.keras import Model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from scipy import odr
import tensorflow_probability as tfp
from tensorflow_probability.python.math import bessel
from tensorflow_probability.python.distributions import noncentral_chi2
from scipy import stats
tfd = tfp.distributions
from tensorflow_probability.python.math import generic as tfp_math



def simulate_data(gwas_ss, n_snps, eqtl_ss_1, med_h2_1, ge_h2_1, n_genes, snps_per_gene, nm_h2):
	sim_alpha_1_sq = med_h2_1/ge_h2_1

	gene_midpoint_lb = int(snps_per_gene/2) +1
	gene_midpoint_ub = n_snps - (int(snps_per_gene/2) +1)

	gene_causal_eqtl_effects = []
	gene_snp_indices_arr = []
	for gene_iter in range(n_genes):

		gene_midpoint = np.random.choice(np.arange(gene_midpoint_lb, gene_midpoint_ub))
		gene_start = gene_midpoint - int(snps_per_gene/2)
		gene_end = gene_midpoint + int(snps_per_gene/2)
		gene_snp_indices = np.arange(gene_start, gene_end)

		gene_snp_causal_eqtl_effects = np.zeros(n_snps)
		gene_snp_causal_eqtl_effects[gene_snp_indices] = np.random.normal(loc=0, scale=np.sqrt(ge_h2_1/snps_per_gene), size=snps_per_gene)

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
	sim_alphas = np.random.normal(loc=0, scale=np.sqrt(sim_alpha_1_sq/n_genes), size=n_genes)
	genetic_trait = np.zeros(gwas_ss)
	for gene_iter in range(n_genes):
		genetic_gene = np.dot(gwas_geno, gene_causal_eqtl_effects[gene_iter, :])
		genetic_trait =genetic_trait + genetic_gene*sim_alphas[gene_iter]
	print(np.var(genetic_trait))
	# Causal nm var-trait effect
	sim_beta = np.random.normal(loc=0, scale=np.sqrt(nm_h2/n_snps), size=n_snps)
	genetic_trait = genetic_trait + np.dot(gwas_geno, sim_beta)
	print(np.var(genetic_trait))

	Y = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0-np.var(genetic_trait)))
	Y = (Y - np.mean(Y))/np.std(Y)

	return sim_alpha_1_sq, Y, Es, gwas_geno, eqtl_geno_1, gene_causal_eqtl_effects, sim_beta, sim_alphas, gene_snp_indices_arr


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





def mediated_effect_est_loss_function(gwas_beta, eqtl_beta, gwas_ss, eqtl_sss, eqtl_beta_deltas_var, alpha_var):
	loss1 = tf.math.square(gwas_beta - tf.linalg.matmul(eqtl_beta + eqtl_beta_deltas_var, alpha_var))*gwas_ss
	loss2 = tf.linalg.matmul(tf.math.square(eqtl_beta_deltas_var), eqtl_sss)

	agg_loss = loss1 + loss2

	return tf.math.reduce_sum(agg_loss)



def mediated_squared_effect_est_loss_function2(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var, alpha_sq_var, eqtl_beta_sq_vars, gene_snp_indices_arr):
	#gwas_pred_chi_sq = (tf.math.reduce_sum(eqtl_beta_sq_var*alpha_sq_var,axis=0)) + gwas_beta_var
	gwas_pred_chi_sq = gwas_beta_var


	log_like2 = 0

	for gene_iter, gene_snp_indices in enumerate(gene_snp_indices_arr):
		gwas_pred_chi_sq = tf.tensor_scatter_nd_update(gwas_pred_chi_sq, tf.expand_dims(gene_snp_indices,axis=1), tf.gather(gwas_pred_chi_sq, indices=(gene_snp_indices)) + eqtl_beta_sq_vars[gene_iter]*alpha_sq_var)


		pred_chi_sq2 = eqtl_beta_sq_vars[gene_iter] + tf.gather(eqtl_beta_var[gene_iter,:], indices=gene_snp_indices)
		tmp_log_like = -tf.math.divide(tf.math.square(tf.gather(eqtl_beta[gene_iter,:],indices=gene_snp_indices)), 2.0*pred_chi_sq2) - (.5*tf.math.log(2.0*pred_chi_sq2))
		log_like2 = log_like2 + tf.math.reduce_sum(tmp_log_like)


	log_like1 = -tf.math.divide(tf.math.square(gwas_beta), 2.0*gwas_pred_chi_sq) - (.5*tf.math.log(2.0*gwas_pred_chi_sq))

	
	agg_loss = -tf.math.reduce_sum(log_like1) - log_like2

	return agg_loss



def mediated_squared_effect_est_loss_function(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var, alpha_sq_var, eqtl_beta_deltas_var):
	pred_chi_sq = (tf.math.reduce_sum(tf.math.square(eqtl_beta + eqtl_beta_deltas_var)*alpha_sq_var,axis=0)) + gwas_beta_var
	log_like = (-.5)*tf.math.log(tf.math.square(gwas_beta)) - tf.math.divide(tf.math.square(gwas_beta), 2.0*pred_chi_sq) - (.5*tf.math.log(2.0*pred_chi_sq))
	
	loss2 = tf.math.reduce_sum(tf.math.square(eqtl_beta_deltas_var)/eqtl_beta_var,axis=0)/2.0
	agg_loss = -log_like + loss2

	return tf.math.reduce_sum(agg_loss)


def mediated_squared_effect_est_loss_function3(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_var, eqtl_beta_sq_var):

	pred_chi_sq = alpha_sq_var*(tf.math.reduce_sum(tf.math.square(eqtl_beta_sq_var)*eqtl_mask,axis=0)) + gwas_beta_var
	log_like = -tf.math.divide(tf.math.square(gwas_beta), 2.0*pred_chi_sq) - (.5*tf.math.log(2.0*pred_chi_sq))
	

	pred_chi_sq2 = (tf.math.square(eqtl_beta_sq_var) + eqtl_beta_var)
	log_like2 = -tf.math.divide(tf.math.square(eqtl_beta), 2.0*pred_chi_sq2) - (.5*tf.math.log(2.0*pred_chi_sq2))

	agg_loss = -tf.reduce_sum(log_like) - tf.reduce_sum(log_like2*eqtl_mask)

	return agg_loss


def mediated_squared_effect_est_loss_function_unobserved_true(gwas_beta, gwas_beta_var, eqtl_beta, alpha_sq_var):
	pred_chi_sq = (tf.math.reduce_sum(tf.math.square(eqtl_beta)*alpha_sq_var,axis=0)) + gwas_beta_var
	log_like = (-.5)*tf.math.log(tf.math.square(gwas_beta)) - tf.math.divide(tf.math.square(gwas_beta), 2.0*pred_chi_sq) - (.5*tf.math.log(2.0*pred_chi_sq))

	return -tf.math.reduce_sum(log_like)


def custom_tf_odr_mediated_squared_effect_est2(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se,gene_snp_indices_arr, max_epochs=1000, conv_thresh=1e-15):
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Initialize variables to optimize over
	alpha_sq_var = tf.Variable(initial_value=0.0,trainable=True, name='alpha_sq')
	#eqtl_beta_sq_var = tf.Variable(initial_value=np.square(eqtl_beta).astype(np.float32),trainable=True, name='eqtl_beta_sqs')
	eqtl_beta_sq_vars = []
	for gene_iter in range(n_genes):
		eqtl_beta_sq_var = tf.Variable(initial_value=np.square(eqtl_beta[gene_iter, gene_snp_indices_arr[gene_iter]]).astype(np.float32),trainable=True, name='eqtl_beta_sqs_' + str(gene_iter))
		eqtl_beta_sq_vars.append(eqtl_beta_sq_var)

	optimizer = tf.keras.optimizers.Adam()

	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)

	# Lopp through windows
	prev_est_alpha_sq = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = mediated_squared_effect_est_loss_function2(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var, alpha_sq_var, eqtl_beta_sq_vars, gene_snp_indices_arr)

		trainable_variables = []
		for gene_iter in range(n_genes):
			trainable_variables.append(eqtl_beta_sq_vars[gene_iter])
		trainable_variables.append(alpha_sq_var)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))

		print(alpha_sq_var)

	pdb.set_trace()

	return np.asarray(alpha_sq_var)[0,0]


def custom_tf_odr_mediated_squared_effect_est6(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, gene_snp_indices_arr, gwas_ss, eqtl_ss, gene_causal_eqtl_effects, nm_var_causal_effects, gene_trait_effects, max_epochs=800, conv_thresh=1e-8):
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)


	# Get variances using (unobserved) simulated data
	squared_trait_causal_effects = np.square(nm_var_causal_effects) + np.dot(np.transpose(np.square(gene_causal_eqtl_effects)), np.square(gene_trait_effects))
	#gwas_variance = 2.0 + (4.0*squared_trait_causal_effects/np.square(gwas_beta_se))
	#eqtl_variance = 2.0 + (4.0*np.square(gene_causal_eqtl_effects)/np.square(eqtl_beta_se))

	# Initialize variables to optimize over
	beta_sq_var = tf.Variable(initial_value=.00001,trainable=True, name='beta_sq')
	alpha_sq_var = tf.Variable(initial_value=.00001,trainable=True, name='alpha_sq')
	eqtl_beta_sq_var = tf.Variable(initial_value=np.square(eqtl_beta),trainable=True, name='eqtl_beta_sq', dtype=tf.float32)
	gwas_precision = tf.Variable(initial_value=1.0,trainable=True, name='gwas_precision', dtype=tf.float32)
	eqtl_precisions = tf.Variable(initial_value=np.ones(n_genes),trainable=True, name='eqtl_precisions',dtype=tf.float32)

	optimizer = tf.keras.optimizers.Adam()

	# Create mask matrix
	eqtl_mask = np.zeros(eqtl_beta.shape)
	for gene_iter, gene_snp_indices in enumerate(gene_snp_indices_arr):
		eqtl_mask[gene_iter, gene_snp_indices] = eqtl_mask[gene_iter, gene_snp_indices] + 1.0 

	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)
	gwas_ss = tf.convert_to_tensor(gwas_ss, dtype=tf.float32)
	eqtl_ss = tf.convert_to_tensor(eqtl_ss, dtype=tf.float32)

	#gwas_variance = tf.convert_to_tensor(gwas_variance, dtype=tf.float32)
	#eqtl_variance = tf.convert_to_tensor(eqtl_variance, dtype=tf.float32)



	# Lopp through windows
	prev_est_alpha_sq = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = mediated_squared_effect_est_loss_function6(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_var, eqtl_beta_sq_var, gwas_ss, eqtl_ss, beta_sq_var, gwas_precision, eqtl_precisions, n_snps)

		trainable_variables = []
		trainable_variables.append(alpha_sq_var)
		trainable_variables.append(eqtl_beta_sq_var)
		trainable_variables.append(beta_sq_var)
		#trainable_variables.append(gwas_precision)
		#trainable_variables.append(eqtl_precisions)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est = np.asarray(alpha_sq_var)*1.0

		diff = np.abs(prev_est_alpha_sq -cur_est)
		if diff < conv_thresh:
			converged = True
			#break

		prev_est_alpha_sq = cur_est

		print(cur_est)
		print(loss_value)

	if converged == False:
		print('did not converge')
		print(diff)

	print(np.mean(np.sum(eqtl_beta_sq_var,axis=1)))
	print(np.mean(np.sum(eqtl_beta_sq_var,axis=1))*cur_est)

	print(loss_value)
	pdb.set_trace()


	return cur_est



def custom_tf_odr_mediated_squared_effect_est5(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, gene_snp_indices_arr, max_epochs=20000, conv_thresh=1e-8):
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Initialize variables to optimize over
	beta_sq_var = tf.Variable(initial_value=0.0000000001,trainable=True, name='beta_sq')
	alpha_sq_var = tf.Variable(initial_value=.0000001,trainable=True, name='alpha_sq')
	eqtl_beta_sq_var = tf.Variable(initial_value=np.square(eqtl_beta),trainable=True, name='eqtl_beta_sq', dtype=tf.float32)

	optimizer = tf.keras.optimizers.Adam()

	# Create mask matrix
	eqtl_mask = np.zeros(eqtl_beta.shape)
	for gene_iter, gene_snp_indices in enumerate(gene_snp_indices_arr):
		eqtl_mask[gene_iter, gene_snp_indices] = eqtl_mask[gene_iter, gene_snp_indices] + 1.0 

	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)



	# Lopp through windows
	prev_est_alpha_sq = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = mediated_squared_effect_est_loss_function5(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_var, eqtl_beta_sq_var, beta_sq_var)

		trainable_variables = []
		trainable_variables.append(alpha_sq_var)
		trainable_variables.append(eqtl_beta_sq_var)
		trainable_variables.append(beta_sq_var)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est = np.asarray(alpha_sq_var)*1.0

		diff = np.abs(prev_est_alpha_sq -cur_est)
		if diff < conv_thresh:
			converged = True
			break

		prev_est_alpha_sq = cur_est

		pdb.set_trace(0)
		print(cur_est)
		print(beta_sq_var)
		print(loss_value)

	if converged == False:
		print('did not converge')
		print(diff)

	pdb.set_trace()


	return cur_est



def custom_tf_odr_mediated_squared_effect_est4(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, gene_snp_indices_arr, max_epochs=20000, conv_thresh=1e-15):
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Initialize variables to optimize over
	beta_sq_var = tf.Variable(initial_value=-12.0,trainable=True, name='beta_sq')
	alpha_sq_var = tf.Variable(initial_value=-7.0,trainable=True, name='alpha_sq')
	eqtl_beta_sq_var = tf.Variable(initial_value=np.asarray(tfp.math.softplus_inverse(np.square(eqtl_beta)/3.0)),trainable=True, name='eqtl_beta_sq', dtype=tf.float32)

	optimizer = tf.keras.optimizers.Adam()

	# Create mask matrix
	eqtl_mask = np.zeros(eqtl_beta.shape)
	for gene_iter, gene_snp_indices in enumerate(gene_snp_indices_arr):
		eqtl_mask[gene_iter, gene_snp_indices] = eqtl_mask[gene_iter, gene_snp_indices] + 1.0 

	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)



	# Lopp through windows
	prev_est_alpha_sq = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = mediated_squared_effect_est_loss_function4(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_var, eqtl_beta_sq_var, beta_sq_var)

		trainable_variables = []
		trainable_variables.append(alpha_sq_var)
		trainable_variables.append(eqtl_beta_sq_var)
		trainable_variables.append(beta_sq_var)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))

		cur_est = np.asarray(tf.math.softplus(alpha_sq_var))*1.0
		nm_est = np.asarray(tf.math.softplus(beta_sq_var))*1.0

		diff = np.abs(prev_est_alpha_sq -cur_est)
		if diff < conv_thresh:
			converged = True
			break

		prev_est_alpha_sq = cur_est

		print(cur_est)
		print(np.sum(np.asarray(tf.math.softplus(eqtl_beta_sq_var)))/n_genes)
		print(nm_est)
		print(loss_value)

	if converged == False:
		print('did not converge')
		print(diff)

	pdb.set_trace()

	return cur_est


def custom_tf_odr_mediated_squared_effect_est3(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, gene_snp_indices_arr, max_epochs=5000, conv_thresh=1e-15):
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Initialize variables to optimize over
	alpha_sq_var = tf.Variable(initial_value=0.0,trainable=True, name='alpha_sq')
	eqtl_beta_sq_var = tf.Variable(initial_value=eqtl_beta.astype(np.float32),trainable=True, name='eqtl_beta_sq')

	optimizer = tf.keras.optimizers.Adam()

	# Create mask matrix
	eqtl_mask = np.zeros(eqtl_beta.shape)
	for gene_iter, gene_snp_indices in enumerate(gene_snp_indices_arr):
		eqtl_mask[gene_iter, gene_snp_indices] = eqtl_mask[gene_iter, gene_snp_indices] + 1.0 

	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)



	# Lopp through windows
	prev_est_alpha_sq = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = mediated_squared_effect_est_loss_function3(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_var, eqtl_beta_sq_var)

		trainable_variables = []
		trainable_variables.append(alpha_sq_var)
		trainable_variables.append(eqtl_beta_sq_var)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))

		print(loss_value)
		print(alpha_sq_var)

	pdb.set_trace()

	return np.asarray(alpha_sq_var)[0,0]


def custom_tf_odr_mediated_squared_effect_est(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, max_epochs=1000, conv_thresh=1e-15):
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Initialize variables to optimize over
	alpha_sq_var = tf.Variable(initial_value=np.zeros((1,1)).astype(np.float32),trainable=True, name='alpha_sq')
	eqtl_beta_deltas_var = tf.Variable(initial_value=np.zeros(eqtl_beta.shape).astype(np.float32),trainable=True, name='eqtl_beta_deltas')

	optimizer = tf.keras.optimizers.Adam()

	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)

	# Lopp through windows
	prev_est_alpha_sq = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = mediated_squared_effect_est_loss_function(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var, alpha_sq_var, eqtl_beta_deltas_var)

		trainable_variables = []
		trainable_variables.append(alpha_sq_var)
		trainable_variables.append(eqtl_beta_deltas_var)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))

		print(np.asarray(alpha_sq_var)[0,0])

	pdb.set_trace()

	return np.asarray(alpha_sq_var)[0,0]


def _temme_expansion_from_scratch(v, x, output_log_space=True):
	"""Compute modified bessel functions using Temme's method."""
	# The implementation of this is based on [1].
	# [1] N. Temme, On the Numerical Evaluation of the Modified Bessel Function
	#   of the Third Kind. Journal of Computational Physics 19, 1975.

	v_less_than_zero = v < 0.
	v = tf.math.abs(v)
	n = tf.math.round(v)
	# Use this to compute Kv(u, x) and Kv(u + 1., x)
	u = v - n
	x_abs = tf.math.abs(x)

	small_x = tf.where(x_abs <= 2., x_abs, 0.1)
	large_x = tf.where(x_abs > 2., x_abs, 1000.)
	temme_kue, temme_kuep1 = bessel._temme_series(u, small_x, output_log_space=output_log_space)
	cf_kue, cf_kuep1 = bessel._continued_fraction_kv(u, large_x, output_log_space=output_log_space)

	kue = tf.where(x_abs <= 2., temme_kue, cf_kue)
	kuep1 = tf.where(x_abs <= 2., temme_kuep1, cf_kuep1)

	# Now use the forward recurrence for modified bessel functions
	# to compute Kv(v, x). That is,
	# K_{v + 1}(z) - (2v / z) K_v(z) - K_{v - 1}(z) = 0.
	# This is known to be forward numerically stable.
	# Note: This recurrence is also satisfied by K_v(z) * exp(z)

	def bessel_recurrence(index, kve, kvep1):
		next_kvep1 = tfp_math.log_add_exp(kvep1 + tf.math.log(u + index) + np.log(2.) - tf.math.log(x_abs), kve)

		kve = tf.where(index > n, kve, kvep1)
		kvep1 = tf.where(index > n, kvep1, next_kvep1)
		return index + 1., kve, kvep1

	_, kve, kvep1 = tf.while_loop(cond=lambda i, *_: tf.reduce_any(i <= n), body=bessel_recurrence, loop_vars=(1.0, kue, kuep1))

	# Finally, it is known that the Wronskian
	# det(I_v * K'_v - K_v * I'_v) = - 1. / x. We can
	# use this to evaluate I_v by taking advantage of identities of Bessel
	# derivatives.

	ive = -tf.math.log(x_abs) - tfp_math.log_add_exp(kve + tf.math.log(bessel.bessel_iv_ratio(v + 1., x)), kvep1)
	# We need to add a correction term for negative v.
	log_ive = ive
	negative_v_correction = kve - 2. * x_abs

	coeff = 2 / np.pi * tf.math.sin(np.pi * u)
	coeff = (1. - 2. * tf.math.mod(n, 2.)) * coeff

	lse, sign = tfp_math.log_sub_exp(log_ive,negative_v_correction + tf.math.log(tf.math.abs(coeff)),return_sign=True)

	sign = tf.where(coeff < 0., sign, 1.)

	log_ive_negative_v = tf.where(coeff < 0.,lse,tfp_math.log_add_exp(log_ive, negative_v_correction + tf.math.log(tf.math.abs(coeff))))

	z = u + tf.math.mod(n, 2.)

	ive = tf.where(v_less_than_zero, log_ive_negative_v, ive)

	ive = tf.where(tf.math.equal(x, 0.),tf.where(tf.math.equal(v, 0.), 0.,-np.inf), ive)

	ive = tf.where(tf.math.equal(x, 0.) & v_less_than_zero,tf.where(tf.math.equal(z, tf.math.floor(z)),ive, np.inf), ive)

	kve = tf.where(tf.math.equal(x, 0.), np.inf, kve)
	ive = tf.where(x < 0., np.nan, ive)
	kve = tf.where(x < 0., np.nan, kve)

	return ive, kve





def bessel_from_scratch(v, z):
	# Taken from https://github.com/tensorflow/probability/blob/v0.23.0/tensorflow_probability/python/math/bessel.py#L992-L1018 (line 879)

	z_abs = tf.math.abs(z)
	# Handle the zero case specially.
	z_abs = tf.where(tf.math.equal(z_abs, 0.), 1.0, z_abs)


	#ive = _temme_expansion_from_scratch(v, z_abs, output_log_space=True)[0]
	ive = bessel._temme_expansion(v, z_abs, output_log_space=True)[0]

	# Handle when z is zero.
	ive = tf.where(tf.math.equal(z, 0.), tf.where(tf.math.equal(v, 0.), 0.0, tf.where(v < 0., np.inf,-np.inf)), ive)

	return ive




def NCX2_log_prob(xx, df, nc):
	'''
	log_pdf_check = stats.ncx2.logpdf(xx, df=df, nc=nc)
	'''
	'''
	# Taken from https://github.com/tensorflow/probability/blob/v0.23.0/tensorflow_probability/python/distributions/noncentral_chi2.py#L424-L666 (line 44)
	df_factor = 0.25 * df - 0.5
	log_pdf = -tf.math.log(2.0)
	#log_pdf = log_pdf - ((xx+nc)/2.0)  # This is what it should be
	log_pdf = log_pdf  - (0.5 * tf.math.square(tf.math.sqrt(xx) - tf.math.sqrt(nc))) # This is used by TF and scipy. The additional 
	log_pdf = log_pdf + tf.math.xlogy(df_factor, xx)
	log_pdf = log_pdf - tf.math.xlogy(df_factor, nc)
	log_pdf = log_pdf + bessel.log_bessel_ive(2. * df_factor, tf.math.sqrt(nc*xx))
	'''


	#ref_bessel = bessel.log_bessel_ive(2. * df_factor, tf.math.sqrt(nc*xx))
	#new_bessel = bessel_from_scratch(2. * df_factor, tf.math.sqrt(nc*xx))


	nc_chi2 = noncentral_chi2.NoncentralChi2(df, nc)
	log_pdf = nc_chi2.log_prob(xx)

	return log_pdf

def NCX2_snp_h2_loss(gwas_beta_sq, gwas_beta_var, beta_sq_var):
	# Prep data for NCX2
	chi_sq_stat = gwas_beta_sq/gwas_beta_var
	df = 1.0
	nc = tf.math.divide(tf.math.softplus(beta_sq_var), gwas_beta_var)/len(gwas_beta_sq)

	ncx2_log_probz = NCX2_log_prob(chi_sq_stat, df, nc)


	return -tf.math.reduce_sum(ncx2_log_probz)


def mediated_squared_effect_est_loss_function5(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_var, eqtl_beta_sq_var, beta_sq_var):

	df = 1.0

	gwas_chi_sq_stat = tf.square(gwas_beta)/gwas_beta_var
	eqtl_chi_sq_stat = tf.square(eqtl_beta)/eqtl_beta_var

	gwas_nc = tf.math.divide(tf.nn.relu(alpha_sq_var), gwas_beta_var)*(tf.math.reduce_sum(tf.nn.relu(eqtl_beta_sq_var)*eqtl_mask,axis=0)) + tf.math.divide(tf.nn.relu(beta_sq_var),gwas_beta_var)

	gwas_ncx2_log_probz = NCX2_log_prob(gwas_chi_sq_stat, df, gwas_nc)

	eqtl_nc = tf.math.divide(tf.nn.relu(eqtl_beta_sq_var), eqtl_beta_var)


	eqtl_ncx2_log_probz = NCX2_log_prob(eqtl_chi_sq_stat[eqtl_mask==1], df, eqtl_nc[eqtl_mask==1])

	agg_loss = -tf.reduce_sum(gwas_ncx2_log_probz) - tf.reduce_sum(eqtl_ncx2_log_probz)

	return agg_loss


def mediated_squared_effect_est_loss_function6(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_var, eqtl_beta_sq_var, gwas_ss, eqtl_ss, beta_sq_var, gwas_precision, eqtl_precisions, n_snps):

	big_eqtl_precisions = tf.tile(tf.reshape(eqtl_precisions, [-1, 1]), [1,n_snps])

	gwas_chi_sq_stat = tf.square(gwas_beta)/gwas_beta_var
	eqtl_chi_sq_stat = tf.square(eqtl_beta)/eqtl_beta_var


	gwas_pred = 1.0 + tf.math.divide(beta_sq_var, gwas_beta_var) + tf.math.divide(alpha_sq_var, gwas_beta_var)*(tf.math.reduce_sum(eqtl_beta_sq_var*eqtl_mask,axis=0))
	#gwas_pred_variance = 2.0 + 4.0*(tf.math.divide(beta_sq_var, gwas_beta_var) + tf.math.divide(alpha_sq_var, gwas_beta_var)*(tf.math.reduce_sum(eqtl_beta_sq_var*eqtl_mask,axis=0)))
	gwas_variance = 2.0 + tf.math.divide((4.0*0.4/50000.0), gwas_beta_var)

	eqtl_pred = 1.0 + tf.math.divide(eqtl_beta_sq_var, eqtl_beta_var)
	eqtl_var = 2.0 + tf.math.divide((4.0*0.1/200.0), eqtl_beta_var)
	#eqtl_pred_variance = 2.0 + 4.0*tf.math.divide(eqtl_beta_sq_var, eqtl_beta_var)


	#gwas_loss = tf.math.square(gwas_chi_sq_stat-gwas_pred)*tf.math.softplus(gwas_precision)
	gwas_loss = tf.math.divide(tf.math.square(gwas_chi_sq_stat-gwas_pred), gwas_variance)


	#eqtl_loss = tf.math.square(eqtl_chi_sq_stat[eqtl_mask==1]-eqtl_pred[eqtl_mask==1])*tf.math.softplus(big_eqtl_precisions)[eqtl_mask==1]
	eqtl_loss = tf.math.divide(tf.math.square(eqtl_chi_sq_stat[eqtl_mask==1]-eqtl_pred[eqtl_mask==1]), eqtl_var[eqtl_mask==1])

	agg_loss = tf.reduce_sum(gwas_loss) + tf.reduce_sum(eqtl_loss)


	#agg_loss = agg_loss + tf.reduce_sum(tf.math.log(1.0/tf.math.softplus(eqtl_precisions))*200) + tf.math.log(1.0/tf.math.softplus(gwas_precision))*n_snps


	return agg_loss





def mediated_squared_effect_est_loss_function4(gwas_beta, gwas_beta_var, eqtl_beta, eqtl_beta_var,eqtl_mask, alpha_sq_var, eqtl_beta_sq_var, beta_sq_var):

	df = 1.0

	gwas_chi_sq_stat = tf.square(gwas_beta)/gwas_beta_var
	eqtl_chi_sq_stat = tf.square(eqtl_beta)/eqtl_beta_var


	gwas_nc = tf.math.divide(tf.math.softplus(alpha_sq_var), gwas_beta_var)*(tf.math.reduce_sum(tf.math.softplus(eqtl_beta_sq_var)*eqtl_mask,axis=0)) + tf.math.divide(tf.math.softplus(beta_sq_var), gwas_beta_var)

	gwas_ncx2_log_probz = NCX2_log_prob(gwas_chi_sq_stat, df, gwas_nc)

	eqtl_nc = tf.math.divide(tf.math.softplus(eqtl_beta_sq_var), eqtl_beta_var)


	eqtl_ncx2_log_probz = NCX2_log_prob(eqtl_chi_sq_stat[eqtl_mask==1], df, eqtl_nc[eqtl_mask==1])

	agg_loss = -tf.reduce_sum(gwas_ncx2_log_probz) - tf.reduce_sum(eqtl_ncx2_log_probz)

	return agg_loss




def estimate_snp_h2_with_non_central_chi_sq2(gwas_beta, gwas_beta_se, max_epochs=10000, conv_thresh=1e-10):
	n_snps = len(gwas_beta)

	gwas_beta_sq = tf.convert_to_tensor(np.square(gwas_beta), dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)


	beta_sq_var = tf.Variable(initial_value=-2.0,trainable=True, name='beta_sq')  # To go through softplus

	#chi_sq_rv = tfd.NoncentralChi2(df=1, noncentrality=tf.Variable(.5, name='nc'), allow_nan_stats=False)

	# NOTE: Lower values may be caused by the fact that variance is not applied to NC
	
	optimizer = tf.keras.optimizers.Adam()

	
	prev_est_h2 = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = NCX2_snp_h2_loss(gwas_beta_sq, gwas_beta_var, beta_sq_var)

		trainable_variables = []
		trainable_variables.append(beta_sq_var)


		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))
		cur_h2 = np.asarray(tf.math.softplus(beta_sq_var))*1.0
		diff = np.abs(cur_h2 - prev_est_h2)
		if diff < conv_thresh:
			converged=True
			break
		prev_est_h2 = cur_h2
		print(cur_h2)

	if converged == False:
		print('did not converge: ' + str(diff))

	return cur_h2



def nll(dist, x_train):
	"""Calculates the negative log-likelihood for a given distribution
	and a data set."""
	return -tf.reduce_mean(dist.log_prob(x_train))


def estimate_snp_h2_with_non_central_chi_sq(gwas_beta, gwas_beta_se, max_epochs=10000, conv_thresh=1e-15):
	n_snps = len(gwas_beta)

	gwas_beta_sq = tf.convert_to_tensor(np.square(gwas_beta), dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)

	chi_sq_rv = tfd.NoncentralChi2(df=1, noncentrality=tf.Variable(.5, name='nc'), allow_nan_stats=False)

	# NOTE: Lower values may be caused by the fact that variance is not applied to NC
	
	optimizer = tf.keras.optimizers.Adam()

	
	prev_est_alpha_sq = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			tape.watch(chi_sq_rv.trainable_variables)
			loss_value = nll(chi_sq_rv, gwas_beta_sq/gwas_beta_var)
			grads = tape.gradient(loss_value, chi_sq_rv.trainable_variables)


		optimizer.apply_gradients(zip(grads, chi_sq_rv.trainable_variables))
		print(loss_value)
		print(chi_sq_rv.trainable_variables)

	pdb.set_trace()

def log_lss_loss(beta_hat_sq, beta_hat_var, beta_sq_var):
	pred_chi_sq = (tf.math.softplus(beta_sq_var)/len(beta_hat_sq)) + beta_hat_var
	log_like = (-.5)*tf.math.log(beta_hat_sq) - tf.math.divide(beta_hat_sq, 2.0*pred_chi_sq) - (.5*tf.math.log(2.0*pred_chi_sq))

	return -tf.math.reduce_sum(log_like)


def estimate_snp_h2_with_log_lss(beta_hat, beta_hat_se, max_epochs=100000, conv_thresh=1e-10):
	n_snps = len(beta_hat)

	# Initialize variables to optimize over
	beta_sq_var = tf.Variable(initial_value=.2,trainable=True, name='beta_sq')

	optimizer = tf.keras.optimizers.Adam()

	beta_hat_sq = tf.convert_to_tensor(np.square(beta_hat), dtype=tf.float32)
	beta_hat_var = tf.convert_to_tensor(np.square(beta_hat_se), dtype=tf.float32)

	# Lopp through windows
	prev_est_h2 = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = log_lss_loss(beta_hat_sq, beta_hat_var, beta_sq_var)

		trainable_variables = []
		trainable_variables.append(beta_sq_var)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))

		cur_h2 = np.asarray(tf.math.softplus(beta_sq_var))*1.0
		diff = np.abs(cur_h2 - prev_est_h2)
		if diff < conv_thresh:
			converged=True
			break
		prev_est_h2 = cur_h2 

	if converged == False:
		print('did not converge: ' + str(diff))

	return cur_h2

	return np.asarray(alpha_sq_var)[0,0]







def custom_tf_odr_mediated_squared_effect_est_unobserved_true(gwas_beta, gwas_beta_se, eqtl_beta, max_epochs=1000, conv_thresh=1e-15):
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)

	# Initialize variables to optimize over
	alpha_sq_var = tf.Variable(initial_value=np.zeros((1,1)).astype(np.float32),trainable=True, name='alpha_sq')

	optimizer = tf.keras.optimizers.Adam()

	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)

	# Lopp through windows
	prev_est_alpha_sq = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = mediated_squared_effect_est_loss_function_unobserved_true(gwas_beta, gwas_beta_var, eqtl_beta, alpha_sq_var)

		trainable_variables = []
		trainable_variables.append(alpha_sq_var)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))





	return np.asarray(alpha_sq_var)[0,0]






#########################
# Command line args
#########################
global_simulation_number = sys.argv[1]
n_sims = int(sys.argv[2])
gwas_ss = int(sys.argv[3])
n_snps = int(sys.argv[4])
eqtl_ss_1 = int(sys.argv[5])
med_h2_1 = float(sys.argv[6])
ge_h2_1 = float(sys.argv[7])
output_root = sys.argv[8]
global_simulation_number = int(sys.argv[9])
n_genes = int(sys.argv[10])
snps_per_gene = int(sys.argv[11])
nm_h2 = float(sys.argv[12])


# Set seed
np.random.seed(global_simulation_number+10)


# Open and print header to output file
output_file = output_root + '_effect_est_res_summary.txt'
t = open(output_file,'w')
t.write('method\tsim\tsim_alpha_1\test_alpha_1_no_noise\test_alpha_1\tloss\n')


# Loop through sims
for sim_iter in range(n_sims):
	print(sim_iter)

	# First simulate data
	sim_alpha_1_sq, Y, Es, Y_geno, E_geno, gene_causal_eqtl_effects, nm_var_causal_effects, gene_trait_effects, gene_snp_indices_arr = simulate_data(gwas_ss, n_snps, eqtl_ss_1, med_h2_1, ge_h2_1, n_genes, snps_per_gene, nm_h2)
	print(sim_alpha_1_sq)

	# Next get summary statistics
	gwas_beta, gwas_beta_se = get_marginal_summary_statistics(Y, Y_geno)
	eqtl_beta, eqtl_beta_se = get_marginal_summary_statistics_across_genes(Es, E_geno, gene_snp_indices_arr)


	#est_squared_alphas_custom_tf_no_noise = custom_tf_odr_mediated_squared_effect_est_unobserved_true(gwas_beta, gwas_beta_se, gene_causal_eqtl_effects)

	#est_squared_alphas_custom_tf = custom_tf_odr_mediated_squared_effect_est_unobserved_true(gwas_beta, gwas_beta_se, eqtl_beta)

	# Quick testing to get gwas snp h2
	'''
	ldsc_snp_h2_est = np.mean(np.square(gwas_beta) - np.square(gwas_beta_se))*n_snps
	print(ldsc_snp_h2_est)
	ncx2_snp_h2_est = estimate_snp_h2_with_non_central_chi_sq2(gwas_beta, gwas_beta_se)
	print(ncx2_snp_h2_est)
	'''
	# Quick testing to get eqtl snp h2
	'''
	gene_index = 16
	ldsc_snp_h2_est = np.mean(np.square(eqtl_beta[gene_index,gene_snp_indices_arr[gene_index]]) - np.square(eqtl_beta_se[gene_index,gene_snp_indices_arr[gene_index]]))*n_snps
	zed = eqtl_beta[gene_index,gene_snp_indices_arr[gene_index]]/eqtl_beta_se[gene_index,gene_snp_indices_arr[gene_index]]
	ldsc_snp_h2_est2 = np.mean(np.square(zed) - 1.0)*n_snps/eqtl_ss_1
	print(ldsc_snp_h2_est)
	print(ldsc_snp_h2_est2)
	loglss_snp_h2_est = estimate_snp_h2_with_log_lss(eqtl_beta[gene_index,gene_snp_indices_arr[gene_index]], eqtl_beta_se[gene_index,gene_snp_indices_arr[gene_index]])
	print(loglss_snp_h2_est)
	ncx2_snp_h2_est = estimate_snp_h2_with_non_central_chi_sq2(eqtl_beta[gene_index,gene_snp_indices_arr[gene_index]], eqtl_beta_se[gene_index,gene_snp_indices_arr[gene_index]])
	print(ncx2_snp_h2_est)
	'''



	#################################
	# Estimate squared alphas using custom tensor flow implementation of ODR
	method = 'custom_tf_odr_squared_effect_est'
	# Optimize fxn
	#est_squared_alphas_custom_tf = custom_tf_odr_mediated_squared_effect_est(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se)
	#est_squared_alphas_custom_tf = custom_tf_odr_mediated_squared_effect_est2(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, gene_snp_indices_arr)
	#est_squared_alphas_custom_tf = custom_tf_odr_mediated_squared_effect_est3(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, gene_snp_indices_arr)
	#est_squared_alphas_custom_tf = custom_tf_odr_mediated_squared_effect_est4(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, gene_snp_indices_arr)
	#est_squared_alphas_custom_tf = custom_tf_odr_mediated_squared_effect_est5(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, gene_snp_indices_arr)
	est_squared_alphas_custom_tf = custom_tf_odr_mediated_squared_effect_est6(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, gene_snp_indices_arr, gwas_ss, eqtl_ss_1, gene_causal_eqtl_effects, nm_var_causal_effects, gene_trait_effects)

	# Print to output
	#t.write(method + '\t' + str(sim_iter) + '\t' + str(sim_alpha_1_sq) + '\t' + str(est_squared_alphas_custom_tf_no_noise) + '\t' + str(est_squared_alphas_custom_tf) + '\n')

t.close()


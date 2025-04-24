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



def simulate_data(sample_size, n_snps, sim_h2, trait_architecture):

	# Simulate GWAS genotype data
	gwas_geno = np.random.normal(loc=0,scale=1.0,size=(sample_size,n_snps))
	for snp_iter in range(n_snps):
		gwas_geno[:,snp_iter] = (gwas_geno[:,snp_iter] - np.mean(gwas_geno[:,snp_iter]))/np.std(gwas_geno[:,snp_iter])



	# Simulate causal trait effects
	if trait_architecture == 'sparse':
		n_causal_snps = int(np.floor(n_snps*.05))
	elif trait_architecture == 'polygenic':
		n_causal_snps = n_snps
	else:
		print('assumption erooror: architecture not yet implemented')
		pdb.set_trace()

	sim_beta = np.zeros(n_snps)
	causal_indices = np.random.choice(np.arange(n_snps), size=n_causal_snps, replace=False)
	sim_beta[causal_indices] = np.random.normal(loc=0, scale=np.sqrt(sim_h2/n_causal_snps), size=n_causal_snps)
	genetic_trait = np.dot(gwas_geno, sim_beta)

	Y = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0-np.var(genetic_trait)))
	Y = (Y - np.mean(Y))/np.std(Y)


	return Y, gwas_geno, sim_beta, np.var(genetic_trait)





def get_marginal_summary_statistics(Y, G):
	n_snps = G.shape[1]
	beta = []
	beta_se = []
	for snp_iter in range(n_snps):
		olser = sm.OLS(Y, sm.add_constant(G[:,snp_iter])).fit()
		beta.append(olser.params[1])
		beta_se.append(olser.bse[1])
	return np.asarray(beta), np.asarray(beta_se)



def compute_univariate_gaussian_log_like(xx, mean_value, variance_vec):
	log_like = -(tf.math.log(variance_vec)/2) - tf.math.divide(tf.square(xx-mean_value), (2.0*variance_vec))
	
	return log_like


def sumstat_reml_no_ld_loss(gwas_beta, gwas_beta_var, h2_raw_variable, num_snps, per_snp_h2_est_lb):

	h2_variable = tf.nn.relu(h2_raw_variable) + (per_snp_h2_est_lb*num_snps)

	log_like = compute_univariate_gaussian_log_like(gwas_beta, 0.0, gwas_beta_var + (h2_variable/num_snps))

	loss = -tf.reduce_sum(log_like)
	
	return loss


def snp_h2_with_sumstat_reml_no_ld(gwas_beta, gwas_beta_se, max_epochs=50000, conv_thresh=1e-15):
	num_snps = len(gwas_beta)


	# Set up optimization function
	optimizer = tf.keras.optimizers.Adam()

	# Convert summary stats to tf tensors
	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)

	# Get low bound on estimate (smallest value estimate can take on without becoming loss becoming undefined)
	# Basically, per snp h2 has to be greater than -1/N.
	per_snp_h2_est_lb = -np.min(gwas_beta_var)*.99


	# Initialize variables to optimize over
	h2_raw_variable = tf.Variable(initial_value=.001 - (per_snp_h2_est_lb*num_snps),trainable=True, name='beta_sq_raw', dtype=tf.float32) 
	#h2_variable = tf.math.softplus(h2_raw_variable) + est_lb



	# Lopp through windows
	prev_est_h2 = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = sumstat_reml_no_ld_loss(gwas_beta, gwas_beta_var, h2_raw_variable, num_snps, per_snp_h2_est_lb)

		trainable_variables = []
		trainable_variables.append(h2_raw_variable)


		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		# Get current estimate
		h2_variable = tf.nn.relu(h2_raw_variable) + (per_snp_h2_est_lb*num_snps)
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



def snp_h2_with_sumstat_ldsc_no_ld(gwas_beta, gwas_beta_se, sample_size):
	num_snps = len(gwas_beta)
	z_scores = gwas_beta/gwas_beta_se

	chi_sq = np.square(z_scores)

	h2_est = np.mean(chi_sq - 1.0)*num_snps/sample_size

	return h2_est




#########################
# Command line args
#########################
n_snps = int(sys.argv[1])
sim_h2 = float(sys.argv[2])
sample_size = int(sys.argv[3])
trait_architecture = sys.argv[4]
output_stem = sys.argv[5]

# Number of unique simulations
n_sims = 100

# Set seed
np.random.seed(1)


# Open and print header to output file
output_file = output_stem + '_est_h2_summary.txt'
t = open(output_file,'w')
t.write('sim_iter\tsim_h2\tsim_h2_v2\tldsc_est\tsumstat_reml_est\n')


# Loop through simulations
for sim_iter in range(n_sims):

	# Simulate the data
	Y, G, sim_causal_effects, obs_sim_h2 = simulate_data(sample_size, n_snps, sim_h2, trait_architecture)
	obs_sim_h2_v2 = np.sum(np.square(sim_causal_effects))


	# Get marginal summary statistics
	gwas_beta, gwas_beta_se = get_marginal_summary_statistics(Y, G)

	# Method 1: Run LDSC
	gwas_snp_h2_ldsc = snp_h2_with_sumstat_ldsc_no_ld(gwas_beta, gwas_beta_se, sample_size)
	# Method 2: estimate gwas snp h2 using sumstat reml
	gwas_snp_h2_sumstat_reml = snp_h2_with_sumstat_reml_no_ld(gwas_beta, gwas_beta_se)


	# Print to output file
	t.write(str(sim_iter) + '\t' + str(obs_sim_h2) + '\t' + str(obs_sim_h2_v2) + '\t' + str(gwas_snp_h2_ldsc) + '\t' + str(gwas_snp_h2_sumstat_reml) + '\n')
	t.flush()


t.close()
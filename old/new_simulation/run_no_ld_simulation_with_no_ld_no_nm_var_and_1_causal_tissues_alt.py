import sys
import numpy as np 
import pandas as pd
import os
import pdb
#import tensorflow as tf
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
from scipy import odr



def simulate_data(gwas_ss, n_snps, eqtl_ss_1, eqtl_ss_2, eqtl_ss_3, med_h2_1, med_h2_2, med_h2_3, ge_h2_1, ge_h2_2, ge_h2_3):
	# Simulate causal gene-trait effects
	sim_alpha_1 = -np.sqrt(med_h2_1)/np.sqrt(ge_h2_1)
	sim_alpha_2 = -np.sqrt(med_h2_2)/np.sqrt(ge_h2_2)
	sim_alpha_3 = -np.sqrt(med_h2_3)/np.sqrt(ge_h2_3)

	# Simulate causal variant_gene effects
	cov_mat = np.ones((3,3))*.75
	cov_mat[0,0] = ge_h2_1/n_snps
	cov_mat[1,1] = ge_h2_2/n_snps
	cov_mat[2,2] = ge_h2_3/n_snps
	corr = 0.7
	cov_mat[0,1] = np.sqrt(ge_h2_1/n_snps)*np.sqrt(ge_h2_2/n_snps)*corr
	cov_mat[1,0] = np.sqrt(ge_h2_1/n_snps)*np.sqrt(ge_h2_2/n_snps)*corr
	cov_mat[0,2] = np.sqrt(ge_h2_1/n_snps)*np.sqrt(ge_h2_3/n_snps)*corr
	cov_mat[2,0] = np.sqrt(ge_h2_1/n_snps)*np.sqrt(ge_h2_3/n_snps)*corr
	cov_mat[1,2] = np.sqrt(ge_h2_2/n_snps)*np.sqrt(ge_h2_3/n_snps)*corr
	cov_mat[2,1] = np.sqrt(ge_h2_2/n_snps)*np.sqrt(ge_h2_3/n_snps)*corr


	gene_causal_eqtl_effects = []
	for snp_iter in range(n_snps):
		gene_snp_causal_eqtl_effects = np.random.multivariate_normal(np.zeros(3), cov_mat)
		gene_causal_eqtl_effects.append(gene_snp_causal_eqtl_effects)
	gene_causal_eqtl_effects = np.asarray(gene_causal_eqtl_effects)


	# Simulate eQTL genotype data
	eqtl_geno_1 = np.random.normal(loc=0,scale=1.0,size=(eqtl_ss_1,n_snps))
	for snp_iter in range(n_snps):
		eqtl_geno_1[:,snp_iter] = (eqtl_geno_1[:,snp_iter] - np.mean(eqtl_geno_1[:,snp_iter]))/np.std(eqtl_geno_1[:,snp_iter])
	eqtl_geno_2 = np.random.normal(loc=0,scale=1.0,size=(eqtl_ss_2,n_snps))
	for snp_iter in range(n_snps):
		eqtl_geno_2[:,snp_iter] = (eqtl_geno_2[:,snp_iter] - np.mean(eqtl_geno_2[:,snp_iter]))/np.std(eqtl_geno_2[:,snp_iter])
	eqtl_geno_3 = np.random.normal(loc=0,scale=1.0,size=(eqtl_ss_3,n_snps))
	for snp_iter in range(n_snps):
		eqtl_geno_3[:,snp_iter] = (eqtl_geno_3[:,snp_iter] - np.mean(eqtl_geno_3[:,snp_iter]))/np.std(eqtl_geno_3[:,snp_iter])


	# Simulate GWAS genotype data
	gwas_geno = np.random.normal(loc=0,scale=1.0,size=(gwas_ss,n_snps))
	for snp_iter in range(n_snps):
		gwas_geno[:,snp_iter] = (gwas_geno[:,snp_iter] - np.mean(gwas_geno[:,snp_iter]))/np.std(gwas_geno[:,snp_iter])

	# Simulate Gene trait value
	genetic_gene_trait_1 = np.dot(eqtl_geno_1, gene_causal_eqtl_effects[:,0])
	E_1 = np.random.normal(loc=genetic_gene_trait_1, scale=np.sqrt(1.0-np.var(genetic_gene_trait_1)))
	E_1 = (E_1 - np.mean(E_1))/np.std(E_1)
	genetic_gene_trait_2 = np.dot(eqtl_geno_2, gene_causal_eqtl_effects[:,1])
	E_2 = np.random.normal(loc=genetic_gene_trait_2, scale=np.sqrt(1.0-np.var(genetic_gene_trait_2)))
	E_2 = (E_2 - np.mean(E_2))/np.std(E_2)
	genetic_gene_trait_3 = np.dot(eqtl_geno_3, gene_causal_eqtl_effects[:,2])
	E_3 = np.random.normal(loc=genetic_gene_trait_3, scale=np.sqrt(1.0-np.var(genetic_gene_trait_3)))
	E_3 = (E_3 - np.mean(E_3))/np.std(E_3)

	# Simulate GWAS trait value
	genetic_gene_1 = np.dot(gwas_geno, gene_causal_eqtl_effects[:,0])
	genetic_gene_2 = np.dot(gwas_geno, gene_causal_eqtl_effects[:,1])
	genetic_gene_3 = np.dot(gwas_geno, gene_causal_eqtl_effects[:,2])

	#std_genetic_gene = genetic_gene/np.std(genetic_gene)
	genetic_trait = genetic_gene_1*sim_alpha_1 + genetic_gene_2*sim_alpha_2 + genetic_gene_3*sim_alpha_3
	Y = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0-np.var(genetic_trait)))
	Y = (Y - np.mean(Y))/np.std(Y)

	return sim_alpha_1, sim_alpha_2, sim_alpha_3, Y, E_1, E_2, E_3, gwas_geno, eqtl_geno_1, eqtl_geno_2, eqtl_geno_3, gene_causal_eqtl_effects


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

def linfit(beta, x):
	return beta[0]*x


def scipy_odr_mediated_effect_est(gwas_beta, gene_betas, eqtl_sss, gwas_ss, sstol_string='default'):
	eqtl_weights = 1.0/(1.0/eqtl_sss)
	gwas_weight =(1.0/(1.0/gwas_ss))

	linear = odr.Model(linfit)

	data = odr.Data(gene_betas, gwas_beta, wd=eqtl_weights,we=gwas_weight)

	if sstol_string == 'default':
		odr_obj = odr.ODR(data, linear, maxit=20000, beta0=[0.0])
	elif sstol_string == 'stringent':
		odr_obj = odr.ODR(data, linear,sstol=1e-19, maxit=20000, beta0=[0.0])
	else:
		print('not currently implemented')
		pdb.set_trace()

	output = odr_obj.run()

	if output.stopreason[0] != 'Sum of squares convergence' and output.stopreason[0] != 'Parameter convergence':
		print('ODR convergence error')
		pdb.set_trace()

	return output.beta, np.transpose(output.delta)


def scipy_odr_mediated_effect_est_per_snp_var(gwas_beta, gene_betas, gene_beta_variances, gwas_beta_variances, sstol_string='default'):
	eqtl_weights = 1.0/gene_beta_variances
	gwas_weights = 1.0/gwas_beta_variances

	linear = odr.Model(linfit)
	data = odr.Data(gene_betas[:,0], gwas_beta, wd=eqtl_weights[:,0],we=gwas_weights)
	if sstol_string == 'default':
		odr_obj = odr.ODR(data, linear, maxit=20000, beta0=[0.0])
	elif sstol_string == 'stringent':
		odr_obj = odr.ODR(data, linear,sstol=1e-19, maxit=20000, beta0=[0.0])
	else:
		print('not currently implemented')
		pdb.set_trace()
	output = odr_obj.run()

	new_betas = gene_betas + output.delta

	if output.stopreason[0] != 'Sum of squares convergence' and output.stopreason[0] != 'Parameter convergence':
		print('ODR convergence error')
		pdb.set_trace()

	return output.beta, np.transpose(output.delta)


def mediated_effect_est_loss_function(gwas_beta, eqtl_beta, gwas_ss, eqtl_sss, eqtl_beta_deltas_var, alpha_var):
	loss1 = tf.math.square(gwas_beta - tf.linalg.matmul(eqtl_beta + eqtl_beta_deltas_var, alpha_var))*gwas_ss
	loss2 = tf.linalg.matmul(tf.math.square(eqtl_beta_deltas_var), eqtl_sss)

	agg_loss = loss1 + loss2

	return tf.math.reduce_sum(agg_loss)


def custom_tf_odr_mediated_effect_est(gwas_beta, eqtl_betas, eqtl_sss, gwas_ss, max_epochs=400000, conv_thresh=1e-15, smart_init=True):
	n_genes = len(eqtl_sss)

	if smart_init:
		olser = sm.OLS(gwas_beta, eqtl_betas).fit()
		alpha_var = tf.Variable(initial_value=olser.params.reshape((n_genes,1)).astype(np.float32),trainable=True, name='alpha')
	else:
		alpha_var = tf.Variable(initial_value=np.zeros((n_genes,1)).astype(np.float32),trainable=True, name='alpha')

	eqtl_beta_deltas_var = tf.Variable(initial_value=np.zeros(eqtl_betas.shape).astype(np.float32),trainable=True, name='eqtl_beta_deltas')

	optimizer = tf.keras.optimizers.Adam()

	gwas_beta = tf.convert_to_tensor(gwas_beta.reshape(len(gwas_beta),1), dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_betas, dtype=tf.float32)
	eqtl_sss = tf.convert_to_tensor(eqtl_sss.reshape(len(eqtl_sss),1), dtype=tf.float32)

	# Lopp through windows
	prev_est_alphas = np.ones(n_genes)*100000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = mediated_effect_est_loss_function(gwas_beta, eqtl_beta, gwas_ss, eqtl_sss, eqtl_beta_deltas_var, alpha_var)
			
		trainable_variables = []
		trainable_variables.append(alpha_var)
		trainable_variables.append(eqtl_beta_deltas_var)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))

		curr_alphas = np.asarray(alpha_var)[:,0]
		diff = np.sum(np.square(prev_est_alphas - curr_alphas))
		prev_est_alphas = curr_alphas

		if diff < conv_thresh:
			converged = True
			break

	if converged == False:
		print(diff)
		print('TF algorithm for loss did not not converge')

	return np.asarray(alpha_var)[:,0], np.asarray(eqtl_beta_deltas_var)

def compute_odr_mediated_effect_loss_np(est_alphas, est_delta, gwas_beta, eqtl_beta,eqtl_sss, gwas_ss):
	loss1 = np.square(gwas_beta - np.matmul(np.transpose(eqtl_beta) + est_delta.reshape((len(est_delta),1)), est_alphas))*gwas_ss
	loss2 = np.square(est_delta)*eqtl_sss[0]

	agg_loss = loss1 + loss2

	return np.sum(agg_loss)

def compute_odr_mediated_effect_loss_per_snp_var_np(est_alphas, est_delta, gwas_beta, eqtl_beta, gene_beta_variances, gwas_beta_variance):
	loss1 = np.square(gwas_beta - np.matmul(np.transpose(eqtl_beta) + est_delta.reshape((len(est_delta),1)), est_alphas))/gwas_beta_variance
	loss2 = np.square(est_delta)/gene_beta_variances[:,0]

	agg_loss = loss1 + loss2

	return np.sum(agg_loss)


#########################
# Command line args
#########################
global_simulation_number = sys.argv[1]
n_sims = int(sys.argv[2])
gwas_ss = int(sys.argv[3])
n_snps = int(sys.argv[4])
eqtl_ss_1 = int(sys.argv[5])
eqtl_ss_2 = int(sys.argv[6])
eqtl_ss_3 = int(sys.argv[7])
med_h2_1 = float(sys.argv[8])
med_h2_2 = float(sys.argv[9])
med_h2_3 = float(sys.argv[10])
ge_h2_1 = float(sys.argv[11])
ge_h2_2 = float(sys.argv[12])
ge_h2_3 = float(sys.argv[13])
output_root = sys.argv[14]
global_simulation_number = int(sys.argv[15])


# Set seed
np.random.seed(global_simulation_number)


# Open and print header to output file
output_file = output_root + '_effect_est_res_summary.txt'
t = open(output_file,'w')
t.write('method\tsim\tsim_alpha_1\test_alpha_1\tloss\n')


# Loop through sims
for sim_iter in range(n_sims):
	print(sim_iter)

	# First simulate data
	sim_alpha_1, sim_alpha_2, sim_alpha_3, Y, E_1, E_2, E_3, Y_geno, E_1_geno, E_2_geno, E_3_geno, gene_causal_eqtl_effects = simulate_data(gwas_ss, n_snps, eqtl_ss_1, eqtl_ss_2, eqtl_ss_3, med_h2_1, med_h2_2, med_h2_3, ge_h2_1, ge_h2_2, ge_h2_3)

	# Next get summary statistics
	gwas_beta, gwas_beta_se = get_marginal_summary_statistics(Y, Y_geno)
	eqtl_1_beta, eqtl_1_beta_se = get_marginal_summary_statistics(E_1, E_1_geno)



	# Organize data
	gene_betas = np.row_stack((eqtl_1_beta))
	eqtl_sss = np.asarray([eqtl_ss_1])
	gene_beta_variances = np.row_stack((np.square(eqtl_1_beta_se)))


	#################################
	# Estimate alphas using scipy ODR per snp var
	method = 'scipy_odr_effect_est_per_snp_var'
	# Optimize fxn
	est_alphas, est_delta = scipy_odr_mediated_effect_est_per_snp_var(gwas_beta, gene_betas, gene_beta_variances, np.square(gwas_beta_se), sstol_string='default')
	new_betas = gene_betas[:,0] + est_delta
	#olser = sm.OLS(gwas_beta, gene_betas[:,0]).fit().params
	new_beta_r_squared = np.square(new_betas)/np.sum(np.square(new_betas))

	new_beta_r_squared2 = []
	pred = np.dot(Y_geno, new_betas)
	for kk in range(Y_geno.shape[1]):
		temp_rsquared = np.square(np.corrcoef(Y_geno[:,kk], pred)[0,1])
		correction_factor = (1.0 - temp_rsquared)/(gwas_ss-2)
		new_beta_r_squared2.append(temp_rsquared)

	new_beta_r_squared2 = np.asarray(new_beta_r_squared2)



	gwas_z = gwas_beta/gwas_beta_se

	olser = sm.OLS(np.square(gwas_z) - 1, sm.add_constant(new_beta_r_squared -.1)).fit()
	olser = sm.OLS(np.square(gwas_z) - 1, sm.add_constant(np.square(gene_causal_eqtl_effects[:,0])/np.sum(np.square(gene_causal_eqtl_effects[:,0])))).fit()

	#olser = sm.OLS(np.square(gwas_z) - 1, np.ones(len(gwas_z))).fit()
	print(olser.params[1]/10000)

	pdb.set_trace()
	'''
	#################################
	# Estimate alphas using linear regression with unobserved (true eqtl effect sizes)
	method = 'linear_reg_unobsered_true_effect_sizes'
	olser = sm.OLS(gwas_beta, gene_causal_eqtl_effects[:,0]).fit()
	t.write(method + '\t' + str(sim_iter) + '\t' + str(sim_alpha_1) + '\t' + str(olser.params[0]) + '\t' + 'NA' + '\n')


	#################################
	# Estimate alphas using scipy ODR
	method = 'scipy_odr_effect_est'
	# Optimize fxn
	est_alphas, est_delta = scipy_odr_mediated_effect_est(gwas_beta, gene_betas[:,0], eqtl_sss[0], gwas_ss, sstol_string='default')
	# Compute loss given optimized parameters
	loss = compute_odr_mediated_effect_loss_np(est_alphas, est_delta, gwas_beta, np.transpose(gene_betas),eqtl_sss, gwas_ss)
	# Print to output
	t.write(method + '\t' + str(sim_iter) + '\t' + str(sim_alpha_1) + '\t' + str(est_alphas[0]) + '\t' + str(loss) + '\n')
	t.flush()

	#################################
	# Estimate alphas using scipy ODR
	method = 'scipy_odr_effect_est_stringent'
	# Optimize fxn
	est_alphas, est_delta = scipy_odr_mediated_effect_est(gwas_beta, gene_betas[:,0], eqtl_sss[0], gwas_ss, sstol_string='stringent')
	# Compute loss given optimized parameters
	loss = compute_odr_mediated_effect_loss_np(est_alphas, est_delta, gwas_beta, np.transpose(gene_betas),eqtl_sss, gwas_ss)
	# Print to output
	t.write(method + '\t' + str(sim_iter) + '\t' + str(sim_alpha_1) + '\t' + str(est_alphas[0]) + '\t' + str(loss) + '\n')
	t.flush()


	#################################
	# Estimate alphas using scipy ODR per snp var
	method = 'scipy_odr_effect_est_per_snp_var'
	# Optimize fxn
	est_alphas, est_delta = scipy_odr_mediated_effect_est_per_snp_var(gwas_beta, gene_betas, gene_beta_variances, np.square(gwas_beta_se), sstol_string='default')
	# Compute loss given optimized parameters
	loss = compute_odr_mediated_effect_loss_per_snp_var_np(est_alphas, est_delta, gwas_beta, np.transpose(gene_betas),gene_beta_variances, np.square(gwas_beta_se))
	# Print to output
	t.write(method + '\t' + str(sim_iter) + '\t' + str(sim_alpha_1) + '\t' + str(est_alphas[0]) + '\t' + str(loss) + '\n')


	#################################
	# Estimate alphas using scipy ODR per snp var
	method = 'scipy_odr_effect_est_stringent_per_snp_var'
	# Optimize fxn
	est_alphas, est_delta = scipy_odr_mediated_effect_est_per_snp_var(gwas_beta, gene_betas, gene_beta_variances, np.square(gwas_beta_se), sstol_string='stringent')
	# Compute loss given optimized parameters
	loss = compute_odr_mediated_effect_loss_per_snp_var_np(est_alphas, est_delta, gwas_beta, np.transpose(gene_betas),gene_beta_variances, np.square(gwas_beta_se))
	# Print to output
	t.write(method + '\t' + str(sim_iter) + '\t' + str(sim_alpha_1) + '\t' + str(est_alphas[0]) + '\t' + str(loss) + '\n')
	'''


	'''
	#################################
	# Estimate alphas using custom tensor flow implementation of ODR
	method = 'custom_tf_odr_effect_est'
	# Optimize fxn
	est_alphas_custom_tf, est_delta_beta_custom_tf = custom_tf_odr_mediated_effect_est(gwas_beta, np.transpose(gene_betas), eqtl_sss, gwas_ss, smart_init=False)
	# Compute loss given optimized parameters
	loss_custom_tf = compute_odr_mediated_effect_loss_np(est_alphas_custom_tf, est_delta_beta_custom_tf, gwas_beta, np.transpose(gene_betas),eqtl_sss, gwas_ss)
	# Print to output
	t.write(method + '\t' + str(sim_iter) + '\t' + str(sim_alpha_1) + '\t' + str(est_alphas_custom_tf[0]) + '\t' + str(sim_alpha_2) + '\t' + str(est_alphas_custom_tf[1]) + '\t' + str(sim_alpha_3) + '\t' + str(est_alphas_custom_tf[2]) + '\t' + str(loss_custom_tf) + '\n')
	'''
	# Print to output

t.close()


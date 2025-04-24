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


def est_reg_coef(xx, yy):
	olser = sm.OLS(yy, sm.add_constant(xx)).fit()
	return olser.params[1]
def est_reg_coef_fast(xx,yy):
	return np.cov(xx,yy)[0,1]/np.var(xx,ddof=1)
def est_squared_reg_coef_fast(xx,yy):
	return np.square(np.cov(xx,yy)[0,1]/np.var(xx,ddof=1))
def print_95_ci(arr):
	meany = np.mean(arr)
	se = np.std(arr)/np.sqrt(len(arr))
	ub = meany + (1.96*se)
	lb = meany - (1.96*se)

	print(str(meany) + ':   ' + '[' + str(lb) + ', ' + str(ub) + ']')
	return


def simulate_gwas_and_trait_summary_statistics(eqtl_ss, gwas_ss, ge_h2, h2, frac_med_h2, n_snps):
	# Mediated h2
	med_h2 = h2*frac_med_h2

	nm_h2 = h2 - med_h2

	# Simulate causal gene-trait effects
	sim_alpha = -np.sqrt(med_h2)
	# Simulate causal variant_gene effects
	gene_causal_eqtl_effects = np.hstack((np.random.normal(loc=0, scale=np.sqrt(ge_h2/(n_snps*.5)), size=int(n_snps*.5)), np.zeros(int(n_snps*.5))))

	# Simulate causal variant_trait effects
	variant_trait_causal_effects = np.random.normal(loc=0, scale=np.sqrt(nm_h2/n_snps), size=n_snps)


	# Simulate eQTL genotype data
	eqtl_geno = np.random.normal(loc=0,scale=1.0,size=(eqtl_ss,n_snps))
	for snp_iter in range(n_snps):
		eqtl_geno[:,snp_iter] = (eqtl_geno[:,snp_iter] - np.mean(eqtl_geno[:,snp_iter]))/np.std(eqtl_geno[:,snp_iter])

	# Simulate GWAS genotype data
	gwas_geno = np.random.normal(loc=0,scale=1.0,size=(gwas_ss,n_snps))
	for snp_iter in range(n_snps):
		gwas_geno[:,snp_iter] = (gwas_geno[:,snp_iter] - np.mean(gwas_geno[:,snp_iter]))/np.std(gwas_geno[:,snp_iter])

	# Simulate Gene trait value
	genetic_gene_trait = np.dot(eqtl_geno, gene_causal_eqtl_effects)
	E = np.random.normal(loc=genetic_gene_trait, scale=np.sqrt(1.0-np.var(genetic_gene_trait)))
	E = (E - np.mean(E))/np.std(E)

	# Simulate GWAS trait value
	genetic_gene = np.dot(gwas_geno, gene_causal_eqtl_effects)
	std_genetic_gene = genetic_gene/np.std(genetic_gene)
	genetic_trait = std_genetic_gene*sim_alpha + np.dot(gwas_geno, variant_trait_causal_effects)
	Y = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0-np.var(genetic_trait)))
	Y = (Y - np.mean(Y))/np.std(Y)

	# Get gwas summary stats
	gene_z = []
	for snp_iter in range(n_snps):
		olser = sm.OLS(E, eqtl_geno[:,snp_iter]).fit()
		beta = olser.params[0]
		beta_se = olser.bse[0]
		zed = beta/beta_se
		gene_z.append(zed)
	gene_z = np.asarray(gene_z)

	# Get gwas summary stats
	gwas_z = []
	for snp_iter in range(n_snps):
		olser = sm.OLS(Y, gwas_geno[:,snp_iter]).fit()
		beta = olser.params[0]
		beta_se = olser.bse[0]
		zed = beta/beta_se
		gwas_z.append(zed)
	gwas_z = np.asarray(gwas_z)


	return np.square(gwas_z), np.square(gene_z)



def run_ld_score_regression(chi_sq, ld_scores, NN, MM, intercept=False):
	if intercept == False:
		olser = sm.OLS(chi_sq -1, NN*ld_scores/MM).fit()
		h2_est = olser.params[0]
	else:
		print('still need to implement with intercept')
		pdb.set_trace()


	return h2_est

def init_linear_mapping_from_genomic_annotations_to_gamma(annotation_data_dimension):
	# Initialize Neural network model
	model = tf.keras.models.Sequential()
	#model.add(tf.keras.layers.Dense(units=1, activation='softplus', input_dim=annotation_data_dimension, bias_initializer=tf.constant_initializer(-10.0)))
	model.add(tf.keras.layers.Dense(units=1, activation='softplus', input_dim=annotation_data_dimension))

	return model


def ldsc_tf_loss_fxn(chi_sq, pred_tau, samp_size, ld_sq, intercept_variable):
	pred_chi_sq = (samp_size*tf.linalg.matmul(ld_sq, pred_tau)) + (tf.math.softplus(intercept_variable))

	log_like = (-.5)*tf.math.log(chi_sq) - tf.math.divide(chi_sq, 2.0*pred_chi_sq) - (.5*tf.math.log(2.0*pred_chi_sq))

	return -tf.math.reduce_sum(log_like), log_like


def ldsc_tf_joint_loss_fxn(gwas_chi_sq, gene_chi_sq, gwas_ss, eqtl_ss, squared_ld, nm_var_trait_h2, log_var_gene_h2, log_gene_trait_h2, log_intercept_variable_gwas, log_intercept_variable_eqtl):
	# Convert variables to be nonnegative with softplus transform
	var_gene_h2 = tf.math.softplus(log_var_gene_h2)
	gene_trait_h2 = tf.math.softplus(log_gene_trait_h2)


	# Get predicted gwas chi squared stats
	#pred_gwas_tau = nm_var_trait_h2 + var_gene_h2*gene_trait_h2
	pred_gwas_tau = nm_var_trait_h2 + (tf.math.divide(var_gene_h2, tf.math.reduce_sum(var_gene_h2))*gene_trait_h2)
	pred_gwas_chi_sq = (gwas_ss*tf.linalg.matmul(squared_ld, pred_gwas_tau)) + (tf.math.softplus(log_intercept_variable_gwas))

	# Get predicted eqtl chi squared stats
	pred_gene_chi_sq = (eqtl_ss*tf.linalg.matmul(squared_ld, var_gene_h2)) + (tf.math.softplus(log_intercept_variable_eqtl))
	#pred_gene_chi_sq = (eqtl_ss*tf.linalg.matmul(squared_ld, nm_var_trait_h2)) + (tf.math.softplus(log_intercept_variable_eqtl))


	log_like = (-.5)*tf.math.log(gwas_chi_sq) - tf.math.divide(gwas_chi_sq, 2.0*pred_gwas_chi_sq) - (.5*tf.math.log(2.0*pred_gwas_chi_sq))
	log_like = log_like + (-.5)*tf.math.log(gene_chi_sq) - tf.math.divide(gene_chi_sq, 2.0*pred_gene_chi_sq) - (.5*tf.math.log(2.0*pred_gene_chi_sq))
	#log_like = (-.5)*tf.math.log(gene_chi_sq) - tf.math.divide(gene_chi_sq, 2.0*pred_gene_chi_sq) - (.5*tf.math.log(2.0*pred_gene_chi_sq))


	return -tf.math.reduce_sum(log_like), log_like




def run_joint_log_lss_regression(gwas_chi_sq, gene_chi_sq, ld_sq, gwas_ss, eqtl_ss, MM, intercept=False, max_epochs=100000, convergence_thresh=1e-8):
	# Initialize mapping from annotations to NM-var per snp heritability
	nm_var_genomic_anno_to_gamma_model = init_linear_mapping_from_genomic_annotations_to_gamma(1)
	# Initialize variant to gene heritabiliites
	log_var_gene_h2 = tf.Variable(initial_value=np.ones((MM,1)).astype(np.float32)*-0.0,trainable=True, name='var_gene_h2')
	# Initialize gene to trait h2
	log_gene_trait_h2 = tf.Variable(initial_value=np.ones((1,1)).astype(np.float32)*0.0,trainable=True, name='gene_trait_h2')

	optimizer = tf.keras.optimizers.Adam()

	# Whether or not to learn intercept in LDSC
	# Initial value is np.log(np.exp(1)-1.0) [which equals 1 when put through softplus activation function]
	if intercept:
		log_intercept_variable_gwas = tf.Variable(initial_value=0.541324854612918,trainable=True, name='intercept_gwas')
		log_intercept_variable_eqtl = tf.Variable(initial_value=0.541324854612918,trainable=True, name='intercept_eqtl')
	elif intercept == False:
		log_intercept_variable_gwas = tf.Variable(initial_value=0.541324854612918,trainable=False, name='intercept_gwas')
		log_intercept_variable_eqtl = tf.Variable(initial_value=0.541324854612918,trainable=False, name='intercept_eqtl')
	else:
		print('assumption error: intercept model called ' + intercept + ' not currently implemented')
		pdb.set_trace()

	gwas_chi_sq = tf.convert_to_tensor(gwas_chi_sq.reshape(len(gwas_chi_sq),1), dtype=tf.float32)
	gene_chi_sq = tf.convert_to_tensor(gene_chi_sq.reshape(len(gene_chi_sq),1), dtype=tf.float32)
	squared_ld = tf.convert_to_tensor(ld_sq, dtype=tf.float32)

	# Lopp through windows
	prev_est_h2 = 1000
	for epoch_iter in range(max_epochs):
		# Keep track of training log likelihoods and weights of each regression snp
		#epoch_training_log_likelihoods = []
		#epoch_training_weights = []

		# Loop through windows
		#print('###################################')
		#print('epoch iter ' + str(epoch_iter))
		#print('###################################')

		window_genomic_anno = np.ones((MM, 1))
		
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			nm_var_trait_h2 = nm_var_genomic_anno_to_gamma_model(window_genomic_anno, training=True)

			loss_value, log_likelihoods = ldsc_tf_joint_loss_fxn(gwas_chi_sq, gene_chi_sq, gwas_ss, eqtl_ss, squared_ld, nm_var_trait_h2, log_var_gene_h2, log_gene_trait_h2, log_intercept_variable_gwas, log_intercept_variable_eqtl)
			
		# Define trainable variables
		trainable_variables = nm_var_genomic_anno_to_gamma_model.trainable_weights
		#trainable_variables = []
		if intercept:
			trainable_variables.append(log_intercept_variable)
		# Add var_gene_h2 and gene_trait_h2
		trainable_variables.append(log_gene_trait_h2)
		trainable_variables.append(log_var_gene_h2)
		# Compute and apply gradients
		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))

		# Extract relevent variables
		pred_gene_trait_var = np.asarray(tf.math.softplus(np.asarray(log_gene_trait_h2)[0,0]))*1.0
		pred_var_gene_var = np.asarray(tf.math.softplus(log_var_gene_h2))[:,0]
		#pred_var_gene_var = np.asarray((log_var_gene_h2))[:,0]

		est_nm_h2 = np.asarray(nm_var_trait_h2[0]*1.0)[0]*MM


		if np.divmod(epoch_iter, 500)[1] == 0:
			print('##############################')
			print('trait nm h2: ' + str(est_nm_h2))
			print('trait med h2: ' + str(pred_gene_trait_var))
			print('gene h2: ' + str(np.sum(pred_var_gene_var)))
			print('')
		if np.divmod(epoch_iter, 20000)[1] == 0 and epoch_iter != 0:
			pdb.set_trace()

def run_joint_log_lss_regression_binned(gwas_chi_sq, gene_chi_sq, ld_sq, gwas_ss, eqtl_ss, MM, intercept=False, max_epochs=100000, convergence_thresh=1e-8, n_bins=21):
	# Initialize mapping from annotations to NM-var per snp heritability
	nm_var_genomic_anno_to_gamma_model = init_linear_mapping_from_genomic_annotations_to_gamma(1)
	# Initialize variant to gene heritabiliites
	log_bin_var_gene_h2 = tf.Variable(initial_value=np.ones((n_bins,1)).astype(np.float32)*-0.0,trainable=True, name='var_gene_h2')
	# Initialize gene to trait h2
	log_gene_trait_h2 = tf.Variable(initial_value=np.ones((1,1)).astype(np.float32)*0.0,trainable=True, name='gene_trait_h2')



	# Split snps into bins
	snp_bin_mat = np.zeros((MM, n_bins))
	bin_splits = np.array_split(np.arange(MM),n_bins)
	for ii, bin_split in enumerate(bin_splits):
		snp_bin_mat[bin_split, ii] = 1.0
	snp_bin_mat = tf.convert_to_tensor(snp_bin_mat, dtype=tf.float32)

	# Initialize optimizer
	optimizer = tf.keras.optimizers.Adam()

	# Whether or not to learn intercept in LDSC
	# Initial value is np.log(np.exp(1)-1.0) [which equals 1 when put through softplus activation function]
	if intercept:
		log_intercept_variable_gwas = tf.Variable(initial_value=0.541324854612918,trainable=True, name='intercept_gwas')
		log_intercept_variable_eqtl = tf.Variable(initial_value=0.541324854612918,trainable=True, name='intercept_eqtl')
	elif intercept == False:
		log_intercept_variable_gwas = tf.Variable(initial_value=0.541324854612918,trainable=False, name='intercept_gwas')
		log_intercept_variable_eqtl = tf.Variable(initial_value=0.541324854612918,trainable=False, name='intercept_eqtl')
	else:
		print('assumption error: intercept model called ' + intercept + ' not currently implemented')
		pdb.set_trace()

	gwas_chi_sq = tf.convert_to_tensor(gwas_chi_sq.reshape(len(gwas_chi_sq),1), dtype=tf.float32)
	gene_chi_sq = tf.convert_to_tensor(gene_chi_sq.reshape(len(gene_chi_sq),1), dtype=tf.float32)
	squared_ld = tf.convert_to_tensor(ld_sq, dtype=tf.float32)

	# Lopp through windows
	prev_est_h2 = 1000
	for epoch_iter in range(max_epochs):
		# Keep track of training log likelihoods and weights of each regression snp
		#epoch_training_log_likelihoods = []
		#epoch_training_weights = []

		# Loop through windows
		#print('###################################')
		#print('epoch iter ' + str(epoch_iter))
		#print('###################################')

		window_genomic_anno = np.ones((MM, 1))
		
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			nm_var_trait_h2 = nm_var_genomic_anno_to_gamma_model(window_genomic_anno, training=True)

			log_var_gene_h2 = tf.linalg.matmul(snp_bin_mat, log_bin_var_gene_h2)

			loss_value, log_likelihoods = ldsc_tf_joint_loss_fxn(gwas_chi_sq, gene_chi_sq, gwas_ss, eqtl_ss, squared_ld, nm_var_trait_h2, log_var_gene_h2, log_gene_trait_h2, log_intercept_variable_gwas, log_intercept_variable_eqtl)
			
		# Define trainable variables
		trainable_variables = nm_var_genomic_anno_to_gamma_model.trainable_weights
		#trainable_variables = []
		if intercept:
			trainable_variables.append(log_intercept_variable)
		# Add var_gene_h2 and gene_trait_h2
		trainable_variables.append(log_gene_trait_h2)
		trainable_variables.append(log_bin_var_gene_h2)
		# Compute and apply gradients
		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))

		# Extract relevent variables
		pred_gene_trait_var = np.asarray(tf.math.softplus(np.asarray(log_gene_trait_h2)[0,0]))*1.0
		pred_var_gene_var = np.asarray(tf.math.softplus(tf.linalg.matmul(snp_bin_mat, log_bin_var_gene_h2)))[:,0]
		#pred_var_gene_var = np.asarray((log_var_gene_h2))[:,0]

		est_nm_h2 = np.asarray(nm_var_trait_h2[0]*1.0)[0]*MM


		if np.divmod(epoch_iter, 500)[1] == 0:
			print('##############################')
			print('trait nm h2: ' + str(est_nm_h2))
			print('trait med h2: ' + str(pred_gene_trait_var))
			print('gene h2: ' + str(np.sum(pred_var_gene_var)))
			print('')
		if np.divmod(epoch_iter, 20000)[1] == 0 and epoch_iter != 0:
			pdb.set_trace()



def run_log_lss_regression(chi_sq, ld_sq, NN, MM, intercept=False, max_epochs=100000, convergence_thresh=1e-8):
	# Initialize mapping from annotations to per snp heritability
	genomic_anno_to_gamma_model = init_linear_mapping_from_genomic_annotations_to_gamma(1)
	optimizer = tf.keras.optimizers.Adam()

	# Whether or not to learn intercept in LDSC
	# Initial value is np.log(np.exp(1)-1.0) [which equals 1 when put through softplus activation function]
	if intercept:
		log_intercept_variable = tf.Variable(initial_value=0.541324854612918,trainable=True, name='intercept')
	elif intercept == False:
		log_intercept_variable = tf.Variable(initial_value=0.541324854612918,trainable=False, name='intercept')
	else:
		print('assumption error: intercept model called ' + intercept + ' not currently implemented')
		pdb.set_trace()


	window_chi_sq = tf.convert_to_tensor(chi_sq.reshape(len(chi_sq),1), dtype=tf.float32)
	squared_ld = tf.convert_to_tensor(ld_sq, dtype=tf.float32)

	# Lopp through windows
	prev_est_h2 = 1000
	for epoch_iter in range(max_epochs):
		# Keep track of training log likelihoods and weights of each regression snp
		#epoch_training_log_likelihoods = []
		#epoch_training_weights = []

		# Loop through windows
		#print('###################################')
		#print('epoch iter ' + str(epoch_iter))
		#print('###################################')

		window_genomic_anno = np.ones((MM, 1))
		
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			window_pred_tau = genomic_anno_to_gamma_model(window_genomic_anno, training=True)

			loss_value, log_likelihoods = ldsc_tf_loss_fxn(window_chi_sq, window_pred_tau, NN, squared_ld, log_intercept_variable)
			
		# Define trainable variables
		trainable_variables = genomic_anno_to_gamma_model.trainable_weights
		if intercept:
			trainable_variables.append(log_intercept_variable)
		# Compute and apply gradients
		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		# Get current h2 esst
		window_pred_tau = genomic_anno_to_gamma_model(window_genomic_anno, training=False)

		est_h2 = np.asarray(window_pred_tau[0]*1.0)[0]*MM
		#print(est_h2)

		diff = np.abs(est_h2-prev_est_h2)
		if diff < convergence_thresh:
			break
		prev_est_h2 = est_h2

	return est_h2



#######################
# Command line args
#######################
simulation_number = sys.argv[1]
eqtl_ss = int(sys.argv[2])
gwas_ss = int(sys.argv[3])
ge_h2 = float(sys.argv[4])
h2 = float(sys.argv[5])
frac_med_h2 = float(sys.argv[6])
simulation_input_dir = sys.argv[7]
simulation_output_dir = sys.argv[8]


# Other parameters
num_snps=1000


# Simulate chi squared stats
gwas_chi_sq, gene_chi_sq = simulate_gwas_and_trait_summary_statistics(eqtl_ss, gwas_ss, ge_h2, h2, frac_med_h2, num_snps)

# Compute LDSC est of gwas variant h2
ldsc_est_gwas_h2 = run_ld_score_regression(gwas_chi_sq, np.ones(len(gwas_chi_sq)), gwas_ss, num_snps, intercept=False)
# Compute LDSC est of eQTL variant h2
ldsc_est_gene_h2 = run_ld_score_regression(gene_chi_sq, np.ones(len(gene_chi_sq)), eqtl_ss, num_snps, intercept=False)
print(ldsc_est_gene_h2)

'''
# Compute log_lss_regression est of gwas variant h2
log_lss_est_gwas_h2 = run_log_lss_regression(gwas_chi_sq, np.eye(len(gwas_chi_sq)), gwas_ss, num_snps, intercept=False)
# Compute log_lss_regression est of eQTL variant h2
log_lss_est_gene_h2 = run_log_lss_regression(gene_chi_sq, np.eye(len(gene_chi_sq)), eqtl_ss, num_snps, intercept=False)
'''

# Joint log_lss_regression est of mediated h2
#log_lss_est_gwas_h2 = run_joint_log_lss_regression(gwas_chi_sq, gene_chi_sq, np.eye(len(gwas_chi_sq)), gwas_ss, eqtl_ss, num_snps, intercept=False)
log_lss_est_gwas_h2 = run_joint_log_lss_regression_binned(gwas_chi_sq, gene_chi_sq, np.eye(len(gwas_chi_sq)), gwas_ss, eqtl_ss, num_snps, intercept=False)


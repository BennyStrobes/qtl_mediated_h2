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



def mediated_squared_effect_est_loss_function_unobserved_true(chi_sq, var_ld_scores, gene_ld_scores, beta_sq_var, alpha_sq_var):
	pred_chi_sq = var_ld_scores*beta_sq_var + gene_ld_scores*alpha_sq_var + 1.0
	log_like = (-.5)*tf.math.log(chi_sq) - tf.math.divide(chi_sq, 2.0*pred_chi_sq) - (.5*tf.math.log(2.0*pred_chi_sq))

	return -tf.math.reduce_sum(log_like)


def custom_tf_odr_mediated_squared_effect_est_unobserved_true(chi_sq_stats, var_ld_scores, gene_ld_scores, max_epochs=30000, conv_thresh=1e-15):
	#n_genes = eqtl_beta.shape[0]
	#n_snps = len(gwas_beta)

	# Initialize variables to optimize over
	alpha_sq_var = tf.Variable(initial_value=14.76,trainable=True, name='alpha_sq')
	beta_sq_var = tf.Variable(initial_value=0.055,trainable=True, name='beta_sq')

	optimizer = tf.keras.optimizers.Adam()

	chi_sq = tf.convert_to_tensor(chi_sq_stats, dtype=tf.float32)


	# Lopp through windows
	prev_est_alpha_sq = 10000000

	converged = False
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			#loss_value = mediated_squared_effect_est_loss_function_unobserved_true(chi_sq, var_ld_scores, gene_ld_scores, beta_sq_var, alpha_sq_var)
			loss_value = NCX2_snp_h2_loss(chi_sq, var_ld_scores, gene_ld_scores, beta_sq_var, alpha_sq_var)

		trainable_variables = []
		trainable_variables.append(alpha_sq_var)
		trainable_variables.append(beta_sq_var)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))

		print(alpha_sq_var)
		print(beta_sq_var)
		print(loss_value)

	pdb.set_trace()


	return np.asarray(alpha_sq_var)[0,0]


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

def NCX2_snp_h2_loss(chi_sq_stat, var_ld_scores, gene_ld_scores, beta_sq_var, alpha_sq_var):
	# Prep data for NCX2
	df = 1.0
	nc = var_ld_scores*beta_sq_var + gene_ld_scores*alpha_sq_var

	ncx2_log_probz = NCX2_log_prob(chi_sq_stat, df, nc)


	return -tf.math.reduce_sum(ncx2_log_probz)

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


gene_ld_scores = np.loadtxt(output_root + '_sum_gene_ld_scores.txt', dtype=str, delimiter='\t').astype(float)
var_ld_scores = np.ones(len(gene_ld_scores))
chi_sq_stats = np.square(np.loadtxt(output_root + '_gwas_z.txt', dtype=str, delimiter='\t').astype(float))

custom_tf_odr_mediated_squared_effect_est_unobserved_true(chi_sq_stats, var_ld_scores, gene_ld_scores)

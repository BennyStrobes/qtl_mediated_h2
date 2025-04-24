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
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != gene_name:
			continue
		sumstats.append(data[4])
		cis_snps.append(data[3])
		sumstat_ses.append(data[5])
		window_name = data[2]
		window_names.append(window_name)
		rsids.append(data[1])
	f.close()

	unique_window_names = np.unique(window_names)
	if len(unique_window_names) != 1:
		print('assumption eroror')
		pdb.set_trace()

	gene_window_name = unique_window_names[0]

	return np.asarray(sumstats).astype(float), np.asarray(sumstat_ses).astype(float), np.asarray(cis_snps).astype(float), gene_window_name, np.asarray(rsids)


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
	beta_se_arr = []
	ldscore_arr = []
	index_arr = []
	cis_snp_index_arr = []
	n_cis_snp_arr = []
	# Now loop through genes
	for gene_name in gene_names:
		# Extract summary stats for specific gene
		eqtl_gene_beta, eqtl_gene_beta_se, gene_cis_snp_indices, gene_window_name, gene_rsids = extract_eqtl_sumstats_for_specific_gene(eqtl_sumstat_file, gene_name)

		# Load in q_mat and w_mat for this window
		ld_file = window_name_to_ld_files[gene_window_name]
		ld_mat = np.load(ld_file)

		squared_ld_mat = np.square(ld_mat)
		squared_adj_ld_mat = squared_ld_mat - ((1.0 - squared_ld_mat)/(100000.0-2.0))
		
		eqtl_ld_scores = np.sum(squared_adj_ld_mat[gene_cis_snp_indices==1,:], axis=0)



		# Get names of snps in gene
		gene_snp_positions = []
		for gene_rsid in gene_rsids:
			gene_snp_positions.append(snp_name_to_position[gene_rsid])
		gene_snp_positions = np.asarray(gene_snp_positions)

		# Add to global array
		gene_arr.append(gene_name)
		beta_arr.append(eqtl_gene_beta)
		beta_se_arr.append(eqtl_gene_beta_se)
		ldscore_arr.append(eqtl_ld_scores)
		index_arr.append(gene_snp_positions)
		cis_snp_index_arr.append(gene_snp_positions[gene_cis_snp_indices==1])
		n_cis_snp_arr.append(np.sum(gene_cis_snp_indices))

	return np.asarray(gene_arr), beta_arr, beta_se_arr, ldscore_arr, index_arr, cis_snp_index_arr, np.asarray(n_cis_snp_arr)


def get_matrix_form_of_eqtl_data(input_data_arr, input_data_positions, null_value, n_snps):
	n_genes = len(input_data_arr)
	output_mat = np.zeros((n_genes, n_snps)) + null_value


	for gene_iter, gene_arr in enumerate(input_data_arr):
		gene_input_position = input_data_positions[gene_iter]

		output_mat[gene_iter, gene_input_position] = gene_arr


	return output_mat


def sumstat_proportional_joint_reml_loss(gwas_beta, gwas_beta_var, gwas_ld_scores, eqtl_beta, eqtl_beta_var, eqtl_ld_scores,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, per_snp_eqtl_h2_est_lb, per_snp_gwas_h2_est_lb, gwas_noise, eqtl_noise):

	# Currently working
	#eqtl_beta_sq_variable_scaled = tf.nn.relu(eqtl_beta_sq_variable) + (per_snp_eqtl_h2_est_lb)
	#big_eqtl_beta_sq_variable = tf.tile(tf.reshape(eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps])

	eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable + (per_snp_eqtl_h2_est_lb)

	big_eqtl_beta_sq_variable = eqtl_ld_scores*(tf.tile(tf.reshape(eqtl_beta_sq_variable_scaled, [-1, 1]), [1,n_snps]))


	pred_per_snp_gwas_h2 = (gwas_ld_scores*beta_sq_variable) + alpha_sq_variable*(tf.math.reduce_sum((big_eqtl_beta_sq_variable)*eqtl_mask,axis=0))

	#thresholded_pred_per_snp_gwas_h2 = tf.maximum(pred_per_snp_gwas_h2, per_snp_gwas_h2_est_lb)
	#pred_per_snp_gwas_h2[pred_per_snp_gwas_h2 < per_snp_gwas_h2_est_lb] = per_snp_gwas_h2_est_lb

	gwas_log_like = compute_univariate_gaussian_log_like(gwas_beta, 0.0, gwas_beta_var + gwas_noise + pred_per_snp_gwas_h2)


	big_eqtl_noise = tf.tile(tf.reshape(eqtl_noise, [-1, 1]), [1,n_snps])

	eqtl_indices = eqtl_mask==1.0
	eqtl_log_like = compute_univariate_gaussian_log_like(eqtl_beta[eqtl_indices], 0.0, eqtl_beta_var[eqtl_indices] + big_eqtl_noise[eqtl_indices] + big_eqtl_beta_sq_variable[eqtl_indices])


	loss = -tf.reduce_sum(gwas_log_like) - tf.reduce_sum(eqtl_log_like)

	return loss

def compute_univariate_gaussian_log_like(xx, mean_value, variance_vec):
	log_like = -(tf.math.log(variance_vec)/2) - tf.math.divide(tf.square(xx-mean_value), (2.0*variance_vec))
	
	return log_like


def med_h2_with_sumstat_reml(gwas_beta, gwas_beta_se, gwas_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ld_scores, eqtl_n_cis_snps, max_epochs=40000, conv_thresh=1e-12, intercept_variables=False, all_intercepts=True):
	# dimensionality of system
	n_genes = eqtl_beta.shape[0]
	n_snps = len(gwas_beta)


	optimizer = tf.keras.optimizers.Adam(learning_rate=7.5e-6)

	# Create mask matrix
	eqtl_mask = 1.0*(eqtl_beta!=0.0)
	snps_per_gene_arr = np.sum(eqtl_mask!=0.0,axis=1)


	# Convert variabels to tf tensors
	gwas_beta = tf.convert_to_tensor(gwas_beta, dtype=tf.float32)
	gwas_beta_var = tf.convert_to_tensor(np.square(gwas_beta_se), dtype=tf.float32)
	gwas_ld_scores = tf.convert_to_tensor(gwas_ld_scores, dtype=tf.float32)
	eqtl_beta = tf.convert_to_tensor(eqtl_beta, dtype=tf.float32)
	eqtl_beta_var = tf.convert_to_tensor(np.square(eqtl_beta_se), dtype=tf.float32)
	eqtl_mask = tf.convert_to_tensor(np.square(eqtl_mask), dtype=tf.float32)
	eqtl_ld_scores = tf.convert_to_tensor(eqtl_ld_scores, dtype=tf.float32)

	per_snp_eqtl_h2_est_lb = -np.min(eqtl_beta_var[eqtl_mask==1])*.99
	per_snp_gwas_h2_est_lb = -np.min(gwas_beta_var)*.99


	# Initialize variables to optimize over
	beta_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='beta_sq')
	alpha_sq_variable = tf.Variable(initial_value=0.0000001,trainable=True, name='alpha_sq')
	eqtl_beta_sq_variable = tf.Variable(initial_value=np.ones(n_genes)*0.0 - per_snp_eqtl_h2_est_lb,trainable=True, name='eqtl_beta_sq', dtype=tf.float32)

	gwas_noise = tf.Variable(initial_value=0.0000000,trainable=True, name='gwas_noise')
	eqtl_noise = tf.Variable(initial_value=np.zeros(n_genes),trainable=True, name='eqtl_noise', dtype=tf.float32)

	converged = False
	prev_est_alpha_sq=10000
	best_loss = 1e10
	for epoch_iter in range(max_epochs):
		# Use tf.gradient tape to compute gradients
		with tf.GradientTape() as tape:
			loss_value = sumstat_proportional_joint_reml_loss(gwas_beta, gwas_beta_var,gwas_ld_scores, eqtl_beta, eqtl_beta_var,eqtl_ld_scores,eqtl_mask, alpha_sq_variable, eqtl_beta_sq_variable, beta_sq_variable, n_snps, per_snp_eqtl_h2_est_lb, per_snp_gwas_h2_est_lb, gwas_noise, eqtl_noise)

		trainable_variables = []
		trainable_variables.append(alpha_sq_variable)
		trainable_variables.append(eqtl_beta_sq_variable)
		trainable_variables.append(beta_sq_variable)
		if intercept_variables:
			trainable_variables.append(gwas_noise)
			trainable_variables.append(eqtl_noise)

		grads = tape.gradient(loss_value, trainable_variables)
		optimizer.apply_gradients(zip(grads, trainable_variables))


		cur_est = np.asarray(alpha_sq_variable)*1.0


		prev_est_alpha_sq = cur_est

		cur_loss = np.asmatrix(loss_value)[0,0]
		if cur_loss < best_loss:
			best_loss = cur_loss
			eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable + (per_snp_eqtl_h2_est_lb)
			best_gene_cis_h2 = np.asarray(eqtl_beta_sq_variable_scaled)*eqtl_n_cis_snps
			best_med_h2 = np.sum(alpha_sq_variable*best_gene_cis_h2)
			best_nm_h2 = np.sum(beta_sq_variable)*n_snps			


		if np.mod(epoch_iter, 100) == 0.0:
			eqtl_beta_sq_variable_scaled = eqtl_beta_sq_variable + (per_snp_eqtl_h2_est_lb)
			gene_cis_h2 = np.asarray(eqtl_beta_sq_variable_scaled)*eqtl_n_cis_snps
			med_h2 = np.sum(alpha_sq_variable*gene_cis_h2)
			nm_h2 = np.sum(beta_sq_variable)*n_snps
			print(epoch_iter)
			print('loss: ' + str(loss_value))
			print('med: ' + str(med_h2))
			print('nm: ' + str(nm_h2))
			print('eqtl: ' + str(np.mean(gene_cis_h2)))
			print('gwas_noise: ' + str(np.asmatrix(gwas_noise)[0,0]))
			print('eqtl noise: ' + str(np.mean(eqtl_noise)))
			print(np.sort(eqtl_noise))

	print('best loss: ' + str(best_loss))




	return best_med_h2, best_nm_h2, best_gene_cis_h2




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
N_gwas = n_gwas_individuals

# Get variant ld scores
#quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_full_space_ld_ld_summary.txt'
quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_summary.txt'
#quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_summary.txt'
gwas_variant_ld_scores, window_to_ld_files = get_gwas_variant_ld_scores(gwas_beta, gwas_rsids, quasi_ld_window_summary_file)

# Create mapping from rsid to position
snp_name_to_position = create_mapping_from_rsid_to_position(gwas_rsids)

# load in eqtl data
# Out of sample eqtl ld
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
#eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_big_window_eqtl_sumstats.txt'
genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, eqtl_cis_snp_position, eqtl_n_cis_snps = load_in_eqtl_data(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, eqtl_ld='out_of_sample')


# Generate matrix form of eqtl data
eqtl_beta_mat = get_matrix_form_of_eqtl_data(eqtl_beta, eqtl_position, 0.0, len(gwas_variant_ld_scores))
eqtl_beta_se_mat = get_matrix_form_of_eqtl_data(eqtl_beta_se, eqtl_position, 1.0/np.sqrt(eqtl_sample_size), len(gwas_variant_ld_scores))
eqtl_ld_score_mat = get_matrix_form_of_eqtl_data(eqtl_ldscore, eqtl_position, 0.0, len(gwas_variant_ld_scores))


output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + str(eqtl_sample_size) + '_joint_reml_w_all_intercepts.txt'
t = open(output_file,'w')
t.write('eQTL_SS\tsim_h2\tsim_med_h2\tsim_nm_h2\test_med_h2_joint_reml\test_nm_h2_joint_reml\test_mean_eqtl_h2_joint_reml\test_h2_ldsc\test_h2_ldsc_constrained_intercept\n')

# Run standard S-LDSC
X = sm.add_constant(gwas_variant_ld_scores)
Y = np.square(gwas_beta/gwas_beta_se)
model = sm.OLS(Y,X).fit()
model_constrained_intercept = sm.OLS(Y-1, gwas_variant_ld_scores).fit()
ldsc_snp_h2 = len(gwas_rsids)*(model.params[1]/N_gwas)
ldsc_constrained_intercept_snp_h2 = len(gwas_rsids)*(model_constrained_intercept.params[0]/N_gwas)



est_med_h2, est_nm_h2, est_gene_cis_h2 = med_h2_with_sumstat_reml(gwas_beta, gwas_beta_se, gwas_variant_ld_scores, eqtl_beta_mat, eqtl_beta_se_mat, eqtl_ld_score_mat, eqtl_n_cis_snps, max_epochs=70000, intercept_variables=True, all_intercepts=True)

t.write(str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(est_med_h2) + '\t' + str(est_nm_h2) + '\t' + str(np.mean(est_gene_cis_h2)) + '\t' + str(ldsc_snp_h2) + '\t' + str(ldsc_constrained_intercept_snp_h2) + '\n')

t.close()
print(output_file)


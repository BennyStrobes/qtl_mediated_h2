import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import bayesian_lmm_rss_med_h2





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


def extract_gwas_info_for_each_window(gwas_beta, gwas_rsids, quasi_ld_window_summary_file):
	window_names = []
	window_info = {}
	pc_sumstats = []
	f = open(quasi_ld_window_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_snp_indices = np.load(data[4])
		window_rsids = np.load(data[5])
		# QUick error checking
		if np.array_equal(window_rsids, gwas_rsids[window_snp_indices]) == False:
			print('assumption eroror')
			pdb.set_trace()
		window_name = data[0]
		q_mat_file = data[7]
		w_premult_mat_file = data[8]
		Q_mat = np.load(q_mat_file)
		w_premult_mat = np.load(w_premult_mat_file)


		window_pc_gwas_beta = np.dot(w_premult_mat, gwas_beta[window_snp_indices])
		pc_sumstats.append(window_pc_gwas_beta)


		# Add to global array
		window_names.append(window_name)
		# Quick error check
		if window_name in window_info:
			print('assumption error: window already seen')
			pdb.set_trace()

		window_info[window_name] = {}
		window_info[window_name]['beta_pc'] = np.copy(window_pc_gwas_beta)
		window_info[window_name]['Q_file'] = q_mat_file
		window_info[window_name]['W_file'] = w_premult_mat_file
		window_info[window_name]['n_snps'] = Q_mat.shape[1]
		window_info[window_name]['rsids'] = window_rsids
		window_info[window_name]['genes'] = []


	f.close()

	return np.asarray(window_names), window_info, np.hstack(pc_sumstats)


def extract_eqtl_sumstats_for_specific_gene(sumstats_file, gene_name):
	f = open(sumstats_file)
	sumstats = []
	cis_snps = []
	window_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != gene_name:
			continue
		sumstats.append(data[4])
		cis_snps.append(data[3])
		window_name = data[2]
		window_names.append(window_name)
	f.close()

	unique_window_names = np.unique(window_names)
	if len(unique_window_names) != 1:
		print('assumption eroror')
		pdb.set_trace()

	gene_window_name = unique_window_names[0]

	return np.asarray(sumstats).astype(float), np.asarray(cis_snps).astype(float), gene_window_name


def load_in_eqtl_data(eqtl_sumstat_file, window_info, window_names):
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


	# Now loop through genes
	gene_info = {}
	for gene_name in gene_names:
		# Extract summary stats for specific gene
		eqtl_gene_beta, gene_cis_snp_indices, gene_window_name = extract_eqtl_sumstats_for_specific_gene(eqtl_sumstat_file, gene_name)

		# Load in q_mat and w_mat for this window
		q_mat_file= window_info[gene_window_name]['Q_file']
		w_mat_file = window_info[gene_window_name]['W_file']
		q_mat = np.load(q_mat_file)
		w_mat = np.load(w_mat_file)

		# Get PC space ss
		gene_pc_beta = np.dot(w_mat, eqtl_gene_beta)

		# Add gene name to window info
		window_info[gene_window_name]['genes'].append(gene_name)

		# Add to gene_info data structure
		if gene_name in gene_info:
			print('assumption eroror')
			pdb.set_trace()
		gene_info[gene_name] = {}
		gene_info[gene_name]['beta_pc'] = np.copy(gene_pc_beta)
		gene_info[gene_name]['cis_snps'] = np.copy(gene_cis_snp_indices) == 1.0
		gene_info[gene_name]['n_cis_snps'] = np.sum(gene_cis_snp_indices)
		gene_info[gene_name]['window_name'] = gene_window_name

		# Error check 
		if q_mat.shape[1] != len(gene_cis_snp_indices):
			print('assumption eroror')
			pdb.set_trace()

	# Organize into np array
	for window_name in window_names:
		window_info[window_name]['genes'] = np.asarray(window_info[window_name]['genes'])


	return window_info, gene_info, gene_names





######################
# Command line args
######################
simulation_number = int(sys.argv[1])
simulation_name_string = sys.argv[2]
simulated_trait_dir = sys.argv[3]
simulated_gwas_dir = sys.argv[4]
simulation_genotype_dir = sys.argv[5]
simulated_learned_gene_models_dir = sys.argv[6]
N_gwas = int(sys.argv[7])
N_eqtl = int(sys.argv[8])
trait_med_h2_inference_dir = sys.argv[9]
resid_var_bool = sys.argv[10]


# Load in true simulated data parameters
genetic_trait_expr_med_file = simulated_trait_dir + simulation_name_string +'_expression_mediated_trait_values.txt'
genetic_trait_nm_file = simulated_trait_dir +simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'
sim_med_h2 = np.var(np.loadtxt(genetic_trait_expr_med_file))
sim_nm_h2 = np.var(np.loadtxt(genetic_trait_nm_file))
sim_h2 = np.var(np.loadtxt(genetic_trait_nm_file) + np.loadtxt(genetic_trait_expr_med_file))


# Load in GWAS summary statistics
'''
gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results.txt'
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(gwas_summary_file)

# Extract information in each window
quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_summary.txt' 
window_names, window_info, gwas_pc_beta = extract_gwas_info_for_each_window(gwas_beta, gwas_rsids, quasi_ld_window_summary_file)



# load in eqtl data
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(N_eqtl) + '_eqtl_sumstats.txt'
window_info, gene_info, gene_names = load_in_eqtl_data(eqtl_sumstat_file, window_info, window_names)

import pickle
f = open('window_info.pkl', 'wb')
pickle.dump(window_info, f)
f.close()
f = open('gene_info.pkl', 'wb')
pickle.dump(gene_info, f)
f.close()
'''
import pickle
f = open('window_info.pkl', 'rb')
window_info = pickle.load(f) 
f.close()
f = open('gene_info.pkl', 'rb')
gene_info = pickle.load(f) 
f.close()

if resid_var_bool == 'True':
	tmp_output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_' + str(N_eqtl) + '_resid_var_variable_rss_tmp_res.txt'
	mod = bayesian_lmm_rss_med_h2.Bayesian_LMM_RSS_med_h2_inference(window_info, gene_info, N_gwas, N_eqtl, tmp_output_file)
	mod.fit(burn_in_iterations=1, total_iterations=10000, update_resid_var_bool=True)
elif resid_var_bool == 'False':
	tmp_output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_' + str(N_eqtl) + '_resid_var_const_rss_tmp_res.txt'
	mod = bayesian_lmm_rss_med_h2.Bayesian_LMM_RSS_med_h2_inference(window_info, gene_info, N_gwas, N_eqtl, tmp_output_file)
	mod.fit(burn_in_iterations=1, total_iterations=10000, update_resid_var_bool=False)




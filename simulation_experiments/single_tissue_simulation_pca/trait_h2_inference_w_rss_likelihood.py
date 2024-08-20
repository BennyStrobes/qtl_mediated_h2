import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import bayesian_lmm_rss_h2


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


		# Add to global array
		window_names.append(window_name)
		# Quick error check
		if window_name in window_info:
			print('assumption error: window already seen')
			pdb.set_trace()

		window_info[window_name] = {}
		window_info[window_name]['beta_pc'] = np.copy(window_pc_gwas_beta)
		window_info[window_name]['Q_file'] = q_mat_file
		window_info[window_name]['n_snps'] = Q_mat.shape[1]

	f.close()

	return np.asarray(window_names), window_info








#####################
# Command line args
#####################
simulation_number = int(sys.argv[1])
simulation_name_string = sys.argv[2]
simulated_trait_dir = sys.argv[3]
simulated_gwas_dir = sys.argv[4]
simulation_genotype_dir = sys.argv[5]
n_gwas_individuals = int(sys.argv[6])
trait_h2_inference_dir = sys.argv[7]




####################
# Load in data
####################
# Load in true simulated data parameters
genetic_trait_expr_med_file = simulated_trait_dir + simulation_name_string +'_expression_mediated_trait_values.txt'
genetic_trait_nm_file = simulated_trait_dir +simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'
sim_med_h2 = np.var(np.loadtxt(genetic_trait_expr_med_file))
sim_nm_h2 = np.var(np.loadtxt(genetic_trait_nm_file))
sim_h2 = np.var(np.loadtxt(genetic_trait_nm_file) + np.loadtxt(genetic_trait_expr_med_file))


# Load in GWAS summary statistics
gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results.txt'
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(gwas_summary_file)
N_gwas = n_gwas_individuals


# Extract information in each window
quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_summary.txt' 

window_names, window_info = extract_gwas_info_for_each_window(gwas_beta, gwas_rsids, quasi_ld_window_summary_file)


mod = bayesian_lmm_rss_h2.Bayesian_LMM_RSS_h2_inference(window_names, window_info, n_gwas_individuals)
mod.fit(burn_in_iterations=1000, total_iterations=2000)

# Print sampling results to output file
output_file = trait_h2_inference_dir + simulation_name_string + '_pca_rss_likelihood_sampling_results.txt'
t = open(output_file,'w')
t.write('sampling_iter\tobs_h2\test_h2\test_resid_var\n')
for ii, sampled_h2 in enumerate(mod.sampled_h2):
	t.write(str(ii) + '\t' + str(sim_h2) + '\t' + str(sampled_h2) + '\t' + str(mod.sampled_resid_var[ii]) + '\n')
t.close()




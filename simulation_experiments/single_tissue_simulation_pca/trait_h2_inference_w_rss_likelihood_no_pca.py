import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import bayesian_lmm_rss_med_h2
import bayesian_lmm_rss_gibbs_h2_no_pca
import bayesian_lmm_rss_VI_h2_no_pca


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

		ld_file = data[3]
		#eqtl_ld_file = ld_file.split('geno_gwas')[0] + 'geno_eqtl_' + str(N_eqtl) + ld_file.split('geno_gwas')[1]



		# Add to global array
		window_names.append(window_name)
		# Quick error check
		if window_name in window_info:
			print('assumption error: window already seen')
			pdb.set_trace()

		window_info[window_name] = {}
		window_info[window_name]['beta'] = np.copy(gwas_beta[window_snp_indices])
		window_info[window_name]['ld_file'] = ld_file
		#window_info[window_name]['eqtl_ld_file'] = eqtl_ld_file
		window_info[window_name]['n_snps'] = len(gwas_beta[window_snp_indices])
		window_info[window_name]['rsids'] = window_rsids
		window_info[window_name]['genes'] = []


	f.close()

	return np.asarray(window_names), window_info


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
	count = 0

	# Extract summary stats for specific gene
	#eqtl_gene_beta, gene_cis_snp_indices, gene_window_name = extract_eqtl_sumstats_for_specific_gene(eqtl_sumstat_file, gene_name)
	f = open(eqtl_sumstat_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count +1
			continue
		gene_name = data[0]
		sumstat = float(data[4])
		sumstat_se = float(data[5])
		cis_snp = float(data[3])
		window_name = data[2]

		if gene_name not in gene_info:
			gene_info[gene_name] = {}
			gene_info[gene_name]['beta'] = []
			gene_info[gene_name]['cis_snps'] = []
			gene_info[gene_name]['n_cis_snps'] = 0
			gene_info[gene_name]['window_name'] = window_name
			gene_info[gene_name]['beta_se'] = []

		gene_info[gene_name]['beta'].append(sumstat)
		gene_info[gene_name]['beta_se'].append(sumstat_se)
		gene_info[gene_name]['cis_snps'].append(cis_snp)
		if gene_info[gene_name]['window_name'] != window_name:
			print('assumption eroror')
			pdb.set_trace()

	f.close()


	for gene in [*gene_info]:
		# Organize arrays
		gene_info[gene]['beta'] = np.asarray(gene_info[gene]['beta'])
		gene_info[gene]['beta_se'] = np.asarray(gene_info[gene]['beta_se'])
		gene_info[gene]['cis_snps'] = np.asarray(gene_info[gene]['cis_snps']) == 1.0
		gene_info[gene]['n_cis_snps'] = np.sum(gene_info[gene]['cis_snps'])

		# add gene to window
		window_info[gene_info[gene]['window_name']]['genes'].append(gene)


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
N_gwas = int(sys.argv[6])
trait_h2_inference_dir = sys.argv[7]
window_version = sys.argv[8]


# Load in true simulated data parameters
genetic_trait_expr_med_file = simulated_trait_dir + simulation_name_string +'_expression_mediated_trait_values.txt'
genetic_trait_nm_file = simulated_trait_dir +simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'
sim_med_h2 = np.var(np.loadtxt(genetic_trait_expr_med_file))
sim_nm_h2 = np.var(np.loadtxt(genetic_trait_nm_file))
sim_h2 = np.var(np.loadtxt(genetic_trait_nm_file) + np.loadtxt(genetic_trait_expr_med_file))

print(sim_nm_h2)
print(sim_med_h2)
# Load in GWAS summary statistics
gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results.txt'
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(gwas_summary_file)

# Extract information in each window
if window_version == 'small':
	quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_summary.txt' 
else:
	quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_summary.txt' 
window_names, window_info = extract_gwas_info_for_each_window(gwas_beta, gwas_rsids, quasi_ld_window_summary_file)


output_file = trait_h2_inference_dir + simulation_name_string + '_' + window_version + '_trait_h2_inference_summary.txt'
print(output_file)
t = open(output_file,'w')
t.write('method\tsim_h2_v1\tsim_h2_v2\th2_est_1\th2_est2\n')

prior_parameters = [0.0, 1e-10, 1e-6, 1e-3]

for prior_parameter in prior_parameters:

	# VI (univariate updates)
	mod = bayesian_lmm_rss_VI_h2_no_pca.Bayesian_LMM_RSS_h2_inference(window_info, N_gwas)
	mod.fit(total_iterations=500, update_resid_var_bool=False, univariate_updates=False,cc=prior_parameter)
	t.write('VI_rss_lmm_multivariate_ig_prior_' + str(prior_parameter) + '\t' + str(sim_h2) + '\t' + str(sim_nm_h2 + sim_med_h2) + '\t' + str(mod.h2_independent) + '\t' + str(mod.h2_full) + '\n')


	# VI (univariate updates)
	mod = bayesian_lmm_rss_VI_h2_no_pca.Bayesian_LMM_RSS_h2_inference(window_info, N_gwas)
	mod.fit(total_iterations=500, update_resid_var_bool=False, univariate_updates=True,cc=prior_parameter)
	t.write('VI_rss_lmm_univariate_ig_prior_' + str(prior_parameter) + '\t' + str(sim_h2) + '\t' + str(sim_nm_h2 + sim_med_h2) + '\t' + str(mod.h2_independent) + '\t' + str(mod.h2_full) + '\n')

	# Gibbs sampler (multivariate updates)
	try:
		mod = bayesian_lmm_rss_gibbs_h2_no_pca.Bayesian_LMM_RSS_h2_inference(window_info, N_gwas)
		mod.fit(burn_in_iterations=250, total_iterations=500, update_resid_var_bool=False, univariate_updates=False,cc=prior_parameter)
		t.write('gibbs_rss_lmm_multivariate_ig_prior_' + str(prior_parameter) + '\t' + str(sim_h2) + '\t' + str(sim_nm_h2 + sim_med_h2) + '\t' + str(np.mean(mod.sampled_nm_h2_independent)) + '\t' + str(np.mean(mod.sampled_nm_h2_full)) + '\n')
	except:
		t.write('gibbs_rss_lmm_multivariate_ig_prior_' + str(prior_parameter) + '\t' + str(sim_h2) + '\t' + str(sim_nm_h2 + sim_med_h2) + '\t' + 'NA' + '\t' + 'NA' + '\n')

	# Gibbs sampler (univariate updates)
	mod = bayesian_lmm_rss_gibbs_h2_no_pca.Bayesian_LMM_RSS_h2_inference(window_info, N_gwas)
	mod.fit(burn_in_iterations=250, total_iterations=500, update_resid_var_bool=False, univariate_updates=True,cc=prior_parameter)
	t.write('gibbs_rss_lmm_univariate_ig_prior_' + str(prior_parameter) + '\t' + str(sim_h2) + '\t' + str(sim_nm_h2 + sim_med_h2) + '\t' + str(np.mean(mod.sampled_nm_h2_independent)) + '\t' + str(np.mean(mod.sampled_nm_h2_full)) + '\n')


	t.flush()

t.close()

import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import bayesian_lmm_ss_h2
import bayesian_lmm_pca_ss_h2_single_region


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


def load_in_variant_ld_scores(var_ld_score_file):
	rsids = []
	ldscores = []

	head_count = 0

	f = open(var_ld_score_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count = head_count + 1
			continue
		rsids.append(data[0])
		ldscores.append(float(data[1]))

	f.close()
	return np.asarray(rsids), np.asarray(ldscores)


def compute_h2_with_ldsc(gwas_beta, ldscores, n_gwas_individuals, n_snps):
	z_scores = gwas_beta/np.sqrt(1/n_gwas_individuals)

	chi_squared = np.square(z_scores)

	olser = sm.OLS(chi_squared-1, ldscores).fit()

	h2_est = n_snps*(olser.params[0]/n_gwas_individuals)

	h2_est_se = n_snps*(olser.bse[0]/n_gwas_individuals)

	return h2_est, h2_est_se

def convert_to_pc_sum_stat_space(gwas_beta, gwas_rsids, quasi_ld_window_summary_file):
	gwas_beta_pc = []
	var_ld_score_pc = []
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
		Q_mat = np.load(data[7])
		w_premult_mat = np.load(data[8])


		window_pc_gwas_beta = np.dot(w_premult_mat, gwas_beta[window_snp_indices])
		window_pc_ld_scores = np.diag(np.dot(Q_mat, np.transpose(Q_mat)))

		gwas_beta_pc.append(window_pc_gwas_beta)
		var_ld_score_pc.append(window_pc_ld_scores)

	gwas_beta_pc = np.hstack(gwas_beta_pc)
	var_ld_score_pc = np.hstack(var_ld_score_pc)


	return gwas_beta_pc, var_ld_score_pc, len(gwas_beta)

def extract_data_from_gene_summary_file(gene_summary_file):
	f = open(gene_summary_file)
	genes = []
	obs_h2s = []
	est_h2s = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		genes.append(data[0])
		obs_h2s.append(float(data[1]))
		est_h2s.append(float(data[2]))

	f.close()

	return np.asarray(genes), np.asarray(obs_h2s), np.asarray(est_h2s)

def extract_gene_betas(gene_pc_sumstats_file, gene_name):
	f = open(gene_pc_sumstats_file)
	betas = []
	snp_names = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] != gene_name:
			continue
		betas.append(float(data[3]))
		snp_names.append(data[1])
		window_name = data[2]
	f.close()
	return np.asarray(betas), np.asarray(snp_names), window_name

def extract_gene_ldscores(gene_pc_ldscores_file, gene_name):
	f = open(gene_pc_ldscores_file)
	betas = []
	snp_names = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] != gene_name:
			continue
		betas.append(float(data[3]))
		snp_names.append(data[1])
	f.close()
	return np.asarray(betas), np.asarray(snp_names)

def extract_gene_cis_snp_indices(sumstats_file, gene_name):
	f = open(sumstats_file)
	sumstats = []
	cis_snps = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != gene_name:
			continue
		sumstats.append(data[4])
		cis_snps.append(data[3])
	f.close()

	return np.asarray(cis_snps).astype(float), np.asarray(sumstats).astype(float)


#####################
# Command line args
#####################
simulation_number = int(sys.argv[1])
simulation_name_string = sys.argv[2]
simulated_learned_gene_models_dir = sys.argv[3]
simulation_genotype_dir = sys.argv[4]
n_eqtl_individuals = int(sys.argv[5])
trait_h2_inference_dir = sys.argv[6]



####################
# Load in data
####################

# Load in genes
gene_summary_file =simulated_learned_gene_models_dir + simulation_name_string + '_' + str(n_eqtl_individuals) + '_h2_est.txt'
gene_names, gene_obs_h2s, gene_est_h2s = extract_data_from_gene_summary_file(gene_summary_file)


gene_pc_sumstats_file =simulated_learned_gene_models_dir + simulation_name_string + '_' + str(n_eqtl_individuals) + '_pc_space_eqtl_sumstats.txt'
gene_pc_ldscores_file =simulated_learned_gene_models_dir + simulation_name_string + '_' + str(n_eqtl_individuals) + '_eqtl_reference_bias_corrected_ld_scores.txt'
sumstats_file =simulated_learned_gene_models_dir + simulation_name_string + '_' + str(n_eqtl_individuals) + '_eqtl_sumstats.txt'

output_file = trait_h2_inference_dir + simulation_name_string + '_eqtl_' + str(n_eqtl_individuals) + '_h2_est.txt'
t = open(output_file,'w')
t.write('gene\tmethod\tobserved_h2\test_h2\test_h2_lb\test_h2_ub\n')

# Loop through genes
for gene_iter, gene_name in enumerate(gene_names):
	gene_obs_h2 = gene_obs_h2s[gene_iter]

	gene_pc_beta, gene_pc_snp_names, window_name = extract_gene_betas(gene_pc_sumstats_file, gene_name)
	gene_pc_ldscores, gene_pc_snp_names2 = extract_gene_ldscores(gene_pc_ldscores_file, gene_name)
	gene_cis_snp_indices, gene_beta = extract_gene_cis_snp_indices(sumstats_file, gene_name)

	# Load in Q-mat
	#q_mat_file = simulation_genotype_dir + 'variant_ref_geno_eqtl_300_quasi_independent_windows_ld_' + window_name + '_Q_mat.npy'
	q_mat_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_' + window_name + '_Q_mat.npy'
	q_mat = np.load(q_mat_file)
	# Load in w-premult mat
	w_mat_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_' + window_name + '_w_premult_mat.npy'
	w_mat = np.load(w_mat_file)
	
	gene_pc_beta = np.dot(w_mat, gene_beta)


	# Quick error checking
	if np.array_equal(gene_pc_snp_names, gene_pc_snp_names2) == False:
		print('assumption error')
		pdb.set_trace()

	# Est h2 with LDSC
	gene_pc_ldscores = np.diag(np.dot(q_mat[:, gene_cis_snp_indices==1], np.transpose(q_mat[:, gene_cis_snp_indices==1])))
	#mod = sm.OLS(np.square(gene_pc_beta) - (1.0/n_eqtl_individuals), gene_pc_ldscores).fit()
	#ldsc_h2_est = mod.params[0]*np.sum(gene_pc_ldscores)


	mod = bayesian_lmm_pca_ss_h2_single_region.Bayesian_LMM(gene_pc_beta, q_mat[:, gene_cis_snp_indices==1], n_eqtl_individuals)
	mod.fit(burn_in_iterations=3000, total_iterations=5000,update_resid_var_bool=True)
	bayesian_lmm_h2_mean = np.mean(mod.sampled_h2)
	bayesian_lmm_h2_lb = np.sort(mod.sampled_h2)[int(np.round(len(mod.sampled_h2)*.025))]
	bayesian_lmm_h2_ub = np.sort(mod.sampled_h2)[int(np.round(len(mod.sampled_h2)*.975))]


	t.write(gene_name + '\t' + 'bayesian_lmm_pca\t' + str(gene_obs_h2) + '\t' + str(bayesian_lmm_h2_mean) + '\t' + str(bayesian_lmm_h2_lb) + '\t' + str(bayesian_lmm_h2_ub) + '\n')


	gene_pc_ldscores = np.diag(np.dot(q_mat[:, gene_cis_snp_indices==1], np.transpose(q_mat[:, gene_cis_snp_indices==1])))
	n_snps = int(np.sum(gene_cis_snp_indices==1))
	mod = bayesian_lmm_ss_h2.Bayesian_LMM_SS_h2_inference(gene_pc_beta, gene_pc_ldscores, n_eqtl_individuals, n_snps)
	mod.fit(burn_in_iterations=5000, total_iterations=10000,update_resid_var_bool=True)
	bayesian_lmm_marginal_h2_mean = np.mean(mod.sampled_h2)
	bayesian_lmm_marginal_h2_lb = np.sort(mod.sampled_h2)[int(np.round(len(mod.sampled_h2)*.025))]
	bayesian_lmm_marginal_h2_ub = np.sort(mod.sampled_h2)[int(np.round(len(mod.sampled_h2)*.975))]

	t.write(gene_name + '\t' + 'bayesian_lmm_marginal_lv_pca\t' + str(gene_obs_h2) + '\t' + str(bayesian_lmm_marginal_h2_mean) + '\t' + str(bayesian_lmm_marginal_h2_lb) + '\t' + str(bayesian_lmm_marginal_h2_ub) + '\n')

	t.flush()
t.close()

'''
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

# Convert gwas_beta to gwas_beta_pc and get var_pc_ld_scores
quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_summary.txt' 
gwas_beta_pc, var_ld_score_pc, n_snps = convert_to_pc_sum_stat_space(gwas_beta, gwas_rsids, quasi_ld_window_summary_file)




####################
# Compute h2 using ldscore regression
####################
ldsc_h2, ldsc_h2_se = compute_h2_with_ldsc(gwas_beta_pc, var_ld_score_pc, N_gwas, n_snps)
ldsc_h2_lb = ldsc_h2 - (1.96*ldsc_h2_se)
ldsc_h2_ub = ldsc_h2 + (1.96*ldsc_h2_se)

####################
# Compute h2 using bayesian lmm (with resid resid var learning)
####################
# Equal weights
mod = bayesian_lmm_ss_h2.Bayesian_LMM_SS_h2_inference(gwas_beta_pc, var_ld_score_pc, N_gwas, n_snps)
mod.fit(burn_in_iterations=20000, total_iterations=30000)
bayesian_lmm_h2_mean = np.mean(mod.sampled_h2)
bayesian_lmm_h2_lb = np.sort(mod.sampled_h2)[int(np.round(len(mod.sampled_h2)*.025))]
bayesian_lmm_h2_ub = np.sort(mod.sampled_h2)[int(np.round(len(mod.sampled_h2)*.975))]
bayesian_lmm_resid_var_mean = np.mean(mod.sampled_resid_vars)
bayesian_lmm_resid_var_lb = np.sort(mod.sampled_resid_vars)[int(np.round(len(mod.sampled_resid_vars)*.025))]
bayesian_lmm_resid_var_ub = np.sort(mod.sampled_resid_vars)[int(np.round(len(mod.sampled_resid_vars)*.975))]


####################
# Compute h2 using bayesian lmm (w/o resid resid var learning)
####################
# Equal weights
mod2 = bayesian_lmm_ss_h2.Bayesian_LMM_SS_h2_inference(gwas_beta_pc, var_ld_score_pc, N_gwas, n_snps)
mod2.fit(burn_in_iterations=20000, total_iterations=30000, update_resid_var_bool=False)
bayesian_lmm_no_resid_h2_mean = np.mean(mod2.sampled_h2)
bayesian_lmm_no_resid_h2_lb = np.sort(mod2.sampled_h2)[int(np.round(len(mod2.sampled_h2)*.025))]
bayesian_lmm_no_resid_h2_ub = np.sort(mod2.sampled_h2)[int(np.round(len(mod2.sampled_h2)*.975))]




####################
# Print results to output file
####################
output_file = trait_h2_inference_dir + simulation_name_string + '_h2_inference_summary.txt'
t = open(output_file,'w')
# header
t.write('method\tsim_h2\test_h2\test_h2_95_lb\test_h2_95_ub\n')
# Results
t.write('LDSC\t' + str(sim_h2) + '\t' + str(ldsc_h2) + '\t' + str(ldsc_h2_lb) + '\t' + str(ldsc_h2_ub) + '\n')
t.write('bayesian_lmm_learn_resid_var\t' + str(sim_h2) + '\t' + str(bayesian_lmm_h2_mean) + '\t' + str(bayesian_lmm_h2_lb) + '\t' + str(bayesian_lmm_h2_ub) + '\n')
t.write('bayesian_lmm_fixed_resid_var\t' + str(sim_h2) + '\t' + str(bayesian_lmm_no_resid_h2_mean) + '\t' + str(bayesian_lmm_no_resid_h2_lb) + '\t' + str(bayesian_lmm_no_resid_h2_ub) + '\n')

t.close()
'''



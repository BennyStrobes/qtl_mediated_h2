import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import bayesian_lmm_ss_h2



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




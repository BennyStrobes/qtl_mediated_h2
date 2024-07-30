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


def compute_h2_with_ldsc(gwas_beta, gwas_beta_se, ldscores, n_gwas_individuals):
	z_scores = gwas_beta/gwas_beta_se

	chi_squared = np.square(z_scores)

	olser = sm.OLS(chi_squared-1, ldscores).fit()

	h2_est = len(chi_squared)*(olser.params[0]/n_gwas_individuals)

	return h2_est


#####################
# Command line args
#####################
simulation_number = int(sys.argv[1])
simulation_name_string = sys.argv[2]
simulated_trait_dir = sys.argv[3]
simulated_gwas_dir = sys.argv[4]
n_gwas_individuals = int(sys.argv[5])
trait_h2_inference_dir = sys.argv[6]

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

# Load in Variant LD-scores
var_ld_score_file = simulated_gwas_dir + simulation_name_string + '_variant_ref_geno_window_2_mb_ld_scores_bias_corrected.txt'
ldscore_rsids, ldscores = load_in_variant_ld_scores(var_ld_score_file)

# Quick error checking
if np.array_equal(gwas_rsids, ldscore_rsids) == False:
	print('rs match assumption eroror')
	pdb.set_trace()


####################
# Compute h2 using ldscore regression
####################
ldsc_h2 = compute_h2_with_ldsc(gwas_beta, gwas_beta_se, ldscores, n_gwas_individuals)


####################
# Compute h2 using bayesian lmm
####################
mod = bayesian_lmm_ss_h2.Bayesian_LMM_SS_h2_inference(gwas_beta, gwas_beta_se, ldscores)
#mod.fit(burn_in_iterations=20000, total_iterations=25000)
mod.fit(burn_in_iterations=15000, total_iterations=20000, gamma_var_update_version='ld_score_weighting')
bayesian_lmm_h2_mean = np.mean(mod.sampled_h2)
bayesian_lmm_h2_lb = np.sort(mod.sampled_h2)[int(np.round(len(mod.sampled_h2)*.025))]
bayesian_lmm_h2_ub = np.sort(mod.sampled_h2)[int(np.round(len(mod.sampled_h2)*.975))]


# Print results to output file
output_file = trait_h2_inference_dir + simulation_name_string + '_h2_inference_summary3.txt'
t = open(output_file,'w')
t.write('sim_h2\tldsc_h2\tbayesian_lmm_h2\tbayesian_lmm_h2_95_lb\tbayesian_lmm_h2_95_ub\n')
t.write(str(sim_h2) + '\t' + str(ldsc_h2) + '\t' + str(bayesian_lmm_h2_mean) + '\t' + str(bayesian_lmm_h2_lb) + '\t' + str(bayesian_lmm_h2_ub) + '\n')
t.close()




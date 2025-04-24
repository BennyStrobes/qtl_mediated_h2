import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import bayesian_lmm_rss_med_h2
import bayesian_lmm_rss_med_h2_no_pca
import bayesian_lmm_ss_h2_single_region
import stan
import bayesian_vi_lmm_ss_h2_single_region


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


def extract_gwas_info_for_each_window(gwas_beta, gwas_rsids, quasi_ld_window_summary_file, N_eqtl):
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
		eqtl_ld_file = ld_file.split('geno_gwas')[0] + 'geno_eqtl_' + str(N_eqtl) + ld_file.split('geno_gwas')[1]



		# Add to global array
		window_names.append(window_name)
		# Quick error check
		if window_name in window_info:
			print('assumption error: window already seen')
			pdb.set_trace()

		window_info[window_name] = {}
		window_info[window_name]['beta'] = np.copy(gwas_beta[window_snp_indices])
		window_info[window_name]['ld_file'] = ld_file
		window_info[window_name]['eqtl_ld_file'] = eqtl_ld_file
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


def load_in_eqtl_data(eqtl_sumstat_file, window_info, window_names, simulated_gene_expression_dir, simulation_name_string):
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

		causal_eqtl_file = simulated_gene_expression_dir + simulation_name_string + '_' + gene + '_causal_eqtl_effects.npy'

		gene_info[gene]['true_beta'] = np.load(causal_eqtl_file)

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

def stan_individual_lmm_h2_inference_resid_var(expr, geno):
	bayes_lmm_code = """
	data {
		int<lower=0> K;         // Number of snps
		int<lower=0> N;         // Sample size
		matrix[N, K] X;         // genotype Matrix
		vector[N] Y;              // expression
	}
	parameters {
		real<lower=0> sigma_b;      // standard deviation causal effcts
		real<lower=0> sigma_e;      // standard deviation on residual variance
		vector[K] beta;          // Causal betas
	}
	model {
		Y ~ normal(X * beta, sigma_e);
		beta ~ normal(0, sigma_b);
	}

	"""
	stan_h2_data = {"K": geno.shape[1],"N": len(expr), "Y": expr, "X": geno}
	posterior = stan.build(bayes_lmm_code, data=stan_h2_data, random_seed=1)
	fit = posterior.sample(num_chains=4, num_samples=1000, num_warmup=3000)

	ind_h2s = (np.square(fit['sigma_b'])*len(beta))[0,:]
	n_samples = fit['beta'].shape[1]
	full_h2s = []
	for sample_iter in range(n_samples):
		beta_slice = fit['beta'][:, sample_iter]
		full_h2 = np.var(np.dot(geno, beta_slice))
		#full_h2 = np.dot(np.dot(beta_slice, eqtl_ld), beta_slice)
		full_h2s.append(full_h2)
	full_h2s = np.asarray(full_h2s)

	return np.mean(full_h2s), np.mean(ind_h2s)



def stan_individual_lmm_h2_inference(expr, geno):
	bayes_lmm_code = """
	data {
		int<lower=0> K;         // Number of snps
		int<lower=0> N;         // Sample size
		matrix[N, K] X;         // genotype Matrix
		vector[N] Y;              // expression
	}
	parameters {
		real<lower=0> sigma_b;      // standard deviation causal effcts
		vector[K] beta;          // Causal betas
	}
	model {
		Y ~ normal(X * beta, 1);
		beta ~ normal(0, sigma_b);
	}

	"""
	stan_h2_data = {"K": geno.shape[1],"N": len(expr), "Y": expr, "X": geno}
	posterior = stan.build(bayes_lmm_code, data=stan_h2_data, random_seed=1)
	fit = posterior.sample(num_chains=4, num_samples=1000, num_warmup=3000)

	ind_h2s = (np.square(fit['sigma_b'])*len(beta))[0,:]
	n_samples = fit['beta'].shape[1]
	full_h2s = []
	for sample_iter in range(n_samples):
		beta_slice = fit['beta'][:, sample_iter]
		full_h2 = np.var(np.dot(geno, beta_slice))
		#full_h2 = np.dot(np.dot(beta_slice, eqtl_ld), beta_slice)
		full_h2s.append(full_h2)
	full_h2s = np.asarray(full_h2s)

	return np.mean(full_h2s), np.mean(ind_h2s)



def stan_individual_lmm_h2_inference_xbeta_prior(expr, geno):
	bayes_lmm_code = """
	data {
		int<lower=0> K;         // Number of snps
		int<lower=0> N;         // Sample size
		matrix[N, K] X;         // genotype Matrix
		vector[N] Y;              // expression
	}
	parameters {
		real<lower=0> sigma_b;      // standard deviation causal effcts
		real<lower=0> sigma_e;
		vector[K] beta;          // Causal betas
	}
	model {
		Y ~ normal(X * beta, sigma_e);
		X * beta ~ normal(0, sigma_b);
	}

	"""
	stan_h2_data = {"K": geno.shape[1],"N": len(expr), "Y": expr, "X": geno}
	posterior = stan.build(bayes_lmm_code, data=stan_h2_data, random_seed=1)
	fit = posterior.sample(num_chains=4, num_samples=1000, num_warmup=5000)

	ind_h2s = (np.square(fit['sigma_b']))[0,:]
	n_samples = fit['beta'].shape[1]
	full_h2s = []
	for sample_iter in range(n_samples):
		beta_slice = fit['beta'][:, sample_iter]
		full_h2 = np.var(np.dot(geno, beta_slice))
		#full_h2 = np.dot(np.dot(beta_slice, eqtl_ld), beta_slice)
		full_h2s.append(full_h2)
	full_h2s = np.asarray(full_h2s)

	return np.mean(full_h2s), np.mean(ind_h2s)

def stan_rss_lmm_h2_inference_halfnormal_prior(beta, eqtl_ld, N_eqtl):

	bayes_lmm_code = """
	data {
		int<lower=0> K;         // Number of snps
		int<lower=0> N;         // Sample size
		matrix[K, K] R;         // LD Matrix
		matrix[K, K] cov;         // LD Matrix
		vector[K] marginal_betas;              // Marginal betas
	}
	parameters {
		real<lower=0> sigma_b;      // standard deviation causal effcts
		vector[K] beta;          // Causal betas
	}
	model {
		marginal_betas ~ multi_normal(R * beta, cov);
		beta ~ normal(0, sigma_b);
		sigma_b ~ student_t(4,0,.1);
	}

	"""
	#stan_h2_data = {"K": len(beta),"N": N_eqtl, "marginal_betas": beta, "R": eqtl_ld, "cov": np.eye(len(beta))/N_eqtl}
	stan_h2_data = {"K": len(beta),"N": N_eqtl, "marginal_betas": beta, "R": eqtl_ld, "cov": (eqtl_ld/N_eqtl) + np.eye(len(beta))*1e-16}
	posterior = stan.build(bayes_lmm_code, data=stan_h2_data, random_seed=1)
	fit = posterior.sample(num_chains=4, num_samples=1000, num_warmup=3000)

	n_samples = fit['beta'].shape[1]
	full_h2s = []
	for sample_iter in range(n_samples):
		beta_slice = fit['beta'][:, sample_iter]
		full_h2 = np.dot(np.dot(beta_slice, eqtl_ld), beta_slice)
		full_h2s.append(full_h2)
	full_h2s = np.asarray(full_h2s)

	#full_h2s = np.diag(np.dot(np.dot(np.transpose(fit['beta']), eqtl_ld), fit['beta']))

	ind_h2s = (np.square(fit['sigma_b'])*len(beta))[0,:]

	return np.mean(full_h2s), np.mean(ind_h2s)


def stan_rss_lmm_h2_inference_IG_prior(beta, eqtl_ld, N_eqtl):

	bayes_lmm_code = """
	data {
		int<lower=0> K;         // Number of snps
		int<lower=0> N;         // Sample size
		matrix[K, K] R;         // LD Matrix
		matrix[K, K] cov;         // LD Matrix
		vector[K] marginal_betas;              // Marginal betas
	}
	parameters {
		real<lower=0> sigma_b;      // standard deviation causal effcts
		vector[K] beta;          // Causal betas
	}
	model {
		marginal_betas ~ multi_normal(R * beta, cov);
		beta ~ normal(0, sqrt(sigma_b));
		sigma_b ~ inv_gamma(0.000001, 0.000001);
	}

	"""
	#stan_h2_data = {"K": len(beta),"N": N_eqtl, "marginal_betas": beta, "R": eqtl_ld, "cov": np.eye(len(beta))/N_eqtl}
	stan_h2_data = {"K": len(beta),"N": N_eqtl, "marginal_betas": beta, "R": eqtl_ld, "cov": (eqtl_ld/N_eqtl) + np.eye(len(beta))*1e-16}
	posterior = stan.build(bayes_lmm_code, data=stan_h2_data, random_seed=1)
	fit = posterior.sample(num_chains=4, num_samples=1000, num_warmup=3000)

	n_samples = fit['beta'].shape[1]
	full_h2s = []
	for sample_iter in range(n_samples):
		beta_slice = fit['beta'][:, sample_iter]
		full_h2 = np.dot(np.dot(beta_slice, eqtl_ld), beta_slice)
		full_h2s.append(full_h2)
	full_h2s = np.asarray(full_h2s)

	#full_h2s = np.diag(np.dot(np.dot(np.transpose(fit['beta']), eqtl_ld), fit['beta']))

	ind_h2s = ((fit['sigma_b'])*len(beta))[0,:]

	return np.mean(full_h2s), np.mean(ind_h2s)

def stan_rss_lmm_h2_inference(beta, eqtl_ld, N_eqtl):

	bayes_lmm_code = """
	data {
		int<lower=0> K;         // Number of snps
		int<lower=0> N;         // Sample size
		matrix[K, K] R;         // LD Matrix
		matrix[K, K] cov;         // LD Matrix
		vector[K] marginal_betas;              // Marginal betas
	}
	parameters {
		real<lower=0> sigma_b;      // standard deviation causal effcts
		vector[K] beta;          // Causal betas
	}
	model {
		marginal_betas ~ multi_normal(R * beta, cov);
		beta ~ normal(0, sigma_b);
	}

	"""
	#stan_h2_data = {"K": len(beta),"N": N_eqtl, "marginal_betas": beta, "R": eqtl_ld, "cov": np.eye(len(beta))/N_eqtl}
	stan_h2_data = {"K": len(beta),"N": N_eqtl, "marginal_betas": beta, "R": eqtl_ld, "cov": (eqtl_ld/N_eqtl) + np.eye(len(beta))*1e-16}
	posterior = stan.build(bayes_lmm_code, data=stan_h2_data, random_seed=1)
	fit = posterior.sample(num_chains=4, num_samples=1000, num_warmup=3000)

	n_samples = fit['beta'].shape[1]
	full_h2s = []
	for sample_iter in range(n_samples):
		beta_slice = fit['beta'][:, sample_iter]
		full_h2 = np.dot(np.dot(beta_slice, eqtl_ld), beta_slice)
		full_h2s.append(full_h2)
	full_h2s = np.asarray(full_h2s)

	#full_h2s = np.diag(np.dot(np.dot(np.transpose(fit['beta']), eqtl_ld), fit['beta']))

	ind_h2s = (np.square(fit['sigma_b'])*len(beta))[0,:]

	return np.mean(full_h2s), np.mean(ind_h2s)

def quick_error_check_for_sum_stat_individual_data_alignment(expr, geno, beta, eqtl_ld):
	# Check LD
	diff = np.corrcoef(np.transpose(geno)) - eqtl_ld
	if np.max(np.abs(diff)) > 1e-10:
		print('assumption error on ld diff')
		pdb.set_trace()

	# Check sum stats
	KK = len(beta)
	for kk in range(KK):
		diff = beta[kk] - sm.OLS(expr, geno[:,kk]).fit().params[0]
		if np.abs(diff) > 1e-7:
			print('assumption eroror')
			pdb.set_trace()

	return 


def multivariate_beta_update(X, y, tau, residual_prec, X_t_X):
	S = np.linalg.inv(np.diag(np.ones(X.shape[1])*tau) + residual_prec*X_t_X)
	mu = residual_prec*np.dot(np.dot(S, np.transpose(X)), y)
	return mu, S

def univariate_beta_update(X, Y, tau, residual_prec, X_t_X, multivariate_mu):
	Y_resid = Y - np.dot(X, multivariate_mu)
	KK = len(multivariate_mu)

	mu_var = np.zeros(KK)

	for kk in np.random.permutation(np.arange(KK)):
		Y_resid = Y_resid + X[:,kk]*multivariate_mu[kk]
		
		posterior_var = 1.0/(X_t_X[kk,kk]*residual_prec + tau)
		posterior_mu = residual_prec*posterior_var*np.dot(X[:,kk], Y_resid)

		mu_var[kk] = posterior_var
		multivariate_mu[kk] = posterior_mu


		Y_resid = Y_resid - X[:,kk]*multivariate_mu[kk]
	cov = np.diag(mu_var)

	return multivariate_mu, cov



def h2_VI_individual_data(Y, X, multivariate_update=True,update_resid_var=False, resid_var_init=1.0, max_iter=2000, cc=0.0, burn_in_iter=0):
	# Initialize parameters
	residual_variance = resid_var_init
	expected_tau = .1


	# Relevent params
	num_snps = X.shape[1]
	num_indi = X.shape[0]

	# Precompute X_t_X
	X_t_X = np.dot(np.transpose(X), X)

	multivariate_mu = np.zeros(num_snps)

	# Iterate
	for global_iter in range(max_iter):
		# Update effects of X on Y
		if multivariate_update:
			multivariate_mu, cov = multivariate_beta_update(X, Y, expected_tau, 1.0/residual_variance, X_t_X)
		else:
			multivariate_mu, cov = univariate_beta_update(X, Y, expected_tau, 1.0/residual_variance, X_t_X, multivariate_mu)

		# Update Tau
		if global_iter >= burn_in_iter:
			E_beta_sq= np.diag(cov) + multivariate_mu*multivariate_mu
			tau_b = 0.5*np.sum(E_beta_sq) + cc
			tau_a = len(E_beta_sq)/2.0 + cc
			expected_tau = tau_a/tau_b

		# Update residual variance
		if update_resid_var:
			Y_resid = Y - np.dot(X, multivariate_mu)
			squared_error = np.sum(np.square(Y_resid)) + np.trace(np.dot(X_t_X, cov))
			#squared_error = np.sum(np.square(Y_resid)) + np.trace(np.dot(np.dot(X, cov), np.transpose(X)))  # Should be the same
			residual_variance = (.5*squared_error + cc)/(.5*num_indi + cc)

		independent_h2 = num_snps*(1.0/expected_tau)
		# Full h2
		e_beta_t_beta = cov + np.dot(multivariate_mu.reshape(num_snps, 1), multivariate_mu.reshape(1, num_snps))
		full_h2 = np.trace(np.dot(X_t_X, e_beta_t_beta))/num_indi

		'''
		if np.mod(global_iter,100) == 0:
			print(str(independent_h2) + '\t' + str(full_h2))
		'''

	return full_h2, independent_h2







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
expr_trait_h2_inference_dir = sys.argv[9]
window_version = sys.argv[10]
delta_updates = sys.argv[11]
simulated_gene_expression_dir = sys.argv[12]


# Load in true simulated data parameters
genetic_trait_expr_med_file = simulated_trait_dir + simulation_name_string +'_expression_mediated_trait_values.txt'
genetic_trait_nm_file = simulated_trait_dir +simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'
sim_med_h2 = np.var(np.loadtxt(genetic_trait_expr_med_file))
sim_nm_h2 = np.var(np.loadtxt(genetic_trait_nm_file))
sim_h2 = np.var(np.loadtxt(genetic_trait_nm_file) + np.loadtxt(genetic_trait_expr_med_file))

# Load in GWAS summary statistics
gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results.txt'
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(gwas_summary_file)

# Extract information in each window
if window_version == 'small':
	quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_summary.txt' 
else:
	quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_summary.txt' 
window_names, window_info = extract_gwas_info_for_each_window(gwas_beta, gwas_rsids, quasi_ld_window_summary_file, N_eqtl)


# load in eqtl data
if window_version == 'small':
	eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(N_eqtl) + '_small_window_eqtl_sumstats.txt'
else:
	eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(N_eqtl) + '_big_window_eqtl_sumstats.txt'


window_info, gene_info, gene_names = load_in_eqtl_data(eqtl_sumstat_file, window_info, window_names, simulated_gene_expression_dir, simulation_name_string)

output_file = expr_trait_h2_inference_dir + simulation_name_string + '_eqtl_' + str(N_eqtl) + '_' + window_version + '_' + delta_updates + '_expression_trait_h2_inference.txt'
print(output_file)
t = open(output_file,'w')
t.write('method\tgene\tsim_h2\th2_est_1\th2_est2\n')


# Loop through genes
sim_gene_h2s = []
ldsc_est_h2s = []
for gene in gene_names:
	# Get name of window gene is in 
	gene_window = gene_info[gene]['window_name']
	# cis snps
	cis_snp_indices = gene_info[gene]['cis_snps']

	# Eqtl LD
	eqtl_ld_file = window_info[gene_window]['eqtl_ld_file']
	eqtl_ld_big = np.load(eqtl_ld_file)
	eqtl_ld = eqtl_ld_big[cis_snp_indices,:][:, cis_snp_indices]

	# True betas
	true_betas = gene_info[gene]['true_beta']

	# Estimated marginal betas
	window_beta =gene_info[gene]['beta']
	beta = window_beta[cis_snp_indices]
	window_beta_se = gene_info[gene]['beta_se']
	beta_se = window_beta_se[cis_snp_indices]
	z_score = beta/beta_se
	
	# Load in individual level data
	expr_file = qtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(N_eqtl) + '_' + window_version + '_' + gene + '_expression.npy'
	geno_file = qtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(N_eqtl) + '_' + window_version + '_' + gene + '_genotype.npy'
	expr = np.load(expr_file)
	geno = np.load(geno_file)
	
	# Quick error check to make individual level data and genotype data align
	#quick_error_check_for_sum_stat_individual_data_alignment(expr, geno, beta, eqtl_ld)


	# Sim gene h2
	sim_gene_h2 = np.dot(np.dot(true_betas, eqtl_ld), true_betas)

	# Need to iterate over:
	## 1. Gibbs and VI
	## 2. hyper-prior parameters
	## 3. Univariate updates vs multivariate updates
	## 4. Individual data vs Summary stat data

	prior_params = [0.0, 1e-10, 1e-6, 1e-3]


	for prior_param in prior_params:

		# VI individual level data
		#full_h2, independent_h2 = h2_VI_individual_data(expr, geno, update_resid_var=False, resid_var_init=1.0, max_iter=15000, cc=prior_param)
		#t.write('VI_individual_const_resid_multivariate_ig_prior_' + str(prior_param) + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(independent_h2) + '\t' + str(full_h2) + '\n')
	
		# VI summary stats (multiviate updates)
		gene_mod = bayesian_vi_lmm_ss_h2_single_region.Bayesian_LMM(beta, eqtl_ld, N_eqtl)
		gene_mod.fit(burn_in_iterations=1, total_iterations=15000,cc=prior_param, univariate_updates=False)
		t.write('VI_rss_lmm_multivariate_ig_prior_' + str(prior_param) + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(gene_mod.independent_h2) + '\t' + str(gene_mod.full_h2) + '\n')

		# VI summary stats (univariate updats)
		gene_mod = bayesian_vi_lmm_ss_h2_single_region.Bayesian_LMM(beta, eqtl_ld, N_eqtl)
		gene_mod.fit(burn_in_iterations=1, total_iterations=15000,cc=prior_param, univariate_updates=True)
		t.write('VI_rss_lmm_univariate_ig_prior_' + str(prior_param) + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(gene_mod.independent_h2) + '\t' + str(gene_mod.full_h2) + '\n')

		# Gibbs summary stats (multivariate updates)
		try:
			gene_mod = bayesian_lmm_ss_h2_single_region.Bayesian_LMM(beta, eqtl_ld, N_eqtl)
			gene_mod.fit(burn_in_iterations=10000, total_iterations=15000,cc=prior_param, univariate_updates=False)
			t.write('gibbs_rss_lmm_multivariate_ig_prior_' + str(prior_param) + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(np.mean(gene_mod.sampled_h2)) + '\t' + str(np.mean(gene_mod.sampled_full_h2s)) + '\n')
		except:
			t.write('gibbs_rss_lmm_multivariate_ig_prior_' + str(prior_param) + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + 'NA' + '\t' + 'NA' + '\n')

		# Gibbs summary stats (univariate updates)
		try:
			gene_mod = bayesian_lmm_ss_h2_single_region.Bayesian_LMM(beta, eqtl_ld, N_eqtl)
			gene_mod.fit(burn_in_iterations=10000, total_iterations=15000,cc=prior_param, univariate_updates=True)
			t.write('gibbs_rss_lmm_univariate_ig_prior_' + str(prior_param) + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(np.mean(gene_mod.sampled_h2)) + '\t' + str(np.mean(gene_mod.sampled_full_h2s)) + '\n')
		except:
			t.write('gibbs_rss_lmm_univariate_ig_prior_' + str(prior_param) + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + 'NA' + '\t' + 'NA' + '\n')


	'''
	sim_gene_h2s.append(sim_gene_h2)

	# LDSC inference (no intercept)
	ldscore = np.diag(np.dot(eqtl_ld, eqtl_ld))
	chi_sq_stats = np.square(z_score)
	mod = sm.OLS(chi_sq_stats-1, ldscore)
	res = mod.fit()
	ldsc_h2 = res.params[0]*len(ldscore)/N_eqtl
	t.write('ldsc_no_intercept' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(ldsc_h2) + '\t' + str(ldsc_h2) + '\n')
	# LDSC inference (intercept)
	ldscore = np.diag(np.dot(eqtl_ld, eqtl_ld))
	chi_sq_stats = np.square(z_score)
	mod = sm.OLS(chi_sq_stats, sm.add_constant(ldscore))
	res = mod.fit()
	ldsc_h2_w_intercept = res.params[1]*len(ldscore)/N_eqtl
	t.write('ldsc_w_intercept' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(ldsc_h2_w_intercept) + '\t' + str(ldsc_h2_w_intercept) + '\n')

	# Run custom gibbs sampling inference
	gene_mod = bayesian_lmm_ss_h2_single_region.Bayesian_LMM(beta, eqtl_ld, N_eqtl)
	gene_mod.fit(burn_in_iterations=15000, total_iterations=25000,cc=1e-14)
	t.write('custom_gibbs_rss_lmm_ig_prior_1e-14' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(np.mean(gene_mod.sampled_h2)) + '\t' + str(np.mean(gene_mod.sampled_full_h2s)) + '\n')


	# Run custom gibbs sampling inference
	gene_mod = bayesian_lmm_ss_h2_single_region.Bayesian_LMM(beta, eqtl_ld, N_eqtl)
	gene_mod.fit(burn_in_iterations=15000, total_iterations=25000,cc=1e-9)
	t.write('custom_gibbs_rss_lmm_ig_prior_1e-9' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(np.mean(gene_mod.sampled_h2)) + '\t' + str(np.mean(gene_mod.sampled_full_h2s)) + '\n')


	# Run custom gibbs sampling inference
	gene_mod = bayesian_lmm_ss_h2_single_region.Bayesian_LMM(beta, eqtl_ld, N_eqtl)
	gene_mod.fit(burn_in_iterations=15000, total_iterations=25000,cc=1e-6)
	t.write('custom_gibbs_rss_lmm_ig_prior_1e-6' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(np.mean(gene_mod.sampled_h2)) + '\t' + str(np.mean(gene_mod.sampled_full_h2s)) + '\n')

	# Run custom gibbs sampling inference
	gene_mod = bayesian_lmm_ss_h2_single_region.Bayesian_LMM(beta, eqtl_ld, N_eqtl)
	gene_mod.fit(burn_in_iterations=15000, total_iterations=25000,cc=1e-3)
	t.write('custom_gibbs_rss_lmm_ig_prior_1e-3' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(np.mean(gene_mod.sampled_h2)) + '\t' + str(np.mean(gene_mod.sampled_full_h2s)) + '\n')
	'''


	#gene_mod = bayesian_vi_lmm_ss_h2_single_region.Bayesian_LMM(beta, eqtl_ld, N_eqtl)
	#gene_mod.fit(burn_in_iterations=1, total_iterations=15000,cc=0.0)


	# Run custom gibbs sampling inference
	#gene_mod = bayesian_lmm_ss_h2_single_region.Bayesian_LMM(beta, eqtl_ld, N_eqtl)
	#gene_mod.fit(burn_in_iterations=15000, total_iterations=25000,cc=0.0)
	#t.write('custom_gibbs_rss_lmm_ig_prior_0' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(np.mean(gene_mod.sampled_h2)) + '\t' + str(np.mean(gene_mod.sampled_full_h2s)) + '\n')


	'''
	full_h2, independent_h2 = h2_VI_individual_data(expr, geno, multivariate_update=False,update_resid_var=False, resid_var_init=1.0, max_iter=50000, burn_in_iter=1000)
	t.write('vi_individual_const_resid_var_prior_0_univariate_update' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(independent_h2) + '\t' + str(full_h2) + '\n')

	print('univariate_vi' + '\t' + str(independent_h2) + '\t' + str(full_h2))


	full_h2, independent_h2 = h2_VI_individual_data(expr, geno, update_resid_var=False, resid_var_init=1.0, max_iter=15000)
	t.write('vi_individual_const_resid_var_prior_0' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(independent_h2) + '\t' + str(full_h2) + '\n')

	print('multivariate_vi' + '\t' + str(independent_h2) + '\t' + str(full_h2))
	'''
	'''
	full_h2, independent_h2 = h2_VI_individual_data(expr, geno, update_resid_var=False, resid_var_init=1.0, max_iter=15000, cc=1e-14)
	t.write('vi_individual_const_resid_var_prior_1e-14' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(independent_h2) + '\t' + str(full_h2) + '\n')
	
	full_h2, independent_h2 = h2_VI_individual_data(expr, geno, update_resid_var=False, resid_var_init=1.0, max_iter=15000, cc=1e-9)
	t.write('vi_individual_const_resid_var_prior_1e-9' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(independent_h2) + '\t' + str(full_h2) + '\n')

	full_h2, independent_h2 = h2_VI_individual_data(expr, geno, update_resid_var=False, resid_var_init=1.0, max_iter=15000, cc=1e-6)
	t.write('vi_individual_const_resid_var_prior_1e-6' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(independent_h2) + '\t' + str(full_h2) + '\n')

	full_h2, independent_h2 = h2_VI_individual_data(expr, geno, update_resid_var=False, resid_var_init=1.0, max_iter=15000, cc=1e-3)
	t.write('vi_individual_const_resid_var_prior_1e-3' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(independent_h2) + '\t' + str(full_h2) + '\n')

	full_h2, independent_h2 = h2_VI_individual_data(expr, geno, update_resid_var=True, resid_var_init=1.0, max_iter=15000)
	t.write('vi_individual_update_resid_var_prior_0' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(independent_h2) + '\t' + str(full_h2) + '\n')

	'''

	# Run STAN RSS inference
	#stan_rss_lmm_full_h2, stan_rss_lmm_ind_h2 = stan_rss_lmm_h2_inference(beta, eqtl_ld, N_eqtl)
	#t.write('stan_rss_lmm_uniform_prior' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(stan_rss_lmm_ind_h2) + '\t' + str(stan_rss_lmm_full_h2) + '\n')

	# Run STAN individual inference
	#stan_individual_lmm_full_h2, stan_individual_lmm_ind_h2 = stan_individual_lmm_h2_inference(expr, geno)
	#t.write('stan_individual_lmm_uniform_prior' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(stan_individual_lmm_ind_h2) + '\t' + str(stan_individual_lmm_full_h2) + '\n')
	
	# Run STAN individual inference with learned residual_variance
	#stan_individual_resid_var_lmm_full_h2, stan_individual_resid_var_lmm_ind_h2 = stan_individual_lmm_h2_inference_resid_var(expr, geno)
	#t.write('stan_individual_resid_var_lmm_uniform_prior' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(stan_individual_resid_var_lmm_ind_h2) + '\t' + str(stan_individual_resid_var_lmm_full_h2) + '\n')
	
	# Run STAN RSS inference (IG prior)
	#stan_rss_lmm_full_h2, stan_rss_lmm_ind_h2 = stan_rss_lmm_h2_inference_IG_prior(beta, eqtl_ld, N_eqtl)
	#t.write('stan_rss_lmm_IG_prior' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(stan_rss_lmm_ind_h2) + '\t' + str(stan_rss_lmm_full_h2) + '\n')

	# Run STAN RSS inference (half-normal prior)
	#stan_rss_lmm_full_h2, stan_rss_lmm_ind_h2 = stan_rss_lmm_h2_inference_halfnormal_prior(beta, eqtl_ld, N_eqtl)
	#t.write('stan_rss_lmm_halfnormal_prior' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(stan_rss_lmm_ind_h2) + '\t' + str(stan_rss_lmm_full_h2) + '\n')

	#est_var = full_inference_fixed_residual_var(np.transpose(geno), expr)
	#print('vi_individual_weird_prior' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(est_var) + '\t' + str(est_var) + '\n')

	# Run STAN individual inference with weird xbeta prior
	#stan_individual_lmm_full_h2_xbeta_prior, stan_individual_lmm_ind_h2_xbeta_prior = stan_individual_lmm_h2_inference_xbeta_prior(expr, geno)
	#print('stan_individual_lmm_uniform_xbeta_prior' + '\t' + gene + '\t' + str(sim_gene_h2) + '\t' + str(stan_individual_lmm_ind_h2_xbeta_prior) + '\t' + str(stan_individual_lmm_full_h2_xbeta_prior) + '\n')
	t.flush()


t.close()



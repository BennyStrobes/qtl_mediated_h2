import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np
import os
import pdb
import statsmodels.api as sm
import pyro
import torch
import pyro.distributions as dist
from dataclasses import dataclass
from pyro import poutine
from pyro.infer.autoguide import AutoDiagonalNormal, AutoGuideList, AutoDelta, AutoMultivariateNormal, init_to_value
from pyro.infer import SVI, Trace_ELBO, RenyiELBO, Predictive
from pyro.distributions import constraints
import math
from sklearn.linear_model import ARDRegression, BayesianRidge, LinearRegression
import linear_regression_mixture_of_gaussian_prior


device = torch.device("cpu")
if str(device) == "cpu": torch.set_num_threads(1)

@dataclass
class Data:
	#A: torch.tensor
	X: torch.tensor
	y: torch.tensor


def simulate_data(gwas_ss, n_snps, eqtl_ss_1, med_h2_1, ge_h2_options, n_genes, snps_per_gene, nm_h2, eqtl_architecture):
	gene_midpoint_lb = int(snps_per_gene/2) +1
	gene_midpoint_ub = n_snps - (int(snps_per_gene/2) +1)

	gene_causal_eqtl_effects = []
	gene_snp_indices_arr = []
	sim_ge_h2s = []
	gene_counts = np.zeros(n_snps)
	for gene_iter in range(n_genes):

		gene_midpoint = np.random.choice(np.arange(gene_midpoint_lb, gene_midpoint_ub))
		gene_start = gene_midpoint - int(snps_per_gene/2)
		gene_end = gene_midpoint + int(snps_per_gene/2)
		gene_snp_indices = np.arange(gene_start, gene_end)
		gene_counts[gene_snp_indices] = gene_counts[gene_snp_indices] + 1.0/len(gene_snp_indices)

		# Choose gene expression h2
		ge_h2 = np.random.choice(ge_h2_options)
		sim_ge_h2s.append(ge_h2)


		if eqtl_architecture == 'polygenic':
			gene_snp_causal_eqtl_effects = np.zeros(n_snps)
			gene_snp_causal_eqtl_effects[gene_snp_indices] = np.random.normal(loc=0, scale=np.sqrt(ge_h2/snps_per_gene), size=snps_per_gene)
		elif eqtl_architecture == 'sparse':
			gene_snp_causal_eqtl_effects = np.zeros(n_snps)
			causal_indices = np.random.choice(gene_snp_indices, size=10, replace=False)
			gene_snp_causal_eqtl_effects[causal_indices] = np.random.normal(loc=0, scale=np.sqrt(ge_h2/10), size=10)
		else:
			print('assumption error: eqtl architecture ' + eqtl_architecture + ' currently not implemented')
			pdb.set_trace()

		gene_causal_eqtl_effects.append(gene_snp_causal_eqtl_effects)
		gene_snp_indices_arr.append(gene_snp_indices)
	gene_causal_eqtl_effects = np.asarray(gene_causal_eqtl_effects)

	# Simulate eQTL genotype data
	eqtl_geno_1 = np.random.normal(loc=0,scale=1.0,size=(eqtl_ss_1,n_snps))
	for snp_iter in range(n_snps):
		eqtl_geno_1[:,snp_iter] = (eqtl_geno_1[:,snp_iter] - np.mean(eqtl_geno_1[:,snp_iter]))/np.std(eqtl_geno_1[:,snp_iter])

	# Simulate GWAS genotype data
	gwas_geno = np.random.normal(loc=0,scale=1.0,size=(gwas_ss,n_snps))
	for snp_iter in range(n_snps):
		gwas_geno[:,snp_iter] = (gwas_geno[:,snp_iter] - np.mean(gwas_geno[:,snp_iter]))/np.std(gwas_geno[:,snp_iter])

	# Simulate Gene trait value
	Es = []
	for gene_iter in range(n_genes):
		genetic_gene_trait_1 = np.dot(eqtl_geno_1, gene_causal_eqtl_effects[gene_iter, :])
		E_1 = np.random.normal(loc=genetic_gene_trait_1, scale=np.sqrt(1.0-np.var(genetic_gene_trait_1)))
		E_1 = (E_1 - np.mean(E_1))/np.std(E_1)
		Es.append(E_1)

	# Simulate gwas trait value
	sim_alphas = np.random.normal(loc=0, scale=np.sqrt(med_h2_1/n_genes), size=n_genes)
	genetic_trait = np.zeros(gwas_ss)
	for gene_iter in range(n_genes):
		genetic_gene = np.dot(gwas_geno, gene_causal_eqtl_effects[gene_iter, :])
		standardized_genetic_gene = genetic_gene/np.std(genetic_gene)
		genetic_trait =genetic_trait + standardized_genetic_gene*sim_alphas[gene_iter]
	# Causal nm var-trait effect
	n_causal_snps = 3000
	sim_beta = np.zeros(n_snps)
	causal_indices = np.random.choice(np.arange(n_snps), size=n_causal_snps, replace=False)
	sim_beta[causal_indices] = np.random.normal(loc=0, scale=np.sqrt(nm_h2/n_causal_snps), size=n_causal_snps)
	genetic_trait = genetic_trait + np.dot(gwas_geno, sim_beta)

	Y = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0-np.var(genetic_trait)))
	Y = (Y - np.mean(Y))/np.std(Y)

	return Y, Es, gwas_geno, eqtl_geno_1, gene_causal_eqtl_effects, sim_beta, sim_alphas, gene_snp_indices_arr, np.asarray(sim_ge_h2s), gene_counts


def get_marginal_summary_statistics_no_intercept(Y, G):
	n_snps = G.shape[1]
	beta = []
	beta_se = []
	for snp_iter in range(n_snps):
		olser = sm.OLS(Y, G[:,snp_iter]).fit()
		beta.append(olser.params[0])
		beta_se.append(olser.bse[0])
	return np.asarray(beta), np.asarray(beta_se)

def get_marginal_summary_statistics(Y, G):
	n_snps = G.shape[1]
	beta = []
	beta_se = []
	for snp_iter in range(n_snps):
		olser = sm.OLS(Y, sm.add_constant(G[:,snp_iter])).fit()
		beta.append(olser.params[1])
		beta_se.append(olser.bse[1])
	return np.asarray(beta), np.asarray(beta_se)

def get_marginal_summary_statistics_across_genes(Es, G, gene_snp_indices_arr):
	n_genes = len(Es)
	n_snps = E_geno.shape[1]

	beta_arr = []
	beta_se_arr = []
	# Loop through genes
	for gene_iter in range(n_genes):
		beta = np.zeros(n_snps)
		beta_se = np.ones(n_snps)*1e-7
		gene_snp_indices = gene_snp_indices_arr[gene_iter]

		for gene_snp_index in gene_snp_indices:
			olser = sm.OLS(Es[gene_iter], sm.add_constant(G[:,gene_snp_index])).fit()
			beta[gene_snp_index] = olser.params[1]
			beta_se[gene_snp_index] = olser.bse[1]

		# Add to global array
		beta_arr.append(beta)
		beta_se_arr.append(beta_se)

	beta_arr = np.asarray(beta_arr)
	beta_se_arr = np.asarray(beta_se_arr)

	return beta_arr, beta_se_arr


def posterior_distribution_on_causal_eqtl_effects_closed_form(expression_vec, genotype_mat, ge_h2, resid_var):
	n_snps = genotype_mat.shape[1]
	#resid_var = 1.0 - ge_h2
	per_snp_h2 = ge_h2/n_snps 

	tmp_prec = (np.eye(n_snps)/per_snp_h2) + (np.dot(np.transpose(genotype_mat), genotype_mat)/(resid_var))
	posterior_var = np.linalg.inv(tmp_prec)

	posterior_mean = (1.0/resid_var)*np.dot(np.dot(posterior_var, np.transpose(genotype_mat)), expression_vec)

	return posterior_mean, posterior_var



# Code from David Knowles
def robust_chol(cov, desired_min_eig = 1e-6):
	try: # faster than the eigendecomposition if already PSD
		chol_cov = torch.linalg.cholesky(cov)
	except:
		L, V = torch.linalg.eigh(cov)
		cov_min_eig = L.min().item()
		if cov_min_eig < desired_min_eig:
			print("Degenerate cov (min eigenvalue=%1.3e)" % cov_min_eig)
			# smallest addition to diagonal to make min(eig) = min_eig
			cov += (desired_min_eig - cov_min_eig) * torch.eye(cov.shape[0], device=cov.device)
		try:
			chol_cov = torch.linalg.cholesky(cov)
		except:
			print("Couldn't fix cov")
			return torch.eye(cov.shape[0], device=cov.device)
	return chol_cov

# Code from David Knowles
def fit_bayesian_lmm(
	data,
	num_samples = 100,
	iterations = 1000,
	brr_init = None):

	X_X_t = data.X @ data.X.T

	device = data.X.device
	type_kwargs = {"device" : device, "dtype" : torch.float}
	print(device)

	N,P = data.X.shape
	one_vec = torch.ones(P, **type_kwargs)
	one = torch.tensor(1., **type_kwargs)

	if brr_init:
		brr_weights_scale = brr_init.lambda_**-0.5
		brr_noise_scale = (1.0 - (brr_init.lambda_**-1.0))**0.5
		# Quick fix if noise is really low
		brr_noise_scale = np.sqrt(.95)
		brr_weights_scale = np.sqrt(.05)

	one = torch.tensor(1., **type_kwargs)
	twenty_five = torch.tensor(1., **type_kwargs)

	#sqrt_phi_prior = dist.HalfCauchy(one)
	#noise_scale_prior = dist.HalfCauchy(one)

	def convertr(hyperparam, name):
		return torch.tensor(hyperparam, device = device) if (
		  type(hyperparam) in [float,np.float32,torch.float]
		) else pyro.sample(name, hyperparam)

	def model(data):
		# could do e.g. dist.InverseGamma(2. * one,2. * one) instead
		sqrt_phi = pyro.sample("sqrt_phi", dist.HalfCauchy(one))
		noise_scale = pyro.sample("noise_scale", dist.HalfCauchy(one))

		# this could be done much more efficiently by precomputing eig(data.X @ data.X.T)
		# and just updating the eigenvalues as sqrt_phi and noise_scale change
		#cov = sqrt_phi**2 * (data.X @ data.X.T) + noise_scale**2 * torch.eye(N, **type_kwargs)
		cov = sqrt_phi**2 * (X_X_t) + noise_scale**2 * torch.eye(N, **type_kwargs)
		L = robust_chol(cov)
		obs = pyro.sample('obs', dist.MultivariateNormal(torch.zeros(N, **type_kwargs), scale_tril = L), obs = data.y)

	init_dic = {
	  "sqrt_psi" : torch.tensor(brr_weights_scale, **type_kwargs),
	  "noise_scale" : torch.tensor( brr_noise_scale, **type_kwargs),
	} if (brr_init != None) else {}

	guide = AutoMultivariateNormal(
		model,
		init_loc_fn = init_to_value(values=init_dic))


	adam = pyro.optim.Adam({"lr": 0.03})
	svi = SVI(model, guide, adam, loss=Trace_ELBO())
	pyro.clear_param_store()
	losses = []
	time_per_iter = []
	for j in range(iterations):
		print(j, end = '\r')
		loss = svi.step(data)
		losses.append(loss)

	samples = Predictive(
		model,
		guide=guide,
		return_sites = ["sqrt_phi", "noise_scale"],
		num_samples=num_samples)(data)

	return losses, samples


def posterior_distribution_on_causal_eqtl_effects_two_versions(expr_vec, eqtl_geno):
	# dimensions
	NN = len(expr_vec)
	PP = eqtl_geno.shape[1]

	# Initialization with sklearn bayesian ridge regression
	brr = BayesianRidge(compute_score=True, n_iter=100000, fit_intercept=False).fit(eqtl_geno, expr_vec)

	# Get data in compact format
	data = Data(
	X=torch.tensor(eqtl_geno, device = device, dtype = torch.float),
	y=torch.tensor(expr_vec, device = device, dtype = torch.float))


	# Run bayesian LMM (done in pyro)
	# Code from David Knowles
	bayesian_losses, bayesian_lmm_samples = fit_bayesian_lmm(data, num_samples=2000, iterations=4000, brr_init=brr)


	# Extract h2 parameters from est
	genetic_var = PP*np.square(bayesian_lmm_samples["sqrt_phi"])
	resid_var = np.square(bayesian_lmm_samples["noise_scale"])
	PVE = genetic_var/(genetic_var + resid_var)
	print(np.mean(np.asarray(genetic_var)))
	print(np.mean(np.asarray(resid_var)))

	# Num samples of variance parameters
	num_samples = len(bayesian_lmm_samples["sqrt_phi"])

	# Calculate Mean and Covariance of BETA (which is a mixture of gaussians)
	# First pass: get mean and covariance of each mixture
	mixture_means = []
	mixture_covs = []
	standardized_mixture_means = []
	standardized_mixture_covs = []

	for i in range(num_samples):
		sqrt_phi = bayesian_lmm_samples["sqrt_phi"][i]
		noise_scale = bayesian_lmm_samples["noise_scale"][i]
		mixture_mean, mixture_cov = posterior_distribution_on_causal_eqtl_effects_closed_form(expr_vec, eqtl_geno, PP*np.square(np.asarray(sqrt_phi)), np.square(np.asarray(noise_scale)))
		mixture_means.append(mixture_mean)
		mixture_covs.append(mixture_cov)

		# standardize
		standardized_mixture_mean = mixture_mean/np.sqrt((np.square(np.asarray(sqrt_phi))*PP))
		standardized_mixture_cov = mixture_cov/(np.square(np.asarray(sqrt_phi))*PP)
		standardized_mixture_means.append(standardized_mixture_mean)
		standardized_mixture_covs.append(standardized_mixture_cov)
	# Second Pass: get global mean and global covariance (mean and covariance across mixture components)
	# See https://math.stackexchange.com/questions/195911/calculation-of-the-covariance-of-gaussian-mixtures for derivation
	# Unstandardized
	global_mixture_mean = np.mean(np.asarray(mixture_means),axis=0)
	global_mixture_cov = np.copy(mixture_cov)*0.0
	for i in range(num_samples):
		term1 = (1.0/num_samples)*(mixture_covs[i])
		diff_vec = mixture_means[i] - global_mixture_mean
		term2 = (1.0/num_samples)*np.dot(diff_vec.reshape(len(diff_vec),1), diff_vec.reshape(1,len(diff_vec)))
		global_mixture_cov = global_mixture_cov + term1 + term2
	normalizer = np.sum(np.square(global_mixture_mean) + np.diag(global_mixture_cov))
	global_mixture_mean = global_mixture_mean/np.sqrt(normalizer)
	global_mixture_cov = global_mixture_cov/normalizer

	# Standardized
	global_standardized_mixture_mean = np.mean(np.asarray(standardized_mixture_means),axis=0)
	global_standardized_mixture_cov = np.copy(standardized_mixture_cov)*0.0
	for i in range(num_samples):
		term1 = (1.0/num_samples)*(standardized_mixture_covs[i])
		diff_vec = standardized_mixture_means[i] - global_standardized_mixture_mean
		term2 = (1.0/num_samples)*np.dot(diff_vec.reshape(len(diff_vec),1), diff_vec.reshape(1,len(diff_vec)))
		global_standardized_mixture_cov = global_standardized_mixture_cov + term1 + term2

	standardized_normalizer = np.sum(np.square(global_standardized_mixture_mean) + np.diag(global_standardized_mixture_cov))
	global_standardized_mixture_mean = global_standardized_mixture_mean/np.sqrt(standardized_normalizer)
	global_standardized_mixture_cov = global_standardized_mixture_cov/standardized_normalizer	

	return global_mixture_mean, global_mixture_cov, global_standardized_mixture_mean, global_standardized_mixture_cov, np.asarray(genetic_var), np.asarray(resid_var)

def posterior_distribution_on_causal_eqtl_effects(expr_vec, eqtl_geno, standardize_genetic_ge=False):
	# dimensions
	NN = len(expr_vec)
	PP = eqtl_geno.shape[1]

	# Initialization with sklearn bayesian ridge regression
	brr = BayesianRidge(compute_score=True, n_iter=100000, fit_intercept=False).fit(eqtl_geno, expr_vec)

	# Get data in compact format
	data = Data(
	X=torch.tensor(eqtl_geno, device = device, dtype = torch.float),
	y=torch.tensor(expr_vec, device = device, dtype = torch.float))


	# Run bayesian LMM (done in pyro)
	# Code from David Knowles
	bayesian_losses, bayesian_lmm_samples = fit_bayesian_lmm(data, num_samples=2000, iterations=4000, brr_init=brr)


	# Extract h2 parameters from est
	genetic_var = PP*np.square(bayesian_lmm_samples["sqrt_phi"])
	resid_var = np.square(bayesian_lmm_samples["noise_scale"])
	PVE = genetic_var/(genetic_var + resid_var)
	print(np.mean(np.asarray(genetic_var)))
	print(np.mean(np.asarray(resid_var)))

	# Num samples of variance parameters
	num_samples = len(bayesian_lmm_samples["sqrt_phi"])

	# Calculate Mean and Covariance of BETA (which is a mixture of gaussians)
	# First pass: get mean and covariance of each mixture
	mixture_means = []
	mixture_covs = []
	for i in range(num_samples):
		sqrt_phi = bayesian_lmm_samples["sqrt_phi"][i]
		noise_scale = bayesian_lmm_samples["noise_scale"][i]
		mixture_mean, mixture_cov = posterior_distribution_on_causal_eqtl_effects_closed_form(expr_vec, eqtl_geno, PP*np.square(np.asarray(sqrt_phi)), np.square(np.asarray(noise_scale)))
		if standardize_genetic_ge:
			mixture_mean = mixture_mean/np.sqrt((np.square(np.asarray(sqrt_phi))*PP))
			mixture_cov = mixture_cov/(np.square(np.asarray(sqrt_phi))*PP)
		mixture_means.append(mixture_mean)
		mixture_covs.append(mixture_cov)
	# Second Pass: get global mean and global covariance (mean and covariance across mixture components)
	# See https://math.stackexchange.com/questions/195911/calculation-of-the-covariance-of-gaussian-mixtures for derivation
	global_mixture_mean = np.mean(np.asarray(mixture_means),axis=0)
	global_mixture_cov = np.copy(mixture_cov)*0.0
	for i in range(num_samples):
		term1 = (1.0/num_samples)*(mixture_covs[i])
		diff_vec = mixture_means[i] - global_mixture_mean
		term2 = (1.0/num_samples)*np.dot(diff_vec.reshape(len(diff_vec),1), diff_vec.reshape(1,len(diff_vec)))
		global_mixture_cov = global_mixture_cov + term1 + term2

	return global_mixture_mean, global_mixture_cov, np.asarray(genetic_var), np.asarray(resid_var)

def get_genes_posterior_on_gene_related_snp_anno(posterior_mean, posterior_var, LD):
	e_delta_delta_t = np.dot(posterior_mean.reshape((len(posterior_mean),1)), posterior_mean.reshape((1,len(posterior_mean)))) + posterior_var
	posterior_gene_related_snp_annotation = np.sum(np.dot(LD, e_delta_delta_t)*LD,axis=1)
	return posterior_gene_related_snp_annotation



def compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, gene_ld_score, variant_ld_scores):
	chi_sq = np.square(gwas_z)
	ld_scores = np.transpose(np.vstack((variant_ld_scores,gene_ld_score)))
	mod = sm.OLS(chi_sq-1, ld_scores)
	res = mod.fit()
	est_per_variant_h2 = (res.params[0]/gwas_ss)
	est_per_gene_h2 = res.params[1]/gwas_ss
	return est_per_variant_h2, est_per_gene_h2


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


# Open and print header to output file
output_file = output_root + '_effect_est_res_summary.txt'
t = open(output_file,'w')
t.write('method\tsim_iter\tsim_nm_var_h2\tsim_med_h2\test_nm_var_h2\test_med_h2\n')


# Open and print header to output file
output_file = output_root + '_gene_var_est.txt'
t_gene = open(output_file,'w')
t_gene.write('sim_iter\tgene_name\tsim_gene_var\tsim_gene_var_obs\tmean_est_gene_var\tmean_est_resid_var\test_gene_var\test_resid_var\n')


# Options for ge_h2
ge_h2_options = [.05,.1, .2]

# Loop through sims
for sim_iter in range(n_sims):
	print(sim_iter)

	# First simulate data
	Y, Es, Y_geno, E_geno, gene_causal_eqtl_effects, nm_var_causal_effects, gene_trait_effects, gene_snp_indices_arr, sim_ge_h2s, gene_window_based_ld_score = simulate_data(gwas_ss, n_snps, eqtl_ss_1, med_h2_1, ge_h2_options, n_genes, snps_per_gene, nm_h2, eqtl_architecture)

	# Sim h2s
	sim_nm_var_h2 = np.sum(np.square(nm_var_causal_effects))
	sim_med_h2 = np.sum(np.square(gene_trait_effects))

	# Next get summary statistics
	gwas_beta, gwas_beta_se = get_marginal_summary_statistics(Y, Y_geno)
	#eqtl_beta, eqtl_beta_se = get_marginal_summary_statistics_across_genes(Es, E_geno, gene_snp_indices_arr)
	# Compute gene level z-scores
	gwas_z = gwas_beta/gwas_beta_se

	# gene ld scores
	sum_gene_ld_scores = np.zeros(n_snps)
	sum_gene_ld_scores2 = np.zeros(n_snps)
	sum_gene_mean_only_ld_scores = np.zeros(n_snps)


	# Fit gene model for each gene
	for gene_iter in range(n_genes):
		expr_vec =Es[gene_iter]
		gene_geno = E_geno[:, gene_snp_indices_arr[gene_iter]]
		n_cis_snps = gene_geno.shape[1]

		mod = linear_regression_mixture_of_gaussian_prior.LR_MoG_Prior()
		mod.fit(expr_vec, gene_geno)
		pdb.set_trace()

		# Compute gene h2 with posterior eqtl est
		eqtl_posterior_mean, eqtl_posterior_var, eqtl_posterior_mean2, eqtl_posterior_var2, est_genetic_varz, est_resid_varz = posterior_distribution_on_causal_eqtl_effects_two_versions(expr_vec, gene_geno)

		# Print estimated genetic varz and resid varz to output file
		t_gene.write(str(sim_iter) + '\t' + str(gene_iter) + '\t' + str(sim_ge_h2s[gene_iter]) + '\t' + str(np.sum(np.square(gene_causal_eqtl_effects[gene_iter,:]))) + '\t' + str(np.mean(est_genetic_varz)) + '\t' + str(np.mean(est_resid_varz)) + '\t')
		t_gene.write(';'.join(est_genetic_varz.astype(str)) + '\t' + ';'.join(est_resid_varz.astype(str)) + '\n')

		# Compute gene ld scores using full posterior distribution
		gene_ld_score = get_genes_posterior_on_gene_related_snp_anno(eqtl_posterior_mean, eqtl_posterior_var, np.eye(n_cis_snps))
		sum_gene_ld_scores[gene_snp_indices_arr[gene_iter]] = sum_gene_ld_scores[gene_snp_indices_arr[gene_iter]] + gene_ld_score

		gene_ld_score2 = get_genes_posterior_on_gene_related_snp_anno(eqtl_posterior_mean2, eqtl_posterior_var2, np.eye(n_cis_snps))
		sum_gene_ld_scores2[gene_snp_indices_arr[gene_iter]] = sum_gene_ld_scores2[gene_snp_indices_arr[gene_iter]] + gene_ld_score2

		# Compute gene ld scores using only mean of posterior distribution
		gene_ld_score_mean_only = get_genes_posterior_on_gene_related_snp_anno(eqtl_posterior_mean/np.sqrt(np.sum(np.square(eqtl_posterior_mean))), np.zeros(eqtl_posterior_var.shape), np.eye(n_cis_snps))
		sum_gene_mean_only_ld_scores[gene_snp_indices_arr[gene_iter]] = sum_gene_mean_only_ld_scores[gene_snp_indices_arr[gene_iter]] + gene_ld_score_mean_only


	var_ld_scores = np.ones(n_snps)
	est_per_variant_h2, est_per_gene_h2 = compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, sum_gene_ld_scores, var_ld_scores)

	est_nm_h2 = est_per_variant_h2*n_snps
	est_med_h2 = est_per_gene_h2*n_genes

	# Print to output
	t.write(str('posterior_blup') + '\t' + str(sim_iter) + '\t' + str(sim_nm_var_h2) + '\t' + str(sim_med_h2) + '\t' + str(est_nm_h2) + '\t' + str(est_med_h2) + '\n')


	var_ld_scores = np.ones(n_snps)
	est_per_variant_h2, est_per_gene_h2 = compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, sum_gene_ld_scores2, var_ld_scores)

	est_nm_h2 = est_per_variant_h2*n_snps
	est_med_h2 = est_per_gene_h2*n_genes

	# Print to output
	t.write(str('posterior_blup2') + '\t' + str(sim_iter) + '\t' + str(sim_nm_var_h2) + '\t' + str(sim_med_h2) + '\t' + str(est_nm_h2) + '\t' + str(est_med_h2) + '\n')

	# Do for PMCES
	var_ld_scores = np.ones(n_snps)
	est_per_variant_h2, est_per_gene_h2 = compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, sum_gene_mean_only_ld_scores, var_ld_scores)

	est_nm_h2 = est_per_variant_h2*n_snps
	est_med_h2 = est_per_gene_h2*n_genes

	# Print to output
	t.write(str('mean_blup') + '\t' + str(sim_iter) + '\t' + str(sim_nm_var_h2) + '\t' + str(sim_med_h2) + '\t' + str(est_nm_h2) + '\t' + str(est_med_h2) + '\n')

	np.savetxt(output_root + '_sum_gene_ld_scores.txt', sum_gene_ld_scores, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_sum_gene_ld_scores2.txt', sum_gene_ld_scores2, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_sum_gene_mean_only_ld_scores.txt', sum_gene_mean_only_ld_scores, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_gwas_z.txt', gwas_z, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_gwas_beta.txt', gwas_beta, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_gwas_beta_se.txt', gwas_beta_se, fmt="%s", delimiter='\t')

	var_ld_scores = np.ones(n_snps)
	est_per_variant_h2, est_per_gene_h2 = compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, gene_window_based_ld_score, var_ld_scores)
	est_nm_h2 = est_per_variant_h2*n_snps
	est_med_h2 = est_per_gene_h2*n_genes
	print(est_nm_h2)
	print(est_med_h2)
	t.write(str('window_based') + '\t' + str(sim_iter) + '\t' + str(sim_nm_var_h2) + '\t' + str(sim_med_h2) + '\t' + str(est_nm_h2) + '\t' + str(est_med_h2) + '\n')

t.close()
t_gene.close()


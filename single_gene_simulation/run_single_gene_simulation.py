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

device = torch.device("cpu")
if str(device) == "cpu": torch.set_num_threads(1)

@dataclass
class Data:
	#A: torch.tensor
	X: torch.tensor
	y: torch.tensor


def simulate_causal_eqtl_effects(n_snps, gene_cis_h2, sparse=False):
	if sparse == False:
		return np.random.normal(loc=0, scale=np.sqrt(gene_cis_h2/n_snps),size=n_snps)
	else:
		n_causal_snps=5
		causal_eqtl_effects = np.zeros(n_snps)
		causal_indices = np.random.choice(np.arange(n_snps), size=n_causal_snps, replace=False)
		causal_eqtl_effects[causal_indices] = np.random.normal(loc=0, scale=np.sqrt(gene_cis_h2/n_causal_snps),size=n_causal_snps)
		return causal_eqtl_effects

def print_95_ci(arr):
	meany = np.mean(arr)
	se = np.std(arr)/np.sqrt(len(arr))
	ub = meany + (1.96*se)
	lb = meany - (1.96*se)

	print(str(meany) + ':   ' + '[' + str(lb) + ', ' + str(ub) + ']')
	return


def simulate_gwas_trait(gwas_geno, causal_eqtl_effects, gene_trait_alpha):
	# Get genetic gene expression
	genetic_ge = np.dot(gwas_geno, causal_eqtl_effects)
	standardized_genetic_ge = genetic_ge/np.std(genetic_ge)

	genetic_trait = gene_trait_alpha*standardized_genetic_ge

	trait = np.random.normal(genetic_trait, scale=np.sqrt(1.0 - np.var(genetic_trait)))


	std_trait = (trait - np.mean(trait))/np.std(trait)


	return std_trait


def compute_gwas_z_scores(trait_vec, geno_mat):
	z_scores = []
	n_snps = geno_mat.shape[1]
	for snp_iter in range(n_snps):
		mod = sm.OLS(trait_vec, sm.add_constant(geno_mat[:, snp_iter]))
		res = mod.fit()
		# Extract results
		effect_size = res.params[1]
		effect_size_se = res.bse[1]
		effect_size_z = effect_size/effect_size_se
		z_scores.append(effect_size_z)

	return np.asarray(z_scores)


def compute_snp_h2_with_ldsc(gwas_z, gwas_ld, gwas_ss):
	n_snps = gwas_ld.shape[0]
	chi_sq = np.square(gwas_z)
	ld_scores = np.sum(np.square(gwas_ld),axis=0)

	mod = sm.OLS(chi_sq-1, ld_scores)
	res = mod.fit()
	est_snp_h2 = n_snps*(res.params[0]/gwas_ss)
	return est_snp_h2

def compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, gene_ld_score, gene_cis_h2):
	n_snps = gwas_ld.shape[0]
	chi_sq = np.square(gwas_z)
	mod = sm.OLS(chi_sq-1, gene_ld_score)
	res = mod.fit()
	est_gene_h2 = (res.params[0]/gwas_ss)*gene_cis_h2
	return est_gene_h2

def simulate_gene_expression(causal_eqtl_effects, eqtl_geno):
	genetic_gene_expression = np.dot(eqtl_geno, causal_eqtl_effects)
	resid_var = 1.0-np.var(genetic_gene_expression)
	gene_expression = np.random.normal(genetic_gene_expression, scale=np.sqrt(resid_var))

	gene_expression = (gene_expression - np.mean(gene_expression))/np.std(gene_expression)

	return gene_expression,genetic_gene_expression, resid_var

def posterior_distribution_on_causal_eqtl_effects_closed_form(expression_vec, genotype_mat, ge_h2, resid_var):
	n_snps = genotype_mat.shape[1]
	#resid_var = 1.0 - ge_h2
	per_snp_h2 = ge_h2/n_snps 

	tmp_prec = (np.eye(n_snps)/per_snp_h2) + (np.dot(np.transpose(genotype_mat), genotype_mat)/(resid_var))
	posterior_var = np.linalg.inv(tmp_prec)

	posterior_mean = (1.0/resid_var)*np.dot(np.dot(posterior_var, np.transpose(genotype_mat)), expression_vec)

	return posterior_mean, posterior_var

def get_genes_posterior_on_gene_related_snp_anno(posterior_mean, posterior_var, LD):
	e_delta_delta_t = np.dot(posterior_mean.reshape((len(posterior_mean),1)), posterior_mean.reshape((1,len(posterior_mean)))) + posterior_var
	posterior_gene_related_snp_annotation = np.sum(np.dot(LD, e_delta_delta_t)*LD,axis=1)
	return posterior_gene_related_snp_annotation

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
		brr_noise_scale = brr_init.alpha_**-0.5

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




def posterior_distribution_on_causal_eqtl_effects(expr_vec, eqtl_geno, standardize_genetic_ge=False):
	# dimensions
	NN = len(expr_vec)
	PP = eqtl_geno.shape[1]

	# Initialization with sklearn bayesian ridge regression
	brr = BayesianRidge(compute_score=True, n_iter=600, fit_intercept=False, alpha_1=1e-12, alpha_2=1e-12, lambda_1=1e-12, lambda_2=1e-12).fit(eqtl_geno, expr_vec)

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
	return global_mixture_mean, global_mixture_cov





#####################
# Command line args
#####################
genotype_dir = sys.argv[1]
output_dir = sys.argv[2]


########################
# Simulation parameters
#######################
gene_cis_h2 = 0.05
gene_trait_alpha = 0.05
n_sims = 100
expected_h2 = np.square(gene_trait_alpha)
eqtl_ss = 1000

########################
# Load in genotype data
#######################
# GWAS genotype
window_names = ['1','4','7', '9', '12', '20','22','24','26','32']
gwas_genos = []
for window_name in window_names:
	gwas_genotype_file = genotype_dir + 'gwas_genotype_1_' + window_name + '.npy'
	gwas_genos.append(np.load(gwas_genotype_file)[:,:20])
gwas_geno = np.hstack(gwas_genos)

# EQTL genotype
eqtl_genotype_file = genotype_dir + 'eqtl_' + str(eqtl_ss) + '_genotype_1.npy'
full_eqtl_geno = np.load(eqtl_genotype_file)
eqtl_genos = []
for window_name in window_names:
	start_index = 175*int(window_name)
	eqtl_genos.append(full_eqtl_geno[:, start_index:(start_index+20)])
eqtl_geno = np.hstack(eqtl_genos)

# GWAS LD
gwas_ld = np.corrcoef(np.transpose(gwas_geno))

# Other rando parameters
gwas_ss = gwas_geno.shape[0]
n_snps = gwas_geno.shape[1]



########################
# Loop through simulations
#######################
est_snp_h2s = []
est_gene_h2s_known_var = []
est_gene_h2s = []

for sim_iter in range(n_sims):
	print('###################')
	print(sim_iter)
	# Simulate causal eqtl effects
	causal_eqtl_effects = simulate_causal_eqtl_effects(n_snps, gene_cis_h2, sparse=True)

	# Simulate gwas trait
	trait_vec = simulate_gwas_trait(gwas_geno, causal_eqtl_effects, gene_trait_alpha)

	# Simulate gene expression
	expr_vec, genetic_gene_expression, resid_var = simulate_gene_expression(causal_eqtl_effects, eqtl_geno)

	# Get gwas z-scores
	gwas_z = compute_gwas_z_scores(trait_vec, gwas_geno)

	# Compute snp h2 with ldsc
	est_snp_h2 = compute_snp_h2_with_ldsc(gwas_z, gwas_ld, gwas_ss)
	est_snp_h2s.append(est_snp_h2)


	# Compute gene h2 with posterior eqtl est (given known prior and known residual variance)
	eqtl_posterior_mean_known_h2_known_resid_var, eqtl_posterior_var_known_h2_known_resid_var = posterior_distribution_on_causal_eqtl_effects_closed_form(expr_vec, eqtl_geno, gene_cis_h2, resid_var)
	gene_ld_score = get_genes_posterior_on_gene_related_snp_anno(eqtl_posterior_mean_known_h2_known_resid_var, eqtl_posterior_var_known_h2_known_resid_var, gwas_ld)
	# Compute gene h2
	est_gene_h2_known_var = compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, gene_ld_score, gene_cis_h2)
	est_gene_h2s_known_var.append(est_gene_h2_known_var)

	# Compute gene h2 with posterior eqtl est
	eqtl_posterior_mean, eqtl_posterior_var = posterior_distribution_on_causal_eqtl_effects(expr_vec, eqtl_geno, standardize_genetic_ge=True)
	#print(np.diag(eqtl_posterior_var)/(np.diag(eqtl_posterior_var) + np.square(eqtl_posterior_mean)))
	gene_ld_score = get_genes_posterior_on_gene_related_snp_anno(eqtl_posterior_mean, eqtl_posterior_var, gwas_ld)
	#est_gene_h2 = compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, gene_ld_score, gene_cis_h2)
	est_gene_h2 = compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, gene_ld_score, 1.0)
	est_gene_h2s.append(est_gene_h2)

	print(est_gene_h2)
	print(est_gene_h2_known_var)



# Convert into nice, clean array
est_snp_h2s = np.asarray(est_snp_h2s)
est_gene_h2s_known_var = np.asarray(est_gene_h2s_known_var)
est_gene_h2s = np.asarray(est_gene_h2s)
print('####################')
print_95_ci(est_snp_h2s)
print_95_ci(est_gene_h2s_known_var)
print_95_ci(est_gene_h2s)


import sys
import numpy as np
import os
import pdb

import jax.numpy as jnp
import jax.random as random
import statsmodels.api as sm

import numpyro
from numpyro.diagnostics import summary
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS
import time




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


# regression model with continuous-valued outputs/responses
def model_normal_likelihood_horseshoe_prior(X, Y):
	D_X = X.shape[1]

	# sample from horseshoe prior
	#lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	tau = numpyro.sample("tau", dist.HalfCauchy(jnp.ones(1)))

	horseshoe_sigma = tau**2*lambdas**2
	Beta = numpyro.sample('beta', dist.Normal(loc=0, scale=(horseshoe_sigma**.5)))
	#prec_obs = numpyro.sample("prec_obs", dist.Gamma(3.0, 1.0))
	#sigma_obs = 1.0 / jnp.sqrt(prec_obs)

	noise_scale = numpyro.sample("noise_scale", dist.HalfCauchy(jnp.ones(1)))


	mu = jnp.dot(X, Beta)
	numpyro.sample('obs', dist.Normal(loc=mu, scale=noise_scale), obs=Y)

	return



# regression model with continuous-valued outputs/responses
def model_normal_likelihood_dirichlet_horseshoe_prior(X, Y):
	D_X = X.shape[1]

	# sample from horseshoe prior
	#lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	lambdas = numpyro.sample("lambdas", dist.Dirichlet(jnp.ones(D_X)))
	tau = numpyro.sample("tau", dist.HalfCauchy(jnp.ones(1)))

	horseshoe_sigma = tau**2*lambdas
	Beta = numpyro.sample('beta', dist.Normal(loc=0, scale=(horseshoe_sigma**.5)))
	#prec_obs = numpyro.sample("prec_obs", dist.Gamma(3.0, 1.0))
	#sigma_obs = 1.0 / jnp.sqrt(prec_obs)

	noise_scale = numpyro.sample("noise_scale", dist.HalfCauchy(jnp.ones(1)))


	mu = jnp.dot(X, Beta)
	numpyro.sample('obs', dist.Normal(loc=mu, scale=noise_scale), obs=Y)

	return


# helper function for HMC inference
def run_inference(model, rng_key, X, Y, num_warmup=3000, num_samples=8000, thinning=5):
	start = time.time()
	kernel = NUTS(model)
	mcmc = MCMC(
		kernel,
		num_warmup=num_warmup,
		num_samples=num_samples,
		thinning=thinning,
		num_chains=1,
		progress_bar=False if "NUMPYRO_SPHINXBUILD" in os.environ else True,
	)

	mcmc.run(rng_key, X, Y)
	mcmc.print_summary(exclude_deterministic=False)

	samples = mcmc.get_samples()
	#summary_dict = summary(samples, group_by_chain=False)

	print("\nMCMC elapsed time:", time.time() - start)

	return samples


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
	sum_gene_ld_scores3 = np.zeros(n_snps)

	sum_gene_mean_only_ld_scores = np.zeros(n_snps)


	# Fit gene model for each gene
	for gene_iter in range(n_genes):
		expr_vec =Es[gene_iter]
		gene_geno = E_geno[:, gene_snp_indices_arr[gene_iter]]
		n_cis_snps = gene_geno.shape[1]


		rng_key, rng_key_predict = random.split(random.PRNGKey(0))
		samples2 = run_inference(model_normal_likelihood_horseshoe_prior, rng_key, gene_geno, expr_vec)
		eqtl_posterior_mean2 = np.mean(samples2['beta'],axis=0)
		eqtl_posterior_var2 = np.cov(np.transpose(samples2['beta']))
		scaler2 = np.sum(np.square(eqtl_posterior_mean2)) + np.sum(np.diag(eqtl_posterior_var2))
		standardized_eqtl_posterior_mean2 = eqtl_posterior_mean2/np.sqrt(scaler2)
		standardized_eqtl_posterior_var2 = eqtl_posterior_var2/scaler2
		#est_genetic_varz = np.var(np.dot(gene_geno, np.transpose(samples2['beta'])), axis=0)
		#est_resid_varz = np.square(samples2['noise_scale'][:,0])
		#posterior_mean2 = np.mean(samples2['lambdas'],axis=0)


		#mod = linear_regression_mixture_of_gaussian_prior.LR_MoG_Prior()
		#mod.fit(expr_vec, gene_geno)
		rng_key, rng_key_predict = random.split(random.PRNGKey(0))
		samples = run_inference(model_normal_likelihood_dirichlet_horseshoe_prior, rng_key, gene_geno, expr_vec)
		eqtl_posterior_mean = np.mean(samples['beta'],axis=0)
		eqtl_posterior_var = np.cov(np.transpose(samples['beta']))
		scaler = np.sum(np.square(eqtl_posterior_mean)) + np.sum(np.diag(eqtl_posterior_var))
		standardized_eqtl_posterior_mean = eqtl_posterior_mean/np.sqrt(scaler)
		standardized_eqtl_posterior_var = eqtl_posterior_var/scaler
		est_genetic_varz = np.var(np.dot(gene_geno, np.transpose(samples['beta'])), axis=0)
		est_resid_varz = np.square(samples['noise_scale'][:,0])
		posterior_mean3 = np.mean(samples['lambdas'],axis=0)

		# Print estimated genetic varz and resid varz to output file
		t_gene.write(str(sim_iter) + '\t' + str(gene_iter) + '\t' + str(sim_ge_h2s[gene_iter]) + '\t' + str(np.sum(np.square(gene_causal_eqtl_effects[gene_iter,:]))) + '\t' + str(np.mean(est_genetic_varz)) + '\t' + str(np.mean(est_resid_varz)) + '\t')
		t_gene.write(';'.join(est_genetic_varz.astype(str)) + '\t' + ';'.join(est_resid_varz.astype(str)) + '\n')

		# Compute gene ld scores using full posterior distribution
		gene_ld_score = get_genes_posterior_on_gene_related_snp_anno(standardized_eqtl_posterior_mean, standardized_eqtl_posterior_var, np.eye(n_cis_snps))
		sum_gene_ld_scores[gene_snp_indices_arr[gene_iter]] = sum_gene_ld_scores[gene_snp_indices_arr[gene_iter]] + gene_ld_score

		gene_ld_score2 = get_genes_posterior_on_gene_related_snp_anno(standardized_eqtl_posterior_mean2, standardized_eqtl_posterior_var2, np.eye(n_cis_snps))
		sum_gene_ld_scores2[gene_snp_indices_arr[gene_iter]] = sum_gene_ld_scores2[gene_snp_indices_arr[gene_iter]] + gene_ld_score2

		gene_ld_score3 = get_genes_posterior_on_gene_related_snp_anno(posterior_mean3, 0.0*standardized_eqtl_posterior_var, np.eye(n_cis_snps))
		sum_gene_ld_scores3[gene_snp_indices_arr[gene_iter]] = sum_gene_ld_scores3[gene_snp_indices_arr[gene_iter]] + gene_ld_score3


	var_ld_scores = np.ones(n_snps)
	est_per_variant_h2, est_per_gene_h2 = compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, sum_gene_ld_scores2, var_ld_scores)

	est_nm_h2 = est_per_variant_h2*n_snps
	est_med_h2 = est_per_gene_h2*n_genes

	# Print to output
	t.write(str('posterior_horseshoe') + '\t' + str(sim_iter) + '\t' + str(sim_nm_var_h2) + '\t' + str(sim_med_h2) + '\t' + str(est_nm_h2) + '\t' + str(est_med_h2) + '\n')


	var_ld_scores = np.ones(n_snps)
	est_per_variant_h2, est_per_gene_h2 = compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, sum_gene_ld_scores, var_ld_scores)

	est_nm_h2 = est_per_variant_h2*n_snps
	est_med_h2 = est_per_gene_h2*n_genes

	# Print to output
	t.write(str('posterior_dirichlet_horseshoe') + '\t' + str(sim_iter) + '\t' + str(sim_nm_var_h2) + '\t' + str(sim_med_h2) + '\t' + str(est_nm_h2) + '\t' + str(est_med_h2) + '\n')


	var_ld_scores = np.ones(n_snps)
	est_per_variant_h2, est_per_gene_h2 = compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, sum_gene_ld_scores3, var_ld_scores)

	est_nm_h2 = est_per_variant_h2*n_snps
	est_med_h2 = est_per_gene_h2*n_genes

	# Print to output
	t.write(str('posterior_dirichlet_horseshoe2') + '\t' + str(sim_iter) + '\t' + str(sim_nm_var_h2) + '\t' + str(sim_med_h2) + '\t' + str(est_nm_h2) + '\t' + str(est_med_h2) + '\n')

	var_ld_scores = np.ones(n_snps)
	est_per_variant_h2, est_per_gene_h2 = compute_gene_h2_with_med_ldsc(gwas_z, gwas_ss, gene_window_based_ld_score, var_ld_scores)
	est_nm_h2 = est_per_variant_h2*n_snps
	est_med_h2 = est_per_gene_h2*n_genes
	t.write(str('window_based') + '\t' + str(sim_iter) + '\t' + str(sim_nm_var_h2) + '\t' + str(sim_med_h2) + '\t' + str(est_nm_h2) + '\t' + str(est_med_h2) + '\n')


	np.savetxt(output_root + '_sum_gene_ld_scores.txt', sum_gene_ld_scores, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_sum_gene_ld_scores2.txt', sum_gene_ld_scores2, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_gwas_z.txt', gwas_z, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_gwas_beta.txt', gwas_beta, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_gwas_beta_se.txt', gwas_beta_se, fmt="%s", delimiter='\t')

	t.flush()


t.close()
t_gene.close()


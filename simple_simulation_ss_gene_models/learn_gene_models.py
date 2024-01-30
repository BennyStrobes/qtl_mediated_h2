import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
import time
import trait_likelihood_gibbs_variant_only
import time


def extract_non_colinear_predictors(corry, abs_correlation_threshold):
	valid_predictors = []
	valid_predictors.append(0)
	for predictor in range(1,corry.shape[0]):
		if np.max(np.abs(corry[predictor, valid_predictors])) < abs_correlation_threshold:
			valid_predictors.append(predictor)
	valid_predictors = np.asarray(valid_predictors)
	return valid_predictors


def estimate_causal_eqtl_effects_for_a_single_gene(expression_vec, genotype_mat):
	n_predictors = genotype_mat.shape[1]
	corry =  np.corrcoef(np.transpose(genotype_mat))
	non_colinear_predictors = extract_non_colinear_predictors(corry, .75)

	ols = sm.OLS(expression_vec, genotype_mat[:, non_colinear_predictors])
	ols_result = ols.fit()

	# Initialize matrix
	beta = np.zeros(n_predictors)
	beta_varcov = np.zeros((n_predictors, n_predictors))
	# Fill in matrix
	beta[non_colinear_predictors] = ols_result.params
	beta_varcov[non_colinear_predictors[:,None], non_colinear_predictors[None,:]] = ols_result.cov_params()

	sampling_variance = (1.0-ols_result.rsquared_adj)/(len(expression_vec)-1)

	marginal_betas = []
	marginal_zeds = []
	for predictor_iter in range(n_predictors):
		marginal_ols = sm.OLS(expression_vec, genotype_mat[:, predictor_iter])
		marginal_ols_result = marginal_ols.fit()
		marginal_betas.append(marginal_ols_result.params[0])

		marginal_zed = marginal_ols_result.params[0]/marginal_ols_result.bse[0]
		marginal_zeds.append(marginal_zed)


	return beta, beta_varcov, ols_result.rsquared_adj,np.asarray(marginal_betas), sampling_variance, marginal_zeds


def estimate_marginal_eqtl_effects_for_a_single_gene(expression_vec, genotype_mat):
	'''
	n_predictors = genotype_mat.shape[1]
	marginal_betas = []
	marginal_beta_variances = []
	for predictor_iter in range(n_predictors):
		marginal_ols = sm.OLS(expression_vec, sm.add_constant(genotype_mat[:, predictor_iter]))
		marginal_ols_result = marginal_ols.fit()
		marginal_betas.append(marginal_ols_result.params[1])
		marginal_beta_variances.append(np.square(marginal_ols_result.bse[1]))
	marginal_betas = np.asarray(marginal_betas)
	marginal_beta_variances = np.asarray(marginal_beta_variances)
	'''
	marginal_betas = np.dot(np.transpose(genotype_mat), expression_vec)/len(expression_vec)
	marginal_beta_variances = np.var(np.transpose((genotype_mat*marginal_betas)) - expression_vec, ddof=2,axis=1)/len(expression_vec)

	return marginal_betas, marginal_beta_variances


def estimate_marginal_eqtl_effects_for_a_single_gene_from_genetic_expression(expression_vec, genotype_mat):
	n_predictors = genotype_mat.shape[1]

	if np.std(expression_vec) == 0.0:
		marginal_betas = np.zeros(n_predictors)
		return marginal_betas

	std_expression_vec = (expression_vec - np.mean(expression_vec))/np.std(expression_vec)

	'''
	marginal_betas = []
	for predictor_iter in range(n_predictors):
		marginal_ols = sm.OLS(expression_vec, genotype_mat[:, predictor_iter])
		marginal_ols_result = marginal_ols.fit()
		marginal_betas.append(marginal_ols_result.params[0])
	marginal_betas = np.asarray(marginal_betas)
	'''

	marginal_betas = np.dot(np.transpose(genotype_mat), expression_vec)/len(expression_vec)



	return marginal_betas





def estimate_heritability_of_single_gene_with_ldsc(expression_vec, genotype_mat, marginal_z):
	n_snps = genotype_mat.shape[1]
	n_indi = genotype_mat.shape[0]
	ld = np.corrcoef(np.transpose(genotype_mat))

	chi_sq = np.square(marginal_z)

	squared_ld = np.square(ld)

	adj_squared_ld = squared_ld - ((1.0 - squared_ld)/(len(expression_vec) - 2))
	#adj_squared_ld = squared_ld

	ld_scores = np.sum(adj_squared_ld,axis=0)

	model = LinearRegression(fit_intercept=True)

	ldsc_fit = model.fit(ld_scores.reshape((len(ld_scores),1)), chi_sq)

	ldsc_gene_h2 = ldsc_fit.coef_[0]*n_snps/n_indi

	return ldsc_gene_h2


def estimate_heritability_of_single_gene_with_rough_ldsc(expression_vec, genotype_mat, marginal_z):
	n_snps = genotype_mat.shape[1]
	n_indi = genotype_mat.shape[0]
	ld = np.corrcoef(np.transpose(genotype_mat))

	chi_sq = np.square(marginal_z)

	squared_ld = np.square(ld)

	adj_squared_ld = squared_ld - ((1.0 - squared_ld)/(len(expression_vec) - 2))
	#adj_squared_ld = squared_ld

	avg_ld_scores = np.mean(np.sum(adj_squared_ld,axis=0))

	rough_ldsc_gene_h2 = (np.mean(chi_sq) - 1.0)*n_snps/n_indi/avg_ld_scores

	#model = LinearRegression(fit_intercept=True, normalize=False)

	#ldsc_fit = model.fit(ld_scores.reshape((len(ld_scores),1)), chi_sq)

	#ldsc_gene_h2 = ldsc_fit.coef_[0]*n_snps/n_indi

	return rough_ldsc_gene_h2




def estimate_causal_eqtl_effect_distribution(expression_vec, genotype_mat):
	variant_only_model = trait_likelihood_gibbs_variant_only.TRAIT_LIKELIHOOD_GIBBS_VARIANT_ONLY(Y=expression_vec, X=genotype_mat, max_iter=3000, burn_in_iter=2000)
	variant_only_model.fit()

	return variant_only_model

def estimate_marginal_eqtl_effects_for_a_single_gene_from_causal_eqtl_effect_distribution(causal_eqtl_effect_distribution, LD, gene_snp_indices):
	beta_mat = np.asarray(causal_eqtl_effect_distribution.sampled_betas)
	mean_causal_effects = np.mean(beta_mat,axis=0)
	var = np.dot(mean_causal_effects, np.dot(LD[:, gene_snp_indices][gene_snp_indices,:], mean_causal_effects))
	marginal_effects = np.dot(LD[:, gene_snp_indices], mean_causal_effects)/np.sqrt(var)

	sampled_varz = np.diag(np.dot(beta_mat, np.dot(LD[:, gene_snp_indices][gene_snp_indices,:], np.transpose(beta_mat))))
	sampled_marginal_effects = np.dot(LD[:, gene_snp_indices], np.transpose(beta_mat))/np.sqrt(sampled_varz)

	return np.mean(sampled_marginal_effects,axis=1), np.square(np.mean(sampled_marginal_effects,axis=1)) - np.var(sampled_marginal_effects,axis=1)


def estimate_marginal_eqtl_effects_for_a_single_gene_from_causal_eqtl_effect_distribution_using_individual_data(causal_eqtl_effect_distribution, genotype, gene_snp_indices):
	beta_mat = np.asarray(causal_eqtl_effect_distribution.sampled_betas)

	predicted_genetic_expressions = np.dot(genotype[:, gene_snp_indices], np.transpose(beta_mat))

	sampled_varz = np.var(predicted_genetic_expressions,axis=0)

	std_predicted_genetic_expressions = predicted_genetic_expressions/np.sqrt(sampled_varz)


	sampled_marginal_effects = np.dot(np.transpose(genotype), std_predicted_genetic_expressions)/(genotype.shape[0])

	return np.mean(sampled_marginal_effects,axis=1), np.square(np.mean(sampled_marginal_effects,axis=1)) - np.var(sampled_marginal_effects,axis=1)


def estimate_rgk_unbiased_from_causal_eqtl_effect_distribution_using_individual_level_data(causal_eqtl_effect_distribution, genotype, gene_snp_indices):
	beta_mat = np.asarray(causal_eqtl_effect_distribution.sampled_betas)
	mean_causal_effects = np.mean(beta_mat,axis=0)

	predicted_genetic_expressions = np.dot(genotype[:, gene_snp_indices], np.transpose(beta_mat))

	sampled_varz = np.var(predicted_genetic_expressions,axis=0)

	std_predicted_genetic_expressions = predicted_genetic_expressions/np.sqrt(sampled_varz)

	sampled_marginal_effects = np.dot(np.transpose(genotype), std_predicted_genetic_expressions)/(genotype.shape[0])


	mean_predicted_genetic_expression = np.dot(genotype[:, gene_snp_indices], mean_causal_effects)
	variance_mean_predicted_genetic_expression = np.var(mean_predicted_genetic_expression)


	attenuation_fct = variance_mean_predicted_genetic_expression/np.var(predicted_genetic_expressions)

	mean_standardized = mean_predicted_genetic_expression/np.sqrt(variance_mean_predicted_genetic_expression)


	mean_marginal_effects = np.dot(np.transpose(genotype), mean_standardized)/genotype.shape[0]

	print(np.sqrt(attenuation_fct))


	return np.mean(sampled_marginal_effects,axis=1)/np.sqrt(attenuation_fct)



def estimate_marginal_z_scores(Y, X):
	n_predictors = X.shape[1]
	marginal_zeds = []
	for predictor_iter in range(n_predictors):
		marginal_ols = sm.OLS(Y, X[:, predictor_iter])
		marginal_ols_result = marginal_ols.fit()
		marginal_zed = marginal_ols_result.params[0]/marginal_ols_result.bse[0]
		marginal_zeds.append(marginal_zed)
	return np.asarray(marginal_zeds)


def estimate_rgk_unbiased(expression_vec, genotype, gene_snp_indices):
	# Get eqtl ss
	curr_eqtl_ss = genotype.shape[0]

	# Estimate marginal z-scores through standard eqtl analysis
	gene_snp_marginal_z = estimate_marginal_z_scores(expression_vec, genotype[:, gene_snp_indices])

	# Estimate expression heritability 
	gene_est_h2_rough_ldsc = estimate_heritability_of_single_gene_with_rough_ldsc(expression_vec, genotype[:, gene_snp_indices], gene_snp_marginal_z)

	# Estimate marginal betas for all snps
	full_marginal_betas = estimate_marginal_eqtl_effects_for_a_single_gene(expression_vec, genotype)

	# Get adjusted r-squared
	rgk_sq_unbiased = (np.square(full_marginal_betas) - ((1.0 - np.square(full_marginal_betas))/(float(curr_eqtl_ss)-2)))/gene_est_h2_rough_ldsc
	
	return rgk_sq_unbiased, gene_est_h2_rough_ldsc

def estimate_causal_eqtl_effect_distribution_ss(marginal_beta, N, gene_genotype_mat, expression_vec, update_beta_variance=False, beta_var_init=1.0):
	local_ld = np.corrcoef(np.transpose(gene_genotype_mat))

	#mod = rss_gibbs_variant_only.RSS_GIBBS_VARIANT_ONLY(LD=local_ld, marginal_beta=marginal_beta, N=N, X=gene_genotype_mat, Y=expression_vec, update_beta_variance=update_beta_variance, beta_var_init=beta_var_init, max_iter=100, burn_in_iter=60)
	mod = rss_gibbs_variant_only.RSS_GIBBS_VARIANT_ONLY(LD=local_ld, marginal_beta=marginal_beta, N=N, X=gene_genotype_mat, Y=expression_vec, update_beta_variance=update_beta_variance, beta_var_init=beta_var_init, max_iter=9000, burn_in_iter=6000)


	mod.fit()
	return mod

def get_rgk_sq_unbiased(causal_eqtl_effects_mean, gene_cov, LD, gene_snp_indices):
	R_gene = LD[gene_snp_indices, :][:, gene_snp_indices]
	# Calculate gene variance
	gene_meanT_mean = np.dot(causal_eqtl_effects_mean.reshape(len(causal_eqtl_effects_mean),1), causal_eqtl_effects_mean.reshape(1,len(causal_eqtl_effects_mean)))


	unbiased_delta_t_delta = gene_meanT_mean - gene_cov
	tot_gene_variance = np.sum(R_gene*unbiased_delta_t_delta)


	n_snps = LD.shape[0]
	sub_ld = LD[:, gene_snp_indices]
	rgk_sq = np.zeros(n_snps)
	for kk in range(n_snps):
		tmp = np.dot(np.dot(sub_ld[kk,:], (unbiased_delta_t_delta/tot_gene_variance)), sub_ld[kk,:])
		rgk_sq[kk] = rgk_sq[kk] + tmp
	
	return rgk_sq

def calculate_gene_ld_score_noise_according_to_gaussian_distr(beta_hat_vec, beta_hat_var_vec, n_samples):
	unbiased_beta_hat_squared = np.square(beta_hat_vec) - beta_hat_var_vec
	unbiased_beta_var = np.copy(beta_hat_var_vec)
	unbiased_beta_var_squared = np.square(unbiased_beta_var)/(1.0 + (2/(n_samples-2.0)))
	unbiased_var_beta_times_beta_sq = unbiased_beta_hat_squared*unbiased_beta_var - ((-2.0*n_samples + 4)*unbiased_beta_var_squared/(np.power(n_samples-2,2)))
	unbiased_var_beta_hat_squared = (2.0*unbiased_beta_var_squared) + (4.0*unbiased_var_beta_times_beta_sq)
	# Get unbiased var(var(\hat{\beta}))
	unbiased_var_var_beta_hat = 2.0*unbiased_beta_var_squared/(n_samples-2)
	# Get total predictor noise
	predictor_noise = unbiased_var_beta_hat_squared +unbiased_var_var_beta_hat
	return predictor_noise


def posterior_distribution_on_gene_related_snp_annotation_known_gene_cis_h2_known_resid_var(LD, full_marginal_betas, gene_snp_indices, ge_h2, eqtl_ss, residual_variance_est):
	# Get residual variance
	if residual_variance_est == 'one_minus_ge_h2':
		resid_var_est = 1.0 - ge_h2
	else:
		print('assumption error: this residual variance estimate method is not currently implemented')

	# Subset LD and marginal betas to gene snp indices
	sub_ld = LD[:, gene_snp_indices][gene_snp_indices,:]
	sub_marginal_betas = full_marginal_betas[gene_snp_indices]
	MM = len(gene_snp_indices)

	tmp_inv_var = np.linalg.inv((sub_ld/resid_var_est) + (np.eye(MM)*MM/(eqtl_ss*ge_h2)))
	posterior_mean = np.dot(tmp_inv_var, sub_marginal_betas)
	posterior_var = tmp_inv_var/eqtl_ss

	e_delta_delta_t = np.dot(posterior_mean.reshape((len(posterior_mean),1)), posterior_mean.reshape((1,len(posterior_mean)))) + posterior_var
	#posterior_gene_related_snp_annotation_slow = np.diag(np.dot(np.dot(LD[:, gene_snp_indices], e_delta_delta_t), LD[gene_snp_indices,:]))
	posterior_gene_related_snp_annotation = np.sum(np.dot(LD[:, gene_snp_indices], e_delta_delta_t)*LD[:, gene_snp_indices],axis=1)

	return posterior_gene_related_snp_annotation


def get_genes_posterior_on_gene_related_snp_anno(posterior_mean, posterior_var, LD, gene_snp_indices):
	e_delta_delta_t = np.dot(posterior_mean.reshape((len(posterior_mean),1)), posterior_mean.reshape((1,len(posterior_mean)))) + posterior_var
	#posterior_gene_related_snp_annotation_slow = np.diag(np.dot(np.dot(LD[:, gene_snp_indices], e_delta_delta_t), LD[gene_snp_indices,:]))
	posterior_gene_related_snp_annotation = np.sum(np.dot(LD[:, gene_snp_indices], e_delta_delta_t)*LD[:, gene_snp_indices],axis=1)

	return posterior_gene_related_snp_annotation

def posterior_distribution_on_causal_eqtl_effects_closed_form_ss(LD, full_marginal_betas, gene_snp_indices, ge_h2, eqtl_ss):
	resid_var_est = 1.0 - ge_h2
	# Subset LD and marginal betas to gene snp indices
	sub_ld = LD[:, gene_snp_indices][gene_snp_indices,:]
	sub_marginal_betas = full_marginal_betas[gene_snp_indices]
	MM = len(gene_snp_indices)

	tmp_inv_var = np.linalg.inv((sub_ld/resid_var_est) + (np.eye(MM)*MM/(eqtl_ss*ge_h2)))
	posterior_mean = np.dot(tmp_inv_var, sub_marginal_betas)
	posterior_var = tmp_inv_var/eqtl_ss
	return posterior_mean, posterior_var

def posterior_distribution_on_causal_eqtl_effects_gibbs(expression_vec, genotype_mat, ge_h2, residual_variance):
	if ge_h2 is not None:
		per_snp_h2 = ge_h2/genotype_mat.shape[1]
	else:
		per_snp_h2 = None
	try: 
		gibbsy = trait_likelihood_gibbs_variant_only.TRAIT_LIKELIHOOD_GIBBS_VARIANT_ONLY(Y=expression_vec, X=genotype_mat, beta_variance=per_snp_h2, residual_variance=residual_variance, max_iter=2000, burn_in_iter=1000, beta_var_prior_scale=1e-3)
		gibbsy.fit()
	except:
		try:
			print('try 2')
			gibbsy = trait_likelihood_gibbs_variant_only.TRAIT_LIKELIHOOD_GIBBS_VARIANT_ONLY(Y=expression_vec, X=genotype_mat, beta_variance=per_snp_h2, residual_variance=residual_variance, max_iter=2000, burn_in_iter=1000, beta_var_prior_scale=1e-3)
			gibbsy.fit()
		except:
			print('try 3')
			gibbsy = trait_likelihood_gibbs_variant_only.TRAIT_LIKELIHOOD_GIBBS_VARIANT_ONLY(Y=expression_vec, X=genotype_mat, beta_variance=per_snp_h2, residual_variance=residual_variance, max_iter=2000, burn_in_iter=1000, beta_var_prior_scale=1e-3)
			gibbsy.fit()			

	print(np.mean(gibbsy.sampled_residual_variances))
	print(np.mean(gibbsy.sampled_beta_variances))
	return np.mean(gibbsy.sampled_betas,axis=0), np.cov(np.transpose(gibbsy.sampled_betas)), np.mean(gibbsy.sampled_beta_variances)*genotype_mat.shape[1]


def rescale_genotype_mat(genotype_mat):
	n_snps = genotype_mat.shape[1]
	# Calculate range of each snp
	ranges = []
	for snp_iter in range(n_snps):
		snp_range = np.max(genotype_mat[:,snp_iter]) - np.min(genotype_mat[:,snp_iter])
		ranges.append(snp_range)
	ranges = np.asarray(ranges)

	# Get biggest range
	max_range = np.max(ranges)

	# Rescale
	rescaled_genotype_mat = 2.0*np.copy(genotype_mat)/max_range

	# Uncenter
	for snp_iter in range(n_snps):
		rescaled_genotype_mat[:, snp_iter] = rescaled_genotype_mat[:, snp_iter] - np.min(rescaled_genotype_mat[:, snp_iter])

	scale_factor = max_range/2.0

	'''
	# Verify that if we mean center and standardize rescaled_genotype_mat we get genotype_mat
	tmps = []
	sdevs = []
	for snp_iter in range(n_snps):
		sdevs.append(np.std(rescaled_genotype_mat[:, snp_iter]))
		tmp = (rescaled_genotype_mat[:, snp_iter] - np.mean(rescaled_genotype_mat[:, snp_iter]))/np.std(rescaled_genotype_mat[:, snp_iter])
		tmps.append(tmp)
	tmps = np.transpose(np.asarray(tmps))
	sdevs = np.asarray(sdevs)
	# Verify that if we mean center and standardize rescaled_genotype_mat we get genotype_mat
	tmps = []
	for snp_iter in range(n_snps):
		tmp = (rescaled_genotype_mat[:, snp_iter] - np.mean(rescaled_genotype_mat[:, snp_iter]))*scale_factor
		tmps.append(tmp)
	tmps = np.transpose(np.asarray(tmps))
	'''
	return rescaled_genotype_mat, scale_factor

def create_bimbam_mean_genotype_file(rescaled_genotype_mat, bimbam_mean_genotype_file):
	t = open(bimbam_mean_genotype_file,'w')
	n_snps = rescaled_genotype_mat.shape[1]
	for snp_iter in range(n_snps):
		snp_vec = rescaled_genotype_mat[:, snp_iter]
		t.write('rs' + str(snp_iter) + ',A,C,' + ','.join(snp_vec.astype(str)) + '\n')
	t.close()
	return

def create_bimbam_phenotype_file(expression_vec, bimbam_phenotype_file):
	t = open(bimbam_phenotype_file,'w')
	for ele in expression_vec:
		t.write(str(ele) + '\n')
	t.close()
	return


def bslmm_posterior_distribution_on_causal_eqtl_effects(expression_vec, genotype_mat, gemma_path, tmp_file_root):
	# Rescale genotype matrix so that each snp has the same variance, however the snp values are bounded between 0 and 2
	# Note that centering each snp of the rescaled_genotype_mat, and multiplying by the genotype_scale_factor will result in the original genotype_mat
	# Where the original_genotype_mat is mean centered and scaled.
	rescaled_genotype_mat, genotype_scale_factor = rescale_genotype_mat(genotype_mat)

	# Create bimbam mean genotype file
	bimbam_mean_genotype_file = tmp_file_root + 'bimbam_mean_genotype_file.txt'
	create_bimbam_mean_genotype_file(rescaled_genotype_mat, bimbam_mean_genotype_file)

	# Create bimbam phenotype file
	bimbam_phenotype_file = tmp_file_root + 'bimbam_phenotype_file.txt'
	create_bimbam_phenotype_file(expression_vec, bimbam_phenotype_file)

	bslmm_output_root = 'bslmm_' 
	bslmm_cmd = gemma_path + ' -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000 -g ' + bimbam_mean_genotype_file + ' -p ' + bimbam_phenotype_file + ' -bslmm 1 -o ' + bslmm_output_root
	os.system(bslmm_cmd)
	return

def get_pve_est_from_blup_log_file(filer):
	pve_est = -2.0
	pve_se_est = -2.0
	f = open(filer)
	for line in f:
		line = line.rstrip()
		if line.startswith('## pve estimate in the null model ='):
			pve_est = float(line.split('## pve estimate in the null model = ')[1])
		if line.startswith('## se(pve) in the null model ='):
			pve_se_est = float(line.split('## se(pve) in the null model = ')[1])
	f.close()
	if pve_est == -2.0 or pve_se_est == -2.0:
		print('assumption erororo')
		pdb.set_trace()
	return pve_est, pve_se_est

def blup_posterior_distribution_on_causal_eqtl_effects(expression_vec, genotype_mat, gemma_path, tmp_file_dir, tmp_file_stem):
	# Rescale genotype matrix so that each snp has the same variance, however the snp values are bounded between 0 and 2
	# Note that centering each snp of the rescaled_genotype_mat, and multiplying by the genotype_scale_factor will result in the original genotype_mat
	# Where the original_genotype_mat is mean centered and scaled.
	rescaled_genotype_mat, genotype_scale_factor = rescale_genotype_mat(genotype_mat)

	# Create bimbam mean genotype file
	bimbam_mean_genotype_file = tmp_file_dir + tmp_file_stem + 'bimbam_mean_genotype_file.txt'
	create_bimbam_mean_genotype_file(rescaled_genotype_mat, bimbam_mean_genotype_file)

	# Create bimbam phenotype file
	bimbam_phenotype_file = tmp_file_dir + tmp_file_stem + 'bimbam_phenotype_file.txt'
	create_bimbam_phenotype_file(expression_vec, bimbam_phenotype_file)

	blup_cmd = gemma_path + ' -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000 -g ' + bimbam_mean_genotype_file + ' -p ' + bimbam_phenotype_file + ' -bslmm 2 -outdir ' + tmp_file_dir + ' -o ' + tmp_file_stem
	os.system(blup_cmd)

	filer = tmp_file_dir + tmp_file_stem + '.log.txt'
	pve_est, pve_se_est = get_pve_est_from_blup_log_file(filer)

	# Remove unncessary files
	#os.system('rm ' + tmp_file_dir + tmp_file_stem + '*')

	'''
	n_snps = rescaled_genotype_mat.shape[1]
	rescaled_centered_genotype_mat = np.copy(rescaled_genotype_mat)
	for snp_iter in range(n_snps):
		rescaled_centered_genotype_mat[:, snp_iter] = rescaled_centered_genotype_mat[:, snp_iter] - np.mean(rescaled_centered_genotype_mat[:, snp_iter])
	pdb.set_trace()
	tmp_prec = (np.eye(n_snps)/(pve_est/n_snps)) + (np.dot(np.transpose(genotype_mat), genotype_mat)/(1.0-pve_est))
	tmp_var = np.linalg.inv(tmp_prec)
	tmp_beta = (1.0/(1.0-pve_est))*np.dot(np.dot(tmp_var, np.transpose(genotype_mat)), expression_vec)
	blup_est = (np.loadtxt('/n/groups/price/ben/code_temp6/output/blup_.param.txt',dtype=str,delimiter='\t')[1:,4]).astype(float)
	'''

	return pve_est

def posterior_distribution_on_causal_eqtl_effects_closed_form(expression_vec, genotype_mat, ge_h2):
	n_snps = genotype_mat.shape[1]
	resid_var = 1.0 - ge_h2
	per_snp_h2 = ge_h2/n_snps 

	tmp_prec = (np.eye(n_snps)/per_snp_h2) + (np.dot(np.transpose(genotype_mat), genotype_mat)/(resid_var))
	posterior_var = np.linalg.inv(tmp_prec)

	posterior_mean = (1.0/resid_var)*np.dot(np.dot(posterior_var, np.transpose(genotype_mat)), expression_vec)

	return posterior_mean, posterior_var




######################
# Command line args
######################
simulation_name_string = sys.argv[1]
processed_genotype_data_dir = sys.argv[2]
simulated_eqtl_data_dir = sys.argv[3]
simulated_expression_data_dir = sys.argv[4]
simulated_gene_models_dir = sys.argv[5]
eqtl_ss = sys.argv[6]
ge_h2 = float(sys.argv[7])


gemma_path="/n/groups/price/tiffany/subpheno/gemma-0.98.1-linux-static"


# Load in genotype matrix
genotype_file = processed_genotype_data_dir + 'eqtl_' + str(eqtl_ss) + '_genotype_1.npy'
genotype = np.load(genotype_file)
LD = np.corrcoef(np.transpose(genotype))

# Load in gene expression
expression_file = simulated_expression_data_dir + simulation_name_string + 'eqtl_ss_' + str(eqtl_ss) + '.npy'
gene_expression = np.load(expression_file)
# load in genetic gene expression
genetic_expression_file = simulated_expression_data_dir + simulation_name_string + 'eqtl_ss_' + str(eqtl_ss) + '_genetic_ge.npy'
genetic_gene_expression = np.load(genetic_expression_file)


# Initialize output vectors
n_snps = LD.shape[0]
noisy_gene_related_snp_anno_vec = np.zeros(n_snps)
noisy_gene_related_snp_anno_noise_vec = np.zeros(n_snps)
marginal_beta_genetic_ge_vec = []

posterior_known_h2_known_resid_var_gene_related_snp_anno_vec = np.zeros(n_snps)
posterior_known_h2_unknown_resid_var_gene_related_snp_anno_vec = np.zeros(n_snps)
posterior_unknown_h2_unknown_resid_var_gene_related_snp_anno_vec = np.zeros(n_snps)
posterior_blup_related_snp_anno_vec = np.zeros(n_snps)


est_gene_h2 = []
blup_pves = []

simulated_causal_eqtl_effect_file = simulated_eqtl_data_dir + simulation_name_string + 'causal_eqtl_effect_summary.txt'
f = open(simulated_causal_eqtl_effect_file)
head_count = 0
gene_index = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract relevent fields
	gene_name = data[0]

	print(gene_name)
	#print(gene_index)
	#print(time.time())
	gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
	gene_snp_causal_effects = np.asarray(data[2].split(',')).astype(float)
	if int(gene_name.split('_')[1]) != gene_index:
		print('assumption eroror')

	expression_vec = gene_expression[gene_index,:]
	genetic_expression_vec = genetic_gene_expression[gene_index,:]

	#################
	# Get gene related snp annotations by knowing true genetic genetic gene expression
	#################
	# Get marginal betas from simulated, unknown genetic gene expression
	full_marginal_betas_genetic_ge = estimate_marginal_eqtl_effects_for_a_single_gene_from_genetic_expression(genetic_expression_vec, genotype)
	marginal_beta_genetic_ge_vec.append(full_marginal_betas_genetic_ge)
	print(np.sum(np.square(gene_snp_causal_effects))/len(gene_snp_causal_effects))
	print(np.var(genetic_expression_vec,ddof=1)/len(gene_snp_causal_effects))

	#################
	# Get gene related snp annotations by getting noisy, unbiased estimates of them
	#################
	# Estimate marginal betas
	full_marginal_betas, full_marginal_beta_variances = estimate_marginal_eqtl_effects_for_a_single_gene(expression_vec, genotype)
	# Estimate gene_related_snp_annotations
	full_adj_r_squared = np.square(full_marginal_betas) - full_marginal_beta_variances
	noisy_gene_related_snp_anno_vec = noisy_gene_related_snp_anno_vec + full_adj_r_squared
	# Estimate noise in gene related snp annotations
	full_adj_r_squared_noise = calculate_gene_ld_score_noise_according_to_gaussian_distr(full_marginal_betas, full_marginal_beta_variances, float(eqtl_ss))
	noisy_gene_related_snp_anno_noise_vec = noisy_gene_related_snp_anno_noise_vec + full_adj_r_squared_noise


	#################
	# Get posterior distribution on gene-related snp-annotation by knowing h2 and residual variance
	#################
	#eqtl_posterior_mean_known_h2_known_resid_var, eqtl_posterior_var_known_h2_known_resid_var = posterior_distribution_on_causal_eqtl_effects_closed_form_ss(LD, full_marginal_betas, gene_snp_indices, ge_h2, float(eqtl_ss))
	#gene_posterior_gene_related_snp_anno_known_h2_known_resid_var = get_genes_posterior_on_gene_related_snp_anno(eqtl_posterior_mean_known_h2_known_resid_var, eqtl_posterior_var_known_h2_known_resid_var, LD, gene_snp_indices)
	#posterior_known_h2_known_resid_var_gene_related_snp_anno_vec = posterior_known_h2_known_resid_var_gene_related_snp_anno_vec + gene_posterior_gene_related_snp_anno_known_h2_known_resid_var
	eqtl_posterior_mean_known_h2_known_resid_var, eqtl_posterior_var_known_h2_known_resid_var = posterior_distribution_on_causal_eqtl_effects_closed_form(expression_vec, genotype[:,gene_snp_indices], ge_h2)
	gene_posterior_gene_related_snp_anno_known_h2_known_resid_var = get_genes_posterior_on_gene_related_snp_anno(eqtl_posterior_mean_known_h2_known_resid_var, eqtl_posterior_var_known_h2_known_resid_var, LD, gene_snp_indices)
	posterior_known_h2_known_resid_var_gene_related_snp_anno_vec = posterior_known_h2_known_resid_var_gene_related_snp_anno_vec + gene_posterior_gene_related_snp_anno_known_h2_known_resid_var

	#################
	# Get BLUP Posterior
	#################
	tmp_file_stem = simulation_name_string + '_' + str(eqtl_ss) + '_tmp_blup_'
	blup_pve = blup_posterior_distribution_on_causal_eqtl_effects(expression_vec, genotype[:,gene_snp_indices], gemma_path, simulated_gene_models_dir, tmp_file_stem)
	eqtl_posterior_mean_blup, eqtl_posterior_var_blup = posterior_distribution_on_causal_eqtl_effects_closed_form(expression_vec, genotype[:,gene_snp_indices], blup_pve)
	gene_posterior_gene_related_snp_anno_blup = get_genes_posterior_on_gene_related_snp_anno(eqtl_posterior_mean_blup, eqtl_posterior_var_blup, LD, gene_snp_indices)
	posterior_blup_related_snp_anno_vec = posterior_blup_related_snp_anno_vec + gene_posterior_gene_related_snp_anno_blup
	blup_pves.append(blup_pve)


	#################
	# Get BSLMM Posterior
	#################
	#tmp_file_root = simulated_gene_models_dir + simulation_name_string + '_' + str(eqtl_ss) + '_' + str(gene_index)
	#bslmm_posterior_mean, bslmm_posterior_var = bslmm_posterior_distribution_on_causal_eqtl_effects(expression_vec, genotype[:,gene_snp_indices], gemma_path, tmp_file_root)
	#bslmm_posterior_distribution_on_causal_eqtl_effects(expression_vec, genotype[:,gene_snp_indices], gemma_path, tmp_file_root)

	'''
	#################
	# Get posterior distribution on gene-related snp-annotation by knowing h2 and not knowing residual variance
	#################
	eqtl_posterior_mean_known_h2_unknown_resid_var, eqtl_posterior_var_known_h2_unknown_resid_var, est_gene_h2_known_h2_unknown_resid_var  = posterior_distribution_on_causal_eqtl_effects_gibbs(expression_vec, genotype[:,gene_snp_indices], ge_h2, None)
	gene_posterior_gene_related_snp_anno_known_h2_unknown_resid_var = get_genes_posterior_on_gene_related_snp_anno(eqtl_posterior_mean_known_h2_unknown_resid_var, eqtl_posterior_var_known_h2_unknown_resid_var, LD, gene_snp_indices)
	posterior_known_h2_unknown_resid_var_gene_related_snp_anno_vec = posterior_known_h2_unknown_resid_var_gene_related_snp_anno_vec + gene_posterior_gene_related_snp_anno_known_h2_unknown_resid_var

	#################
	# Get posterior distribution on gene-related snp-annotation by not knowning h2 and not knowing residual variance
	#################
	eqtl_posterior_mean_unknown_h2_unknown_resid_var, eqtl_posterior_var_unknown_h2_unknown_resid_var, est_gene_h2_unknown_h2_unknown_resid_var  = posterior_distribution_on_causal_eqtl_effects_gibbs(expression_vec, genotype[:,gene_snp_indices], None, None)
	gene_posterior_gene_related_snp_anno_unknown_h2_unknown_resid_var = get_genes_posterior_on_gene_related_snp_anno(eqtl_posterior_mean_unknown_h2_unknown_resid_var, eqtl_posterior_var_unknown_h2_unknown_resid_var, LD, gene_snp_indices)
	posterior_unknown_h2_unknown_resid_var_gene_related_snp_anno_vec = posterior_unknown_h2_unknown_resid_var_gene_related_snp_anno_vec + gene_posterior_gene_related_snp_anno_unknown_h2_unknown_resid_var
	est_gene_h2.append(est_gene_h2_unknown_h2_unknown_resid_var)
	'''

	gene_index = gene_index + 1



f.close()

# Covert list of vectors to matrix
marginal_beta_genetic_ge_vec = np.asarray(marginal_beta_genetic_ge_vec)

#est_gene_h2 = np.asarray(est_gene_h2)
blup_pves = np.asarray(blup_pves)


# Save genetic gene expression marginal betas
genetic_ge_marginal_beta_file = simulated_gene_models_dir + simulation_name_string + 'estimated_genetic_ge_eqtl_marginal_beta_' + str(eqtl_ss) + '.npy'
np.save(genetic_ge_marginal_beta_file, marginal_beta_genetic_ge_vec)

# Save estimated noisy gene related snp annotations
est_noisy_gene_related_snp_anno_file = simulated_gene_models_dir + simulation_name_string + 'estimated_noisy_gene_related_snp_annotations_' + str(eqtl_ss) + '.npy'
np.save(est_noisy_gene_related_snp_anno_file, noisy_gene_related_snp_anno_vec)

# Save estimated noisy gene related snp annotation noises
est_noisy_gene_related_snp_anno_noise_file = simulated_gene_models_dir + simulation_name_string + 'estimated_noisy_gene_related_snp_annotation_noises_' + str(eqtl_ss) + '.npy'
np.save(est_noisy_gene_related_snp_anno_noise_file, noisy_gene_related_snp_anno_noise_vec)

# Save posterior estimated gene related snp annotation noises (known gene h2, known resid var)
posterior_known_h2_known_resid_var_gene_related_snp_anno_file = simulated_gene_models_dir + simulation_name_string + 'posterior_known_h2_known_resid_var_gene_related_snp_annotations_' + str(eqtl_ss) + '.npy'
np.save(posterior_known_h2_known_resid_var_gene_related_snp_anno_file, posterior_known_h2_known_resid_var_gene_related_snp_anno_vec)

'''
# Save posterior estimated gene related snp annotation noises (known gene h2, unknown resid var)
posterior_known_h2_unknown_resid_var_gene_related_snp_anno_file = simulated_gene_models_dir + simulation_name_string + 'posterior_known_h2_unknown_resid_var_gene_related_snp_annotations_' + str(eqtl_ss) + '.npy'
np.save(posterior_known_h2_unknown_resid_var_gene_related_snp_anno_file, posterior_known_h2_unknown_resid_var_gene_related_snp_anno_vec)

# Save posterior estimated gene related snp annotation noises (unknown gene h2, unknown resid var)
posterior_unknown_h2_unknown_resid_var_gene_related_snp_anno_file = simulated_gene_models_dir + simulation_name_string + 'posterior_unknown_h2_unknown_resid_var_gene_related_snp_annotations_' + str(eqtl_ss) + '.npy'
np.save(posterior_unknown_h2_unknown_resid_var_gene_related_snp_anno_file, posterior_unknown_h2_unknown_resid_var_gene_related_snp_anno_vec)
'''

# Save posterior estimated gene related snp annotation noises (BLUP)
posterior_blup_gene_related_snp_anno_file = simulated_gene_models_dir + simulation_name_string + 'posterior_blup_gene_related_snp_annotations_' + str(eqtl_ss) + '.npy'
np.save(posterior_blup_gene_related_snp_anno_file, posterior_blup_related_snp_anno_vec)

# Save estimated gene h2 (unknown gene h2, unknown resid var)
est_gene_h2_unknown_h2_unknown_resid_var_file = simulated_gene_models_dir + simulation_name_string + 'est_gene_h2_unknown_h2_unknown_resid_var_' + str(eqtl_ss) + '.npy'
np.save(est_gene_h2_unknown_h2_unknown_resid_var_file, est_gene_h2)


# Save estimated gene h2 (BLUP PVES)
est_gene_h2_blup_file = simulated_gene_models_dir + simulation_name_string + 'est_gene_h2_blup_' + str(eqtl_ss) + '.npy'
np.save(est_gene_h2_blup_file, blup_pves)


print("DONE")





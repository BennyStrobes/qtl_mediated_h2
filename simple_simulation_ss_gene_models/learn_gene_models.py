import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
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
	n_predictors = genotype_mat.shape[1]

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


######################
# Command line args
######################
simulation_name_string = sys.argv[1]
processed_genotype_data_dir = sys.argv[2]
simulated_eqtl_data_dir = sys.argv[3]
simulated_expression_data_dir = sys.argv[4]
simulated_gene_models_dir = sys.argv[5]
eqtl_ss = sys.argv[6]

n_bootstraps = 10


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


aa = []
bb = []
cc = []
# Open output file handle
output_file = simulated_gene_models_dir + simulation_name_string + 'model_summaries_' + str(eqtl_ss) + '.txt'
t = open(output_file,'w')
t.write('gene_name\tgene_snp_indices\tgene_snp_causal_eqtl_effects\tmarginal_beta_file\tadj_gene_r_squared_file\tgenetic_ge_marginal_beta_file\n')
# Loop through genes
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
	#print(gene_index)
	#print(time.time())
	gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
	gene_snp_causal_effects = np.asarray(data[2].split(',')).astype(float)
	if int(gene_name.split('_')[1]) != gene_index:
		print('assumption eroror')

	expression_vec = gene_expression[gene_index,:]
	genetic_expression_vec = genetic_gene_expression[gene_index,:]

	#beta, beta_varcov, gene_h2_est, marginal_betas, gene_sampling_var, marginal_z = estimate_causal_eqtl_effects_for_a_single_gene(expression_vec, genotype[:, gene_snp_indices])
	full_marginal_betas = estimate_marginal_eqtl_effects_for_a_single_gene(expression_vec, genotype)
	full_marginal_betas_genetic_ge = estimate_marginal_eqtl_effects_for_a_single_gene_from_genetic_expression(genetic_expression_vec, genotype)
	full_adj_r_squared = np.square(full_marginal_betas) - ((1.0 - np.square(full_marginal_betas))/(float(eqtl_ss)-2))

	'''
	#gene_est_h2_ldsc = estimate_heritability_of_single_gene_with_ldsc(expression_vec, genotype[:, gene_snp_indices], marginal_z)
	marginal_betas = full_marginal_betas[gene_snp_indices]

	gene_genotype_mat = genotype[:, gene_snp_indices]

	causal_eqtl_effect_distribution = estimate_causal_eqtl_effect_distribution_ss(marginal_betas, gene_genotype_mat.shape[0], gene_genotype_mat, expression_vec)
	sampled_betas = np.asarray(causal_eqtl_effect_distribution.sampled_betas)
	causal_eqtl_effects_mean = np.mean(sampled_betas,axis=0)
	causal_eqtl_effects_cov = np.cov(np.transpose(sampled_betas))


	pred_expr = np.dot(gene_genotype_mat, causal_eqtl_effects_mean)


	R_gene = np.corrcoef(np.transpose(gene_genotype_mat))

	pred_expr_var = np.sum(R_gene*np.dot(causal_eqtl_effects_mean.reshape(len(causal_eqtl_effects_mean),1), causal_eqtl_effects_mean.reshape(1,len(causal_eqtl_effects_mean))))
	est_noise_var = np.sum(R_gene*causal_eqtl_effects_cov)


	#est_noise_var = np.mean(np.var(pred_expr_distr,axis=1,ddof=1))
	diff = (genetic_expression_vec - np.dot(gene_genotype_mat, causal_eqtl_effects_mean))
	true_noise_var = np.var(diff,ddof=1)

	est_noise_ratio = est_noise_var/pred_expr_var
	true_noise_ratio = true_noise_var/pred_expr_var



	blup_causal_eqtl_effect_distribution = estimate_causal_eqtl_effect_distribution_ss(marginal_betas, gene_genotype_mat.shape[0], gene_genotype_mat, expression_vec, update_beta_variance=True)
	blup_sampled_betas = np.asarray(blup_causal_eqtl_effect_distribution.sampled_betas)
	blup_causal_eqtl_effects_mean = np.mean(blup_sampled_betas,axis=0)
	blup_pred_expr = np.dot(gene_genotype_mat, blup_causal_eqtl_effects_mean)
	blup_beta_variance = np.mean(np.mean(blup_causal_eqtl_effect_distribution.sampled_beta_variances))


	ridge_regr_causal_eqtl_effect_distribution = estimate_causal_eqtl_effect_distribution_ss(marginal_betas, gene_genotype_mat.shape[0], gene_genotype_mat, expression_vec, update_beta_variance=False, beta_var_init=1.0/len(marginal_betas))
	ridge_regr_sampled_betas = np.asarray(ridge_regr_causal_eqtl_effect_distribution.sampled_betas)
	ridge_regr_causal_eqtl_effects_mean = np.mean(ridge_regr_sampled_betas,axis=0)
	ridge_regr_pred_expr = np.dot(gene_genotype_mat, ridge_regr_causal_eqtl_effects_mean)


	ridge_regr_true_var_causal_eqtl_effect_distribution = estimate_causal_eqtl_effect_distribution_ss(marginal_betas, gene_genotype_mat.shape[0], gene_genotype_mat, expression_vec, update_beta_variance=False, beta_var_init=0.05/len(marginal_betas))
	ridge_regr_true_var_sampled_betas = np.asarray(ridge_regr_true_var_causal_eqtl_effect_distribution.sampled_betas)
	ridge_regr_true_var_causal_eqtl_effects_mean = np.mean(ridge_regr_true_var_sampled_betas,axis=0)
	ridge_regr_true_var_pred_expr = np.dot(gene_genotype_mat, ridge_regr_true_var_causal_eqtl_effects_mean)


	'''

	# Save marginal beta file
	marginal_beta_file = simulated_gene_models_dir + simulation_name_string + 'estimated_eqtl_marginal_beta_' + str(eqtl_ss) + '_' + gene_name + '.npy'
	np.save(marginal_beta_file, full_marginal_betas)
	# Save marginal beta file from genetic ge
	genetic_ge_marginal_beta_file = simulated_gene_models_dir + simulation_name_string + 'estimated_genetic_ge_eqtl_marginal_beta_' + str(eqtl_ss) + '_' + gene_name + '.npy'
	np.save(genetic_ge_marginal_beta_file, full_marginal_betas_genetic_ge)
	# Save adjusted r-squared file
	adj_r_squared_file = simulated_gene_models_dir + simulation_name_string + 'estimated_adj_r_squared_' + str(eqtl_ss) + '_' + gene_name + '.npy'
	np.save(adj_r_squared_file, full_adj_r_squared)
	'''
	# Save estimated causal eqtl effects
	estimated_causal_eqtl_effects_file = simulated_gene_models_dir + simulation_name_string + 'estimated_causal_eqtl_effects_mean_' + str(eqtl_ss) + '_' + gene_name + '.npy'
	np.save(estimated_causal_eqtl_effects_file, causal_eqtl_effects_mean)
	# Save estimated causal eqtl effects
	ridge_regr_estimated_causal_eqtl_effects_file = simulated_gene_models_dir + simulation_name_string + 'ridge_regr_estimated_causal_eqtl_effects_mean_' + str(eqtl_ss) + '_' + gene_name + '.npy'
	np.save(ridge_regr_estimated_causal_eqtl_effects_file, ridge_regr_causal_eqtl_effects_mean)
	# Save estimated causal eqtl effects
	ridge_regr_true_var_estimated_causal_eqtl_effects_file = simulated_gene_models_dir + simulation_name_string + 'ridge_regr_true_var_estimated_causal_eqtl_effects_mean_' + str(eqtl_ss) + '_' + gene_name + '.npy'
	np.save(ridge_regr_true_var_estimated_causal_eqtl_effects_file, ridge_regr_true_var_causal_eqtl_effects_mean)
	# Save estimated causal eqtl effects
	blup_estimated_causal_eqtl_effects_file = simulated_gene_models_dir + simulation_name_string + 'blup_estimated_causal_eqtl_effects_mean_' + str(eqtl_ss) + '_' + gene_name + '.npy'
	np.save(blup_estimated_causal_eqtl_effects_file, blup_causal_eqtl_effects_mean)
	'''

	t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + str(marginal_beta_file) + '\t' + adj_r_squared_file + '\t' + genetic_ge_marginal_beta_file + '\n')

	gene_index = gene_index + 1

f.close()
t.close()

print(output_file)






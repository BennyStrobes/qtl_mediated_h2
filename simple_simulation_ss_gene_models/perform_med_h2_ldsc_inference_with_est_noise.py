import numpy as np 
import os
import sys
import pdb
from sklearn.linear_model import LinearRegression




def load_in_gwas_z_scores(gwas_sum_stats_file):
	f = open(gwas_sum_stats_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(float(data[3]))
	f.close()
	arr = np.asarray(arr)
	return arr

def load_in_gwas_se(gwas_sum_stats_file):
	f = open(gwas_sum_stats_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(float(data[2]))
	f.close()
	arr = np.asarray(arr)
	return arr

def load_in_trait_file(trait_file):
	f = open(trait_file)
	arr = []
	for line in f:
		line = line.rstrip()
		arr.append(float(line))
	f.close()
	return np.asarray(arr)

def extract_causal_eqtl_effects(gene_summary_file):
	causal_eqtl_indices = []
	true_causal_eqtl_effects = []
	head_count = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
		true_causal_eqtl_effect = np.asarray(data[2].split(',')).astype(float)

		causal_eqtl_indices.append(gene_snp_indices)
		true_causal_eqtl_effects.append(true_causal_eqtl_effect)
	f.close()
	return causal_eqtl_indices, true_causal_eqtl_effects



def extract_previously_estimsted_causal_eqtl_effects(gene_summary_file):
	est_causal_eqtl_effects = []
	causal_eqtl_indices = []
	true_causal_eqtl_effects = []
	blup_causal_eqtl_effects = []
	ridge_regr_causal_eqtl_effects = []
	ridge_regr_true_var_causal_eqtl_effects = []
	head_count = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
		true_causal_eqtl_effect = np.asarray(data[2].split(',')).astype(float)
		
		estimated_causal_eqtl_file = data[5]
		estimated_causal_eqtl_effect = np.load(estimated_causal_eqtl_file)
		
		ridge_regr_causal_eqtl_file = data[6]
		ridge_regr_causal_eqtl_effect = np.load(ridge_regr_causal_eqtl_file)

		blup_causal_eqtl_effect_file = data[7]
		blup_causal_eqtl_effect = np.load(blup_causal_eqtl_effect_file)

		ridge_regr_true_var_causal_eqtl_file = data[8]
		ridge_regr_true_var_causal_eqtl_effect = np.load(ridge_regr_true_var_causal_eqtl_file)


		est_causal_eqtl_effects.append(estimated_causal_eqtl_effect)
		causal_eqtl_indices.append(gene_snp_indices)
		true_causal_eqtl_effects.append(true_causal_eqtl_effect)
		blup_causal_eqtl_effects.append(blup_causal_eqtl_effect)
		ridge_regr_causal_eqtl_effects.append(ridge_regr_causal_eqtl_effect)
		ridge_regr_true_var_causal_eqtl_effects.append(ridge_regr_true_var_causal_eqtl_effect)

	f.close()
	return causal_eqtl_indices, true_causal_eqtl_effects, est_causal_eqtl_effects, blup_causal_eqtl_effects, ridge_regr_causal_eqtl_effects, ridge_regr_true_var_causal_eqtl_effects


def sldsc_with_variants_and_distribution_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_effect_covs, causal_eqtl_indices):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	LD_sq = np.square(LD)

	#####################################################
	# Standardize genetic distribution of gene expression
	#####################################################
	std_causal_eqtl_effects = []
	std_causal_eqtl_covs = []
	for gg in range(n_genes):
		# LD specifically for gene
		R_gene = LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]
		# Calculate gene variance
		causal_eqtl_effects_mean = causal_eqtl_effects[gg]
		#gene_meanT_mean = np.dot(causal_eqtl_effects_mean.reshape(len(causal_eqtl_effects_mean),1), causal_eqtl_effects_mean.reshape(1,len(causal_eqtl_effects_mean)))
		#gene_mean_variance = np.sum(R_gene*gene_meanT_mean)
		gene_mean_variance = np.dot(np.dot(causal_eqtl_effects[gg], R_gene), causal_eqtl_effects[gg])
		gene_cov_variance = np.sum(R_gene*causal_eqtl_effect_covs[gg])
		tot_gene_variance = gene_mean_variance + gene_cov_variance

		# Standardize
		std_causal_eqtl_effects_mean = causal_eqtl_effects_mean/np.sqrt(tot_gene_variance)
		std_causal_eqtl_effect_cov = causal_eqtl_effect_covs[gg]/tot_gene_variance


		# Check that it sums to one
		#std_gene_mean_variance = np.dot(np.dot(std_causal_eqtl_effects_mean, R_gene), std_causal_eqtl_effects_mean)
		#std_gene_cov_variance = np.sum(R_gene*std_causal_eqtl_effect_cov)

		std_causal_eqtl_effects.append(std_causal_eqtl_effects_mean)
		std_causal_eqtl_covs.append(std_causal_eqtl_effect_cov)

	#####################################################
	# Create gene LD scores
	#####################################################
	gene_ld_scores = np.zeros(n_snps)
	for gg in range(n_genes):
		causal_eqtl_effects_mean = std_causal_eqtl_effects[gg]
		gene_meanT_mean = np.dot(causal_eqtl_effects_mean.reshape(len(causal_eqtl_effects_mean),1), causal_eqtl_effects_mean.reshape(1,len(causal_eqtl_effects_mean)))
		
		sub_ld = LD[:, causal_eqtl_indices[gg]]
		for kk in range(n_snps):
			gene_ld_scores[kk] = gene_ld_scores[kk] + np.dot(np.dot(sub_ld[kk,:], gene_meanT_mean + std_causal_eqtl_covs[gg]), sub_ld[kk,:])

		#gene_score = np.diag(np.dot(np.dot(LD[:, causal_eqtl_indices[gg]], (gene_meanT_mean + std_causal_eqtl_covs[gg])), LD[causal_eqtl_indices[gg], :]))

	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores)))
	var_ld_scores = np.hstack((np.sum(LD_sq,axis=0)))
	joint_ld_scores = np.transpose(np.vstack((var_ld_scores,gene_ld_scores)))

	############################
	# Run SLDSC
	############################
	model = LinearRegression(fit_intercept=False)
	ldsc_fit = model.fit(joint_ld_scores, chi_sq-1)

	nm_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas
	med_h2 = ldsc_fit.coef_[1]*n_genes/N_gwas

	return nm_h2, med_h2

def sldsc_with_variants_and_marginal_gene_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, causal_eqtl_indices):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	LD_sq = np.square(LD)
	
	#####################################################
	# Standardize genetic distribution of gene expression
	#####################################################
	unbiased_delta_t_deltas = []
	unbiased_gene_variances = []		
	for gg in range(n_genes):
		# LD specifically for gene
		R_gene = LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]
		# Get number of snps for gene
		n_gene_snps = len(causal_eqtl_indices[gg])

		unbiased_delta_t_delta = np.eye(n_gene_snps)/n_gene_snps


		tot_gene_variance = np.sum(R_gene*unbiased_delta_t_delta)
		unbiased_delta_t_deltas.append(unbiased_delta_t_delta)
		unbiased_gene_variances.append(tot_gene_variance)

	#####################################################
	# Create gene LD scores
	#####################################################
	gene_ld_scores = np.zeros(n_snps)
	for gg in range(n_genes):
		unbiased_gene_variance = unbiased_gene_variances[gg]
		unbiased_delta_t_delta = unbiased_delta_t_deltas[gg]


		sub_ld = LD[:, causal_eqtl_indices[gg]]
		for kk in range(n_snps):
			tmp = np.dot(np.dot(sub_ld[kk,:], (unbiased_delta_t_delta/unbiased_gene_variance)), sub_ld[kk,:])
			gene_ld_scores[kk] = gene_ld_scores[kk] + tmp

	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores)))
	var_ld_scores = np.hstack((np.sum(LD_sq,axis=0)))
	joint_ld_scores = np.transpose(np.vstack((var_ld_scores,gene_ld_scores)))

	############################
	# Run SLDSC
	############################
	model = LinearRegression(fit_intercept=False)
	ldsc_fit = model.fit(joint_ld_scores, chi_sq-1)

	nm_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas
	med_h2 = ldsc_fit.coef_[1]*n_genes/N_gwas


	return nm_h2, med_h2




def sldsc_with_variants_and_unbiased_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_effect_covs, causal_eqtl_indices):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	LD_sq = np.square(LD)

	#####################################################
	# Standardize genetic distribution of gene expression
	#####################################################
	unbiased_delta_t_deltas = []
	unbiased_gene_variances = []
	for gg in range(n_genes):
		# LD specifically for gene
		R_gene = LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]
		# Calculate gene variance
		causal_eqtl_effects_mean = causal_eqtl_effects[gg]
		gene_meanT_mean = np.dot(causal_eqtl_effects_mean.reshape(len(causal_eqtl_effects_mean),1), causal_eqtl_effects_mean.reshape(1,len(causal_eqtl_effects_mean)))
		gene_cov = causal_eqtl_effect_covs[gg]

		unbiased_delta_t_delta = gene_meanT_mean - gene_cov

		tot_gene_variance = np.sum(R_gene*unbiased_delta_t_delta)
		unbiased_delta_t_deltas.append(unbiased_delta_t_delta)
		unbiased_gene_variances.append(tot_gene_variance)


	#####################################################
	# Create gene LD scores
	#####################################################
	gene_ld_scores = np.zeros(n_snps)
	for gg in range(n_genes):
		unbiased_gene_variance = unbiased_gene_variances[gg]
		unbiased_delta_t_delta = unbiased_delta_t_deltas[gg]


		sub_ld = LD[:, causal_eqtl_indices[gg]]
		for kk in range(n_snps):
			tmp = np.dot(np.dot(sub_ld[kk,:], (unbiased_delta_t_delta/unbiased_gene_variance)), sub_ld[kk,:])
			gene_ld_scores[kk] = gene_ld_scores[kk] + tmp

		#gene_score = np.diag(np.dot(np.dot(LD[:, causal_eqtl_indices[gg]], (gene_meanT_mean + std_causal_eqtl_covs[gg])), LD[causal_eqtl_indices[gg], :]))

	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores)))
	var_ld_scores = np.hstack((np.sum(LD_sq,axis=0)))
	joint_ld_scores = np.transpose(np.vstack((var_ld_scores,gene_ld_scores)))

	############################
	# Run SLDSC
	############################
	model = LinearRegression(fit_intercept=False)
	ldsc_fit = model.fit(joint_ld_scores, chi_sq-1)

	nm_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas
	med_h2 = ldsc_fit.coef_[1]*n_genes/N_gwas


	return nm_h2, med_h2


def mesc_he_with_variants_and_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, LD_LD_t):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	n_snps = len(gwas_z_scores)


	z_mat = np.dot(gwas_z_scores.reshape(len(gwas_z_scores),1), gwas_z_scores.reshape(1,len(gwas_z_scores)))

	ld_scores = np.copy(LD_LD_t)
	gene_ld_scores = np.copy(z_mat)*0.0

	gene_variances = []

	for gg in range(n_genes):
		# Calculate gene variance
		gene_variance = np.dot(np.dot(causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), causal_eqtl_effects[gg])

		gene_variances.append(gene_variance)
		# Standardize gene effects
		#std_causal_eqtl_effect = causal_eqtl_effects[gg]/np.sqrt(gene_variance)
		std_causal_eqtl_effect = np.copy(causal_eqtl_effects[gg])

		'''
		delta_delta_t = np.dot(std_causal_eqtl_effect.reshape(len(std_causal_eqtl_effect),1), std_causal_eqtl_effect.reshape(1,len(std_causal_eqtl_effect)))
		sub_ld = LD[:, causal_eqtl_indices[gg]]
		gene_ld_scores = gene_ld_scores + np.dot(np.dot(sub_ld, delta_delta_t), np.transpose(sub_ld))
		'''
		tmp =np.dot(LD[:, causal_eqtl_indices[gg]], std_causal_eqtl_effect)
		gene_ld_scores = gene_ld_scores + np.dot(tmp.reshape(len(tmp),1),tmp.reshape(1,len(tmp)))

	model = LinearRegression(fit_intercept=False, normalize=False)
	z_mat = z_mat - LD

	joint_ld_scores = np.transpose(np.vstack((np.matrix.flatten(ld_scores),np.matrix.flatten(gene_ld_scores))))
	mpldsc_fit = model.fit(joint_ld_scores, np.matrix.flatten(z_mat))


	mpldsc_var_h2 = mpldsc_fit.coef_[0]*n_snps/N_gwas
	mpldsc_gene_h2 = mpldsc_fit.coef_[1]*n_genes*np.mean(gene_variances)/N_gwas

	return mpldsc_var_h2, mpldsc_gene_h2


def sldsc_he_with_variants_and_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, LD_LD_t):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	n_snps = len(gwas_z_scores)


	z_mat = np.dot(gwas_z_scores.reshape(len(gwas_z_scores),1), gwas_z_scores.reshape(1,len(gwas_z_scores)))

	ld_scores = np.copy(LD_LD_t)
	gene_ld_scores = np.copy(z_mat)*0.0

	for gg in range(n_genes):
		# Calculate gene variance
		gene_variance = np.dot(np.dot(causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), causal_eqtl_effects[gg])

		# Standardize gene effects
		std_causal_eqtl_effect = causal_eqtl_effects[gg]/np.sqrt(gene_variance)
		#std_causal_eqtl_effect = np.copy(causal_eqtl_effects[gg])

		'''
		delta_delta_t = np.dot(std_causal_eqtl_effect.reshape(len(std_causal_eqtl_effect),1), std_causal_eqtl_effect.reshape(1,len(std_causal_eqtl_effect)))
		sub_ld = LD[:, causal_eqtl_indices[gg]]
		gene_ld_scores = gene_ld_scores + np.dot(np.dot(sub_ld, delta_delta_t), np.transpose(sub_ld))
		'''
		tmp =np.dot(LD[:, causal_eqtl_indices[gg]], std_causal_eqtl_effect)
		gene_ld_scores = gene_ld_scores + np.dot(tmp.reshape(len(tmp),1),tmp.reshape(1,len(tmp)))

	model = LinearRegression(fit_intercept=False, normalize=False)
	z_mat = z_mat - LD

	joint_ld_scores = np.transpose(np.vstack((np.matrix.flatten(ld_scores),np.matrix.flatten(gene_ld_scores))))
	mpldsc_fit = model.fit(joint_ld_scores, np.matrix.flatten(z_mat))

	mpldsc_var_h2 = mpldsc_fit.coef_[0]*n_snps/N_gwas
	mpldsc_gene_h2 = mpldsc_fit.coef_[1]*n_genes/N_gwas


	return mpldsc_var_h2, mpldsc_gene_h2


def mesc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq_v3(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, est_gene_ld_scores, avg_prediction_noise, fit_intercept_bool=False):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	LD_sq = np.square(LD)

	#####################################################
	# Compute gene variances
	#####################################################
	gene_variances = []
	for gg in range(n_genes):
		# Calculate gene variance
		gene_variance = np.dot(np.dot(causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), causal_eqtl_effects[gg])

		gene_variances.append(gene_variance)
	gene_variances = np.asarray(gene_variances)


	##########################
	# Create gene-variant LD mat
	##########################
	gene_variant_ld = np.zeros((n_genes, n_snps))
	for gg in range(n_genes):
		gene_variant_ld[gg,:] = np.dot(LD[:, causal_eqtl_indices[gg]], causal_eqtl_effects[gg])
	# Square LD
	gene_variant_ld_sq = np.square(gene_variant_ld)

	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores)))
	var_ld_scores = np.hstack((np.sum(LD_sq,axis=0)))
	true_gene_ld_scores = np.hstack((np.sum(gene_variant_ld_sq,axis=0)))
	
	############################
	# Run SLDSC
	############################
	if fit_intercept_bool == False:
		n_samples = len(est_gene_ld_scores)
		covs = np.transpose(np.vstack((var_ld_scores,est_gene_ld_scores)))
		S = np.zeros((2,2))
		true_noise = np.var(est_gene_ld_scores,ddof=1) - np.var(true_gene_ld_scores,ddof=1)
		S[1,1] = avg_prediction_noise
		adj_coefs = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(covs), covs) - S*n_samples), np.transpose(covs)), chi_sq-1)
	
		nm_h2_adj = adj_coefs[0]*n_snps/N_gwas
		med_h2_adj = adj_coefs[1]*n_genes*np.mean(gene_variances)/N_gwas
	else: 		
		n_samples = len(est_gene_ld_scores)
		covs = np.transpose(np.vstack((np.ones(n_samples),var_ld_scores,est_gene_ld_scores)))
		S = np.zeros((3,3))
		true_noise = np.var(est_gene_ld_scores,ddof=1) - np.var(true_gene_ld_scores,ddof=1)
		S[2,2] = avg_prediction_noise
		adj_coefs = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(covs), covs) - S*n_samples), np.transpose(covs)), chi_sq-1)
	
		nm_h2_adj = adj_coefs[1]*n_snps/N_gwas
		med_h2_adj = adj_coefs[2]*n_genes*np.mean(gene_variances)/N_gwas


	return med_h2_adj, true_noise

def mesc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq_v4(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, est_gene_ld_scores, avg_prediction_noise, fit_intercept_bool=False):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	LD_sq = np.square(LD)

	#####################################################
	# Compute gene variances
	#####################################################
	gene_variances = []
	for gg in range(n_genes):
		# Calculate gene variance
		gene_variance = np.dot(np.dot(causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), causal_eqtl_effects[gg])

		gene_variances.append(gene_variance)
	gene_variances = np.asarray(gene_variances)


	##########################
	# Create gene-variant LD mat
	##########################
	gene_variant_ld = np.zeros((n_genes, n_snps))
	for gg in range(n_genes):
		gene_variant_ld[gg,:] = np.dot(LD[:, causal_eqtl_indices[gg]], causal_eqtl_effects[gg])
	# Square LD
	gene_variant_ld_sq = np.square(gene_variant_ld)

	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores)))
	var_ld_scores = np.hstack((np.sum(LD_sq,axis=0)))
	true_gene_ld_scores = np.hstack((np.sum(gene_variant_ld_sq,axis=0)))
	
	############################
	# Run SLDSC
	############################
	joint_ld_scores = np.transpose(np.vstack((var_ld_scores,est_gene_ld_scores)))

	model = LinearRegression(fit_intercept=fit_intercept_bool)
	ldsc_fit = model.fit(joint_ld_scores, chi_sq-1)

	nm_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas
	med_h2 = ldsc_fit.coef_[1]*n_genes*np.mean(gene_variances)/N_gwas


	######################
	# Disattenuate
	######################
	true_joint_ld_scores = np.transpose(np.vstack((var_ld_scores,true_gene_ld_scores)))

	VV = np.cov(np.transpose(joint_ld_scores))
	VV_prime = np.cov(np.transpose(joint_ld_scores))
	VV[1,1] = VV[1,1] - avg_prediction_noise


	lambda_term = np.dot(np.linalg.inv(VV_prime), VV)

	adjusted_coef = np.linalg.solve(lambda_term, ldsc_fit.coef_)
	
	nm_h2_adj = adjusted_coef[0]*n_snps/N_gwas
	med_h2_adj = adjusted_coef[1]*n_genes*np.mean(gene_variances)/N_gwas

	return med_h2_adj



def mesc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq_v2(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, est_gene_ld_scores,fit_intercept_bool=False):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	LD_sq = np.square(LD)

	#####################################################
	# Compute gene variances
	#####################################################
	gene_variances = []
	for gg in range(n_genes):
		# Calculate gene variance
		gene_variance = np.dot(np.dot(causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), causal_eqtl_effects[gg])

		gene_variances.append(gene_variance)
	gene_variances = np.asarray(gene_variances)


	##########################
	# Create gene-variant LD mat
	##########################
	gene_variant_ld = np.zeros((n_genes, n_snps))
	for gg in range(n_genes):
		gene_variant_ld[gg,:] = np.dot(LD[:, causal_eqtl_indices[gg]], causal_eqtl_effects[gg])
	# Square LD
	gene_variant_ld_sq = np.square(gene_variant_ld)

	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores)))
	var_ld_scores = np.hstack((np.sum(LD_sq,axis=0)))
	true_gene_ld_scores = np.hstack((np.sum(gene_variant_ld_sq,axis=0)))
	
	############################
	# Run SLDSC
	############################
	joint_ld_scores = np.transpose(np.vstack((var_ld_scores,est_gene_ld_scores)))

	model = LinearRegression(fit_intercept=fit_intercept_bool)
	ldsc_fit = model.fit(joint_ld_scores, chi_sq-1)

	nm_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas
	med_h2 = ldsc_fit.coef_[1]*n_genes*np.mean(gene_variances)/N_gwas


	######################
	# Disattenuate
	######################
	true_joint_ld_scores = np.transpose(np.vstack((var_ld_scores,true_gene_ld_scores)))

	VV = np.cov(np.transpose(true_joint_ld_scores))
	VV_prime = np.cov(np.transpose(joint_ld_scores))

	lambda_term = np.dot(np.linalg.inv(VV_prime), VV)

	adjusted_coef = np.linalg.solve(lambda_term, ldsc_fit.coef_)
	
	nm_h2_adj = adjusted_coef[0]*n_snps/N_gwas
	med_h2_adj = adjusted_coef[1]*n_genes*np.mean(gene_variances)/N_gwas

	return med_h2, med_h2_adj




def mesc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, est_gene_ld_scores,fit_intercept_bool=False):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	LD_sq = np.square(LD)

	#####################################################
	# Compute gene variances
	#####################################################
	gene_variances = []
	for gg in range(n_genes):
		# Calculate gene variance
		gene_variance = np.dot(np.dot(causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), causal_eqtl_effects[gg])

		gene_variances.append(gene_variance)
	gene_variances = np.asarray(gene_variances)


	##########################
	# Create gene-variant LD mat
	##########################
	gene_variant_ld = np.zeros((n_genes, n_snps))
	for gg in range(n_genes):
		gene_variant_ld[gg,:] = np.dot(LD[:, causal_eqtl_indices[gg]], causal_eqtl_effects[gg])
	# Square LD
	gene_variant_ld_sq = np.square(gene_variant_ld)

	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores)))
	var_ld_scores = np.hstack((np.sum(LD_sq,axis=0)))
	true_gene_ld_scores = np.hstack((np.sum(gene_variant_ld_sq,axis=0)))

	############################
	# Run SLDSC
	############################
	joint_ld_scores = np.transpose(np.vstack((var_ld_scores,est_gene_ld_scores)))

	model = LinearRegression(fit_intercept=fit_intercept_bool)
	ldsc_fit = model.fit(joint_ld_scores, chi_sq-1)

	nm_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas
	med_h2 = ldsc_fit.coef_[1]*n_genes*np.mean(gene_variances)/N_gwas


	######################
	# Disattenuate
	######################
	shrink_factor = np.var(true_gene_ld_scores)/np.var(est_gene_ld_scores)
	r_sq = np.square(np.corrcoef(var_ld_scores, est_gene_ld_scores)[0,1])
	new_shrink_factor = (shrink_factor - r_sq)/(1.0 - r_sq)
	med_h2_adj = med_h2/new_shrink_factor

	return med_h2, med_h2_adj




def sldsc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, noise_ratios):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	#LD_sq = np.square(LD)

	#####################################################
	# Standardize genetically predicted gene expression
	#####################################################
	std_causal_eqtl_effects = []
	for gg in range(n_genes):
		# Calculate gene variance
		gene_variance = np.dot(np.dot(causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), causal_eqtl_effects[gg])

		# Standardize gene effects
		std_causal_eqtl_effect = causal_eqtl_effects[gg]/np.sqrt(gene_variance)

		# Append to array
		std_causal_eqtl_effects.append(std_causal_eqtl_effect)


	##########################
	# Create gene-variant LD mat
	##########################
	gene_variant_ld = np.zeros((n_genes, n_snps))
	for gg in range(n_genes):
		gene_variant_ld[gg,:] = np.dot(LD[:, causal_eqtl_indices[gg]], std_causal_eqtl_effects[gg])
	# Square LD
	gene_variant_ld_sq = np.square(gene_variant_ld)

	##########################
	# Create variant ld scores
	##########################
	# Get variant ld scores 
	var_ld_scores = []
	for snp_iter in range(n_snps):
		var_ld_scores.append(np.sum(np.square(LD[snp_iter,:])))


	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores)))
	var_ld_scores = np.asarray(var_ld_scores)
	gene_ld_scores = np.hstack((np.sum(gene_variant_ld_sq,axis=0)))
	joint_ld_scores = np.transpose(np.vstack((var_ld_scores,gene_ld_scores)))
	
	############################
	# Run SLDSC
	############################
	model = LinearRegression(fit_intercept=False)
	ldsc_fit = model.fit(joint_ld_scores, chi_sq-1)

	nm_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas
	med_h2 = ldsc_fit.coef_[1]*n_genes/N_gwas

	return nm_h2, med_h2





def sldsc_with_variants_and_noisy_genetic_gene_expression(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, noise_ratios):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	LD_sq = np.square(LD)

	#####################################################
	# Standardize genetically predicted gene expression
	#####################################################
	std_causal_eqtl_effects = []
	for gg in range(n_genes):
		# Calculate gene variance
		gene_variance = np.dot(np.dot(causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), causal_eqtl_effects[gg])

		# Standardize gene effects
		std_causal_eqtl_effect = causal_eqtl_effects[gg]/np.sqrt(gene_variance)

		# Append to array
		std_causal_eqtl_effects.append(std_causal_eqtl_effect)

	##########################
	# Compute TWAS Z scores
	##########################
	twas_z_scores = []
	for gg in range(n_genes):
		tmp_numerator = np.dot(std_causal_eqtl_effects[gg], gwas_z_scores[causal_eqtl_indices[gg]])
		# NOTE: Denomenator should always be equal to one cause we standardized. but for completeness..
		tmp_denomenator = np.sqrt(np.dot(np.dot(std_causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), std_causal_eqtl_effects[gg])) 
		twas_z = tmp_numerator/tmp_denomenator
		# Add to global array
		twas_z_scores.append(twas_z)
	twas_z_scores = np.asarray(twas_z_scores)


	##########################
	# Create gene-gene LD mat
	##########################
	gene_gene_ld = np.zeros((n_genes, n_genes))
	for gg in range(n_genes):
		for kk in range(n_genes):
			corry = np.dot(np.dot(std_causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg],:][:,causal_eqtl_indices[kk]]), std_causal_eqtl_effects[kk])
			gene_gene_ld[gg,kk] = corry
			gene_gene_ld[kk,gg] = corry


	##########################
	# Correct gene-gene LD mat for gene noise
	##########################
	noise_ratios = np.asarray(noise_ratios)
	# Set noise ratios larger than 1 to 1
	noise_ratios[noise_ratios > 1.0] = 1.0
	# Correction
	gene_gene_ld = gene_gene_ld - np.diag(noise_ratios)
	# Squared LD
	gene_gene_ld_sq = np.square(gene_gene_ld)

	##########################
	# Create gene-variant LD mat
	##########################
	gene_variant_ld = np.zeros((n_genes, n_snps))
	for gg in range(n_genes):
		gene_variant_ld[gg,:] = np.dot(LD[:, causal_eqtl_indices[gg]], std_causal_eqtl_effects[gg])
	# Square LD
	gene_variant_ld_sq = np.square(gene_variant_ld)

	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores), np.square(twas_z_scores)))
	var_ld_scores = np.hstack((np.sum(LD_sq,axis=0), np.sum(gene_variant_ld_sq,axis=1)))
	gene_ld_scores = np.hstack((np.sum(gene_variant_ld_sq,axis=0), np.sum(gene_gene_ld_sq,axis=0)))
	joint_ld_scores = np.transpose(np.vstack((var_ld_scores,gene_ld_scores)))
	
	############################
	# Run SLDSC
	############################
	model = LinearRegression(fit_intercept=False)
	ldsc_fit = model.fit(joint_ld_scores, chi_sq-1)

	nm_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas
	med_h2 = ldsc_fit.coef_[1]*n_genes/N_gwas

	return nm_h2, med_h2

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


def extract_gene_ld_scores(simulated_gene_models_dir, eqtl_ss, simulation_name_string):
	gene_model_summary_file = simulated_gene_models_dir + simulation_name_string + 'model_summaries_' + str(eqtl_ss) + '.txt'
	gene_ld_scores = []
	gene_ld_score_noises = []
	f = open(gene_model_summary_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		counter = counter + 1
		gene_ld_score_file = data[5]
		gene_ld_score = np.load(gene_ld_score_file)
		gene_ld_scores.append(gene_ld_score)

		gene_beta_file = data[3]
		gene_beta = np.load(gene_beta_file)
		gene_beta_var_file = data[4]
		gene_beta_var = np.load(gene_beta_var_file)
		gene_ld_score_noise = calculate_gene_ld_score_noise_according_to_gaussian_distr(gene_beta, gene_beta_var, float(eqtl_ss))
		gene_ld_score_noises.append(gene_ld_score_noise)
	f.close()
	gene_ld_scores = np.asarray(gene_ld_scores)
	gene_ld_score_noises = np.asarray(gene_ld_score_noises)

	return np.sum(gene_ld_scores,axis=0), np.mean(np.sum(gene_ld_score_noises,axis=0))




#########################
# Command line args
#########################
simulation_name_string = sys.argv[1]
simulated_gwas_data_dir = sys.argv[2]
mediated_h2_results_dir = sys.argv[3]
processed_genotype_data_dir = sys.argv[4]
simulated_eqtl_data_dir = sys.argv[5]
selection = sys.argv[6]
eqtl_ss = sys.argv[7]
simulated_gene_models_dir = sys.argv[8]



# File summarizing inferred gene models
gene_summary_file = simulated_eqtl_data_dir + simulation_name_string + 'causal_eqtl_effect_summary.txt'


# Extract previously estimated causal eqtl effects
causal_eqtl_indices, true_causal_eqtl_effects = extract_causal_eqtl_effects(gene_summary_file)


# File containing gwas summary statistics
gwas_sum_stats_file = simulated_gwas_data_dir + simulation_name_string + selection + '_selection_' + 'simulated_gwas_summary_stats.txt'
gwas_z_scores = load_in_gwas_z_scores(gwas_sum_stats_file)
gwas_se = load_in_gwas_se(gwas_sum_stats_file)
gwas_beta = gwas_z_scores*gwas_se

# Simulated total genetic vars
simulated_genetic_var_file = simulated_gwas_data_dir + simulation_name_string+ selection + '_selection_' + 'simulated_genetic_var.npy'
sim_genetic_var = np.load(simulated_genetic_var_file) + 0.0
simulated_nm_genetic_var_file = simulated_gwas_data_dir + simulation_name_string + selection + '_selection_'+ 'simulated_nm_genetic_var.npy'
sim_nm_genetic_var = np.load(simulated_nm_genetic_var_file) + 0.0
simulated_mediated_genetic_var_file = simulated_gwas_data_dir + simulation_name_string+ selection + '_selection_' + 'simulated_mediated_genetic_var.npy'
sim_med_genetic_var = np.load(simulated_mediated_genetic_var_file) + 0.0

# Print what simulated total genetic vars are
print(sim_genetic_var)
print(sim_nm_genetic_var)
print(sim_med_genetic_var)



# Load in Genotype file
# LD file
ld_file = processed_genotype_data_dir + 'gwas_genotype_LD_1.npy'
LD = np.load(ld_file)

ld_ld_t_file = processed_genotype_data_dir + 'gwas_genotype_LD_LD_t_1.npy'
#LD_LD_t = np.load(ld_ld_t_file)

# Genotype file
geno_file = processed_genotype_data_dir + 'gwas_genotype_1.npy'
#genotype_mat = np.load(geno_file)


# Simulated Trait
trait_file = simulated_gwas_data_dir + simulation_name_string + selection + '_selection_' + 'simulated_trait.txt'
trait_vec = load_in_trait_file(trait_file)

# Gwas sample size
N_gwas = len(trait_vec)

# Number of genes
n_genes = len(causal_eqtl_indices)


# Extract gene LD Scores
est_gene_ld_scores, avg_prediction_noise = extract_gene_ld_scores(simulated_gene_models_dir, eqtl_ss, simulation_name_string)

# Print results to output
output_file = mediated_h2_results_dir + simulation_name_string + selection + '_selection_with_est_noise_eqtl_ss_' + str(eqtl_ss)+ 'med_h2_ldsc_style_summary.txt'
t = open(output_file,'w')
# Header
t.write('sim_h2\tsim_nm_h2\tsim_med_h2')
t.write('\tsim_noise\tpred_noise\tmed_sldsc_est_noisy_gene\tmed_sldsc_est_w_intercept_noisy_gene_adj_v3\tmed_sldsc_est_w_intercept_noisy_gene_adj_v4')
t.write('\n')

t.write(str(sim_genetic_var) + '\t' + str(sim_nm_genetic_var) + '\t' + str(sim_med_genetic_var))


# Run SLDSC with variants and true genetic gene expression (accounting for zero noise) and without modeling twas chi-squared stats
med_sldsc_est_noisy_gene, med_sldsc_est_noisy_gene_adj = mesc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, true_causal_eqtl_effects, causal_eqtl_indices, est_gene_ld_scores)
#med_sldsc_est_w_intercept_noisy_gene, med_sldsc_est_w_intercept_noisy_gene_adj = mesc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, true_causal_eqtl_effects, causal_eqtl_indices, est_gene_ld_scores, fit_intercept_bool=True)

#med_sldsc_est_noisy_gene_v2, med_sldsc_est_noisy_gene_adj_v2 = mesc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq_v2(gwas_z_scores, N_gwas, LD, true_causal_eqtl_effects, causal_eqtl_indices, est_gene_ld_scores, fit_intercept_bool=False)
#med_sldsc_est_w_intercept_noisy_gene_v2, med_sldsc_est_w_intercept_noisy_gene_adj_v2 = mesc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq_v2(gwas_z_scores, N_gwas, LD, true_causal_eqtl_effects, causal_eqtl_indices, est_gene_ld_scores, fit_intercept_bool=True)

med_sldsc_est_w_intercept_noisy_gene_adj_v3, true_noise = mesc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq_v3(gwas_z_scores, N_gwas, LD, true_causal_eqtl_effects, causal_eqtl_indices, est_gene_ld_scores, avg_prediction_noise, fit_intercept_bool=True)
med_sldsc_est_w_intercept_noisy_gene_adj_v4 = mesc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq_v4(gwas_z_scores, N_gwas, LD, true_causal_eqtl_effects, causal_eqtl_indices, est_gene_ld_scores, avg_prediction_noise, fit_intercept_bool=True)

t.write('\t' + str(true_noise) + '\t' + str(avg_prediction_noise))
t.write('\t' + str(med_sldsc_est_noisy_gene))
t.write('\t' + str(med_sldsc_est_w_intercept_noisy_gene_adj_v3) + '\t' + str(med_sldsc_est_w_intercept_noisy_gene_adj_v4))
t.write('\n')


t.close()


print(output_file)












import numpy as np
import os
import sys
import pdb
import time
import statsmodels.api as sm


def linear_regression(XX, YY, intercept=False):
	if intercept:
		print('not yet implemented. add column of 1s to X')
		pdb.set_trace()

	return np.dot(np.dot(np.linalg.pinv(np.dot(np.transpose(XX), XX)), np.transpose(XX)), YY)

def fixed_variance_update_eqtl_resid_vars(eqtl_ss, gene_info, genes, eqtl_category_to_index):
	n_cat = len(eqtl_category_to_index)
	sums = np.zeros(n_cat)
	counts = np.zeros(n_cat)
	for gene in genes:
		info = gene_info[gene]
		for ii, category_name in enumerate(info['eqtl_category_names']):
			cat_index = eqtl_category_to_index[category_name]
			sums[cat_index] = sums[cat_index] + np.sum(info['squared_sumstats'][:,ii])
			counts[cat_index] = counts[cat_index] + len(info['squared_sumstats'][:,ii])
	tmp_E_beta_sq = sums/counts
	resid_vars = (4.0*(tmp_E_beta_sq)/eqtl_ss) + 2.0*np.square(1/eqtl_ss)
	return resid_vars

def initialize_variant_to_gene_effect_size_variances(genes, gene_info, n_reg_snps, n_eqtl_classes, eqtl_category_to_index):
	# Initialize output dictionary
	delta_sq = {}
	# Keep track of estimated heritability
	gene_cis_h2 = {}
	# Initialize matrix of gene LD scores
	gene_ld_scores = np.zeros((n_reg_snps, n_eqtl_classes))
	# Initialize vector to keep track of residual varainces
	resid_var_sum_sq = np.zeros(n_eqtl_classes)
	resid_var_count = np.zeros(n_eqtl_classes)

	# loop through genes
	for gene in genes:
		info = gene_info[gene]
		# Load in squared ld 
		squared_ld = np.load(info['squared_ld_file'])

		# Initialize gene delta_sq mat
		gene_delta_sq = np.zeros((len(info['n_low_dimensional_snps']), len(info['eqtl_category_names'])))

		# Initialize gene cis_h2 vec
		gene_cis_h2_vec = np.zeros(len(info['eqtl_category_names']))

		# Load in eqtl sumstats
		eqtl_sq_sumstats = info['squared_sumstats']

		# Load in indices of regression snp
		gene_regression_snp_indices = info['regression_snp_indices']

		# Load in gene eqtl categories
		gene_eqtl_categories = info['eqtl_category_names']

		# Get number of eqtl categories for this gene
		n_eqtl_categories_per_gene = eqtl_sq_sumstats.shape[1]

		# Quick error check
		if n_eqtl_categories_per_gene != gene_delta_sq.shape[1]:
			print('assumption error')
			pdb.set_trace()

		# Loop through eQTL categories
		# Initialize seperately for each category
		for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):
			# Fit model
			#model = sm.OLS(eqtl_sq_sumstats[:, eqtl_category_iter], squared_ld).fit()
			#weights2 = model.params
			weights = linear_regression(squared_ld, eqtl_sq_sumstats[:, eqtl_category_iter])

			# update gene_delta_sq mat
			gene_delta_sq[:, eqtl_category_iter] = weights

			# keep track of estimated h2s
			gene_cis_h2_vec[eqtl_category_iter] = np.sum(weights*info['n_low_dimensional_snps'])

			# Update gene LD scores
			global_category_index = eqtl_category_to_index[eqtl_category_name]
			gene_ld_scores[gene_regression_snp_indices, global_category_index] = gene_ld_scores[gene_regression_snp_indices, global_category_index] + np.dot(squared_ld, weights)

			# Update resid var counts
			resid = eqtl_sq_sumstats[:, eqtl_category_iter] - np.dot(squared_ld, weights)
			resid_var_sum_sq[global_category_index] = resid_var_sum_sq[global_category_index] + np.sum(np.square(resid))
			resid_var_count[global_category_index] = resid_var_count[global_category_index] + len(resid)

		# Add to dictionary
		delta_sq[gene] = gene_delta_sq
		gene_cis_h2[gene] = gene_cis_h2_vec

	# Update resid variances
	eqtl_resid_vars = resid_var_sum_sq/resid_var_count


	return delta_sq, gene_cis_h2, gene_ld_scores, eqtl_resid_vars




class med_h2(object):
	def __init__(self, version='default',max_iter=600, convergence_thresh=1e-10, fixed_variance=False):
		self.max_iter = max_iter  # Total iterations
		self.convergence_thresh = convergence_thresh
		self.version = version
		self.fixed_variance = fixed_variance
	def fit(self, genes, gene_info, gwas_variant_ld_scores, gwas_E_beta_sq, eqtl_dataset_names, n_ref_snps):
		""" Fit the model.
			Args:
			tgfm_data_obj
		"""
		if self.version == 'two_step':
			print('###############################')
			print(' Two-step LDSC   ')
			print('###############################')
		else:
			print('###############################')
			print(' Joint LDSC   ')
			print('###############################')			


		# Initialize variables
		self.initialize_variables(genes, gene_info, gwas_variant_ld_scores, gwas_E_beta_sq, eqtl_dataset_names, n_ref_snps)
		print('NM h2: ' + str(self.nm_h2))
		print('Total med h2: ' + str(np.sum(self.med_h2)))
		print('med h2:')
		print(self.med_h2)
		print('eqtl h2:')
		print(self.avg_eqtl_h2)
		print('Resid vars')
		print(self.gwas_resid_var)
		print(self.eqtl_resid_vars)
		if self.version == 'two_step':
			return

		# Used to assess convergence
		cur_med_h2 = np.copy(self.med_h2)
		converged = False

		# Begin iterative optimization
		for itera in range(self.max_iter):
			t1 = time.time()
			gene_ldscores = self.update_variant_to_gene_effect_size_variances(genes, gene_info)

			# Update non-mediated per snp h2 (gamma_sq) and gene_trait varaince (alpha_sq)
			mesc_style_params = linear_regression(np.hstack((gwas_variant_ld_scores, gene_ldscores)), gwas_E_beta_sq)
			self.gamma_sq = mesc_style_params[:self.n_non_med_categories]
			self.alpha_sq = mesc_style_params[self.n_non_med_categories:]


			# Compute current estimates of non-mediated and mediated heritabilities
			self.update_nm_h2()
			self.update_med_h2(genes, gene_info)
			self.update_average_eqtl_h2(genes)
			# Compute correct estimate of gwas resid var
			if self.fixed_variance == False:
				self.update_gwas_resid_var()
			t2 = time.time()

			print('################')
			print('Iteration: ' + str(itera))
			print('NM h2: ' + str(self.nm_h2))
			print('Total med h2: ' + str(np.sum(self.med_h2)))
			print('med h2:')
			print(self.med_h2)
			print('eqtl h2:')
			print(self.avg_eqtl_h2)
			print('Resid vars')
			print(self.gwas_resid_var)
			print(self.eqtl_resid_vars)
			print('Time')
			print(t2-t1)
			# Assess convergence
			max_abs_diff = np.max(np.abs((cur_med_h2 - self.med_h2)))
			cur_med_h2 = np.copy(self.med_h2)
			print('Diff:')
			print(max_abs_diff)
			print(' ')
			if max_abs_diff <= self.convergence_thresh:
				converged = True
				break

		if converged == False:
			print('Did not converge after ' + str(self.max_iter) + ' iterations')



		return


	def update_variant_to_gene_effect_size_variances(self, genes, gene_info):
		# Initialize matrix keeping track of gene ld scores
		gene_ld_scores = np.zeros((self.n_reg_snps, self.n_eqtl_classes))
		# Initialize vector to keep track of residual varainces
		resid_var_sum_sq = np.zeros(self.n_eqtl_classes)
		resid_var_count = np.zeros(self.n_eqtl_classes)


		# loop through genes
		for gene in np.random.permutation(genes):
			info = gene_info[gene]
			# Load in squared ld 
			squared_ld = np.load(info['squared_ld_file'])

			# Initialize gene delta_sq mat
			gene_delta_sq = np.copy(self.delta_sq[gene])

			# get number of low dimensional snps for gene
			n_low_dim_snps = squared_ld.shape[1]

			# Load in eqtl sumstats
			eqtl_sq_sumstats = info['squared_sumstats']

			# Load in indices of regression snp
			gene_regression_snp_indices = info['regression_snp_indices']

			# Load in gene eqtl categories
			gene_eqtl_categories = info['eqtl_category_names']

			# Get number of eqtl categories for this gene
			n_eqtl_categories_per_gene = eqtl_sq_sumstats.shape[1]

			if n_eqtl_categories_per_gene == 0:
				continue

			# Get alpha vec corresponding to this gene
			tmp_gene_alpha_sq = np.zeros(n_eqtl_categories_per_gene)
			for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):
				tmp_gene_alpha_sq[eqtl_category_iter] = self.alpha_sq[self.eqtl_category_to_index[eqtl_category_name]]

			# Un-regress out gene effects
			pred_gene_effects = np.dot(np.dot(squared_ld, self.delta_sq[gene]), tmp_gene_alpha_sq)
			self.gwas_resid_E_beta_sq[gene_regression_snp_indices] = self.gwas_resid_E_beta_sq[gene_regression_snp_indices] + pred_gene_effects


			precomp_eqtl_xt_x = np.dot(np.transpose(squared_ld), squared_ld)
			xt_y_terms = np.zeros(len(gene_eqtl_categories)*n_low_dim_snps)
			xt_x_terms = np.zeros((len(gene_eqtl_categories)*n_low_dim_snps, len(gene_eqtl_categories)*n_low_dim_snps))
			gwas_X = []
			for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):
				# Get global eqtl category index
				global_eqtl_category_index = self.eqtl_category_to_index[eqtl_category_name]
				# Get indices for gene
				gene_indices = np.arange(eqtl_category_iter*n_low_dim_snps, (eqtl_category_iter+1)*n_low_dim_snps)

				# Compute xt_y for this gene
				tmp_xt_y = np.dot(np.transpose(squared_ld), eqtl_sq_sumstats[:, eqtl_category_iter])/(self.eqtl_resid_vars[global_eqtl_category_index])
				tmp_xt_y = tmp_xt_y + np.dot(np.transpose(squared_ld)*self.alpha_sq[global_eqtl_category_index], self.gwas_resid_E_beta_sq[gene_regression_snp_indices])/(self.gwas_resid_var)
				xt_y_terms[gene_indices] = tmp_xt_y

				# add eqtl xt_x
				xt_x_terms[np.ix_(gene_indices, gene_indices)] = precomp_eqtl_xt_x/(self.eqtl_resid_vars[global_eqtl_category_index])

				gwas_X.append(squared_ld*self.alpha_sq[global_eqtl_category_index])
			gwas_X = np.hstack(gwas_X)
			xt_x_terms = xt_x_terms + (np.dot(np.transpose(gwas_X), gwas_X)/self.gwas_resid_var)

			# Run regression
			S = np.linalg.pinv(xt_x_terms)
			mean = np.dot(S, xt_y_terms)

			tmp_gene_cis_h2_vec = np.zeros(len(info['eqtl_category_names']))
			# Update delta_sqs
			for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):
				# Get indices for gene
				gene_indices = np.arange(eqtl_category_iter*n_low_dim_snps, (eqtl_category_iter+1)*n_low_dim_snps)
				gene_delta_sq[:, eqtl_category_iter] = mean[gene_indices]
				tmp_gene_cis_h2_vec[eqtl_category_iter] = np.sum(mean[gene_indices]*info['n_low_dimensional_snps'])

				# Get global eqtl category index
				global_eqtl_category_index = self.eqtl_category_to_index[eqtl_category_name]
				gene_ld_scores[gene_regression_snp_indices, global_eqtl_category_index] = gene_ld_scores[gene_regression_snp_indices, global_eqtl_category_index] + np.dot(squared_ld, mean[gene_indices])


				# Update resid var counts
				resid = eqtl_sq_sumstats[:, eqtl_category_iter] - np.dot(squared_ld, mean[gene_indices])
				resid_var_sum_sq[global_eqtl_category_index] = resid_var_sum_sq[global_eqtl_category_index] + np.sum(np.square(resid))
				resid_var_count[global_eqtl_category_index] = resid_var_count[global_eqtl_category_index] + len(resid)


			self.delta_sq[gene] = gene_delta_sq
			self.gene_cis_h2[gene] = tmp_gene_cis_h2_vec
			# Re-regress out updated gene effects
			pred_gene_effects = np.dot(np.dot(squared_ld, self.delta_sq[gene]), tmp_gene_alpha_sq)
			self.gwas_resid_E_beta_sq[gene_regression_snp_indices] = self.gwas_resid_E_beta_sq[gene_regression_snp_indices] - pred_gene_effects

		if self.fixed_variance == False:
			self.eqtl_resid_vars = resid_var_sum_sq/resid_var_count

		return gene_ld_scores


	def update_nm_h2(self):
		self.nm_h2 = np.sum((self.gamma_sq)*(self.n_ref_snps))
		return

	def update_med_h2(self, genes, gene_info):
		self.med_h2 = np.zeros(self.n_eqtl_classes)
		for gene in genes:
			tmp_gene_cis_h2s = self.gene_cis_h2[gene]
			gene_eqtl_categories = gene_info[gene]['eqtl_category_names']
			for category_iter, category_name in enumerate(gene_eqtl_categories):
				global_category_index = self.eqtl_category_to_index[category_name]
				self.med_h2[global_category_index] = self.med_h2[global_category_index] + ((tmp_gene_cis_h2s[category_iter])*(self.alpha_sq[global_category_index]))
		return

	def update_average_eqtl_h2(self, genes):
		tmp_arr = []
		for gene in genes:
			tmp_arr.append(self.gene_cis_h2[gene])
		tmp_arr = np.hstack(tmp_arr)
		self.avg_eqtl_h2 = np.mean(tmp_arr)
		return

	def update_gwas_resid_var(self):
		self.gwas_resid_var = np.sum(np.square(self.gwas_resid_E_beta_sq))/len(self.gwas_resid_E_beta_sq)
		return

	def initialize_variables(self, genes, gene_info, gwas_variant_ld_scores, gwas_E_beta_sq, eqtl_dataset_names, n_ref_snps):
		# Number of genes
		self.n_genes = len(genes)
		self.n_reg_snps = gwas_variant_ld_scores.shape[0]
		self.n_ref_snps = n_ref_snps
		self.n_eqtl_classes = len(eqtl_dataset_names)
		self.gwas_E_beta_sq = gwas_E_beta_sq
		self.n_non_med_categories = gwas_variant_ld_scores.shape[1]

		# Create dictionary mapping eqtl_class_category to index
		self.eqtl_class_categories = eqtl_dataset_names
		self.eqtl_category_to_index = {}
		for ii, val in enumerate(self.eqtl_class_categories):
			self.eqtl_category_to_index[val] = ii

		# Initialize causal variant to get effect size variances (done in smart way)
		self.delta_sq, self.gene_cis_h2, gene_ld_scores, self.eqtl_resid_vars = initialize_variant_to_gene_effect_size_variances(genes, gene_info, self.n_reg_snps, self.n_eqtl_classes, self.eqtl_category_to_index)


		# Initialize non-mediated per snp h2 (gamma_sq) and gene_trait varaince (alpha_sq)
		mesc_style_params = linear_regression(np.hstack((gwas_variant_ld_scores, gene_ld_scores)), gwas_E_beta_sq)

		self.gamma_sq = mesc_style_params[:self.n_non_med_categories]
		self.alpha_sq = mesc_style_params[self.n_non_med_categories:]

		# Compute current estimates of non-mediated and mediated heritabilities
		self.update_nm_h2()
		self.update_med_h2(genes, gene_info)
		self.update_average_eqtl_h2(genes)

		# Get residual gwas_E_beta_sq
		self.gwas_resid_E_beta_sq = np.copy(gwas_E_beta_sq) - np.dot(gwas_variant_ld_scores, self.gamma_sq) - np.dot(gene_ld_scores, self.alpha_sq)

		self.update_gwas_resid_var()
		'''
		if self.fixed_variance:
			self.gwas_resid_var = (4.0*(np.mean(self.gwas_E_beta_sq))/100000.0) + 2.0*np.square(1/100000)
			self.eqtl_resid_vars = fixed_variance_update_eqtl_resid_vars(200, gene_info, genes, self.eqtl_category_to_index)

		else:
			# Compute correct estimate of gwas resid var
			self.update_gwas_resid_var()
		'''

		return


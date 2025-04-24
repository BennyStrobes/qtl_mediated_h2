import numpy as np 
import os
import sys
import pdb
from scipy.stats import invgamma








class RSS_GIBBS(object):
	def __init__(self, LD=None, marginal_beta=None, N=None, delta=None,delta_indices=None, X=None, Y=None):
		self.LD = LD
		self.marginal_beta = marginal_beta
		self.X = X
		self.Y = Y
		self.N = N
		self.delta = delta
		self.delta_indices = delta_indices
		print('##############################################################################')
		print('Non-mediated variant and expression mediated heritability analysis by optimizing RSS likelihood with Gibbs sampling')
		print('##############################################################################')
	def fit(self, update_residual_variance=True, max_iter=3000, burn_in_iter=2000, standardize_expression=True):
		# Add model fitting parameters
		self.max_iter = max_iter
		self.burn_in_iter = burn_in_iter
		self.update_residual_variance = update_residual_variance
		self.standardize_expression = standardize_expression

		# Quick error check
		if self.LD is None or self.marginal_beta is None or self.N is None:
			raise ValueError('RSS_VI_VARIANT_ONLY requires LD, marginal betas, and GWAS sample sizes')

		# If we don't have reference X, we cannot update residual variance
		if self.X is None or self.Y is None:
			print('Cannot update residual variance without individual level data')
			self.update_residual_variance = False

		# Initiailze model parameters
		self.initialize_parameters()

		# Perform iterative optimzation
		for itera in range(self.max_iter):
			# Updated mediated effect size of each gene
			self.alpha_updates(mean_field=True)

			# Update non-mediated effect size of each snp
			self.beta_updates(mean_field=True)

			# Update genetic variance parameters
			self.beta_variance_updates()
			self.alpha_variance_updates()

			# Update residual variance parameters
			if self.update_residual_variance:
				self.residual_variance_updates()

			local_h2 = np.var(self.pred_y)
			local_med_h2 = np.var(self.pred_y_med)
			local_nm_h2 = np.var(self.pred_y_nm)

			print(str(local_med_h2) + '\t' + str(local_nm_h2))
			print(local_h2)

			if itera > self.burn_in_iter:
				self.sampled_betas.append(self.beta)
				self.sampled_beta_variances.append(self.beta_variance)
				self.sampled_residual_variances.append(self.residual_variance)
				self.sampled_local_h2s.append(local_h2)
				self.sampled_local_med_h2s.append(local_med_h2)
				self.sampled_local_nm_h2s.append(local_nm_h2)
				
		return

	def beta_variance_updates(self):
		temp_beta_variance_b = np.sum(np.square(self.beta))/2
		temp_beta_variance_a = len(self.beta)/2

		self.beta_variance = 1.0/np.random.gamma(shape=temp_beta_variance_a, scale=1.0/temp_beta_variance_b)
		return

	def alpha_variance_updates(self):
		temp_alpha_variance_b = np.sum(np.square(self.alpha))/2
		temp_alpha_variance_a = len(self.alpha)/2

		self.alpha_variance = 1.0/np.random.gamma(shape=temp_alpha_variance_a, scale=1.0/temp_alpha_variance_b)
		return

	def residual_variance_updates(self):
		# Compute predicted y (non-mediated)
		self.pred_y_nm = np.dot(self.X,self.beta)
		# Compute predicted y (mediated)
		agg_gene_effects = np.copy(self.beta)*0.0
		for gg in range(self.G):
			agg_gene_effects[self.delta_indices[gg]] = agg_gene_effects[self.delta_indices[gg]] + self.delta[gg]*self.alpha[gg]
		self.pred_y_med = np.dot(self.X, agg_gene_effects)
		# Compute predicted y (total)
		self.pred_y = self.pred_y_med + self.pred_y_nm


		# Compute conditional sampling distribution
		temp_residual_variance_b = np.sum(np.square((self.Y - self.pred_y)))/2.0
		temp_residual_variance_a = len(self.Y)/2.0

		# Sample
		self.residual_variance = 1.0/np.random.gamma(shape=temp_residual_variance_a, scale=1.0/temp_residual_variance_b)
		
		return



	def beta_updates(self, mean_field=True):
		# Loop through snps
		if mean_field:
			for kk in np.random.permutation(range(self.K)):
				# De-residual current snp effect
				self.residual_marginal_beta = self.residual_marginal_beta + self.LD[kk,:]*self.beta[kk]

				# Calculate variance of snp effect
				temp_beta_var = 1.0/(((self.N)/self.residual_variance) + (1.0/self.beta_variance))
				# Calculate mean of snp effect
				temp_beta_mu = (temp_beta_var/self.residual_variance)*self.residual_marginal_beta[kk]*self.N

				# Need to sample beta
				self.beta[kk] = np.random.normal(loc=temp_beta_mu, scale=np.sqrt(temp_beta_var))

				# Re-residualize current snp effect
				self.residual_marginal_beta = self.residual_marginal_beta - self.LD[kk,:]*self.beta[kk]
		return

	def alpha_updates(self, mean_field=True):
		# Loop through genes
		if mean_field:
			for gg in np.random.permutation(range(self.G)):
				# De-residualize current gene effect
				predicted_gene_effect = self.delta[gg]*self.alpha[gg]
				self.residual_marginal_beta = self.residual_marginal_beta + np.dot(self.LD[:, self.delta_indices[gg]], predicted_gene_effect)
				

				# Calculate variance of gene effect
				temp_alpha_var = 1.0/(((self.precomputed_gene_variance_terms[gg])/self.residual_variance) + (1.0/self.alpha_variance))
				# Calculate mean of gene effect
				temp_alpha_mu = (temp_alpha_var/self.residual_variance)*self.N*np.dot(self.delta[gg], self.residual_marginal_beta[self.delta_indices[gg]])

				# Need to sample alpha
				self.alpha[gg] = np.random.normal(loc=temp_alpha_mu, scale=np.sqrt(temp_alpha_var))

				# re-residualize current gene effect
				predicted_gene_effect = self.delta[gg]*self.alpha[gg]
				self.residual_marginal_beta = self.residual_marginal_beta - np.dot(self.LD[:, self.delta_indices[gg]], predicted_gene_effect)

		return



	def initialize_parameters(self):
		# Number of snps
		self.K = len(self.marginal_beta)

		# Number of genes
		self.G = len(self.delta)

		# Variance of each beta predictor
		self.predictor_variances = np.ones(self.K)

		# Intialize residual vector of marginal betas(to setting where all effects are zero)
		self.residual_marginal_beta = np.copy(self.marginal_beta)

		# Initialize value on beta
		self.beta = np.zeros(self.K)

		# Initialize value on alpha
		self.alpha = np.zeros(self.G)

		# Intialize genetic variance parameters
		self.beta_variance = 1.0/1.0
		self.alpha_variance = 1.0/1.0

		# Initialize residual variance parameters
		self.residual_variance = 1.0/1.0

		# Standardize the genetic component of gene expression for each gene
		if self.standardize_expression:
			for g_iter in range(self.G):
				# compute gene variance
				gene_variance = np.dot(np.dot(self.delta[g_iter], self.LD[self.delta_indices[g_iter],:][:, self.delta_indices[g_iter]]), self.delta[g_iter])
				
				# Standardize expression
				self.delta[g_iter] = self.delta[g_iter]/np.sqrt(gene_variance)

		# Precompute gene-variance quantities
		self.precomputed_gene_variance_terms = np.zeros(self.G)
		for g_iter in range(self.G):
			gene_xt_x = self.N*self.LD[self.delta_indices[g_iter],:][:, self.delta_indices[g_iter]]
			var_term = np.dot(np.dot(self.delta[g_iter], gene_xt_x), self.delta[g_iter])
			self.precomputed_gene_variance_terms[g_iter] = var_term

		# Keep track of post burn in parameters
		self.sampled_betas = []
		self.alphas = []
		self.sampled_beta_variances = []
		self.sampled_residual_variances = []
		self.sampled_local_med_h2s = []
		self.sampled_local_nm_h2s = []
		self.sampled_local_h2s = []



		if np.array_equal(self.predictor_variances, np.ones(self.K)) == False:
			print('error: Method currently not implemented for non-1 predictor variance')
			pdb.set_trace()
		return


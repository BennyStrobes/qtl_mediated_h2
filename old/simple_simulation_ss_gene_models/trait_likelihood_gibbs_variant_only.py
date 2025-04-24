import numpy as np 
import os
import sys
import pdb
from scipy.stats import invgamma






class TRAIT_LIKELIHOOD_GIBBS_VARIANT_ONLY(object):
	def __init__(self, Y=None, X=None, beta_variance=None, residual_variance=None, max_iter=3000, burn_in_iter=2000, expected_gene_h2_for_beta_var_prior=.2, beta_var_prior_scale=1e-3, residual_variance_prior_scale=0.0):
		self.Y = Y
		self.X = X
		self.max_iter = max_iter
		self.burn_in_iter = burn_in_iter
		if beta_variance is None:
			self.update_beta_variance = True
			self.beta_variance = 1.0
		else:
			self.update_beta_variance = False
			self.beta_variance = beta_variance

		if residual_variance is None:
			self.update_residual_variance = True
			self.residual_variance = 1.0
		else:
			self.update_residual_variance = False
			self.residual_variance = residual_variance

		# Prior on beta variance temp_beta_variance_b/(temp_beta_variance_a - 1.0)
		n_snps = X.shape[1]
		self.beta_variance_beta_prior = beta_var_prior_scale
		expected_per_snp_h2 = expected_gene_h2_for_beta_var_prior/n_snps
		self.beta_variance_alpha_prior = (self.beta_variance_beta_prior/expected_per_snp_h2) + 1.0
		self.residual_variance_prior_scale = residual_variance_prior_scale


	def fit(self):
		# Quick error check
		if self.Y is None or self.X is None:
			raise ValueError('TRAIT_LIKELIHOOD_GIBBS_VARIANT_ONLY requires LD, z-scores, and summary statistic standard errorrs')

		# Initiailze model parameters
		self.initialize_parameters()

		# Perform iterative optimzation
		for itera in range(self.max_iter):
			# Update effect size and effect size variance of each snp
			self.beta_updates(mean_field=False)

			# Update genetic variance parameter
			if self.update_beta_variance:
				self.beta_variance_updates()

			# Update residual variance parameter
			if self.update_residual_variance:
				self.residual_variance_updates()


			# Record post-burn in iteration parameters
			if itera > self.burn_in_iter:
				self.sampled_betas.append(self.beta)
				self.sampled_beta_variances.append(self.beta_variance)
				self.sampled_residual_variances.append(self.residual_variance)

				# Sample local h2
				local_h2 = np.var(np.dot(self.X, self.beta))
				self.sampled_local_h2s.append(local_h2)

		# Organize post-burn in iteration parameters into neat matrices
		self.sampled_betas = np.asarray(self.sampled_betas)
		self.sampled_beta_variances = np.asarray(self.sampled_beta_variances)
		self.sampled_residual_variances = np.asarray(self.sampled_residual_variances)
		self.sampled_local_h2s = np.asarray(self.sampled_local_h2s)
		return


	def beta_updates(self, mean_field=False):
		if mean_field == True:
			print('assumption error: not currently implemented for mean field == true')
			pdb.set_trace()

		# Compute conditional sampling distribution
		beta_cov = np.linalg.inv((np.eye(self.K)*(1.0/self.beta_variance)) + (self.X_transpose_X/self.residual_variance))
		beta_mean = (1.0/self.residual_variance)*np.dot(beta_cov, np.dot(np.transpose(self.X), self.Y))

		# Sample
		self.beta = np.random.multivariate_normal(beta_mean, beta_cov)
		return

	def beta_variance_updates(self):
		# Compute conditional sampling distribution
		temp_beta_variance_b = (np.sum(np.square(self.beta))/2) + self.beta_variance_beta_prior
		temp_beta_variance_a = (len(self.beta)/2) + self.beta_variance_alpha_prior

		# Sample
		# This samples from an inverse gamma
		# Mean is temp_beta_variance_b/(temp_beta_variance_a - 1.0)
		self.beta_variance = 1.0/np.random.gamma(shape=temp_beta_variance_a, scale=1.0/temp_beta_variance_b)
		return

	def residual_variance_updates(self):
		# Compute predicted y
		pred_y = np.dot(self.X,self.beta)

		# Compute conditional sampling distribution
		temp_residual_variance_b = (np.sum(np.square((self.Y - pred_y)))/2.0) + self.residual_variance_prior_scale
		temp_residual_variance_a = (len(self.Y)/2.0) + self.residual_variance_prior_scale

		# Sample
		# This samples from an inverse gamma
		self.residual_variance = 1.0/np.random.gamma(shape=temp_residual_variance_a, scale=1.0/temp_residual_variance_b)
		return


	def initialize_parameters(self):
		# Number of snps
		self.K = self.X.shape[1]
		
		#self.N = np.mean(1.0/np.square(self.se))
		self.N = self.X.shape[0]

		# Compute X_transose X
		self.X_transpose_X = np.dot(np.transpose(self.X), self.X)

		# Intialize residual vector of z-scores (to setting where all effects are zero)
		self.residual_Y = np.copy(self.Y)

		# Initialize variational distribution on beta
		self.beta = np.zeros(self.K)

		# Intialize genetic variance parameters ( already done)
		#self.beta_variance = 1.0/1.0

		# Initialize sample Y
		self.predicted_Y = np.zeros(self.N)

		# Initialize residual variance (alread done)
		#self.residual_variance = 1.0/1.0

		# Keep track of post burn in parameters
		self.sampled_betas = []
		self.sampled_beta_variances = []
		self.sampled_residual_variances = []
		self.sampled_local_h2s = []

		return
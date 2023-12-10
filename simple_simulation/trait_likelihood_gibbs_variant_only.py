import numpy as np 
import os
import sys
import pdb
from scipy.stats import invgamma






class TRAIT_LIKELIHOOD_GIBBS_VARIANT_ONLY(object):
	def __init__(self, Y=None, X=None, update_beta_var=False, max_iter=3000, burn_in_iter=2000):
		self.Y = Y
		self.X = X
		self.max_iter = max_iter
		self.burn_in_iter = burn_in_iter
		self.update_beta_var = update_beta_var
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
			if self.update_beta_var:
				self.beta_variance_updates()

			# Update residual variance parameter
			self.update_residual_variance()


			# Record post-burn in iteration parameters
			if itera > self.burn_in_iter:
				self.sampled_betas.append(self.beta)
				self.sampled_beta_variances.append(self.beta_variance)
				self.sampled_residual_variances.append(self.residual_variance)

				# Sample local h2
				local_h2 = np.var(np.dot(self.X, self.beta))
				self.sampled_local_h2s.append(local_h2)

		return


	def beta_updates(self, mean_field=False):
		# Compute conditional sampling distribution
		beta_cov = np.linalg.inv((np.eye(self.K)*(1.0/self.beta_variance)) + (self.X_transpose_X/self.residual_variance))
		beta_mean = (1.0/self.residual_variance)*np.dot(beta_cov, np.dot(np.transpose(self.X), self.Y))

		# Sample
		self.beta = np.random.multivariate_normal(beta_mean, beta_cov)


		return

	def beta_variance_updates(self):
		# Compute conditional sampling distribution
		temp_beta_variance_b = np.sum(np.square(self.beta))/2
		temp_beta_variance_a = len(self.beta)/2

		# Sample
		self.beta_variance = 1.0/np.random.gamma(shape=temp_beta_variance_a, scale=1.0/temp_beta_variance_b)
		return

	def update_residual_variance(self):
		# Compute predicted y
		pred_y = np.dot(self.X,self.beta)

		# Compute conditional sampling distribution
		temp_residual_variance_b = np.sum(np.square((self.Y - pred_y)))/2.0
		temp_residual_variance_a = len(self.Y)/2.0

		# Sample
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

		# Intialize genetic variance parameters
		self.beta_variance = 1.0/1.0

		# Initialize sample Y
		self.predicted_Y = np.zeros(self.N)

		# Initialize residual variance
		self.residual_variance = 1.0/1.0

		# Keep track of post burn in parameters
		self.sampled_betas = []
		self.sampled_beta_variances = []
		self.sampled_residual_variances = []
		self.sampled_local_h2s = []

		return

import numpy as np 
import os
import sys
import pdb
from scipy.stats import invgamma






class TRAIT_LIKELIHOOD_GIBBS_VARIANT_ONLY(object):
	def __init__(self, Y=None, X=None, max_iter=10000, burn_in_iter=2000):
		self.Y = Y
		self.X = X
		self.max_iter = max_iter
		self.burn_in_iter = burn_in_iter
		n_snps = X.shape[1]



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

			# Update residual variance parameter
			self.residual_variance_updates()

			# Update genetic variance parameter
			self.beta_variance_updates()


			# Record post-burn in iteration parameters
			if itera > self.burn_in_iter:
				self.sampled_betas.append(self.beta)
				self.sampled_beta_variances.append(self.beta_variance)
				self.sampled_residual_variances.append(self.residual_variance)
				self.sampled_pves.append(self.pve)

				# Sample local h2
				#local_h2 = np.var(np.dot(self.X, self.beta))
				#self.sampled_local_h2s.append(local_h2)

		# Organize post-burn in iteration parameters into neat matrices
		self.sampled_betas = np.asarray(self.sampled_betas)
		self.sampled_beta_variances = np.asarray(self.sampled_beta_variances)
		self.sampled_residual_variances = np.asarray(self.sampled_residual_variances)
		self.sampled_pves = np.asarray(self.sampled_pves)
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
		pred_y = np.dot(self.X,self.beta)
		var_pred_Y = np.var(pred_y)

		cur_pve = var_pred_Y/(var_pred_Y + self.residual_variance)

		new_pve = cur_pve + np.random.uniform(-.01,.01)

		if new_pve > 1.0:
			diff = new_pve - 1.0
			new_pve = 1.0 - diff
		elif new_pve < 0.0:
			new_pve = np.abs(new_pve)

		self.pve = new_pve

		self.beta_variance = self.pve*self.residual_variance/(self.K*(1.0 -self.pve))
		return

	def residual_variance_updates(self):
		# Compute predicted y
		pred_y = np.dot(self.X,self.beta)

		# Compute conditional sampling distribution
		temp_residual_variance_b = (np.sum(np.square((self.Y - pred_y)))/2.0) 
		temp_residual_variance_a = (len(self.Y)/2.0)

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
		self.beta_variance = 1.0/self.K

		# Initialize sample Y
		self.predicted_Y = np.zeros(self.N)

		# Initialize residual variance (alread done)
		self.residual_variance = 1.0/1.0

		self.pve = 1.0

		# Keep track of post burn in parameters
		self.sampled_betas = []
		self.sampled_beta_variances = []
		self.sampled_residual_variances = []
		self.sampled_local_h2s = []
		self.sampled_pves = []

		return
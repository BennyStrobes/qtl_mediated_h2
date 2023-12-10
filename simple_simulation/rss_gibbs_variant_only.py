import numpy as np 
import os
import sys
import pdb
from scipy.stats import invgamma








class RSS_GIBBS_VARIANT_ONLY(object):
	def __init__(self, LD=None, marginal_beta=None, N=None, X=None, Y=None, update_residual_variance=True, max_iter=3000, burn_in_iter=2000):
		self.LD = LD
		self.marginal_beta = marginal_beta
		self.X = X
		self.Y = Y
		self.N = N
		self.max_iter = max_iter
		self.burn_in_iter = burn_in_iter
		self.update_residual_variance = update_residual_variance
		print('##############################################################################')
		print('Variant heritability analysis by optimizing RSS likelihood with Gibbs sampling')
		print('##############################################################################')
	def fit(self):
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
			self.itera = itera
			# Update effect size and effect size variance of each snp
			self.beta_updates(mean_field=True)

			# Update genetic variance parameters
			#self.beta_variance_updates()

			# Update residual variance parameters
			if self.update_residual_variance:
				self.residual_variance_updates()

			local_h2 = np.var(self.pred_y)

			if itera > self.burn_in_iter:
				self.sampled_betas.append(np.copy(self.beta))
				self.sampled_beta_variances.append(self.beta_variance)
				self.sampled_residual_variances.append(self.residual_variance)
				self.sampled_local_h2s.append(local_h2)
		return

	def beta_variance_updates(self):
		temp_beta_variance_b = np.sum(np.square(self.beta))/2
		temp_beta_variance_a = len(self.beta)/2

		self.beta_variance = 1.0/np.random.gamma(shape=temp_beta_variance_a, scale=1.0/temp_beta_variance_b)
		return

	def residual_variance_updates(self):
		# Compute predicted y
		self.pred_y = np.dot(self.X,self.beta)

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

				# Update variance of snp effect
				temp_beta_var = 1.0/(((self.N*self.predictor_variances[kk])/self.residual_variance) + (1.0/self.beta_variance))
				temp_beta_mu = (temp_beta_var/self.residual_variance)*self.residual_marginal_beta[kk]*self.N*self.predictor_variances[kk]


				# Need to sample beta_mu
				self.beta[kk] = np.random.normal(loc=temp_beta_mu, scale=np.sqrt(temp_beta_var))

				# Re-residualize current snp effect
				self.residual_marginal_beta = self.residual_marginal_beta - self.LD[kk,:]*self.beta[kk]
		return



	def initialize_parameters(self):
		# Number of snps
		self.K = len(self.marginal_beta)

		# Variance of each predictor
		self.predictor_variances = np.ones(self.K)

		# Intialize residual vector of marginal betas(to setting where all effects are zero)
		self.residual_marginal_beta = np.copy(self.marginal_beta)

		# Initialize variational distribution on beta
		self.beta = np.zeros(self.K)

		# Intialize genetic variance parameters
		self.beta_variance = 1.0/1.0

		# Initialize residual variance parameters
		self.residual_variance = 1.0/1.0


		# Keep track of post burn in parameters
		self.sampled_betas = []
		self.sampled_beta_variances = []
		self.sampled_residual_variances = []
		self.sampled_local_h2s = []

		return


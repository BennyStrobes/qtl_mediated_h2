import numpy as np 
import os
import sys
import pdb








class RSS_VI_VARIANT_ONLY(object):
	def __init__(self, LD=None, z=None, se=None, max_iter=2000):
		self.LD = LD
		self.z = z
		self.se = se
		self.max_iter = max_iter
	def fit(self):
		# Quick error check
		if self.LD is None or self.z is None or self.se is None:
			raise ValueError('RSS_VI_VARIANT_ONLY requires LD, z-scores, and summary statistic standard errorrs')

		# Initiailze model parameters
		self.initialize_parameters()

		# Perform iterative optimzation
		for itera in range(self.max_iter):
			# Update effect size and effect size variance of each snp
			self.beta_updates(mean_field=True)

			# Update genetic variance parameters
			self.beta_variance_updates()

			print(self.beta_variance)
		pdb.set_trace()

	def beta_variance_updates(self):
		self.beta_variance_b = np.sum(np.square(self.beta_mu)) + np.sum(self.beta_var)
		self.beta_variance_a = len(self.beta_mu)

		self.beta_variance = self.beta_variance_b/self.beta_variance_a


	def beta_updates(self, mean_field=True):
		# Loop through snps
		if mean_field:
			for kk in np.random.permutation(range(self.K)):
				# De-residual current snp effect 
				self.residual_z = self.residual_z + (self.LD[kk,:]*self.beta_mu[kk]/self.se[kk])

				# Update variance of snp effect
				self.beta_var[kk] = 1.0/((1.0/np.square(self.se[kk])) + (1.0/self.beta_variance))
				self.beta_mu[kk] = (self.beta_var[kk])*self.residual_z[kk]/self.se[kk]


				# Re-residualize current snp effect
				self.residual_z = self.residual_z - (self.LD[kk,:]*self.beta_mu[kk]/self.se[kk])
		else:
			residual_prec = 1.428
			tau = 1.0/self.beta_variance
			pdb.set_trace()
			S = np.linalg.inv((residual_prec)*(self.N*self.LD) + np.diag(np.ones(self.LD.shape[0])*tau))
			mu = residual_prec*np.dot(S, self.z*np.sqrt(self.N))

			self.beta_mu = np.copy(mu)
			self.beta_var = np.diag(S)
			pdb.set_trace()


	def initialize_parameters(self):
		# Number of snps
		self.K = len(self.z)
		
		#self.N = np.mean(1.0/np.square(self.se))
		self.N = 50000

		# Intialize residual vector of z-scores (to setting where all effects are zero)
		self.residual_z = np.copy(self.z)

		# Initialize variational distribution on beta
		self.beta_mu = np.zeros(self.K)
		self.beta_var = np.ones(self.K)

		# Intialize genetic variance parameters
		self.beta_variance_a = 1.0
		self.beta_variance_b = 1.0
		self.beta_variance = self.beta_variance_b/self.beta_variance_a


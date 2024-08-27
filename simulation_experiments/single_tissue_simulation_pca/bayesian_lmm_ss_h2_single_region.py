import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm








class Bayesian_LMM(object):
	def __init__(self, beta, ld, NN):
		self.gwas_beta = beta
		self.gwas_beta_var = 1.0/NN
		self.LD = ld

		# Number of components
		self.KK = ld.shape[1]

		return


	def fit(self, total_iterations=15000, burn_in_iterations=10000, update_resid_var_bool=False):
		""" Fit the model.
		"""
		# Initialize model params
		self.initialize_variables()

		# Keep track of iterations
		self.itera = 0

		# Iterative Gibbs sampling algorithm
		for itera in range(total_iterations):
			# Update gamma
			self.update_gamma()
	
			# Update gamma_var
			self.update_gamma_var()

			if update_resid_var_bool:
				print('not yet implemented!')
				pdb.set_trace()
				if itera > 10:
					self.update_resid_var()

			# Update iteration number
			self.itera = self.itera + 1
	
			'''
			if np.mod(itera, 50) == 0.0:
				print(self.resid_var)
				print(self.gamma_var)
				print(itera)
			'''

			if itera > burn_in_iterations:
				self.sampled_gamma_vars.append(self.gamma_var)
				self.sampled_resid_vars.append(self.resid_var)
				self.sampled_gammas.append(np.copy(self.gamma))

		self.sampled_gamma_vars = np.asarray(self.sampled_gamma_vars)
		self.sampled_h2 = self.KK*self.sampled_gamma_vars
		self.sampled_resid_vars = np.asarray(self.sampled_resid_vars)
		self.sampled_gammas = np.asarray(self.sampled_gammas)

		return

	def update_resid_var(self, v0=0.0, s_sq=0.0):
		weights = np.ones(self.QQ)
		vv = np.sum(weights) + v0

		tau_sq = np.sum(weights*np.square(self.gwas_pca_beta_resid)/self.gwas_beta_var) + s_sq

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma(vv/2 + 1e-5, scale=tau_sq/2 + 1e-5)
		# Sample from it
		self.resid_var = invgamma_dist.rvs(size=1)[0]
		
		return


	def update_gamma_var(self, v0=0.0, s_sq=0.0):
		weights = np.ones(self.KK)
		vv = np.sum(weights) + v0
		tau_sq = np.sum(weights*np.square(self.gamma)) + s_sq			

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma(vv/2 + 1e-8, scale=tau_sq/2 + 1e-8)
		# Sample from it
		self.gamma_var = invgamma_dist.rvs(size=1)[0]

		#self.gamma_var = tau_sq/vv

		return



	def update_gamma(self, univariate_updates=True):

		if univariate_updates:
			for kk in np.random.permutation(range(self.KK)):
				# Re include current effects
				self.gwas_beta_resid = self.gwas_beta_resid + self.LD[:, kk]*self.gamma[kk]


				# Get posterior distributions
				posterior_var = 1.0/((1.0/(self.resid_var*self.gwas_beta_var)) + (1.0/self.gamma_var))
				posterior_mean = ((self.gwas_beta_resid[kk])/(self.resid_var*self.gwas_beta_var))*posterior_var

				# Sample from posterior distribution
				self.gamma[kk] = np.random.normal(loc=posterior_mean, scale=np.sqrt(posterior_var))


				# Remove updated current effects
				self.gwas_beta_resid = self.gwas_beta_resid - self.LD[:, kk]*self.gamma[kk]

		elif univariate_updates == False:
			print('not yet implemented')
			pdb.set_trace()
		return



	def initialize_variables(self):
		# Initialize causal effcts
		self.gamma = np.zeros(self.KK)

		# Initialize variance parameters
		self.gamma_var = 1.0
		self.resid_var = 1.0

		# Remove causal effects from gwas beta
		self.gwas_beta_resid = np.copy(self.gwas_beta) - np.dot(self.LD, self.gamma)


		# Keep track of sampled gamma_vars
		self.sampled_gamma_vars = []
		self.sampled_resid_vars = []
		self.sampled_gammas = []

		return
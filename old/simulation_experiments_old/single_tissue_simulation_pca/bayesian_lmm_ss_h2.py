import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm








class Bayesian_LMM_SS_h2_inference(object):
	def __init__(self, gwas_beta, variant_ld_scores, NN, n_snps):
		self.gwas_beta = gwas_beta
		self.gwas_beta_var = 1.0/NN
		self.var_ldscores = variant_ld_scores

		# Number of components
		self.KK = len(gwas_beta)

		# Number of snps
		self.n_snps = n_snps


	def fit(self, total_iterations=15000, burn_in_iterations=10000, update_resid_var_bool=True):
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
				if itera > 100:
					self.update_resid_var()

			# Update iteration number
			self.itera = self.itera + 1
	
			'''
			if np.mod(itera, 1000) == 0.0:
				print(self.resid_var)
				print(self.gamma_var)
				print(itera)
			'''

			if itera > burn_in_iterations:
				self.sampled_gamma_vars.append(self.gamma_var)
				self.sampled_resid_vars.append(self.resid_var)

		self.sampled_gamma_vars = np.asarray(self.sampled_gamma_vars)
		self.sampled_h2 = self.n_snps*self.sampled_gamma_vars
		self.sampled_resid_vars = np.asarray(self.sampled_resid_vars)

		return


	def update_resid_var(self, v0=0.0, s_sq=0.0):
		weights = np.ones(self.KK)
		vv = np.sum(weights) + v0

		tau_sq = np.sum(weights*np.square(self.gwas_beta_resid)/self.gwas_beta_var) + s_sq

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma((vv/2) + 1e-5, scale=(tau_sq/2) + 1e-5)
		# Sample from it
		self.resid_var = invgamma_dist.rvs(size=1)[0]
		
		return


	def update_gamma_var(self, v0=0.0, s_sq=0.0):
		weights = np.ones(self.KK)
		vv = np.sum(weights) + v0
		tau_sq = np.sum(weights*np.square(self.gamma)/self.var_ldscores) + s_sq			

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma((vv/2) + 1e-5, scale=(tau_sq/2)+1e-5)
		# Sample from it
		self.gamma_var = invgamma_dist.rvs(size=1)[0]

		return

	def update_gamma(self):
		# Re-include current effects
		self.gwas_beta_resid = self.gwas_beta_resid + self.gamma

		# Compute posterior distribution
		marginal_gamma_vars = self.gamma_var*self.var_ldscores
		posterior_vars = 1.0/((1.0/(self.resid_var*self.gwas_beta_var)) + (1.0/marginal_gamma_vars))
		posterior_means = (self.gwas_beta_resid/(self.resid_var*self.gwas_beta_var))*posterior_vars

		# Sample from posterior distribution
		self.gamma = np.random.normal(loc=posterior_means, scale=np.sqrt(posterior_vars))

		# Remove current effects
		self.gwas_beta_resid = self.gwas_beta_resid - self.gamma

		return



	def initialize_variables(self):
		# Initialize causal effcts
		self.gamma = np.zeros(self.KK)

		# Initialize variance parameters
		self.gamma_var = 1.0
		self.resid_var = 1.0

		# Remove causal effects from gwas beta
		self.gwas_beta_resid = np.copy(self.gwas_beta) - self.gamma

		# Keep track of sampled gamma_vars
		self.sampled_gamma_vars = []
		self.sampled_resid_vars = []
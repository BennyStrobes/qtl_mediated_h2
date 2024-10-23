import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import time
import bayesian_lmm_ss_h2_single_region




class Bayesian_LMM_RSS_h2_inference(object):
	def __init__(self, window_info, N_gwas):
		self.N_gwas = N_gwas
		self.window_info = window_info

		self.windows = np.sort([*self.window_info])

		# Number of snps
		self.KK = 0.0
		for window in self.windows:
			self.KK = self.KK + self.window_info[window]['n_snps']

		return



	def fit(self, total_iterations=15000, burn_in_iterations=10000, update_resid_var_bool=False, univariate_updates=True,cc=1e-6):
		""" Fit the model.
		"""
		self.univariate_updates = univariate_updates
		self.cc=cc
		# Initialize model params
		self.initialize_variables()


		# Keep track of iterations
		self.itera = 0

		# Iterative Gibbs sampling algorithm
		for itera in range(total_iterations):
			# Loop through windows
			for window_name in self.windows:
				# Update gamma, delta, and alpha in each window seperately
				self.update_gamma_in_single_window(window_name)

			# Update gamma_var
			self.update_gamma_var(cc=self.cc)


			# Update resid var
			if update_resid_var_bool:
				print('not currently updated')
				pdb.set_trace()
				#gwas_resid_vars = self.update_gwas_resid_var()
				#self.update_eqtl_resid_vars()

			# Update iteration number
			self.itera = self.itera + 1

	
			if itera > burn_in_iterations:
				nm_h2 = self.gamma_var*self.KK
				nm_h2_ld_depen = self.get_ld_dependent_h2s()
				self.sampled_nm_h2_full.append(nm_h2_ld_depen)
				self.sampled_nm_h2_independent.append(nm_h2)


		self.sampled_nm_h2_full = np.asarray(self.sampled_nm_h2_full)
		self.sampled_nm_h2_independent = np.asarray(self.sampled_nm_h2_independent)


		return

	def get_ld_dependent_h2s(self):
		nm = 0.0
		for window in self.windows:
			nm = nm + self.window_nm_h2[window]
		return nm



	def update_gamma_var(self, v0=0.0, s_sq=0.0, cc=1e-6):
		# First get middle gammas
		all_gammas = []
		for window_name in self.windows:
			window_gammas = self.gamma[window_name]
			all_gammas.append(window_gammas)
		all_gammas = np.hstack(all_gammas)

		vv = len(all_gammas) + v0
		tau_sq = np.sum(np.square(all_gammas)) + s_sq

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
		# Sample from it
		self.gamma_var = invgamma_dist.rvs(size=1)[0]

		return

	def update_gamma_in_single_window(self, window_name):
		# Load in relevent quantities for this window
		window_ld = np.load(self.window_info[window_name]['ld_file'])

		#order_bool = np.random.choice([0,1])

		# Update gamma and alpha in single window
		if self.univariate_updates:
			self.univariate_update_gamma_in_single_window(window_name, window_ld)
		else:
			self.multivariate_update_gamma_in_single_window(window_name, window_ld)


		self.update_window_genetic_var(window_name, window_ld)

		return

	def update_window_genetic_var(self, window_name, LD_mat):
		window_nm_causal_effects = np.copy(self.gamma[window_name])
		self.window_nm_h2[window_name] = np.dot(np.dot(window_nm_causal_effects, LD_mat), window_nm_causal_effects)
		return



	def univariate_update_gamma_in_single_window(self, window_name, LD):
		# Extract relevent quantities
		window_gamma = self.gamma[window_name]
		window_gwas_beta_resid = self.gwas_beta_resid[window_name]
		window_gwas_resid_var = self.gwas_resid_var[window_name]

		# Precompute posterior variance (same for all snps)
		snp_posterior_var = 1.0/((self.N_gwas/window_gwas_resid_var) + (1.0/self.gamma_var))

		# Loop through genetic elements in random order
		window_genetic_elements = self.window_info[window_name]['genetic_elements']
		n_genetic_elements = len(window_genetic_elements)
		ordering = np.random.permutation(np.arange(n_genetic_elements))
		for genetic_element_index in ordering:
			# Name of genetic element
			genetic_element_class, genetic_element_name = window_genetic_elements[genetic_element_index]

			# Updates for snps
			if genetic_element_class == 'snp':
				window_gamma, window_gwas_beta_resid = self.gamma_update_for_single_snp(window_gamma, window_gwas_beta_resid, LD, genetic_element_name, snp_posterior_var, window_gwas_resid_var)



		# Re update quantities
		self.gamma[window_name] = window_gamma
		self.gwas_beta_resid[window_name] = window_gwas_beta_resid

		return


	def multivariate_update_gamma_in_single_window(self, window_name, LD):
		# Extract relevent quantities
		window_gamma = self.gamma[window_name]
		window_gwas_beta_resid = self.gwas_beta_resid[window_name]
		window_gwas_resid_var = self.gwas_resid_var[window_name]
		KK = len(window_gamma)

		inv_cov = (LD/(window_gwas_resid_var/self.N_gwas)) + (np.eye(KK)/self.gamma_var)
		cov = np.linalg.inv(inv_cov)
		posterior_mean = (1.0/(window_gwas_resid_var/self.N_gwas))*np.dot(cov,window_gwas_beta_resid)

		self.gamma[window_name] = np.random.multivariate_normal(mean=posterior_mean, cov=cov)

		return




	def gamma_update_for_single_snp(self, window_gamma, window_gwas_beta_resid, LD, snp_index, snp_posterior_var, window_gwas_resid_var):

		# Re include current effect
		window_gwas_beta_resid = window_gwas_beta_resid + (LD[:, snp_index]*window_gamma[snp_index])

		# Compute posterior mean for this snp
		snp_posterior_mean = snp_posterior_var*(self.N_gwas/window_gwas_resid_var)*window_gwas_beta_resid[snp_index]

		# Sample from posterior distribution
		window_gamma[snp_index] = np.random.normal(loc=snp_posterior_mean, scale=np.sqrt(snp_posterior_var))

		# Remove updated effect
		window_gwas_beta_resid = window_gwas_beta_resid - (LD[:, snp_index]*window_gamma[snp_index])

		return window_gamma, window_gwas_beta_resid




	def initialize_variables(self):
		# Initialize causal variant effcts
		self.gamma = {}
		self.gwas_beta_resid = {}
		self.gwas_resid_var = {}
		self.window_med_h2 = {}
		self.window_nm_h2 = {}
		for window_name in self.windows:
			self.gamma[window_name] = np.zeros(self.window_info[window_name]['n_snps'])
			self.gwas_beta_resid[window_name] = np.copy(self.window_info[window_name]['beta'])

			# Add sum of square terms (precomputed)
			self.gwas_resid_var[window_name] = 1.0

			self.window_med_h2[window_name] = 0.0
			self.window_nm_h2[window_name] = 0.0


		# Create list of genetic elements in each window
		for window_name in self.windows:
			window_snps = self.window_info[window_name]['rsids']
			window_genetic_elements = []
			for ii in range(len(window_snps)):
				window_genetic_elements.append(('snp', ii))
			# Add to window info
			self.window_info[window_name]['genetic_elements'] = window_genetic_elements




		# Initialize variance parameters
		self.gamma_var = 1e-5


		# Keep track of sampled gamma_vars
		self.sampled_nm_h2_full = []
		self.sampled_nm_h2_independent = []

		return





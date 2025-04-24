import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import time

def update_gamma_from_single_window(LD, gamma_vec, gwas_beta_resid, gamma_var, N_gwas, univariate=True):
	if univariate == False:
		KK = len(gamma_vec)
		SS_inv = ((LD*N_gwas) + np.eye(KK)/gamma_var)
		SS= np.linalg.inv(SS_inv)

		posterior_mean = N_gwas*np.dot(SS,gwas_beta_resid)
		gamma_vec = np.random.multivariate_normal(mean=posterior_mean, cov=SS)



	else:
		# N-snps in window
		KK = len(gamma_vec)

		# Precompute posterior variance (same for all snps)
		posterior_var = 1.0/(N_gwas + (1.0/gamma_var))

		# Update each snp in turn
		for snp_index in np.random.permutation(range(KK)):
			# Re include current effect
			gwas_beta_resid = gwas_beta_resid + (LD[snp_index,:]*gamma_vec[snp_index])


			# Compute posterior mean for this snp
			posterior_mean = posterior_var*N_gwas*gwas_beta_resid[snp_index]

			# Sample from posterior distribution
			gamma_vec[snp_index] = np.random.normal(loc=posterior_mean, scale=np.sqrt(posterior_var))


			# Remove updated effect
			gwas_beta_resid = gwas_beta_resid - (LD[snp_index,:]*gamma_vec[snp_index])

	return gamma_vec, gwas_beta_resid





class Bayesian_LMM_RSS_h2_inference(object):
	def __init__(self, gwas_beta, n_gwas_individuals, window_names, window_info, temp_output_file):
		self.gwas_beta = gwas_beta
		self.N_gwas = n_gwas_individuals
		self.window_names = window_names
		self.window_info = window_info

		self.temp_output_file = temp_output_file

		# Number of snps
		self.KK = len(self.gwas_beta)


	def fit(self, total_iterations=15000, burn_in_iterations=10000):
		""" Fit the model.
		"""
		# Initialize model params
		self.initialize_variables()

		# Keep track of iterations
		self.itera = 0

		# Iterative Gibbs sampling algorithm
		for itera in range(total_iterations):
			t1 = time.time()
			# Update gamma
			self.update_gamma()
	
			# Update gamma_var
			self.update_gamma_var()

			# Update iteration number
			self.itera = self.itera + 1
	

			if itera > burn_in_iterations:
				self.sampled_gamma_vars.append(self.gamma_var)

				if np.mod(itera, 5) == 0:
					t = open(self.temp_output_file,'w')
					t.write('sample_iter\tsampled_h2\n')
					for ii,sample_gamma_var in enumerate(self.sampled_gamma_vars):
						t.write(str(ii) + '\t' + str(sample_gamma_var*self.KK) + '\n')
					t.close()

			t2 = time.time()
			print('Iteration ' + str(itera) + ' completed in ' + str(t2-t1) + ' seconds')


		return



	def update_gamma_var(self, v0=2.0, s_sq=0.0):
		# First get middle gammas
		middle_gammas = []
		for window_name in self.window_names:
			window_middle_gammas = self.window_to_gamma[window_name][self.window_info[window_name]['middle_indices']]
			middle_gammas.append(window_middle_gammas)
		middle_gammas = np.hstack(middle_gammas)

		vv = len(middle_gammas) + v0
		tau_sq = np.sum(np.square(middle_gammas)) + s_sq

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma(vv/2, scale=tau_sq/2)
		# Sample from it
		self.gamma_var = invgamma_dist.rvs(size=1)[0]

		return

	def update_gamma(self):
		for window_name in self.window_names:


			gamma_vec, gwas_beta_resid_vec = update_gamma_from_single_window(np.load(self.window_info[window_name]['ld_file']), self.window_to_gamma[window_name], self.window_to_gwas_beta_resid[window_name], self.gamma_var, self.N_gwas)
			self.window_to_gamma[window_name] = gamma_vec
			self.window_to_gwas_beta_resid[window_name] = gwas_beta_resid_vec

		return



	def initialize_variables(self):
		# Initialize causal effcts
		self.window_to_gamma = {}
		self.window_to_gwas_beta_resid = {}


		for window_name in self.window_names:
			self.window_to_gamma[window_name] = np.zeros(len(self.window_info[window_name]['positions']))
			self.window_to_gwas_beta_resid[window_name] = np.copy(self.gwas_beta[self.window_info[window_name]['positions']])

		# Initialize variance parameters
		self.gamma_var = 1e-6

		# Keep track of sampled gamma_vars
		self.sampled_gamma_vars = []
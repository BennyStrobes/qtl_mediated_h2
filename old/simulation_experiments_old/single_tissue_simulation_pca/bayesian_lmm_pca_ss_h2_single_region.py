import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm








class Bayesian_LMM(object):
	def __init__(self, gwas_pca_beta, Q_mat, NN):
		self.gwas_pca_beta = gwas_pca_beta
		self.gwas_beta_var = 1.0/NN
		self.Q_mat = Q_mat

		# Number of components
		self.QQ = len(gwas_pca_beta)
		self.KK = Q_mat.shape[1]

		self.QQ_sq = np.sum(np.square(self.Q_mat), axis=0)

		return


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
		#invgamma_dist = invgamma(vv/2 + 1e-8, scale=tau_sq/2 + 1e-8)
		# Sample from it
		#self.gamma_var = invgamma_dist.rvs(size=1)[0]

		self.gamma_var = tau_sq/vv

		return



	def update_gamma(self, univariate_updates=True):

		if univariate_updates:
			for kk in np.random.permutation(range(self.KK)):
				# Re include current effects
				self.gwas_pca_beta_resid = self.gwas_pca_beta_resid + self.Q_mat[:, kk]*self.gamma[kk]


				# Get posterior distributions
				posterior_var = 1.0/((self.QQ_sq[kk]/(self.resid_var*self.gwas_beta_var)) + (1.0/self.gamma_var))
				posterior_mean = (np.dot(self.gwas_pca_beta_resid,self.Q_mat[:,kk])/(self.resid_var*self.gwas_beta_var))*posterior_var

				# Sample from posterior distribution
				self.gamma[kk] = np.random.normal(loc=posterior_mean, scale=np.sqrt(posterior_var))


				# Remove updated current effects
				self.gwas_pca_beta_resid = self.gwas_pca_beta_resid - self.Q_mat[:, kk]*self.gamma[kk]
		elif univariate_updates == False:
			cov_inv = (np.dot(np.transpose(self.Q_mat), self.Q_mat)/(self.resid_var*self.gwas_beta_var)) + (np.eye(len(self.gamma))/self.gamma_var)

			posterior_cov = np.linalg.inv(cov_inv)
			posterior_mean = (1.0/(self.resid_var*self.gwas_beta_var))*np.dot(posterior_cov, np.dot(np.transpose(self.Q_mat), self.gwas_pca_beta_resid))
			self.gamma = np.random.multivariate_normal(posterior_mean, posterior_cov)
			#posterior_mean = (np.dot(self.gwas_pca_beta_resid,self.Q_mat[:,kk])/(self.resid_var*self.gwas_beta_var))*posterior_var

		return



	def initialize_variables(self):
		# Initialize causal effcts
		self.gamma = np.zeros(self.KK)

		# Initialize variance parameters
		self.gamma_var = 1.0
		self.resid_var = 1.0

		# Remove causal effects from gwas beta
		self.gwas_pca_beta_resid = np.copy(self.gwas_pca_beta) - np.dot(self.Q_mat, self.gamma)


		# Keep track of sampled gamma_vars
		self.sampled_gamma_vars = []
		self.sampled_resid_vars = []
		self.sampled_gammas = []

		return
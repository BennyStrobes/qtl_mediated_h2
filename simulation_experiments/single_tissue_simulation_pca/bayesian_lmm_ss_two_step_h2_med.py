import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import time







class Bayesian_LMM_SS_h2_med_inference(object):
	def __init__(self, gwas_beta, variant_ld_scores, N_gwas, n_gwas_snps, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_sample_size, n_eqtl_snps):
		# Load in gwas data
		self.gwas_beta = gwas_beta
		self.gwas_beta_var = 1.0/N_gwas
		self.var_ldscores = variant_ld_scores

		# Load in eqtl data
		self.eqtl_beta = eqtl_beta
		self.eqtl_ldscore = eqtl_ldscore
		self.eqtl_position = eqtl_position
		self.eqtl_beta_var = 1.0/(eqtl_sample_size)

		# Number of components
		self.KK = len(self.gwas_beta)
		# Number of genes
		self.GG = len(self.eqtl_beta)

		# Get cis snps
		self.n_gwas_snps = n_gwas_snps
		self.n_eqtl_snps = n_eqtl_snps

		return


	def fit(self, total_iterations=15000, burn_in_iterations=10000, update_gwas_resid_var=True, update_eqtl_resid_var=True, v0=2.0, s_sq=0.0, cc=0.0):
		""" Fit the model.
		"""
		# Initialize model params
		self.total_iterations = total_iterations
		self.burn_in_iterations = burn_in_iterations
		genome_delta_sq, gene_h2s = self.fit_gene_models()

		model = sm.OLS(np.square(self.gwas_beta) - (self.gwas_beta_var), np.transpose(np.vstack((genome_delta_sq,self.var_ldscores)))).fit()

		self.est_nm_h2 = model.params[1]*np.sum(self.var_ldscores)
		self.est_med_h2 = model.params[0]*self.GG*np.mean(gene_h2s)
		self.eqtl_h2 = np.mean(gene_h2s)

		return



	def compute_overlapping_med_h2(self):
		genome_delta_alpha = np.zeros(self.KK)
		eqtl_alt = []
		genome_delta_alpha_sq = np.zeros(self.KK)
		#genome_delta_alpha2 = np.zeros(self.KK)
		for gg in range(self.GG):
			genome_delta_alpha[self.eqtl_position[gg]] = genome_delta_alpha[self.eqtl_position[gg]] + (self.deltas[gg]*self.alpha[gg])
			eqtl_alt.append(np.sum(np.square(self.deltas[gg])))

			if gg not in self.eqtl_tracker:
				self.eqtl_tracker[gg] = []
			self.eqtl_tracker[gg].append(self.deltas[gg])
			#genome_delta_alpha_sq[self.eqtl_position[gg]] = genome_delta_alpha_sq[self.eqtl_position[gg]] + np.square(self.deltas[gg])

		if self.itera == 15000:
			for gg in range(self.GG):
				eqtl_distr = np.asarray(self.eqtl_tracker[gg])
				eqtl_mean = np.mean(eqtl_distr,axis=0)
				eqtl_var = np.var(eqtl_distr,axis=0)
				genome_delta_alpha_sq[self.eqtl_position[gg]] = genome_delta_alpha_sq[self.eqtl_position[gg]] + np.square(eqtl_mean) + eqtl_var

			pdb.set_trace()
			model = sm.OLS(np.square(self.gwas_beta) - (self.gwas_beta_resid), np.transpose(np.vstack((genome_delta_alpha_sq,self.var_ldscores)))).fit()


		med_alt = np.sum(np.square(genome_delta_alpha))

		total = np.sum(np.square(genome_delta_alpha + self.gamma))

		nm_alt = np.sum(np.square(self.gamma))
		
		return med_alt, nm_alt, total, np.asarray(eqtl_alt)


	

	def fit_gene_models(self, v0=0.0, s_sq=0.0, cc=0.0):
		genome_delta_sq = np.zeros(self.KK)
		gene_h2s = []
		for gg in range(self.GG):
			# precompute quantities for gene
			gamma_var = .1/np.sum(self.eqtl_ldscore[gg])
			gamma = np.zeros(len(self.eqtl_ldscore[gg]))
			eqtl_resid_var = 1.0
			weights = np.ones(len(self.eqtl_ldscore[gg]))
			sampled_gammas = []

			# Iterative algorithm
			for itera in range(self.total_iterations):
				# Update gamma
				marginal_gamma_vars = gamma_var*self.eqtl_ldscore[gg]
				posterior_vars = 1.0/((1.0/(eqtl_resid_var*self.eqtl_beta_var)) + (1.0/marginal_gamma_vars))
				posterior_means = (self.eqtl_beta[gg]/(eqtl_resid_var*self.eqtl_beta_var))*posterior_vars
				gamma = np.random.normal(loc=posterior_means, scale=np.sqrt(posterior_vars))

				# Update gamma var
				vv = np.sum(weights) + v0
				tau_sq = np.sum(weights*np.square(gamma)/self.eqtl_ldscore[gg]) + s_sq
				gamma_var = (1.0/np.random.gamma(shape=(vv/2) + cc, scale=1.0/((tau_sq/2) + cc),size=1))[0]

				if itera > self.burn_in_iterations:
					sampled_gammas.append(gamma)

			sampled_gammas = np.asarray(sampled_gammas)

			E_gamma_sq = np.square(np.mean(sampled_gammas,axis=0)) + np.var(sampled_gammas,axis=0)
			
			genome_delta_sq[self.eqtl_position[gg]] = genome_delta_sq[self.eqtl_position[gg]] + E_gamma_sq

			gene_h2s.append(np.mean(np.sum(np.square(sampled_gammas),axis=1)))

		return genome_delta_sq, gene_h2s
				


	def update_gamma(self):
		# Re-include current effects
		self.gwas_beta_resid = self.gwas_beta_resid + self.gamma

		# Compute posterior distribution
		marginal_gamma_vars = self.gamma_var*self.var_ldscores

		posterior_vars = 1.0/((1.0/(self.gwas_resid_var*self.gwas_beta_var)) + (1.0/marginal_gamma_vars))
		posterior_means = (self.gwas_beta_resid/(self.gwas_resid_var*self.gwas_beta_var))*posterior_vars

		# Sample from posterior distribution
		self.gamma = np.random.normal(loc=posterior_means, scale=np.sqrt(posterior_vars))

		# Remove current effects
		self.gwas_beta_resid = self.gwas_beta_resid - self.gamma

		return







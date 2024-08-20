import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm










class Bayesian_LMM_SS_h2_med_inference(object):
	def __init__(self, gwas_beta, gwas_beta_se, variant_ld_scores, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, cis_snps):
		# Load in gwas data
		self.gwas_beta = gwas_beta
		self.gwas_beta_var = np.square(gwas_beta_se)
		self.var_ldscores = variant_ld_scores

		# Load in eqtl data
		self.eqtl_beta = eqtl_beta
		self.eqtl_ldscore = eqtl_ldscore
		self.eqtl_position = eqtl_position
		self.eqtl_beta_var = []
		for beta_se_tmp in eqtl_beta_se:
			self.eqtl_beta_var.append(np.square(beta_se_tmp))

		# Number of components
		self.KK = len(self.gwas_beta)
		# Number of genes
		self.GG = len(self.eqtl_beta)

		# Get cis snps
		self.cis_snps = cis_snps

		avg_eqtl_ld_scores = []
		for gg in range(self.GG):
			avg_eqtl_ld_scores.append(np.mean(self.eqtl_ldscore[gg]))
		self.avg_eqtl_ld_scores = np.asarray(avg_eqtl_ld_scores)


	def fit(self, total_iterations=15000, burn_in_iterations=10000, gamma_var_update_version='ld_score_weighting', update_resid_var=False):
		""" Fit the model.
		"""
		# Initialize model params
		self.initialize_variables()

		# Keep track of iterations
		self.itera = 0

		# Iterative Gibbs sampling algorithm
		for itera in range(total_iterations):
			# Update alpha
			self.update_alpha()

			# Update alpha
			self.update_deltas()

			# Update gamma
			self.update_gamma()
	
			# Update gamma_var
			self.update_gamma_var(version=gamma_var_update_version)

			# Update alpha var
			self.update_alpha_var()

			# Update delta var
			self.update_delta_vars(version=gamma_var_update_version)

			if update_resid_var:
				self.update_resid_var()

			# Update iteration number
			self.itera = self.itera + 1
	

			if np.mod(itera, 50) == 0.0:
				print(itera)

			if itera > burn_in_iterations:
				self.sampled_gamma_vars.append(self.gamma_var)

				tmp_eqtl_h2 = []
				for gene_iter in range(self.GG):
					tmp_eqtl_h2.append(self.delta_vars[gene_iter]*np.sum(self.cis_snps[gene_iter]))
				tmp_eqtl_h2 = np.asarray(tmp_eqtl_h2)
				self.sampled_eqtl_h2s.append(tmp_eqtl_h2)
				self.sampled_alpha_vars.append(self.alpha_var)
				self.sampled_med_h2.append(np.sum(tmp_eqtl_h2*self.alpha_var))

				if np.mod(itera, 10) == 0.0:
					med_h2 = np.sum(tmp_eqtl_h2*self.alpha_var)
					nm_h2 = self.gamma_var*self.KK
					eqtl_h2 = np.mean(tmp_eqtl_h2)

					print('med: ' + str(med_h2))
					print('nm: ' + str(nm_h2))
					print('eqtl: ' + str(eqtl_h2))
					print('resid_var: ' + str(self.resid_var))

		self.sampled_gamma_vars = np.asarray(self.sampled_gamma_vars)
		self.sampled_nm_h2 = self.KK*self.sampled_gamma_vars
		self.sampled_eqtl_h2s = np.asarray(self.sampled_eqtl_h2s)
		self.sampled_alpha_vars = np.asarray(self.sampled_alpha_vars)
		self.sampled_med_h2 = np.asarray(self.sampled_med_h2)

		return

	def update_resid_var(self):
		self.resid_var = np.sum(np.square(self.gwas_beta_resid)/self.gwas_beta_var)/len(self.gwas_beta_resid)
		return

	def update_delta_vars(self, version='ld_score_weighting', v0=2.0, s_sq=0.0):
		# Loop through genes
		for gg in np.random.permutation(range(self.GG)):

			'''
			if version == 'ld_score_weighting':
				weights = 1.0/self.eqtl_ldscore[gg]
			elif version == 'equal_weights':
				weights = np.ones(len(self.eqtl_position[gg]))
			else:
				print('asumption error: weighting not yet implemented')
				pdb.set_trace()
			vv = np.sum(weights) + v0
			tau_sq = np.sum(weights*np.square(self.deltas[gg])/self.eqtl_ldscore[gg]) + s_sq			

			tmp_indices = self.eqtl_ldscore[gg] > 1.0

			vv = np.sum(tmp_indices) + v0
			tau_sq = np.sum(np.square(self.deltas[gg][tmp_indices])/self.eqtl_ldscore[gg][tmp_indices]) + s_sq
			'''
			#tmp_indices = self.cis_snps[gg] == 1
			tmp_indices = self.cis_snps[gg] < 10

			if version == 'equal_weights':
				weights = np.ones(np.sum(tmp_indices))
			elif version == 'ld_score_weighting':
				weights = 1.0/self.eqtl_ldscore[gg][tmp_indices]
			else:
				print('asumption error: weighting not yet implemented')
				pdb.set_trace()			

			vv = np.sum(weights) + v0
			tau_sq = np.sum(weights*np.square(self.deltas[gg][tmp_indices])/self.eqtl_ldscore[gg][tmp_indices]) + s_sq

			# Initialize inverse gamma distribution
			#invgamma_dist = invgamma(vv/2, scale=tau_sq/2)
			# Sample from it
			#self.delta_vars[gg] = invgamma_dist.rvs(size=1)[0]
			self.delta_vars[gg] = np.sum(weights*np.square(self.deltas[gg][tmp_indices])/self.eqtl_ldscore[gg][tmp_indices])/np.sum(weights)

		return	

	def update_deltas(self):
		# Loop through genes
		for gg in np.random.permutation(range(self.GG)):

			# Re include effects of current gene
			cis_gwas_beta = self.gwas_beta_resid[self.eqtl_position[gg]] + (self.deltas[gg]*self.alpha[gg])
			cis_gwas_beta_var = self.gwas_beta_var[self.eqtl_position[gg]]

			# Compute posterior distribution
			marginal_delta_vars = self.delta_vars[gg]*self.eqtl_ldscore[gg]

			#gwas_gene_weights = 1.0/np.sqrt(self.var_ldscores[self.eqtl_position[gg]])
			#eqtl_gene_weights = 1.0/np.sqrt(self.eqtl_ldscore[gg])

			posterior_vars = 1.0/((np.square(self.alpha[gg])/(self.resid_var*cis_gwas_beta_var)) + (1.0/self.eqtl_beta_var[gg]) + (1.0/marginal_delta_vars))
			posterior_means = ((cis_gwas_beta*self.alpha[gg]/(self.resid_var*cis_gwas_beta_var)) + (self.eqtl_beta[gg]/self.eqtl_beta_var[gg]))*posterior_vars

			# Sample from posterior
			self.deltas[gg] = np.random.normal(loc=posterior_means, scale=np.sqrt(posterior_vars))

			# Remove updated effects of this gene
			self.gwas_beta_resid[self.eqtl_position[gg]] = cis_gwas_beta - (self.deltas[gg]*self.alpha[gg])

		return

	def update_alpha(self):
		# Loop through genes
		for gg in np.random.permutation(range(self.GG)):
			# Re include effects of current gene
			tmp_indices = self.cis_snps[gg] == 1

			cis_gwas_beta = self.gwas_beta_resid[self.eqtl_position[gg]] + self.deltas[gg]*self.alpha[gg]
			cis_gwas_beta_var = self.gwas_beta_var[self.eqtl_position[gg]]

			# Compute posterior distribution
			# Consider weighting by gene ld scores here
			#weight = 1.0/np.sqrt(self.avg_eqtl_ld_scores[gg])
			posterior_var = 1.0/(np.sum(np.square(self.deltas[gg])/(self.resid_var*cis_gwas_beta_var)) + (1.0/self.alpha_var))
			posterior_mean = np.sum(cis_gwas_beta*self.deltas[gg]/(self.resid_var*cis_gwas_beta_var))*posterior_var

			# Sample
			self.alpha[gg] = np.random.normal(loc=posterior_mean, scale=np.sqrt(posterior_var))

			# Remove updated effects of this gene
			self.gwas_beta_resid[self.eqtl_position[gg]] = cis_gwas_beta - (self.deltas[gg]*self.alpha[gg])
		
		return


	def update_gamma_var(self, version='ld_score_weighting', v0=2.0, s_sq=0.0):
		if version == 'ld_score_weighting':
			weights = 1.0/self.var_ldscores
			vv = np.sum(weights) + v0
			tau_sq = np.sum(weights*np.square(self.gamma)/self.var_ldscores) + s_sq
		elif version == 'equal_weights':
			weights = np.ones(self.KK)
			vv = np.sum(weights) + v0
			tau_sq = np.sum(weights*np.square(self.gamma)/self.var_ldscores) + s_sq			
		else:
			print('asssumption erorr: version not implimented yet')
			pdb.set_trace()

		# Initialize inverse gamma distribution
		#invgamma_dist = invgamma(vv/2, scale=tau_sq/2)
		# Sample from it
		#self.gamma_var = invgamma_dist.rvs(size=1)[0]
		self.gamma_var = np.sum(weights*np.square(self.gamma)/self.var_ldscores)/np.sum(weights)

		return

	def update_alpha_var(self, v0=2.0, s_sq=0.0):
		weights = np.ones(self.GG)
		#weights = 1.0/self.avg_eqtl_ld_scores
		vv = np.sum(weights) + v0
		
		tau_sq = np.sum(weights*np.square(self.alpha)) + s_sq

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma(vv/2, scale=tau_sq/2)
		# Sample from it
		#self.alpha_var = invgamma_dist.rvs(size=1)[0]
		self.alpha_var = np.sum(weights*np.square(self.alpha))/np.sum(weights)

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
		self.gamma = np.zeros(self.KK)  # Non-mediated variant (marginal) effects
		self.alpha = np.zeros(self.GG)  # Gene-trait effects
		self.deltas = []  # Variant to gene (marginal) effects for each gene
		for gene_iter in range(self.GG):
			self.deltas.append(self.eqtl_beta[gene_iter])

		# Initialize variance parameters
		self.gamma_var = 1e-6
		self.alpha_var = 1e-6
		self.delta_vars = np.ones(self.GG)*1e-3

		self.resid_var = 1.0

		# Remove causal effects from gwas beta
		self.gwas_beta_resid = np.copy(self.gwas_beta) - self.gamma  # Remove non-mediated variant effects
		# Remove mediated effects of each gene
		for gene_iter in range(self.GG):
			self.gwas_beta_resid[self.eqtl_position[gene_iter]] = self.gwas_beta_resid[self.eqtl_position[gene_iter]] - (self.alpha[gene_iter]*self.deltas[gene_iter])

		# Keep track of sampled gamma_vars
		self.sampled_gamma_vars = []
		self.sampled_alpha_vars = []
		self.sampled_delta_vars = []
		self.sampled_eqtl_h2s = []
		self.sampled_med_h2 = []



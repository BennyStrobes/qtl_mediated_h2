import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm










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


	def fit(self, total_iterations=15000, burn_in_iterations=10000, update_resid_var=True, v0=2.0, s_sq=0.0, cc=0.0, gene_gamma_shape=1.0):
		""" Fit the model.
		"""
		# Initialize model params
		self.gene_gamma_shape = gene_gamma_shape
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
			self.update_gamma_var(v0=v0, s_sq=s_sq, cc=cc)

			# Update alpha var
			self.update_alpha_var(v0=v0, s_sq=s_sq, cc=cc)

			# Update delta var
			self.update_delta_vars(v0=v0, s_sq=s_sq, cc=cc)

			# Update gene gamma shape
			self.update_gene_gamma_scale_param()

			if update_resid_var:
				if itera >= 20:
					self.update_gwas_resid_var(v0=v0, s_sq=s_sq, cc=cc)
					self.update_eqtl_resid_var(v0=v0, s_sq=s_sq, cc=cc)

			# Update iteration number
			self.itera = self.itera + 1
	


			if itera > burn_in_iterations:

				nm_h2 = self.gamma_var*self.n_gwas_snps

				tmp_eqtl_h2 = []
				for gene_iter in range(self.GG):
					tmp_eqtl_h2.append(self.delta_vars[gene_iter]*self.n_eqtl_snps[gene_iter])
				eqtl_h2s = np.asarray(tmp_eqtl_h2)
				eqtl_h2 = np.mean(eqtl_h2s)
				med_h2 = np.sum(eqtl_h2s*self.alpha_var)
				alt_med_h2, alt2_med_h2, total_h2, alt_nm_h2 = self.compute_overlapping_med_h2()

				self.sampled_iters.append(itera)
				self.sampled_nm_h2_v1.append(nm_h2)
				self.sampled_nm_h2_v2.append(alt_nm_h2)
				self.sampled_med_h2_v1.append(med_h2)
				self.sampled_med_h2_v2.append(alt_med_h2)
				self.sampled_med_h2_v3.append(alt2_med_h2)
				self.sampled_total_h2.append(total_h2)
				self.sampled_eqtl_h2s.append(eqtl_h2)


				if np.mod(itera, 50) == 0.0:

					print('ITERA ' + str(itera))
					print('med: ' + str(med_h2))
					print('alt_med: ' + str(alt_med_h2))
					print('alt_med2: ' + str(alt2_med_h2))
					print('nm: ' + str(nm_h2))
					print('alt_nm: ' + str(alt_nm_h2))
					print('total: ' + str(total_h2))
					print('eqtl: ' + str(eqtl_h2))
					print('resid_var: ' + str(self.gwas_resid_var))
					print('eqtl_resid_var: ' + str(np.mean(self.eqtl_resid_vars)))
					print(print(self.gene_gamma_scale/(self.gene_gamma_shape-1.0)))
					print(np.sort(eqtl_h2s))


		self.sampled_iters = np.asarray(self.sampled_iters)
		self.sampled_nm_h2_v1 = np.asarray(self.sampled_nm_h2_v1)
		self.sampled_nm_h2_v2 = np.asarray(self.sampled_nm_h2_v2)
		self.sampled_med_h2_v1 = np.asarray(self.sampled_med_h2_v1)
		self.sampled_med_h2_v2 = np.asarray(self.sampled_med_h2_v2)
		self.sampled_med_h2_v3 = np.asarray(self.sampled_med_h2_v3)
		self.sampled_total_h2 = np.asarray(self.sampled_total_h2)
		self.sampled_eqtl_h2s = np.asarray(self.sampled_eqtl_h2s)
		
		return

	def update_gwas_resid_var(self, v0=0.0, s_sq=0.0, cc=1e-6):
		vv = len(self.gwas_beta_resid) + v0
		tau_sq = np.sum(np.square(self.gwas_beta_resid)/self.gwas_beta_var) + s_sq

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
		# Sample from it
		self.gwas_resid_var = invgamma_dist.rvs(size=1)[0]

		return

	def update_delta_vars(self, v0=0.0, s_sq=0.0, cc=1e-6):
		# Loop through genes
		for gg in np.random.permutation(range(self.GG)):


			vv = len(self.deltas[gg]) + v0
			tau_sq = np.sum(np.square(self.deltas[gg])/self.eqtl_ldscore[gg]) + s_sq


			#indices = self.eqtl_ldscore[gg] > .1		
			#vv = np.sum(indices)
			#tau_sq = np.sum(np.square(self.deltas[gg][indices])/self.eqtl_ldscore[gg][indices]) + s_sq

			# Initialize inverse gamma distribution
			# Shape (alpha) stays same
			# Scale (beta) gets update cross genes
			invgamma_dist = invgamma((vv/2) + self.gene_gamma_shape, scale=(tau_sq/2) + self.gene_gamma_scale)
			# Sample from it
			self.delta_vars[gg] = invgamma_dist.rvs(size=1)[0]
			#self.delta_vars[gg] = np.sum(weights*np.square(self.deltas[gg][tmp_indices])/self.eqtl_ldscore[gg][tmp_indices])/np.sum(weights)
		return

	def update_gene_gamma_scale_param(self, beta0=0.0, alpha0=0.0):

		gamma_b =beta0 + np.sum(1.0/self.delta_vars)
		gamma_a = alpha0 + len(self.delta_vars)*self.gene_gamma_shape


		self.gene_gamma_scale = np.random.gamma(gamma_a, scale=1/gamma_b,size=1)[0]

		return

	def update_eqtl_resid_var(self, v0=0.0, s_sq=0.0, cc=1e-6):
		for gg in range(self.GG):
			resid_vec = self.eqtl_beta[gg] - self.deltas[gg]

			vv = len(resid_vec) + v0
			tau_sq = np.sum(np.square(resid_vec)/self.eqtl_beta_var) + s_sq

			# Initialize inverse gamma distribution
			invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
			# Sample from it
			self.eqtl_resid_vars[gg] = invgamma_dist.rvs(size=1)[0]
		
		return

	def compute_overlapping_med_h2(self):
		genome_delta_alpha = np.zeros(self.KK)
		eqtl_ld_scores = np.zeros(self.KK)
		v1 = 0
		total_eqtl_snps = 0
		for gg in range(self.GG):
			genome_delta_alpha[self.eqtl_position[gg]] = genome_delta_alpha[self.eqtl_position[gg]] + (self.deltas[gg]*self.alpha[gg])
			v1 = v1 + np.sum(np.square((self.deltas[gg]*self.alpha[gg])))
			eqtl_ld_scores[self.eqtl_position[gg]] = eqtl_ld_scores[self.eqtl_position[gg]] + self.eqtl_ldscore[gg]
			total_eqtl_snps = total_eqtl_snps + self.n_eqtl_snps[gg]


		v2 = np.sum(np.square(genome_delta_alpha))
		weighted = np.square(genome_delta_alpha)/eqtl_ld_scores
		v3 = np.mean(weighted[np.isnan(weighted) == False])*total_eqtl_snps

		total = np.sum(np.square(genome_delta_alpha + self.gamma))

		nm_alt = np.sum(np.square(self.gamma))
		
		return v2, v3, total, nm_alt

	def update_deltas(self):
		# Loop through genes
		for gg in np.random.permutation(range(self.GG)):

			# Re include effects of current gene
			cis_gwas_beta = self.gwas_beta_resid[self.eqtl_position[gg]] + (self.deltas[gg]*self.alpha[gg])
			cis_gwas_beta_var = self.gwas_beta_var

			# Compute posterior distribution
			marginal_delta_vars = self.delta_vars[gg]*self.eqtl_ldscore[gg]

			posterior_vars = 1.0/((np.square(self.alpha[gg])/(self.gwas_resid_var*cis_gwas_beta_var)) + (1.0/(self.eqtl_beta_var*self.eqtl_resid_vars[gg])) + (1.0/marginal_delta_vars))
			posterior_means = ((cis_gwas_beta*self.alpha[gg]/(self.gwas_resid_var*cis_gwas_beta_var)) + (self.eqtl_beta[gg]/(self.eqtl_beta_var*self.eqtl_resid_vars[gg])))*posterior_vars

			# Sample from posterior
			self.deltas[gg] = np.random.normal(loc=posterior_means, scale=np.sqrt(posterior_vars))

			# Remove updated effects of this gene
			self.gwas_beta_resid[self.eqtl_position[gg]] = cis_gwas_beta - (self.deltas[gg]*self.alpha[gg])

		return

	def update_alpha(self):
		# Loop through genes
		for gg in np.random.permutation(range(self.GG)):
			# Re include effects of current gene
			cis_gwas_beta = self.gwas_beta_resid[self.eqtl_position[gg]] + self.deltas[gg]*self.alpha[gg]
			cis_gwas_beta_var = self.gwas_beta_var

			# Compute posterior distribution
			# Consider weighting by gene ld scores here
			#weight = 1.0/np.sqrt(self.avg_eqtl_ld_scores[gg])
			posterior_var = 1.0/(np.sum(np.square(self.deltas[gg])/(self.gwas_resid_var*cis_gwas_beta_var)) + (1.0/self.alpha_var))
			posterior_mean = np.sum(cis_gwas_beta*self.deltas[gg]/(self.gwas_resid_var*cis_gwas_beta_var))*posterior_var

			# Sample
			self.alpha[gg] = np.random.normal(loc=posterior_mean, scale=np.sqrt(posterior_var))

			# Remove updated effects of this gene
			self.gwas_beta_resid[self.eqtl_position[gg]] = cis_gwas_beta - (self.deltas[gg]*self.alpha[gg])
		
		return


	def update_gamma_var(self, version='ld_score_weighting', v0=0.0, s_sq=0.0, cc=1e-6):
		weights = np.ones(self.KK)
		vv = np.sum(weights) + v0
		tau_sq = np.sum(weights*np.square(self.gamma)/self.var_ldscores) + s_sq

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
		# Sample from it
		self.gamma_var = invgamma_dist.rvs(size=1)[0]

		return

	def update_alpha_var(self, v0=0.0, s_sq=0.0, cc=1e-6):
		weights = np.ones(self.GG)
		#weights = 1.0/self.avg_eqtl_ld_scores
		vv = np.sum(weights) + v0
		
		tau_sq = np.sum(weights*np.square(self.alpha)) + s_sq

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2)+cc)
		# Sample from it
		self.alpha_var = invgamma_dist.rvs(size=1)[0]
		#self.alpha_var = np.sum(weights*np.square(self.alpha))/np.sum(weights)

		return

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



	def initialize_variables(self):
		# Initialize causal effcts
		self.gamma = np.zeros(self.KK)  # Non-mediated variant (marginal) effects
		self.alpha = np.zeros(self.GG)  # Gene-trait effects
		self.deltas = []  # Variant to gene (marginal) effects for each gene
		for gene_iter in range(self.GG):
			self.deltas.append(1.0*self.eqtl_beta[gene_iter])

		# Initialize variance parameters
		self.gamma_var = 1e-6
		self.alpha_var = 1e-3
		self.delta_vars = np.ones(self.GG)*1e-3

		self.gene_gamma_scale = 2.0
		self.gwas_resid_var = 1.0
		self.eqtl_resid_vars = []
		for gene_iter in range(self.GG):
			self.eqtl_resid_vars.append(1.0)

		# Remove causal effects from gwas beta
		self.gwas_beta_resid = np.copy(self.gwas_beta) - self.gamma  # Remove non-mediated variant effects
		# Remove mediated effects of each gene
		for gene_iter in range(self.GG):
			self.gwas_beta_resid[self.eqtl_position[gene_iter]] = self.gwas_beta_resid[self.eqtl_position[gene_iter]] - (self.alpha[gene_iter]*self.deltas[gene_iter])

		# Keep track of sampled gamma_vars
		self.sampled_nm_h2_v1 = []
		self.sampled_nm_h2_v2 = []
		self.sampled_med_h2_v1 = []
		self.sampled_med_h2_v2 = []
		self.sampled_med_h2_v3 = []
		self.sampled_total_h2 = []
		self.sampled_eqtl_h2s = []
		self.sampled_iters = []

		return



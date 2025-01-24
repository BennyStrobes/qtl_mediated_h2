import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import time







class Bayesian_LMM_SS_h2_med_inference(object):
	def __init__(self, gwas_beta, variant_ld_scores, N_gwas, n_gwas_snps, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_sample_size, n_eqtl_snps, gene_classes):
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

		# Number of gene classes
		self.gene_classes = gene_classes
		self.unique_gene_classes = np.unique(self.gene_classes)

		# Get cis snps
		self.n_gwas_snps = n_gwas_snps
		self.n_eqtl_snps = n_eqtl_snps

		return


	def fit(self, total_iterations=15000, burn_in_iterations=10000, update_gwas_resid_var=True, update_eqtl_resid_var=True, v0=2.0, s_sq=0.0, cc=0.0):
		""" Fit the model.
		"""
		# Initialize model params
		self.initialize_variables()

		# Keep track of iterations
		self.itera = 0

		# Iterative Gibbs sampling algorithm
		for itera in range(total_iterations):
			# Update alpha
			t1 = time.time()
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
			self.update_delta_vars_fast(v0=v0, s_sq=s_sq, cc=cc)

			# Update residual variances
			if update_gwas_resid_var and itera >= 20:
				self.update_gwas_resid_var(v0=v0, s_sq=s_sq, cc=cc)
			if update_eqtl_resid_var and itera >= 20:
				self.update_eqtl_resid_var_fast(v0=v0, s_sq=s_sq, cc=cc)


			# Update iteration number
			self.itera = self.itera + 1
	


			if itera > burn_in_iterations:


				alt_med_h2, agg_alt_med_h2, alt_nm_h2, total_h2, alt_eqtl_h2s = self.compute_overlapping_med_h2()


				self.sampled_iters.append(itera)
				self.sampled_nm_h2_v2.append(alt_nm_h2)
				self.sampled_med_h2_v2.append(alt_med_h2)
				self.sampled_agg_med_h2_v2.append(agg_alt_med_h2)
				self.sampled_total_h2.append(total_h2)
				self.sampled_eqtl_h2s_v2.append(np.mean(alt_eqtl_h2s))


				if np.mod(itera, 500) == 0.0:

					print('#######################')
					print('ITERA ' + str(itera))
					print('Mediated')
					print(alt_med_h2)
					print(agg_alt_med_h2)
					print('Non-mediated')
					print(alt_nm_h2)
					print('Total')
					print(total_h2)
					print('eqtl h2')
					print(np.sort(alt_eqtl_h2s))
			t2 = time.time()
			if itera > burn_in_iterations and np.mod(itera, 500) == 0.0:
				print(t2-t1)



		self.sampled_iters = np.asarray(self.sampled_iters)
		self.sampled_nm_h2_v2 = np.asarray(self.sampled_nm_h2_v2)
		self.sampled_med_h2_v2 = np.asarray(self.sampled_med_h2_v2)
		self.sampled_total_h2 = np.asarray(self.sampled_total_h2)
		self.sampled_eqtl_h2s_v2 = np.asarray(self.sampled_eqtl_h2s_v2)
		self.sampled_agg_med_h2_v2 = np.asarray(self.sampled_agg_med_h2_v2)

		return

	def update_gwas_resid_var(self, v0=0.0, s_sq=0.0, cc=1e-6):
		vv = len(self.gwas_beta_resid) + v0
		tau_sq = np.sum(np.square(self.gwas_beta_resid)/self.gwas_beta_var) + s_sq

		# Initialize inverse gamma distribution
		#invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
		# Sample from it
		#self.gwas_resid_var = invgamma_dist.rvs(size=1)[0]

		self.gwas_resid_var = (1.0/np.random.gamma(shape=(vv/2) + cc, scale=1.0/((tau_sq/2) + cc),size=1))[0]

		return

	def update_delta_vars(self, v0=0.0, s_sq=0.0, cc=1e-6, weighted=False):
		# Loop through genes
		for gg in np.random.permutation(range(self.GG)):

			if weighted == False:
				vv = len(self.deltas[gg]) + v0
				tau_sq = np.sum(np.square(self.deltas[gg])/self.eqtl_ldscore[gg]) + s_sq
			else:
				vv = np.sum(self.eqtl_ldscore[gg]) + v0
				tau_sq = np.sum(np.square(self.deltas[gg])) + s_sq



			#indices = self.eqtl_ldscore[gg] > .1		
			#vv = np.sum(indices)
			#tau_sq = np.sum(np.square(self.deltas[gg][indices])/self.eqtl_ldscore[gg][indices]) + s_sq

			# Initialize inverse gamma distribution
			#invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2)+cc)
			# Sample from it
			#self.delta_vars[gg] = invgamma_dist.rvs(size=1)[0]
			
			# Sample from inverse gamma using numpy
			self.delta_vars[gg] = 1.0/np.random.gamma(shape=(vv/2) + cc, scale=1.0/((tau_sq/2) + cc),size=1)
		return

	def update_delta_vars_fast(self, v0=0.0, s_sq=0.0, cc=1e-6, weighted=False):
		# Loop through genes
		param1 = []
		param2 = []
		for gg in range(self.GG):

			if weighted == False:
				vv = len(self.deltas[gg]) + v0
				tau_sq = np.sum(np.square(self.deltas[gg])/self.eqtl_ldscore[gg]) + s_sq
			else:
				vv = np.sum(self.eqtl_ldscore[gg]) + v0
				tau_sq = np.sum(np.square(self.deltas[gg])) + s_sq

			param1.append((vv/2) + cc)
			param2.append(1.0/((tau_sq/2) + cc))

			#indices = self.eqtl_ldscore[gg] > .1		
			#vv = np.sum(indices)
			#tau_sq = np.sum(np.square(self.deltas[gg][indices])/self.eqtl_ldscore[gg][indices]) + s_sq

			# Initialize inverse gamma distribution
			#invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2)+cc)
			# Sample from it
			#self.delta_vars[gg] = invgamma_dist.rvs(size=1)[0]
			
			# Sample from inverse gamma using numpy
			#self.delta_vars[gg] = 1.0/np.random.gamma(shape=(vv/2) + cc, scale=1.0/((tau_sq/2) + cc),size=1)
		param1 = np.asarray(param1)
		param2 = np.asarray(param2)

		self.delta_vars = 1.0/np.random.gamma(shape=param1, scale=param2)
		return


	def update_eqtl_resid_var(self, v0=0.0, s_sq=0.0, cc=1e-6):
		for gg in range(self.GG):
			resid_vec = self.eqtl_beta[gg] - self.deltas[gg]

			vv = len(resid_vec) + v0
			tau_sq = np.sum(np.square(resid_vec)/self.eqtl_beta_var) + s_sq

			# Initialize inverse gamma distribution
			#invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
			# Sample from it
			#self.eqtl_resid_vars[gg] = invgamma_dist.rvs(size=1)[0]
			
			# Sample from inverse gamma using numpy
			self.eqtl_resid_vars[gg] = 1.0/np.random.gamma(shape=(vv/2) + cc, scale=1.0/((tau_sq/2) + cc),size=1)
		return

	def update_eqtl_resid_var_fast(self, v0=0.0, s_sq=0.0, cc=1e-6):
		param1 = []
		param2 = []
		for gg in range(self.GG):
			resid_vec = self.eqtl_beta[gg] - self.deltas[gg]

			vv = len(resid_vec) + v0
			tau_sq = np.sum(np.square(resid_vec)/self.eqtl_beta_var) + s_sq

			param1.append((vv/2) + cc)
			param2.append(1.0/((tau_sq/2) + cc))

			# Initialize inverse gamma distribution
			#invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
			# Sample from it
			#self.eqtl_resid_vars[gg] = invgamma_dist.rvs(size=1)[0]
			
			# Sample from inverse gamma using numpy
		self.eqtl_resid_vars = 1.0/np.random.gamma(shape=param1, scale=param2)
		return


	def compute_overlapping_med_h2(self):
		genome_delta_alpha = np.zeros((len(self.unique_gene_classes), self.KK))
		eqtl_alt = []
		#genome_delta_alpha2 = np.zeros(self.KK)
		for gg in range(self.GG):
			gene_class_assignment = self.gene_class_assignments[gg]

			genome_delta_alpha[gene_class_assignment, self.eqtl_position[gg]] = genome_delta_alpha[gene_class_assignment, self.eqtl_position[gg]] + (self.deltas[gg]*self.alpha[gg])
			eqtl_alt.append(np.sum(np.square(self.deltas[gg])))


		med_alt = np.sum(np.square(genome_delta_alpha),axis=1)


		genome_delta_alpha_agg = np.sum(genome_delta_alpha,axis=0)

		agg_med_alt = np.sum(np.square(genome_delta_alpha_agg))


		total = np.sum(np.square(genome_delta_alpha_agg + self.gamma))

		nm_alt = np.sum(np.square(self.gamma))
		
		return med_alt, agg_med_alt, nm_alt, total, np.asarray(eqtl_alt)

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
			posterior_var = 1.0/(np.sum(np.square(self.deltas[gg])/(self.gwas_resid_var*cis_gwas_beta_var)) + (1.0/self.alpha_vars[self.gene_class_assignments[gg]]))
			posterior_mean = np.sum(cis_gwas_beta*self.deltas[gg]/(self.gwas_resid_var*cis_gwas_beta_var))*posterior_var

			# Sample
			self.alpha[gg] = np.random.normal(loc=posterior_mean, scale=np.sqrt(posterior_var))

			# Remove updated effects of this gene
			self.gwas_beta_resid[self.eqtl_position[gg]] = cis_gwas_beta - (self.deltas[gg]*self.alpha[gg])
		
		return


	def update_gamma_var(self, v0=0.0, s_sq=0.0, cc=1e-6, weighted=False):
		if weighted == False:
			weights = np.ones(self.KK)
		else:
			weights = np.copy(self.var_ldscores)
		vv = np.sum(weights) + v0
		tau_sq = np.sum(weights*np.square(self.gamma)/self.var_ldscores) + s_sq

		# Initialize inverse gamma distribution
		#invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
		# Sample from it
		#self.gamma_var = invgamma_dist.rvs(size=1)[0]

		# Inverse gamma distribution in numpy
		self.gamma_var = (1.0/np.random.gamma(shape=(vv/2) + cc, scale=1.0/((tau_sq/2) + cc),size=1))[0]

		return

	def update_alpha_var(self, v0=0.0, s_sq=0.0, cc=1e-6):
		vv = np.zeros(len(self.unique_gene_classes))
		tau_sq = np.zeros(len(self.unique_gene_classes))
		
		for gg in range(self.GG):
			gene_class_assignment = self.gene_class_assignments[gg]
			vv[gene_class_assignment] = vv[gene_class_assignment] + 1
			tau_sq[gene_class_assignment] = tau_sq[gene_class_assignment] + np.square(self.alpha[gg])

		if self.itera < 500:
			agg_vv = np.sum(vv) + v0
			agg_tau_sq = np.sum(tau_sq) + s_sq

			alpha_var = (1.0/np.random.gamma(shape=(agg_vv/2) + cc, scale=1.0/((agg_tau_sq/2) + cc),size=1))[0]
			self.alpha_vars = (self.alpha_vars*0.0) + alpha_var
		else:
			vv = vv + v0
			tau_sq = tau_sq + s_sq

			param1 = (vv/2) + cc
			param2 = 1.0/((tau_sq/2) + cc)

			self.alpha_vars = 1.0/np.random.gamma(shape=param1, scale=param2)

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
		self.alpha_vars = np.ones(len(self.unique_gene_classes))*1e-3
		self.delta_vars = np.ones(self.GG)*1e-3

		print(self.unique_gene_classes)
		# Gene class to alpha var index
		self.gene_class_mapping = {}
		for ii,gene_class in enumerate(self.unique_gene_classes):
			self.gene_class_mapping[gene_class] = ii
		# Gene class assignments
		self.gene_class_assignments = []
		for gene_class in self.gene_classes:
			self.gene_class_assignments.append(self.gene_class_mapping[gene_class])
		self.gene_class_assignments = np.asarray(self.gene_class_assignments)


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
		self.sampled_eqtl_h2s_v1 = []
		self.sampled_eqtl_h2s_v2 = []
		self.sampled_agg_med_h2_v2 = []
		self.sampled_iters = []

		return



import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import time
import bayesian_lmm_pca_ss_h2_single_region


def q_delta_sq_sum_debug(small_Q, delta_mu, delta_mu_var, Q_delta_sq_sum, n_iter=1000000):
	valz = []
	for itera in range(n_iter):
		delta = np.random.normal(loc=delta_mu, scale=np.sqrt(delta_mu_var))
		val = np.sum(np.square(np.dot(small_Q, delta)))
		valz.append(val)
	valz = np.asarray(valz)

	return

class Bayesian_LMM_RSS_med_h2_inference(object):
	def __init__(self, window_info, gene_info, N_gwas, N_eqtl, tmp_output_file):
		self.N_gwas = N_gwas
		self.N_eqtl = N_eqtl
		self.window_info = window_info
		self.gene_info = gene_info

		self.genes = np.sort([*self.gene_info])
		self.GG = len(self.genes)
		self.windows = np.sort([*self.window_info])

		# Number of snps
		self.KK = 0.0
		for window in self.windows:
			self.KK = self.KK + self.window_info[window]['n_snps']


		self.tt = open(tmp_output_file,'w')
		self.tt.write('Iteration\tnm_h2\tmed_h2\talt_med_h2\teqtl_h2\tresid_var\teqtl_resid_var\n')
		

		return



	def fit(self, total_iterations=15000, burn_in_iterations=10000, update_resid_var_bool=True):
		""" Fit the model.
		"""
		# Initialize model params
		self.initialize_variables()

		# Keep track of iterations
		self.itera = 0

		# Iterative Gibbs sampling algorithm
		print('start')
		for itera in range(total_iterations):
			# Loop through windows
			for window_name in self.windows:
				# Update gamma, delta, and alpha in each window seperately
				self.update_gamma_delta_and_alpha_in_single_window(window_name)

			# Update gamma_var
			self.update_gamma_var()

			# Update delta var
			self.update_delta_vars()

			# Update alpha var
			self.update_alpha_var()

			# Update resid var
			if update_resid_var_bool:
				gwas_resid_vars = self.update_gwas_resid_var()
				self.update_eqtl_resid_vars()

			# Update iteration number
			self.itera = self.itera + 1

	

			if itera > burn_in_iterations:
				nm_h2 = self.gamma_var*self.KK
				eqtl_h2s = self.compute_eqtl_h2()
				#med_h2_alt, med_h2_alt2 = self.compute_med_h2_alt()
				med_h2 = np.sum(eqtl_h2s*self.alpha_var)
				avg_eqtl_h2 = np.mean(eqtl_h2s)
				#eqtl_resid_vars = self.get_eqtl_resid_vars()

				#nm_h2_ld_depen, med_h2_ld_depen = self.get_ld_dependent_h2s()


				self.sampled_nm_h2.append(nm_h2)
				self.sampled_med_h2.append(med_h2)
				self.sampled_eqtl_h2.append(avg_eqtl_h2)
				self.sampled_alpha_var.append(self.alpha_var)
				#self.sampled_gwas_resid_var.append(np.mean(gwas_resid_vars))
				#self.sampled_eqtl_resid_var.append(np.mean(eqtl_resid_vars))

				print('NM: ' + str(nm_h2))
				#print('NM2: ' + str(nm_h2_ld_depen))
				print('MED: ' + str(med_h2))
				print('eQTL: ' + str(avg_eqtl_h2))
				#print('MED2: ' + str(med_h2_alt))
				#print('MED3 ' + str(med_h2_alt2))
				#print('MED4 ' + str(med_h2_ld_depen))


				#self.tt.write(str(itera) + '\t' + str(nm_h2) + '\t' + str(nm_h2_ld_depen) + '\t' + str(med_h2) + '\t' + str(med_h2_alt) + '\t' + str(med_h2_ld_depen) + '\t' + str(avg_eqtl_h2) + '\t' + str('NA') + '\n')
				#self.tt.flush()


		self.sampled_nm_h2 = np.asarray(self.sampled_nm_h2)
		self.sampled_med_h2 = np.asarray(self.sampled_med_h2)
		self.sampled_eqtl_h2 = np.asarray(self.sampled_eqtl_h2)
		self.sampled_alpha_var = np.asarray(self.sampled_alpha_var)
		self.sampled_gwas_resid_var = np.asarray(self.sampled_gwas_resid_var)
		self.tt.close()

		return

	def get_ld_dependent_h2s(self):
		med = 0.0
		nm = 0.0
		for window in self.windows:
			med = med + self.window_med_h2[window]
			nm = nm + self.window_nm_h2[window]
		return nm, med

	def check_resid_vars(self):
		for window in self.windows:
			window_Q = np.load(self.window_info[window]['Q_file'])
			window_causal_effects = np.copy(self.gamma[window])
			for gene in self.window_info[window]['genes']:
				cis_indices = self.gene_info[gene]['cis_snps']
				window_causal_effects[cis_indices] = window_causal_effects[cis_indices] + self.deltas[gene]*self.alpha[gene]

				pred_gene = np.dot(window_Q[:, cis_indices], self.deltas[gene])
				pred_resid = self.gene_info[gene]['beta_pc'] - pred_gene

				diff = self.eqtl_beta_resid[gene] - pred_resid
				max_diff = np.max(np.abs(diff))
				if max_diff > 1e-13:
					print('assumption eroror')
					pdb.set_trace()


			pred_resid = self.window_info[window]['beta_pc'] - np.dot(window_Q, window_causal_effects)
			diff = pred_resid - self.gwas_beta_resid[window]
			max_diff = np.max(np.abs(diff))
			if max_diff > 1e-13:
				print('assumption eroror')
				pdb.set_trace()

	def get_eqtl_resid_vars(self):
		resid_vars = []
		for gene in self.genes:
			resid_vars.append(self.eqtl_resid_vars[gene])
		return np.asarray(resid_vars)

	def compute_eqtl_h2(self):
		tmp_eqtl_h2s = []
		for gene in self.genes:
			tmp_eqtl_h2s.append(self.delta_vars[gene]*self.gene_info[gene]['n_cis_snps'])
		return np.asarray(tmp_eqtl_h2s)

	def compute_med_h2_alt(self):
		delta_alphas = []
		for gene in self.genes:
			delta_alphas.append(self.deltas[gene]*self.alpha[gene])
		delta_alphas = np.hstack(delta_alphas)
		sum_delta_alphas = []
		for window in self.windows:
			window_causal_effects = np.zeros(int(self.window_info[window]['n_snps']))
			for gene in self.window_info[window]['genes']:
				cis_indices = self.gene_info[gene]['cis_snps']

				window_causal_effects[cis_indices] = window_causal_effects[cis_indices] + self.deltas[gene]*self.alpha[gene]
			sum_delta_alphas.append(window_causal_effects)

		sum_delta_alphas = np.hstack(sum_delta_alphas)


		return np.sum(np.square(delta_alphas)), np.sum(np.square(sum_delta_alphas))


	def update_alpha_var(self, v0=0.0, s_sq=0.0, cc=0.0):
		# Get vector of alphas (of length number of genes)
		alpha_sq_vec = []
		for gene in self.genes:
			alpha_sq_vec.append(np.square(self.alpha_mu[gene]) + self.alpha_mu_var[gene])
		alpha_sq_vec = np.asarray(alpha_sq_vec)

		vv = len(alpha_sq_vec) + v0
		tau_sq = np.sum(alpha_sq_vec) + s_sq

		self.alpha_var = (.5*tau_sq + cc)/(.5*vv + cc)

		# Initialize inverse gamma distribution
		#invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
		# Sample from it
		#self.alpha_var = invgamma_dist.rvs(size=1)[0]

		return


	def update_delta_vars(self, v0=0.0, s_sq=0.0, cc=0.0):
		# Loop through genes
		for gene in self.genes:
		
			delta_sq = np.square(self.deltas_mu[gene]) + self.deltas_mu_var[gene]
			vv = len(delta_sq) + v0
			tau_sq = np.sum(delta_sq) + s_sq

			self.delta_vars[gene] = (.5*tau_sq + cc)/(.5*vv + cc)

			# Initialize inverse gamma distribution
			#invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
			# Sample from it
			#self.delta_vars[gene] = invgamma_dist.rvs(size=1)[0]
		
		return

	def update_eqtl_resid_vars(self, v0=0.0, s_sq=0.0, cc=1e-6):
		# Loop through genes
		for gene in self.genes:
			gene_resids = self.eqtl_beta_resid[gene]

			vv = len(gene_resids) + v0
			tau_sq = np.sum(np.square(gene_resids)*self.N_eqtl) + s_sq

			# Initialize inverse gamma distribution
			invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
			# Sample from it
			self.eqtl_resid_vars[gene] = invgamma_dist.rvs(size=1)[0]

		return

	def update_gwas_resid_var(self, v0=0.0, s_sq=0.0, cc=1e-6, window_specific_variance=False):
		gwas_resid_vars = []
		if window_specific_variance == False:
			all_resids = []
			for window_name in self.windows:
				window_resid = self.gwas_beta_resid[window_name]
				all_resids.append(window_resid)
			all_resids = np.hstack(all_resids)

			vv = len(all_resids) + v0
			tau_sq = np.sum(np.square(all_resids)*self.N_gwas) + s_sq

			# Initialize inverse gamma distribution
			invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
			# Sample from it
			tmp_gwas_resid_var = invgamma_dist.rvs(size=1)[0]
			for window_name in self.windows:
				self.gwas_resid_var[window_name] = tmp_gwas_resid_var
				gwas_resid_vars.append(tmp_gwas_resid_var)
		elif window_specific_variance:
			for window_name in self.windows:
				window_resid = self.gwas_beta_resid[window_name]
				vv = len(window_resid) + v0
				tau_sq = np.sum(np.square(window_resid)*self.N_gwas) + s_sq
				# Initialize inverse gamma distribution
				invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
				# Sample from it
				tmp_gwas_resid_var = invgamma_dist.rvs(size=1)[0]
				gwas_resid_vars.append(tmp_gwas_resid_var)
				self.gwas_resid_var[window_name] = tmp_gwas_resid_var

		return np.asarray(gwas_resid_vars)


	def update_gamma_var(self, v0=0.0, s_sq=0.0, cc=0.0):
		# First get middle gammas
		all_gammas = []
		for window_name in self.windows:
			window_gammas = np.square(self.gamma_mu[window_name]) + self.gamma_mu_var[window_name]
			all_gammas.append(window_gammas)
		all_gammas = np.hstack(all_gammas)

		vv = len(all_gammas) + v0
		tau_sq = np.sum(all_gammas) + s_sq
	
		self.gamma_var = (.5*tau_sq + cc)/(.5*vv + cc)

		# Initialize inverse gamma distribution
		#invgamma_dist = invgamma((vv/2) + cc, scale=(tau_sq/2) + cc)
		# Sample from it
		#self.gamma_var = invgamma_dist.rvs(size=1)[0]

		return

	def update_gamma_delta_and_alpha_in_single_window(self, window_name):
		# Load in relevent quantities for this window
		window_Q = np.load(self.window_info[window_name]['Q_file'])

		order_bool = np.random.choice([0,1])

		if order_bool == 0:
			# Update gamma and alpha in single window
			self.update_gamma_and_alpha_in_single_window(window_name, window_Q)
			# Update delta in single window
			self.update_delta_in_single_window(window_name, window_Q)
		elif order_bool == 1:
			# Update delta in single window
			self.update_delta_in_single_window(window_name, window_Q)	
			# Update gamma and alpha in single window
			self.update_gamma_and_alpha_in_single_window(window_name, window_Q)		


		window_ld = np.load(self.window_info[window_name]['LD_file'])
		#self.update_window_genetic_var(window_name, window_ld)

		return


	def update_window_genetic_var(self, window_name, LD_mat):

		window_genes = self.window_info[window_name]['genes']

		window_med_causal_effects = np.zeros(self.window_info[window_name]['n_snps'])
		window_nm_causal_effects = np.copy(self.gamma[window_name])

		for gene in window_genes:
			window_cis_indices = self.gene_info[gene]['cis_snps']
			window_med_causal_effects[window_cis_indices] = window_med_causal_effects[window_cis_indices] + (self.deltas[gene]*self.alpha[gene])

		self.window_med_h2[window_name] = np.dot(np.dot(window_med_causal_effects, LD_mat), window_med_causal_effects)
		self.window_nm_h2[window_name] = np.dot(np.dot(window_nm_causal_effects, LD_mat), window_nm_causal_effects)

		return



	def update_delta_in_single_window(self, window_name, window_Q):
		# Get list of genes in this window
		window_genes = self.window_info[window_name]['genes']
		window_gwas_resid_var = self.gwas_resid_var[window_name]

		window_gwas_beta_resid = self.gwas_beta_resid[window_name]


		# Loop through genes 
		for gene in np.random.permutation(window_genes):
			# Update delta corresponding to this gene
			delta_mu = self.deltas_mu[gene]
			delta_mu_var = self.deltas_mu_var[gene]
			window_cis_indices = self.gene_info[gene]['cis_snps']

			# Get Q matrix corresponding to only cis snps of the gene
			gene_Q = window_Q[:, window_cis_indices]


			# Extract eqtl summary stats for gene
			eqtl_gene_beta_resid = self.eqtl_beta_resid[gene]

			# Compute update for single gene
			delta_mu, delta_mu_var, window_gwas_beta_resid, eqtl_gene_beta_resid = self.update_delta_for_single_gene(gene_Q, delta_mu, delta_mu_var, window_gwas_beta_resid, eqtl_gene_beta_resid, self.alpha_mu[gene], self.alpha_mu_var[gene], self.delta_vars[gene], self.eqtl_resid_vars[gene], window_gwas_resid_var)

			self.eqtl_beta_resid[gene] = eqtl_gene_beta_resid
			self.deltas_mu[gene] = delta_mu
			self.deltas_mu_var[gene] = delta_mu_var


		self.gwas_beta_resid[window_name] = window_gwas_beta_resid

		return

	def update_delta_for_single_gene(self, gene_Q, delta_mu, delta_mu_var, window_gwas_beta_resid, eqtl_gene_beta_resid, gene_alpha_mu, gene_alpha_mu_var, gene_delta_var, gene_resid_var, gwas_resid_var):
		# N-snps in the the gene
		gene_K = len(delta_mu)


		# Precompute snp posterior variances
		gene_Q_sum_sq = np.sum(np.square(gene_Q),axis=0)
		snp_posterior_vars = 1.0/((gene_Q_sum_sq*(np.square(gene_alpha_mu) + gene_alpha_mu_var)*self.N_gwas/gwas_resid_var) + (gene_Q_sum_sq*self.N_eqtl/gene_resid_var) + (1.0/gene_delta_var))


		for kk in np.random.permutation(range(gene_K)):

			# Re include current effect
			window_gwas_beta_resid = window_gwas_beta_resid + (gene_Q[:, kk]*delta_mu[kk]*gene_alpha_mu)
			eqtl_gene_beta_resid = eqtl_gene_beta_resid + (gene_Q[:, kk]*delta_mu[kk])

			# Compute posterior distribution
			gwas_term = (self.N_gwas/gwas_resid_var)*gene_alpha_mu*np.dot(window_gwas_beta_resid, gene_Q[:, kk])
			eqtl_term = (self.N_eqtl/gene_resid_var)*np.dot(eqtl_gene_beta_resid, gene_Q[:, kk])
			snp_posterior_mean = snp_posterior_vars[kk]*(gwas_term + eqtl_term)

			# Sample from posterior distribution
			delta_mu[kk] = snp_posterior_mean
			delta_mu_var[kk] = snp_posterior_vars[kk]

			# Remove Updated effect
			window_gwas_beta_resid = window_gwas_beta_resid - (gene_Q[:, kk]*delta_mu[kk]*gene_alpha_mu)
			eqtl_gene_beta_resid = eqtl_gene_beta_resid - (gene_Q[:, kk]*delta_mu[kk])


		return delta_mu, delta_mu_var, window_gwas_beta_resid, eqtl_gene_beta_resid

	def update_gamma_and_alpha_in_single_window(self, window_name, Q_mat):
		# Extract relevent quantities
		window_gamma_mu = self.gamma_mu[window_name]
		window_gamma_mu_var = self.gamma_mu_var[window_name]
		window_gwas_beta_resid = self.gwas_beta_resid[window_name]
		window_gwas_resid_var = self.gwas_resid_var[window_name]

		# Precompute posterior variance (same for all snps)
		snp_posterior_var = 1.0/((self.window_info[window_name]['Q_sum_sq']*self.N_gwas/window_gwas_resid_var) + (1.0/self.gamma_var))

		# Loop through genetic elements in random order
		window_genetic_elements = self.window_info[window_name]['genetic_elements']
		n_genetic_elements = len(window_genetic_elements)
		ordering = np.random.permutation(np.arange(n_genetic_elements))
		for genetic_element_index in ordering:
			# Name of genetic element
			genetic_element_class, genetic_element_name = window_genetic_elements[genetic_element_index]

			# Updates for snps
			if genetic_element_class == 'snp':
				window_gamma_mu, window_gamma_mu_var, window_gwas_beta_resid = self.gamma_update_for_single_snp(window_gamma_mu, window_gamma_mu_var, window_gwas_beta_resid, Q_mat, genetic_element_name, snp_posterior_var, window_gwas_resid_var)
			#Updates for genes
			elif genetic_element_class == 'gene':
				window_gwas_beta_resid = self.alpha_update_for_single_gene(window_gwas_beta_resid, Q_mat, genetic_element_name, window_gwas_resid_var)
			else:
				print('asssumptino eroror: should not have this class of genetic element')
				pdb.set_trace()

		# Re update quantities
		self.gamma_mu[window_name] = window_gamma_mu
		self.gamma_mu_var[window_name] = window_gamma_mu_var
		self.gwas_beta_resid[window_name] = window_gwas_beta_resid

		return


	def alpha_update_for_single_gene(self, window_gwas_beta_resid, Q_mat, gene_name, window_gwas_resid_var):
		delta_mu = self.deltas_mu[gene_name]
		delta_mu_var = self.deltas_mu_var[gene_name]
		window_cis_indices = self.gene_info[gene_name]['cis_snps']

		# Precompute Q-delta
		small_Q = Q_mat[:, window_cis_indices]
		Q_delta = np.dot(small_Q, delta_mu)

		n_cis_snps = len(delta_mu)
		E_delta_delta_t = np.diag(delta_mu_var) + np.dot(delta_mu.reshape((n_cis_snps,1)), delta_mu.reshape((1, n_cis_snps)))
		Q_delta_sq_sum = np.trace(np.dot(np.dot(np.transpose(small_Q), small_Q), E_delta_delta_t))
		#q_delta_sq_sum_debug(small_Q, delta_mu, delta_mu_var, Q_delta_sq_sum)

		# Re include current effect
		window_gwas_beta_resid = window_gwas_beta_resid + (Q_delta*self.alpha_mu[gene_name])

		# Compute posterior distribution
		gene_posterior_var = 1.0/((Q_delta_sq_sum*self.N_gwas/window_gwas_resid_var) + (1.0/self.alpha_var))
		gene_posterior_mean = gene_posterior_var*(self.N_gwas/window_gwas_resid_var)*np.dot(window_gwas_beta_resid, Q_delta)

		# Draw from posterior distribution
		self.alpha_mu[gene_name] = gene_posterior_mean
		self.alpha_mu_var[gene_name] = gene_posterior_var

		# Remove updated effect
		window_gwas_beta_resid = window_gwas_beta_resid - (Q_delta*self.alpha_mu[gene_name])

		return window_gwas_beta_resid

	def gamma_update_for_single_snp(self, window_gamma_mu, window_gamma_mu_var, window_gwas_beta_resid, Q_mat, snp_index, snp_posterior_var, window_gwas_resid_var):

		# Re include current effect
		window_gwas_beta_resid = window_gwas_beta_resid + (Q_mat[:, snp_index]*window_gamma_mu[snp_index])

		# Compute posterior mean for this snp
		snp_posterior_mean = snp_posterior_var[snp_index]*(self.N_gwas/window_gwas_resid_var)*np.dot(window_gwas_beta_resid,Q_mat[:,snp_index])

		# Sample from posterior distribution
		window_gamma_mu[snp_index] = snp_posterior_mean
		window_gamma_mu_var[snp_index] = snp_posterior_var[snp_index]

		# Remove updated effect
		window_gwas_beta_resid = window_gwas_beta_resid - (Q_mat[:, snp_index]*window_gamma_mu[snp_index])

		return window_gamma_mu, window_gamma_mu_var, window_gwas_beta_resid




	def initialize_variables(self):
		# Initialize causal variant effcts
		self.gamma_mu = {}
		self.gamma_mu_var = {}
		self.gwas_beta_resid = {}
		self.gwas_resid_var = {}
		self.window_med_h2 = {}
		self.window_nm_h2 = {}
		for window_name in self.windows:
			self.gamma_mu[window_name] = np.zeros(self.window_info[window_name]['n_snps'])
			self.gamma_mu_var[window_name] = np.ones(self.window_info[window_name]['n_snps'])
			self.gwas_beta_resid[window_name] = np.copy(self.window_info[window_name]['beta_pc'])

			# Add sum of square terms (precomputed)
			window_Q = np.load(self.window_info[window_name]['Q_file'])
			self.window_info[window_name]['Q_sum_sq'] = np.sum(np.square(window_Q),axis=0)

			self.gwas_resid_var[window_name] = 1.0

			self.window_med_h2[window_name] = 0.0
			self.window_nm_h2[window_name] = 0.0




		# Create list of genetic elements in each window
		for window_name in self.windows:
			window_snps = self.window_info[window_name]['rsids']
			window_genes = self.window_info[window_name]['genes']
			window_genetic_elements = []
			for ii in range(len(window_snps)):
				window_genetic_elements.append(('snp', ii))
			for ii, gene_name in enumerate(window_genes):
				window_genetic_elements.append(('gene', gene_name))
			# Add to window info
			self.window_info[window_name]['genetic_elements'] = window_genetic_elements


		# Initialize causal gene effects
		self.alpha_mu = {}
		self.alpha_mu_var = {}
		self.deltas_mu = {}
		self.deltas_mu_var = {}
		self.delta_vars = {}
		self.eqtl_resid_vars = {}
		self.eqtl_beta_resid = {}
		for gene in self.genes:
			# Initialize causal gene-trait effects to zero
			self.alpha_mu[gene] = 0.0
			self.alpha_mu_var[gene] = 0.0

			# Initialize variant-gene effects with best fit
			gene_beta_pc = self.gene_info[gene]['beta_pc']
			gene_window = self.gene_info[gene]['window_name']
			gene_cis_snp_indices = self.gene_info[gene]['cis_snps']
			gene_Q = np.load(self.window_info[gene_window]['Q_file'])[:, gene_cis_snp_indices]
			gene_mod = bayesian_lmm_pca_ss_h2_single_region.Bayesian_LMM(gene_beta_pc, gene_Q, self.N_eqtl)
			gene_mod.fit(burn_in_iterations=2, total_iterations=5)

			# Initialize
			self.deltas_mu[gene] = np.mean(gene_mod.sampled_gammas,axis=0)
			self.deltas_mu_var[gene] = np.mean(gene_mod.sampled_gammas,axis=0)*0.0 + 0.0
			self.delta_vars[gene] = np.mean(gene_mod.sampled_gamma_vars)
			#self.eqtl_resid_vars[gene] = np.mean(gene_mod.sampled_resid_vars)
			self.eqtl_resid_vars[gene] = 1.0
			self.eqtl_beta_resid[gene] = gene_beta_pc - np.dot(gene_Q, self.deltas_mu[gene])


		# Initialize variance parameters
		self.gamma_var = 1e-5
		self.alpha_var = 1e-2


		# Keep track of sampled gamma_vars
		self.sampled_nm_h2 = []
		self.sampled_med_h2 = []
		self.sampled_eqtl_h2 = []
		self.sampled_gamma_var = []
		self.sampled_alpha_var = []
		self.sampled_delta_vars = []
		self.sampled_gwas_resid_var = []
		self.sampled_eqtl_resid_var = []


		return





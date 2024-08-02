import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import time

def update_gamma_from_single_window(LD, gamma_vec, gwas_beta_resid, gamma_var, N_gwas):
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





class Bayesian_LMM_RSS_h2_med_inference(object):
	def __init__(self, gwas_beta, n_gwas_individuals, window_names, window_info, gene_info, temp_output_file):
		self.gwas_beta = gwas_beta
		self.N_gwas = n_gwas_individuals
		self.window_names = window_names
		self.window_info = window_info
		self.gene_info = gene_info

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
			# Update gamma
			self.update_gamma_and_alpha()

			# Update deltas
			self.update_deltas()

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

	def update_gamma_and_delta_in_single_window(self, window_name, LD):
		# Precompute posterior variance (same for all snps)
		snp_posterior_var = 1.0/(self.N_gwas + (1.0/self.gamma_var))

		# Update each element in turn
		for element_name in np.random.permutation(self.window_info[window_name]['elements']):
			if element_name.startswith('ENSG'):
				# Relevent fields
				gene_window_indices = self.window_info[window_name]['gene_to_window_indices'][element_name]
				gene_index = self.window_info[window_name]['gene_name_to_gene_index'][element_name]
				# Re-include current effect
				genetic_pred_expr = np.dot(LD[:, gene_window_indices], self.gene_to_deltas[element_name])
				self.window_to_gwas_beta_resid[window_name] = self.window_to_gwas_beta_resid[window_name] + genetic_pred_expr*self.window_to_alpha[window_name][gene_index]

				# Compute posterior variance for this gene
				gene_variance = np.dot(np.dot(self.gene_to_deltas[element_name], LD[:, gene_window_indices][gene_window_indices,:]), self.gene_to_deltas[element_name])
				gene_posterior_var = 1.0/((self.N_gwas*gene_variance) + (1.0/self.alpha_var))
				# Compute posterior mean for this gene
				gene_posterior_mean = gene_posterior_var*self.N_gwas*np.dot(self.gene_to_deltas[element_name], self.window_to_gwas_beta_resid[window_name][gene_window_indices])

				self.window_to_alpha[window_name][gene_index] = np.random.normal(loc=gene_posterior_mean, scale=np.sqrt(gene_posterior_var))

				# Remove updated effect
				self.window_to_gwas_beta_resid[window_name] = self.window_to_gwas_beta_resid[window_name] - genetic_pred_expr*self.window_to_alpha[window_name][gene_index]
			else:
				# SNP
				snp_index = int(element_name)
				# Re-include current effect
				self.window_to_gwas_beta_resid[window_name] = self.window_to_gwas_beta_resid[window_name] + (LD[snp_index, :]*self.window_to_gamma[window_name][snp_index])
				
				# Compute posterior mean for this snp
				snp_posterior_mean = snp_posterior_var*self.N_gwas*self.window_to_gwas_beta_resid[window_name][snp_index]

				# Sample from posterior distribution
				self.window_to_gamma[window_name][snp_index] = np.random.normal(loc=snp_posterior_mean, scale=np.sqrt(snp_posterior_var))

				# Remove updated effect
				self.window_to_gwas_beta_resid[window_name] = self.window_to_gwas_beta_resid[window_name] - (LD[snp_index, :]*self.window_to_gamma[window_name][snp_index])
		return


	def update_gamma_and_alpha(self):
		for window_name in self.window_names:
			window_ld = np.load(self.window_info[window_name]['ld_file'])
			self.update_gamma_and_delta_in_single_window(window_name, window_ld)

		return


	def initialize_gene_delta_from_only_eqtl_data(self, gene_name, window_name, window_ld):
		gene_window_indices = self.window_info[window_name]['gene_to_window_indices'][gene_name]
		eqtl_beta = self.gene_info[gene_name]['eqtl_beta']
		N_eqtl = self.gene_info[gene_name]['N']

		gene_window_ld = window_ld[gene_window_indices, :][:, gene_window_indices]
		delta_var = self.gene_to_delta_var[gene_name]


		S_inv = N_eqtl*gene_window_ld + (np.eye(len(eqtl_beta))/delta_var)
		SS = np.linalg.inv(S_inv)

		self.gene_to_deltas[gene_name] = N_eqtl*np.dot(SS, eqtl_beta)

		return

	def update_delta_for_specific_gene(self, gene_name, window_name, LD):
		gene_window_indices = self.window_info[window_name]['gene_to_window_indices'][gene_name]
		eqtl_beta = self.gene_info[gene_name]['eqtl_beta']
		N_eqtl = self.gene_info[gene_name]['N']	
		
		gene_window_ld = LD[gene_window_indices, :][:, gene_window_indices]
		delta_var = self.gene_to_delta_var[gene_name]
		
		gene_index = self.window_info[window_name]['gene_name_to_gene_index'][gene_name]
		# Re-include current effect
		genetic_pred_expr = np.dot(LD[:, gene_window_indices], self.gene_to_deltas[gene_name])
		self.window_to_gwas_beta_resid[window_name] = self.window_to_gwas_beta_resid[window_name] + genetic_pred_expr*self.window_to_alpha[window_name][gene_index]

		pdb.set_trace()


	def update_deltas(self):
		used_genes = {}
		for window_name in self.window_names:
			window_ld = np.load(self.window_info[window_name]['ld_file'])
			for gene_name in self.window_info[window_name]['middle_genes']:
				used_genes[gene_name] = 1
				self.update_delta_for_specific_gene(gene_name, window_name, window_ld)
		for window_name in self.window_names:
			window_ld = np.load(self.window_info[window_name]['ld_file'])
			for gene_name in self.window_info[window_name]['all_genes']:
				if gene_name in used_genes:
					continue
				used_genes[gene_name] = 1
				self.update_delta_for_specific_gene(gene_name, window_name, window_ld)
		return


	def initialize_variables(self):
		# Initialize causal effcts
		self.window_to_gamma = {}
		self.window_to_gwas_beta_resid = {}
		self.window_to_alpha = {}
		self.gene_to_deltas = {}
		self.genes = []
		self.gene_to_delta_var = {}

		used_genes = {}
		for window_name in self.window_names:
			self.window_to_gamma[window_name] = np.zeros(len(self.window_info[window_name]['positions']))
			self.window_to_gwas_beta_resid[window_name] = np.copy(self.gwas_beta[self.window_info[window_name]['positions']])
			self.window_to_alpha[window_name] = np.zeros(len(self.window_info[window_name]['all_genes']))

			window_ld = np.load(self.window_info[window_name]['ld_file'])
			for gene_name in self.window_info[window_name]['middle_genes']:
				used_genes[gene_name] = 1
				self.genes.append(gene_name)
				self.gene_to_deltas[gene_name] = np.zeros(len(self.gene_info[gene_name]['eqtl_beta']))
				if len(self.gene_info[gene_name]['eqtl_beta']) != len(self.window_info[window_name]['gene_to_window_indices'][gene_name]):
					print('assumption eroror')
					pdb.ste_trace()
				n_cis_snps = len(self.gene_info[gene_name]['eqtl_beta'])

				# Initialize gene quantities
				self.gene_to_delta_var[gene_name] = 0.2/n_cis_snps
				
				self.initialize_gene_delta_from_only_eqtl_data(gene_name, window_name, window_ld)
				#self.gene_to_deltas[gene_name] = np.copy(self.gene_info[gene_name]['eqtl_beta'])
		for window_name in self.window_names:
			window_ld = np.load(self.window_info[window_name]['ld_file'])
			for gene_name in self.window_info[window_name]['all_genes']:
				if gene_name in used_genes:
					continue
				used_genes[gene_name] = 1
				self.genes.append(gene_name)
				self.gene_to_deltas[gene_name] = np.zeros(len(self.gene_info[gene_name]['eqtl_beta']))
				if len(self.gene_info[gene_name]['eqtl_beta']) != len(self.window_info[window_name]['gene_to_window_indices'][gene_name]):
					print('assumption eroror')
					pdb.ste_trace()
				n_cis_snps = len(self.gene_info[gene_name]['eqtl_beta'])

				# Initialize gene quantities
				self.gene_to_delta_var[gene_name] = 0.2/n_cis_snps
				
				self.initialize_gene_delta_from_only_eqtl_data(gene_name, window_name, window_ld)
				#self.gene_to_deltas[gene_name] = np.copy(self.gene_info[gene_name]['eqtl_beta'])



		self.genes = np.asarray(self.genes)


		# Initialize variance parameters
		self.gamma_var = 1e-6
		self.alpha_var = 1e-6

		# Keep track of sampled gamma_vars
		self.sampled_gamma_vars = []
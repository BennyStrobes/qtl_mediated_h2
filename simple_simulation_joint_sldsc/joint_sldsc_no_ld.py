import numpy as np
import os
import pdb
import statsmodels.api as sm





def compute_emperical_variance(pred_value, true_value):
	return np.sum(np.square(true_value - pred_value))/(len(true_value)-2)


class JOINT_SLDSC(object):
	def __init__(self, max_iter=100):
		self.max_iter = max_iter
	def fit(self, trait_chi_sq, trait_ss, eqtl_info):
		self.trait_chi_sq = trait_chi_sq
		self.trait_ss = trait_ss
		self.eqtl_info = eqtl_info
		self.initialize()

		# Now run iterative algorithm until convergence
		for itera in range(self.max_iter):
			# First update alpha_sq and beta_sq
			self.update_gwas_h2_parameters()

			# Update eqtl variance
			self.update_eqtl_variances()

			# Now update delta_sq
			self.update_delta_sq_parameters()

			#print(self.alpha_sq*60*.2)
			#print(self.beta_sq*self.n_snps)
		self.gene_med_h2 = 0
		self.gene_h2s = []
		for gene_name in self.gene_names:
			self.gene_med_h2 = self.gene_med_h2 + np.sum(self.eqtl_info[gene_name]['delta_sq']*self.alpha_sq)
			gene_h2 = np.sum(self.eqtl_info[gene_name]['delta_sq'])
			self.gene_h2s.append(gene_h2)
		self.gene_h2s = np.asarray(self.gene_h2s)
		self.nm_h2 = self.beta_sq*self.n_snps
		return
		


	def update_eqtl_variances(self):
		# update each gene
		variances = []
		for gene_name in self.gene_names:
			eqtl_chi_sq = self.eqtl_info[gene_name]['chi_sq']
			eqtl_ss = self.eqtl_info[gene_name]['ss']

			variance = compute_emperical_variance(self.eqtl_info[gene_name]['delta_sq']*eqtl_ss, eqtl_chi_sq-1.0)
			self.eqtl_info[gene_name]['eqtl_var'] = variance
			variances.append(variance)

		self.eqtl_var = np.mean(variances)

		return



	def update_delta_sq_parameters(self):
		# First get predicted triat chi-sq
		pred_trait = self.trait_ss*self.beta_sq*np.ones(self.n_snps)
		for gene_name in self.gene_names:
			gene_indices = self.eqtl_info[gene_name]['indices']
			pred_trait[gene_indices] = pred_trait[gene_indices] + self.trait_ss*self.eqtl_info[gene_name]['delta_sq']*self.alpha_sq
	
		# update each gene
		for gene_name in self.gene_names:
			gene_indices = self.eqtl_info[gene_name]['indices']
			# Remove current genes effect from pred trait
			pred_trait[gene_indices] = pred_trait[gene_indices] - self.trait_ss*self.eqtl_info[gene_name]['delta_sq']*self.alpha_sq

			eqtl_chi_sq = self.eqtl_info[gene_name]['chi_sq']
			eqtl_ss = self.eqtl_info[gene_name]['ss']

			for eqtl_snp_iter, gwas_snp_iter in enumerate(gene_indices):
				output = np.zeros(2)
				output[0] = (eqtl_chi_sq[eqtl_snp_iter] -1)
				output[1] = self.trait_chi_sq[gwas_snp_iter] - 1.0 - pred_trait[gwas_snp_iter]
				input_vec = np.zeros(2)
				input_vec[0] = eqtl_ss
				input_vec[1] = self.trait_ss*self.alpha_sq
				# Compute variances
				eqtl_var = np.square((1.0 + (eqtl_ss*self.eqtl_info[gene_name]['delta_sq'][eqtl_snp_iter])))
				gwas_var = np.square((1.0 + (self.trait_ss*(self.beta_sq + (self.alpha_sq*self.eqtl_info[gene_name]['delta_sq'][eqtl_snp_iter])))))
				# Consider weighting by pred variance (see LOGLSS paper or supp of sldsc)
				weights = np.zeros(2)
				weights[0] = 1.0/eqtl_var
				weights[1] = 1.0/gwas_var
				small_reg = sm.WLS(output, input_vec, weights=weights).fit()
				#small_reg = sm.OLS(output, input_vec).fit()
				self.eqtl_info[gene_name]['delta_sq'][eqtl_snp_iter] = small_reg.params[0]
				#self.eqtl_info[gene_name]['delta_sq'][eqtl_snp_iter] = (eqtl_chi_sq[eqtl_snp_iter] -1)/eqtl_ss

			# Add back in current gene into pred trait
			pred_trait[gene_indices] = pred_trait[gene_indices] + self.trait_ss*self.eqtl_info[gene_name]['delta_sq']*self.alpha_sq
		return


	def update_gwas_h2_parameters(self):
		# Get LD scores
		var_ld_scores = np.ones(self.n_snps)*self.trait_ss
		# get gene ld scores
		gene_ld_scores = np.zeros(self.n_snps)
		for gene_name in self.gene_names:
			gene_indices = self.eqtl_info[gene_name]['indices']
			gene_ld_scores[gene_indices] = gene_ld_scores[gene_indices] + self.trait_ss*self.eqtl_info[gene_name]['delta_sq']

		# Put together in a matrix
		joint_ld_scores = np.transpose(np.vstack((var_ld_scores,gene_ld_scores)))

		# Run regression
		ldsc_trait_reg = sm.OLS(self.trait_chi_sq - 1.0, joint_ld_scores).fit()

		# Update gwas h2 params
		self.beta_sq = ldsc_trait_reg.params[0]
		self.alpha_sq = ldsc_trait_reg.params[1]

		variance = compute_emperical_variance(np.dot(joint_ld_scores, ldsc_trait_reg.params), self.trait_chi_sq-1.0)
		self.gwas_var = variance

		return

	def initialize(self):
		# Get gene names
		self.gene_names = np.sort([*self.eqtl_info])

		# Initialize gene related annotations
		for gene_name in self.gene_names:
			self.eqtl_info[gene_name]['delta_sq'] = np.copy(self.eqtl_info[gene_name]['init'])
			self.eqtl_info[gene_name]['eqtl_var'] = 1.0

		# Get number of snps
		self.n_snps = len(self.trait_chi_sq)

		# Initialize alpha_sq and beta_sq
		# Doesn't really matter as it will get updated first
		self.alpha_sq = .1
		self.beta_sq = .1

		# Variances
		self.gwas_var = 1.0


import numpy as np
import pdb
import statsmodels.api as sm
from sklearn.linear_model import ARDRegression, BayesianRidge, LinearRegression
from scipy.stats import invgamma
import time


class LR_SaS_Prior(object):
	def __init__(self, K=2):
		# Number of prior mixture components
		self.K = K
	def fit(self, Y, G, total_iterations=10000, burn_in_iterations=5000):
		""" Fit the model.
			Args:
			G: A genotype matrix of floats with shape [num_samples, num_snps].
			Y: A trait vector vector of floats of length num_samples.
		"""
		# Add data to object
		self.Y = Y
		self.G = G

		#####################
		# Initialize variables
		print('###############################')
		print('Initialize variables')
		print('###############################')
		self.initialize_variables()

		#####################
		# Initialize variables
		print('###############################')
		print('Sampling')
		print('###############################')
		for itera in range(total_iterations):
			if np.mod(itera,1000) == 0:
				print(itera)
			# First update causal variant-trait effects
			self.update_causal_variant_trait_effects_shell()
			# Updates pis (prior probabilities on class assignments)
			self.update_pi()
			# Update beta variance
			self.update_beta_variance()
			# Update residual variance
			self.update_resid_variance()
			if itera > burn_in_iterations:
				if np.mod(itera, 10) == 0:
					self.sampled_betas.append(np.copy(self.beta))
					self.sampled_h2.append(np.var(np.dot(self.G, self.beta)))
					self.sampled_pis.append(np.copy(self.pis))
					self.sampled_beta_vars.append(self.beta_var)
					self.sampled_resid_vars.append(self.resid_var)
		self.sampled_pis = np.asarray(self.sampled_pis)
		self.sampled_betas = np.asarray(self.sampled_betas)
		self.sampled_h2 = np.asarray(self.sampled_h2)
		self.sampled_beta_vars = np.asarray(self.sampled_beta_vars)
		self.sampled_resid_vars = np.asarray(self.sampled_resid_vars)
		return

	def update_resid_variance(self):
		# Sum of sequared errors
		SSE = np.dot(self.Y_resid, self.Y_resid)
		# Scaled inverse chi-squared params
		#v = self.v0 + self.N
		#tau_sq = (self.v0*self.s2_0 + SSE)/(self.v0 + self.N)
		v = self.N
		tau_sq = SSE/self.N

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma(v/2 + 1e-5, scale=v*tau_sq/2 + 1e-5)
		# Sample from it
		self.resid_var = invgamma_dist.rvs(size=1)[0]

		return

	def update_beta_variance(self):
		# Count up number of snps with non-zero effect
		qq = np.sum(self.class_membership == 1)

		# Scaled inverse chi-squared params
		#v = self.v0 + qq
		#tau_sq = (self.v0*self.s2_0 + np.sum(np.square(self.beta)))/(self.v0 + qq)

		if qq == 0.0:
			v = self.prior_v
			tau_sq = self.prior_tau_sq
			print('zero')
		else:
			v = qq
			tau_sq = np.sum(np.square(self.beta))/(qq)


		# Initialize inverse gamma distribution
		invgamma_dist = invgamma(v/2 + 1e-5, scale=v*tau_sq/2 + 1e-5)
		# Sample from it
		self.beta_var = invgamma_dist.rvs(size=1)[0]

		self.prior_v = v
		self.prior_tau_sq = tau_sq
		
		return


	def update_pi(self):
		# Count up number of variants in each class
		counts = np.asarray([np.sum(self.class_membership == 0.0), np.sum(self.class_membership == 1.0)])
		# Randomly sample pi from dirichlet
		self.pis = np.random.dirichlet(counts + self.alpha_0)
		return



	def update_causal_variant_trait_effects_shell(self):
		# Randomly permute snp update ordering
		original_snp_ordering = np.arange(self.P)
		permuted_snp_ordering = np.random.permutation(original_snp_ordering)

		# Loop through snps
		for ppp in permuted_snp_ordering:

			# Keep track of current beta
			old_beta = self.beta[ppp]

			# Compute x^ty with Y_resid (after re-including current beta)
			r_j = np.dot(self.G[:,ppp], self.Y_resid) + (self.GtG[ppp]*self.beta[ppp])

			# Compute some temporary terms
			lhs_1 = self.GtG[ppp] + (self.resid_var/self.beta_var)

			# Compute log likelihoods of both classes
			log_like_0 = np.log(self.pis[0]) 
			log_like_1 = np.log(self.pis[1]) - (.5*np.log(self.beta_var*lhs_1/self.resid_var)) + (.5*(np.square(r_j)/(self.resid_var*lhs_1)))

			# Get probabilities of each
			p0 = 1.0/(np.exp(log_like_1 - log_like_0) + np.exp(0))
			p1= 1.0/(np.exp(log_like_0 - log_like_1) + np.exp(0))

			# Sample class membership
			class_p = np.random.choice([0,1], p=[p0,p1])

			# Update class membership
			self.class_membership[ppp] = class_p

			# Update betas
			if class_p == 0:  # The spike
				self.beta[ppp] = 0.0
			else:  # The slab
				effect_mean = r_j/lhs_1
				effect_var = self.resid_var/lhs_1
				self.beta[ppp] = np.random.normal(loc=effect_mean, scale=np.sqrt(effect_var))

			# Re-include current snp effect in residual Y
			self.Y_resid = self.Y_resid - self.G[:,ppp]*(self.beta[ppp] - old_beta)

		return



	def initialize_variables(self):
		# Get sample size and num snps
		self.N = self.G.shape[0]
		self.P = self.G.shape[1]

		# Run bayesian ridge regression for initialization
		brr = BayesianRidge(compute_score=True, n_iter=100000, fit_intercept=False).fit(self.G, self.Y)

		# Initialize causal variant-trait effects
		self.beta = np.copy(brr.coef_)

		# Initialize pis
		self.pis = np.ones(self.K)/self.K

		# Initialize hyperparameter on pis
		self.alpha_0 = np.ones(self.K)

		# Initialize snp assignments
		self.class_membership = np.zeros(self.P)

		# Initialize causal effect variances
		pred_genetic_component = np.dot(self.G, self.beta)
		self.beta_var = np.var(pred_genetic_component)
		if self.beta_var > 0.7:
			self.beta_var = 0.7

		# Initialize residual variance
		self.resid_var = 1.0 - self.beta_var

		# Variance hyperparams
		self.v0 = -2.0
		self.s2_0 = 0.0

		# Initialize residual Y
		self.Y_resid = self.Y - pred_genetic_component

		# Precompute GtG
		self.GtG = np.zeros(self.P)
		for ppp in range(self.P):
			self.GtG[ppp] = np.dot(self.G[:, ppp], self.G[:,ppp])

		# Keep track of betas
		self.sampled_betas = []
		self.sampled_h2 = []
		self.sampled_pis = []
		self.sampled_beta_vars = []
		self.sampled_resid_vars = []



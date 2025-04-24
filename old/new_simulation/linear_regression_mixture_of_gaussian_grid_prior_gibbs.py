import numpy as np
import pdb
import statsmodels.api as sm
from sklearn.linear_model import ARDRegression, BayesianRidge, LinearRegression
from scipy.stats import invgamma
import time
import scipy.special


class LR_MoG_Grid_Prior(object):
	def __init__(self, variance_grid_start=1e-7, variance_grid_end=0.5, standard_dev_multiplier=np.sqrt(2)):
		# Variance Grid
		self.beta_variance_grid = []
		cur_variance =variance_grid_start
		while cur_variance < variance_grid_end:
			self.beta_variance_grid.append(cur_variance)
			cur_variance = np.square(np.sqrt(cur_variance)*standard_dev_multiplier)
		self.beta_variance_grid.append(cur_variance)
		self.beta_variance_grid = np.asarray(self.beta_variance_grid)

		# Number of components
		self.K = len(self.beta_variance_grid)


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
			# Update residual variance
			self.update_resid_variance()
			if itera > burn_in_iterations:
				if np.mod(itera, 10) == 0:
					self.sampled_betas.append(np.copy(self.beta))
					self.sampled_h2.append(np.var(np.dot(self.G, self.beta)))
					self.sampled_pis.append(np.copy(self.pis))
					self.sampled_resid_vars.append(self.resid_var)
		self.sampled_pis = np.asarray(self.sampled_pis)
		self.sampled_betas = np.asarray(self.sampled_betas)
		self.sampled_h2 = np.asarray(self.sampled_h2)
		self.sampled_resid_vars = np.asarray(self.sampled_resid_vars)
		return

	def update_resid_variance(self):
		# Sum of sequared errors
		SSE = np.dot(self.Y_resid, self.Y_resid)
		# Scaled inverse chi-squared params
		v = self.v0 + self.N
		tau_sq = (self.v0*self.s2_0 + SSE)/(self.v0 + self.N)
		#v = self.N
		#tau_sq = SSE/self.N

		# Initialize inverse gamma distribution
		invgamma_dist = invgamma(v/2, scale=v*tau_sq/2)
		# Sample from it
		self.resid_var = invgamma_dist.rvs(size=1)[0]

		return


	def update_pi(self):
		# Count up number of variants in each class
		counts = []
		for kk in range(self.K):
			counts.append(np.sum(self.class_membership == kk))
		counts = np.asarray(counts)
		# Randomly sample pi from dirichlet
		self.pis = np.random.dirichlet(counts + self.alpha_0)

		# Set values of zero to really small number
		self.pis[self.pis==0.0]=1e-30

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

			# Temporary terms
			lhs = self.GtG[ppp] + (self.resid_var/self.beta_variance_grid)

			# Compute log likelihoods of each class classes
			log_like = np.log(self.pis) - (.5*np.log(self.beta_variance_grid*lhs/self.resid_var)) + (.5*(np.square(r_j)/(self.resid_var*lhs)))

			# Get probability of each class
			temp = scipy.special.logsumexp(log_like - log_like[:, np.newaxis],axis=1)
			if np.sum(temp >600):
				temp[temp>600] = 600
			probs = 1.0/np.exp(temp)
			#probs = 1.0/np.sum(np.exp(log_like - log_like[:, np.newaxis]),axis=1)

			# Sample class membership
			class_p = np.random.choice(self.classes, p=probs)

			# Update class membership
			self.class_membership[ppp] = class_p

			# Update betas
			effect_mean = r_j/(lhs[class_p])
			effect_var = self.resid_var/(lhs[class_p])
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

		# Initialize residual variance
		pred_genetic_component = np.dot(self.G, self.beta)
		self.resid_var = 1.0 - np.var(pred_genetic_component)

		# Variance hyperparams
		self.v0 = -2.0
		self.s2_0 = 0.0

		# Initialize residual Y
		self.Y_resid = self.Y - pred_genetic_component

		# Precompute GtG
		self.GtG = np.zeros(self.P)
		for ppp in range(self.P):
			self.GtG[ppp] = np.dot(self.G[:, ppp], self.G[:,ppp])

		# Self class names
		self.classes = np.arange(self.K)

		# Keep track of betas
		self.sampled_betas = []
		self.sampled_h2 = []
		self.sampled_pis = []
		self.sampled_resid_vars = []



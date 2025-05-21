import numpy as np
import os
import sys
import pdb
import scipy.special
import time


def get_sigma_sq_grid(m_factor, min_sigma, max_sigma):
	sigma_grid = []
	cur_sigma = min_sigma
	while cur_sigma < max_sigma:
		sigma_grid.append(cur_sigma)
		cur_sigma = cur_sigma*m_factor
	sigma_grid = np.asarray(sigma_grid)
	sigma_sq_grid = np.square(sigma_grid)

	return sigma_sq_grid

def sample_rows(probs):
	"""
	probs: array of shape (N, K), each row sums to 1
	returns: integer array of shape (N,) with values in 0..K-1
	"""
	# 1. cumulative sums along each row → shape (N, K)
	cum = np.cumsum(probs, axis=1)
	# 2. one uniform per row → shape (N,)
	r = np.random.rand(probs.shape[0])
	# 3. for each row i, find first index j where cum[i,j] > r[i]
	#    this gives the sampled column
	return np.argmax(cum > r[:, None], axis=1)


class model(object):
	def __init__(self, max_iter=10, burn_in_iter=5, m_factor=np.sqrt(2), alpha_0=1e-10, variance_grid_lb_divisor=10):
		self.max_iter = max_iter  # Total iterations
		self.burn_in_iter = burn_in_iter
		self.m_factor = m_factor
		self.alpha_0 = alpha_0
		self.variance_grid_lb_divisor = variance_grid_lb_divisor

	def fit(self, beta_hats, beta_hat_ses):
		self.beta_hats = np.copy(beta_hats)
		self.beta_hat_ses = np.copy(beta_hat_ses)
		self.beta_hat_vars = np.square(self.beta_hat_ses)

		self.initialize_model_parameters()

		self.sampled_betas = []
		self.sampled_pis = []
		self.n_samples = 0
		self.sampled_beta_sum_mean = np.zeros(len(beta_hats))
		self.sampled_beta_sum_sq_mean = np.zeros(len(beta_hats))
		self.sampled_pi_sum_mean = np.zeros(self.C)


		for itera in range(self.max_iter):
			print(itera)
			self.beta_update()
			self.pi_update()

			if itera > self.burn_in_iter:
				#self.sampled_betas.append(self.betas)
				#self.sampled_pis.append(self.pis)
				self.sampled_beta_sum_mean = self.sampled_beta_sum_mean + np.copy(self.betas)
				self.sampled_beta_sum_sq_mean = self.sampled_beta_sum_sq_mean + np.square(np.copy(self.betas))
				self.sampled_pi_sum_mean = self.sampled_pi_sum_mean + np.copy(self.pis)
				self.n_samples = self.n_samples + 1

		#self.sampled_betas = np.asarray(self.sampled_betas)
		#self.sampled_pis = np.asarray(self.sampled_pis)

		self.sampled_beta_mean = (self.sampled_beta_sum_mean/self.n_samples)
		self.sampled_beta_var = (self.sampled_beta_sum_sq_mean/self.n_samples) - np.square(self.sampled_beta_mean)
		self.sampled_pis_mean = (self.sampled_pi_sum_mean/self.n_samples)

		return



	def pi_update(self):
		counts = np.zeros(self.C)
		for class_name in self.classes:
			counts[class_name] = counts[class_name] + np.sum(self.class_memberships ==class_name)

		# Randomly sample pi from dirichlet
		self.pis = np.random.dirichlet(counts + self.alpha_0)

		# Set values of zero to really small number
		self.pis[self.pis==0.0]=1e-80

		return

	def beta_update(self):
		'''
		# Compute posterior variance grid
		num = np.outer(self.beta_hat_vars, self.sigma_sq_grid)   # shape (K, C)
		den = self.beta_hat_vars[:, None] + self.sigma_sq_grid[None, :]
		var_post_grid = num / den
		# Compute posterior mean grid
		mu_post_grid = (var_post_grid / self.beta_hat_vars[:, np.newaxis])*self.beta_hats[:, np.newaxis] 
		'''
		# Compute log likelihoods of each class classes
		log_like = np.log(self.pis) + self.part_log_like


		lse = scipy.special.logsumexp(log_like, axis=1)            # shape (N,)
		# 2) subtract each element in the row:
		temp = lse[:, None] - log_like   

		#if np.sum(temp >600):
			#temp[temp>600] = 600
		probs = 1.0/np.exp(temp)

		# Randomly sample class memberships
		self.class_memberships = sample_rows(probs)

		selected_mu_post = self.mu_post_grid[np.arange(self.mu_post_grid.shape[0]), self.class_memberships]
		selected_var_post = self.var_post_grid[np.arange(self.var_post_grid.shape[0]), self.class_memberships]
		self.betas = np.random.normal(selected_mu_post, np.sqrt(selected_var_post))

		return

	def initialize_model_parameters(self):
		# Get number of snps
		self.KK = len(self.beta_hats)
		if len(self.beta_hats) != len(self.beta_hat_vars):
			print('assumpriton oeroorr')
			pdb.set_trace()

		# Initialize true betas
		self.betas = np.copy(self.beta_hats)

		# Set up sigma_sq grid
		self.sigma_sq_grid = get_sigma_sq_grid(self.m_factor, np.min(self.beta_hat_ses)/self.variance_grid_lb_divisor, 2.0*np.sqrt(np.max(np.square(self.beta_hats) - self.beta_hat_vars)))

		# Number of categories
		self.C = len(self.sigma_sq_grid)

		# Self class names
		self.classes = np.arange(self.C)

		# Initialize pis
		self.pis = np.ones(self.C)/self.C

		# Initialize hyperparameter on pis
		self.alpha_0 = np.ones(self.C)*self.alpha_0

		# Initialize class memberships
		self.class_memberships = np.random.choice(self.classes,size=self.KK)


		# Precompute some terms
		num = np.outer(self.beta_hat_vars, self.sigma_sq_grid)   # shape (K, C)
		den = self.beta_hat_vars[:, None] + self.sigma_sq_grid[None, :]
		self.var_post_grid = num / den
		# Compute posterior mean grid
		self.mu_post_grid = (self.var_post_grid / self.beta_hat_vars[:, np.newaxis])*self.beta_hats[:, np.newaxis] 
		# Compute log likelihoods of each class classes
		self.part_log_like =  - (.5*np.log(self.sigma_sq_grid/self.var_post_grid)) + (.5*(np.square(self.mu_post_grid)/(self.var_post_grid)))





		return

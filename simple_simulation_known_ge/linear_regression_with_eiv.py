import numpy as np 
import os
import sys
import pdb







class REGRESSION_EIV(object):
	def __init__(self, y=None, X=None, psi=None, max_iter=25000, burn_in_iter=20000):
		self.y = y
		self.X = X
		self.psi = psi
		self.max_iter = max_iter
		self.burn_in_iter = burn_in_iter

	def fit(self):
		self.initialize_parameters()
		self.sampled_coefs = []
		for itera in range(self.max_iter):
			self.update_coef()
			self.update_residual_variance()
			self.update_X_star()
			if itera > self.burn_in_iter:
				self.sampled_coefs.append(self.coef)
		self.sampled_coefs = np.asarray(self.sampled_coefs)


	def update_X_star(self):
		n_col = self.X_star.shape[1]

		for col_iter in range(n_col):
			valid_indices = self.psi[:,col_iter] != 0.0
			if np.sum(valid_indices) == 0.0:
				continue
			if np.sum(valid_indices) != len(self.X_star[:,col_iter]):
				print('assumption error: still need to implement')
				pdb.set_trace()
			# De-residualize current column
			self.resid = self.resid + self.X_star[:,col_iter]*self.coef[col_iter]

			varz = 1.0/((np.square(self.coef[col_iter])/self.residual_var) + (1.0/self.psi[:,col_iter]))
			pdb.set_trace()
			mus = varz*((self.resid*self.coef[col_iter]/self.residual_var) + (self.X[:,col_iter]/self.psi[:,col_iter]))


			self.X_star[:, col_iter] = np.random.normal(loc=mus,scale=np.sqrt(varz))

			# Re-residualize current column
			self.resid = self.resid - self.X_star[:,col_iter]*self.coef[col_iter]

		return
			

	def update_coef(self):
		S = np.linalg.inv(np.dot(np.transpose(self.X_star), self.X_star)/self.residual_var)
		mu = np.dot(np.dot(S, np.transpose(self.X_star)), self.y)/(self.residual_var)
		self.coef = np.random.multivariate_normal(mean=mu, cov=S)
		return

	def update_residual_variance(self):
		# Compute predicted y
		pred_y = np.dot(self.X_star,self.coef)

		self.resid = self.y-pred_y

		# Compute conditional sampling distribution
		temp_residual_variance_b = np.sum(np.square((self.y - pred_y)))/2.0
		temp_residual_variance_a = len(self.y)/2.0

		# Sample
		self.residual_var = 1.0/np.random.gamma(shape=temp_residual_variance_a, scale=1.0/temp_residual_variance_b)
		return


	def initialize_parameters(self):
		# X_star is the latent, true value of X
		self.X_star = np.copy(self.X)

		# coef is coefficients vector
		self.coef = np.zeros(self.X.shape[1])

		# Residual variance
		self.residual_var = 1.0


		return

import numpy as np
import pdb
import statsmodels.api as sm
from sklearn.linear_model import ARDRegression, BayesianRidge, LinearRegression
from scipy.stats import invgamma
import time
import scipy.special
import pyro
import torch
import pyro.distributions as dist
from dataclasses import dataclass
from pyro import poutine
from pyro.infer.autoguide import AutoDiagonalNormal, AutoGuideList, AutoDelta, AutoMultivariateNormal, init_to_value
from pyro.infer import SVI, Trace_ELBO, RenyiELBO, Predictive
from pyro.distributions import constraints

@dataclass
class Data:
	#A: torch.tensor
	X: torch.tensor
	y: torch.tensor

# Code from David Knowles
def robust_chol(cov, desired_min_eig = 1e-6):
	try: # faster than the eigendecomposition if already PSD
		chol_cov = torch.linalg.cholesky(cov)
	except:
		L, V = torch.linalg.eigh(cov)
		cov_min_eig = L.min().item()
		if cov_min_eig < desired_min_eig:
			print("Degenerate cov (min eigenvalue=%1.3e)" % cov_min_eig)
			# smallest addition to diagonal to make min(eig) = min_eig
			cov += (desired_min_eig - cov_min_eig) * torch.eye(cov.shape[0], device=cov.device)
		try:
			chol_cov = torch.linalg.cholesky(cov)
		except:
			print("Couldn't fix cov")
			return torch.eye(cov.shape[0], device=cov.device)
	return chol_cov





class LR_horseshoe_Prior(object):

	def __init__(self):
		return

	def fit(self, Y, G,num_samples=2000, iterations=4000):
		""" Fit the model.
			Args:
			G: A genotype matrix of floats with shape [num_samples, num_snps].
			Y: A trait vector vector of floats of length num_samples.
		"""
		device = torch.device("cpu")
		if str(device) == "cpu": torch.set_num_threads(1)


		# Add data to object
		# Initialization with sklearn bayesian ridge regression
		brr_init = BayesianRidge(compute_score=True, n_iter=100000, fit_intercept=False).fit(G, Y)

		# Get data in compact format
		data = Data(
		X=torch.tensor(G, device = device, dtype = torch.float),
		y=torch.tensor(Y, device = device, dtype = torch.float))



		device = data.X.device
		type_kwargs = {"device" : device, "dtype" : torch.float}
		print(device)

		N,P = data.X.shape
		one_vec = torch.ones(P, **type_kwargs)
		one = torch.tensor(1., **type_kwargs)

		# Quick fix if noise is really low
		brr_noise_scale = np.sqrt(.95)
		brr_weights_scale = np.sqrt(.05)

		one = torch.tensor(1., **type_kwargs)


		def convertr(hyperparam, name):
			return torch.tensor(hyperparam, device = device) if (
			  type(hyperparam) in [float,np.float32,torch.float]
			) else pyro.sample(name, hyperparam)

		def model(data):
			# could do e.g. dist.InverseGamma(2. * one,2. * one) instead
			sqrt_phi = pyro.sample("sqrt_phi", dist.HalfCauchy(one))
			noise_scale = pyro.sample("noise_scale", dist.HalfCauchy(one))
			sqrt_lambdas = pyro.sample("sqrt_lambdas", dist.HalfCauchy(one_vec).to_event(1))

			horseshoe_sigma = noise_scale**2*sqrt_lambdas**2
			Beta = pyro.sample('beta', dist.Normal(loc=0, scale=horseshoe_sigma).to_event(1))

			obs = pyro.sample('obs', dist.Normal(loc=(data.X @ Beta), scale=noise_scale**2).to_event(1), obs=data.y)
			# this could be done much more efficiently by precomputing eig(data.X @ data.X.T)
			# and just updating the eigenvalues as sqrt_phi and noise_scale change
			#cov = sqrt_phi**2 * (data.X @ torch.diag(sqrt_lambdas**2) @ data.X.T) + noise_scale**2 * torch.eye(N, **type_kwargs)
			#cov = sqrt_phi**2 * (X_X_t) + noise_scale**2 * torch.eye(N, **type_kwargs)
			#L = robust_chol(cov)
			#obs = pyro.sample('obs', dist.MultivariateNormal(torch.zeros(N, **type_kwargs), scale_tril = L), obs = data.y)

		init_dic = {
	  	"sqrt_psi" : torch.tensor(brr_weights_scale, **type_kwargs),
	  	"noise_scale" : torch.tensor( brr_noise_scale, **type_kwargs),
		} if (brr_init != None) else {}

		guide = AutoMultivariateNormal(
			model,
			init_loc_fn = init_to_value(values=init_dic))


		adam = pyro.optim.Adam({"lr": 0.03})
		svi = SVI(model, guide, adam, loss=Trace_ELBO())
		pyro.clear_param_store()
		losses = []
		time_per_iter = []
		for j in range(iterations):
			print(j, end = '\r')
			loss = svi.step(data)
			losses.append(loss)

		pdb.set_trace()

		samples = Predictive(
			model,
			guide=guide,
			return_sites = ["sqrt_phi", "noise_scale"],
			num_samples=num_samples)(data)






		pdb.set_trace()
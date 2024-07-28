import sys
import pdb
import numpy as np
import os
from scipy.special import expit

import jax.numpy as jnp
import jax.random as random

import numpyro
from numpyro.diagnostics import summary
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS
import time
import statsmodels.api as sm


def model_normal_likelihood_gaussian_prior(X, Y):
	D_X = X.shape[1]

	# sample from horseshoe prior
	#lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	tau = numpyro.sample("tau", dist.HalfCauchy(jnp.ones(1)))


	with numpyro.plate("plate_i", D_X):
		Beta = numpyro.sample('beta', dist.Normal(loc=0, scale=tau))

	noise_scale = numpyro.sample("noise_scale", dist.HalfCauchy(jnp.ones(1)))


	mu = jnp.dot(X, Beta)
	numpyro.sample('obs', dist.Normal(loc=mu, scale=noise_scale), obs=Y)

	return

# regression model with continuous-valued outputs/responses
def model_normal_likelihood_horseshoe_prior(X, Y):
	D_X = X.shape[1]

	# sample from horseshoe prior
	#lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	tau = numpyro.sample("tau", dist.HalfCauchy(jnp.ones(1)))

	horseshoe_sigma = tau**2*lambdas**2
	Beta = numpyro.sample('beta', dist.Normal(loc=0, scale=(horseshoe_sigma**.5)))
	#prec_obs = numpyro.sample("prec_obs", dist.Gamma(3.0, 1.0))
	#sigma_obs = 1.0 / jnp.sqrt(prec_obs)

	noise_scale = numpyro.sample("noise_scale", dist.HalfCauchy(jnp.ones(1)))


	mu = jnp.dot(X, Beta)
	numpyro.sample('obs', dist.Normal(loc=mu, scale=noise_scale), obs=Y)

	return

def model_normal_likelihood_horseshoe_prior_fast(X,Y):
	D_X = X.shape[1]

	# sample from horseshoe prior
	#lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	tau = numpyro.sample("tau", dist.HalfCauchy(jnp.ones(1)))

	unscaled_betas = numpyro.sample("unscaled_betas", dist.Normal(0.0, jnp.ones(D_X)))

	scaled_betas = numpyro.deterministic("beta", tau * lambdas * unscaled_betas)

	noise_scale = numpyro.sample("noise_scale", dist.HalfCauchy(jnp.ones(1)))

	mu = jnp.dot(X, scaled_betas)
	numpyro.sample('obs', dist.Normal(loc=mu, scale=noise_scale), obs=Y)

	return

def model_normal_likelihood_regularized_horseshoe_prior_fast(X,Y):
	D_X = X.shape[1]
	NN = X.shape[0]


	noise_scale = numpyro.sample("noise_scale", dist.HalfCauchy(jnp.ones(1)))

	# sample from horseshoe prior
	#lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	tau = numpyro.sample("tau", dist.HalfCauchy(jnp.ones(1)))
	#tau = numpyro.sample("tau", dist.HalfNormal(0.05*noise_scale/(NN**.5)))

	c2 = numpyro.sample("c2", dist.InverseGamma(jnp.ones(1), jnp.ones(1)))

	lambdas_tilde = lambdas*((c2/(c2 + (tau**2 * lambdas**2)))**.5)




	unscaled_betas = numpyro.sample("unscaled_betas", dist.Normal(0.0, jnp.ones(D_X)))

	scaled_betas = numpyro.deterministic("beta", tau * lambdas_tilde * unscaled_betas)


	mu = jnp.dot(X, scaled_betas)
	numpyro.sample('obs', dist.Normal(loc=mu, scale=noise_scale), obs=Y)

	return


# regression model with continuous-valued outputs/responses
def model_normal_likelihood_dirichlet_horseshoe_prior(X, Y):
	D_X = X.shape[1]

	# sample from horseshoe prior
	#lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	beta_0 = numpyro.sample('beta_0', dist.Normal(loc=0, scale=1000))
	lambdas = numpyro.sample("lambdas", dist.Dirichlet(100*jnp.ones(D_X)))
	tau = numpyro.sample("tau", dist.HalfCauchy(jnp.ones(1)))

	horseshoe_sigma = tau**2*lambdas
	Beta = numpyro.sample('beta', dist.Normal(loc=0, scale=(horseshoe_sigma**.5)))
	#prec_obs = numpyro.sample("prec_obs", dist.Gamma(3.0, 1.0))
	#sigma_obs = 1.0 / jnp.sqrt(prec_obs)

	noise_scale = numpyro.sample("noise_scale", dist.HalfCauchy(jnp.ones(1)))
	#noise_scale = pyro.sample("noise_scale", dist.HalfNormal(scale=10.0))


	mu = jnp.dot(X, Beta) + beta_0
	numpyro.sample('obs', dist.Normal(loc=mu, scale=noise_scale), obs=Y)

	return


# regression model with continuous-valued outputs/responses
def model_normal_likelihood_dirichlet_horseshoe_prior_fast(X, Y):
	D_X = X.shape[1]

	# sample from horseshoe prior
	#lambdas = numpyro.sample("lambdas", dist.HalfCauchy(jnp.ones(D_X)))
	beta_0 = numpyro.sample('beta_0', dist.Normal(loc=0, scale=100000))

	lambdas_sq = numpyro.sample("lambdas", dist.Dirichlet(1.0*jnp.ones(D_X)))
	tau = numpyro.sample("tau", dist.HalfCauchy(jnp.ones(1)))
	noise_scale = numpyro.sample("noise_scale", dist.HalfCauchy(jnp.ones(1)))

	#horseshoe_sigma = tau**2*lambdas
	#Beta = numpyro.sample('beta', dist.Normal(loc=0, scale=(horseshoe_sigma**.5)))
	#prec_obs = numpyro.sample("prec_obs", dist.Gamma(3.0, 1.0))
	#sigma_obs = 1.0 / jnp.sqrt(prec_obs)
	unscaled_betas = numpyro.sample("unscaled_betas", dist.Normal(0.0, jnp.ones(D_X)))

	scaled_betas = numpyro.deterministic("beta", tau * (D_X**.5) * (lambdas_sq**.5) * unscaled_betas)



	mu = jnp.dot(X, scaled_betas) + beta_0
	numpyro.sample('obs', dist.Normal(loc=mu, scale=noise_scale), obs=Y)

	return



# helper function for HMC inference
def run_inference(model, rng_key, X, Y, num_warmup=3000, num_samples=8000, thinning=5):
	start = time.time()
	kernel = NUTS(model)
	mcmc = MCMC(
		kernel,
		num_warmup=num_warmup,
		num_samples=num_samples,
		thinning=thinning,
		num_chains=1,
		progress_bar=False if "NUMPYRO_SPHINXBUILD" in os.environ else True,
	)

	mcmc.run(rng_key, X, Y)
	mcmc.print_summary(exclude_deterministic=False)

	samples = mcmc.get_samples()
	#summary_dict = summary(samples, group_by_chain=False)

	print("\nMCMC elapsed time:", time.time() - start)

	return samples



eqtl_arch = sys.argv[1]
calibration_output_file = sys.argv[2]

# Numpyro sampling params
num_warmup=2000
num_samples=3000
thinning=5

t = open(calibration_output_file,'w')
t.write('sim_iter\tmethod_name\tcalibration\tsquared_calibration\n')


# Sim data 
NN = 300
PP = 100
NN2 = 10000
n_causal = PP*0.05
expr_h2 = 0.6

n_sims = 100

models = ['numpyro_lmm']
models = ['numpyro_horseshoe']
models = ['numpyro_regularized_horseshoe']
models = ['numpyro_dirichlet_horseshoe']


for sim_iter in range(n_sims):
	# Sim geno 
	G = []
	for pp in range(PP):
		geno = np.random.normal(loc=0,scale=1,size=NN)
		G.append((geno-np.mean(geno))/np.std(geno))
	G = np.transpose(np.asarray(G))

	# Sim geno 
	G2 = []
	for pp in range(PP):
		geno = np.random.normal(loc=0,scale=1,size=NN2)
		G2.append((geno-np.mean(geno))/np.std(geno))
	G2 = np.transpose(np.asarray(G2))


	beta_sim = np.zeros(PP)
	if eqtl_arch == 'sparse':
		causal_snp_indices = np.random.choice(np.arange(PP), size=int(n_causal), replace=False)
		beta_sim[causal_snp_indices] = np.random.normal(loc=0, scale=np.sqrt(expr_h2/n_causal), size=int(n_causal))
	elif eqtl_arch == 'infinitesimal':
		beta_sim = np.random.normal(loc=0, scale=np.sqrt(expr_h2/PP), size=int(PP))



	genetic_trait = np.dot(G, beta_sim)	
	vary = np.var(genetic_trait)
	if vary > 1:
		vary = .95
	Y = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0 - vary))
	genetic_trait2 = np.dot(G2, beta_sim)
	vary2 = np.var(genetic_trait2)
	if vary2 > 1:
		vary2 = .95
	Y2 = np.random.normal(loc=genetic_trait2, scale=np.sqrt(1.0 - vary2))

	for model_name in models:
		try:
			if model_name == 'numpyro_lmm':
				rng_key, rng_key_predict = random.split(random.PRNGKey(0))
				samples = run_inference(model_normal_likelihood_gaussian_prior, rng_key, G, Y, num_warmup=num_warmup, num_samples=num_samples, thinning=thinning)
				beta_mean = np.mean(samples['beta'],axis=0)
				beta_var = np.var(samples['beta'],axis=0)

			elif model_name == 'numpyro_horseshoe':
				rng_key, rng_key_predict = random.split(random.PRNGKey(0))
				samples = run_inference(model_normal_likelihood_horseshoe_prior_fast, rng_key, G, Y, num_warmup=num_warmup, num_samples=num_samples, thinning=thinning)
				#samples2 = run_inference(model_normal_likelihood_horseshoe_prior, rng_key, G, Y, num_warmup=num_warmup, num_samples=num_samples, thinning=thinning)

				beta_mean = np.asarray(np.mean(samples['beta'],axis=0))
				beta_var = np.asarray(np.var(samples['beta'],axis=0))
				'''
				beta_mean2 = np.asarray(np.mean(samples2['beta'],axis=0))
				beta_var2 = np.asarray(np.var(samples2['beta'],axis=0))

				aa = np.corrcoef(beta_mean, beta_mean2)
				bb = np.corrcoef(np.square(beta_mean) + beta_var, np.square(beta_mean2) + beta_var2)
				cc = sm.OLS(beta_mean, beta_mean2).fit().params
				dd = sm.OLS(np.square(beta_mean) + beta_var, np.square(beta_mean2) + beta_var2).fit().params

				pdb.set_trace()
				'''
			elif model_name == 'numpyro_regularized_horseshoe':
				rng_key, rng_key_predict = random.split(random.PRNGKey(0))
				samples = run_inference(model_normal_likelihood_regularized_horseshoe_prior_fast, rng_key, G, Y, num_warmup=num_warmup, num_samples=num_samples, thinning=thinning)
				#samples2 = run_inference(model_normal_likelihood_horseshoe_prior, rng_key, G, Y, num_warmup=num_warmup, num_samples=num_samples, thinning=thinning)

				beta_mean = np.asarray(np.mean(samples['beta'],axis=0))
				beta_var = np.asarray(np.var(samples['beta'],axis=0))
			elif model_name == 'numpyro_dirichlet_horseshoe':
				rng_key, rng_key_predict = random.split(random.PRNGKey(0))
				samples = run_inference(model_normal_likelihood_dirichlet_horseshoe_prior_fast, rng_key, G, Y, num_warmup=num_warmup, num_samples=num_samples, thinning=thinning)
				#samples2 = run_inference(model_normal_likelihood_horseshoe_prior, rng_key, G, Y, num_warmup=num_warmup, num_samples=num_samples, thinning=thinning)

				beta_mean = np.asarray(np.mean(samples['beta'],axis=0))
				beta_var = np.asarray(np.var(samples['beta'],axis=0))


				posterior_squared = np.asarray(np.mean(samples['lambdas'],axis=0))

				n_samp = samples['lambdas'].shape[0]
				calibration_distr = []
				for samp_iter in range(n_samp):
					olser = sm.OLS(np.square(beta_sim)/np.sum(np.square(beta_sim)), np.asarray(samples['lambdas'][samp_iter,:])).fit()
					calibration_distr.append(olser.params[0])
				calibration_distr = np.asarray(calibration_distr)

				
				lb_iter = int(np.floor(n_samp*.05))
				ub_iter = int(np.floor(n_samp*.95))
				print(str(np.mean(calibration_distr)) + ' : [' + str(np.sort(calibration_distr)[lb_iter]) + ', ' + str(np.sort(calibration_distr)[ub_iter]) + ']')

				pdb.set_trace()



			# Assess calibration
			#olser = sm.OLS(Y2, sm.add_constant(np.dot(G2, beta_mean))).fit()
			#calibration_param = olser.params[1]
			olser = sm.OLS(beta_sim, np.asarray(beta_mean)).fit()
			calibration_param = olser.params[0]

			# Assess squared calibration
			#posterior_squared = np.square(np.asarray(beta_mean)) + np.asarray(beta_var)
			#posterior_squared = posterior_squared/np.sum(posterior_squared)
			#print(np.sort(posterior_squared))
			olser = sm.OLS(np.square(beta_sim)/np.sum(np.square(beta_sim)), posterior_squared).fit()
			squared_calibration_param = olser.params[0]

			# Assess squared calibration

			# Print to output file
			t.write(str(sim_iter) + '\t' + str(model_name) + '\t' + str(calibration_param) + '\t' + str(squared_calibration_param) + '\n')
			print(str(calibration_param) + '\t' + str(squared_calibration_param))

			t.flush()

		except:
			print('ERROR: skip sim-model pair')
			samples = run_inference(model_normal_likelihood_dirichlet_horseshoe_prior_fast, rng_key, G, Y, num_warmup=num_warmup, num_samples=num_samples, thinning=thinning)






'''
# do inference
rng_key, rng_key_predict = random.split(random.PRNGKey(0))
summary, samples = run_inference(model_normal_likelihood_horseshoe_prior, rng_key, G, Y)

beta_mean = np.mean(samples['beta'],axis=0)
beta_var = np.var(samples['beta'],axis=0)

pdb.set_trace()

olser = sm.OLS(Y2, sm.add_constant(np.dot(G2, beta_mean))).fit()
print(olser.params[1])

#posterior_squared = np.square(beta_mean) + beta_var

#pi1 = posterior_squared/np.sum(posterior_squared)
#pi2 = np.mean(samples['lambdas'],axis=0)

mod = linear_regression_mixture_of_gaussian_grid_prior_gibbs.LR_MoG_Grid_Prior()
mod.fit(Y, G)


mod2 = linear_regression_spike_and_slab_prior_gibbs.LR_SaS_Prior()
mod2.fit(Y, G)

pdb.set_trace()
'''

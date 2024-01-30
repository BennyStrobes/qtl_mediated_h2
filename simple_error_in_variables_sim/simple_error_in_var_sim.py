import numpy as np 
import os
import sys
import pdb
import stan
from sklearn.linear_model import LinearRegression


eiv_reg_code = """
data {
  int<lower=0> N;      // number of cases
  array[N] real x_meas; 
  real<lower=0> tau;   // measurement noise
  array[N] real y;
}
parameters {
  array[N] real x;
  real beta;           // slope
  real<lower=0> sigma; // outcome noise
  real<lower=0> sigma_x;
}
model {
  x ~ normal(0, sigma_x);            // prior on x
  for (i in 1:N){
  	x_meas[i] ~ normal(x[i], tau);
  	y[i] ~ normal(beta*x[i], sigma);
  }
}
"""








# Simulation parameters
num_samples=300
n_sims = 40
coef_effect_size = .45
predictor_variance = 1.0
residual_variance = 1.0



learned_coefs = []
learned_noisy_coefs = []
learned_adj_noisy_coefs = []
for sim_iter in range(n_sims):
	# Simulate input
	X = np.random.normal(size=num_samples)
	X_hat = np.random.normal(X, scale=np.sqrt(predictor_variance))
	Y = np.random.normal(X*coef_effect_size, scale=np.sqrt(residual_variance))

	# Run regressions true, unobserved predictors
	reg = LinearRegression().fit(X.reshape(-1,1),Y)
	learned_coefs.append(reg.coef_)

	# Run regressions with noisy predictors 
	reg = LinearRegression().fit(X_hat.reshape(-1,1),Y)
	learned_noisy_coefs.append(reg.coef_)
	print(reg.coef_)

	# STAN EIV MODEL
	stan_data = {"N": len(Y), "y": Y, "x_meas": X_hat, "tau": predictor_variance}
	posterior = stan.build(eiv_reg_code, data=stan_data)
	fit = posterior.sample(num_chains=4, num_samples=5000)
	df = fit.to_frame() 
	sampled_betas = fit['beta']
	learned_adj_noisy_coefs.append(np.mean(sampled_betas))
	print(np.mean(sampled_betas))




learned_coefs = np.asarray(learned_coefs)
learned_noisy_coefs = np.asarray(learned_noisy_coefs)
learned_adj_noisy_coefs = np.asarray(learned_adj_noisy_coefs)

print(np.mean(learned_coefs))
print(np.mean(learned_noisy_coefs))
print(np.mean(learned_adj_noisy_coefs))

pdb.set_trace()

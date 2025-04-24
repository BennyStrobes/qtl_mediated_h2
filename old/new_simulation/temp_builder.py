import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np
import os
import pdb
import statsmodels.api as sm
from sklearn.linear_model import ARDRegression, BayesianRidge, LinearRegression
import linear_regression_spike_and_slab_prior_gibbs
import linear_regression_mixture_of_gaussian_grid_prior_gibbs
import linear_regression_horseshoe_prior_pyro











# Sim data 
NN = 300
PP = 400
n_causal = PP*0.025

per_snp_h2 = 0.005


# Sim geno 
G = []
for pp in range(PP):
	geno = np.random.normal(loc=0,scale=1,size=NN)
	G.append((geno-np.mean(geno))/np.std(geno))
G = np.transpose(np.asarray(G))

beta_sim = np.zeros(PP)
causal_snp_indices = np.random.choice(np.arange(PP), size=int(n_causal), replace=False)
beta_sim[causal_snp_indices] = np.random.normal(loc=0, scale=np.sqrt(per_snp_h2), size=int(n_causal))

genetic_trait = np.dot(G, beta_sim)
print(np.var(genetic_trait))
Y = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0 - np.var(genetic_trait)))
Y = (Y - np.mean(Y))/np.std(Y)


mod = linear_regression_horseshoe_prior_pyro.LR_horseshoe_Prior()
mod.fit(Y, G)


'''
mod = linear_regression_mixture_of_gaussian_grid_prior_gibbs.LR_MoG_Grid_Prior()
mod.fit(Y, G)


mod2 = linear_regression_spike_and_slab_prior_gibbs.LR_SaS_Prior()
mod2.fit(Y, G)
'''

pdb.set_trace()

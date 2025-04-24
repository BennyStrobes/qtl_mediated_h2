import numpy as np
import os
import sys
import pdb
import statsmodels.api as sm







def simulate_trait(geno, snp_heritability, n_causal_snps):
	n_snps = geno.shape[0]

	beta = np.zeros(n_snps)


	causal_snp_indices = np.random.choice(np.arange(n_snps), size=n_causal_snps, replace=False)

	beta[causal_snp_indices] = np.random.normal(loc=0, scale=np.sqrt(snp_heritability/n_causal_snps), size=n_causal_snps)

	genetic_trait = np.dot(beta, geno)

	trait_h2 = np.var(genetic_trait)

	print(trait_h2)

	trait = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0-trait_h2))

	trait = (trait - np.mean(trait))/np.std(trait)
	return trait


def get_marginal_gwas_z_scores(trait, geno):
	n_snps = geno.shape[0]
	zeds = []
	for snp_iter in range(n_snps):

		model = sm.OLS(trait,sm.add_constant(geno[snp_iter,:])).fit()

		zeds.append(model.params[1]/model.bse[1])

	return np.asarray(zeds)


###########
# start

snp_heritability = .3
n_causal_snps = 500


geno = np.load('genos.npy')
ld = np.load('ld.npy')

for itera in range(10):
	trait = simulate_trait(geno, snp_heritability, n_causal_snps)

	zeds = get_marginal_gwas_z_scores(trait, geno)
	#ld = np.corrcoef(geno)

	ldscores = np.sum(np.square(ld),axis=0)

	model = sm.OLS(np.square(zeds),sm.add_constant(ldscores)).fit()

	print(model.params)


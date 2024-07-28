import sys
import numpy as np 
import pandas as pd
import os
import pdb
import tensorflow as tf
import gzip
import time
#from scipy.stats import gamma
#os.environ['OPENBLAS_NUM_THREADS'] = '20'
#os.environ['MKL_NUM_THREADS'] = '20'
#from joblib.externals.loky import set_loky_pickler
#from joblib import wrap_non_picklable_objects
from tensorflow.keras.layers import Conv1D, Dense, BatchNormalization, Activation, GlobalAveragePooling1D, Add, Input, Flatten
from tensorflow.keras import Model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from scipy import odr
import tensorflow_probability as tfp
from tensorflow_probability.python.math import bessel
from tensorflow_probability.python.distributions import noncentral_chi2
from scipy import stats
tfd = tfp.distributions
from tensorflow_probability.python.math import generic as tfp_math



def simulate_data(gwas_ss, n_snps, eqtl_ss_1, med_h2_1, ge_h2_1, n_genes, snps_per_gene):
	sim_alpha_1_sq = med_h2_1/ge_h2_1

	gene_midpoint_lb = int(snps_per_gene/2) +1
	gene_midpoint_ub = n_snps - (int(snps_per_gene/2) +1)

	gene_causal_eqtl_effects = []
	gene_snp_indices_arr = []
	for gene_iter in range(n_genes):

		gene_midpoint = np.random.choice(np.arange(gene_midpoint_lb, gene_midpoint_ub))
		gene_start = gene_midpoint - int(snps_per_gene/2)
		gene_end = gene_midpoint + int(snps_per_gene/2)
		gene_snp_indices = np.arange(gene_start, gene_end)

		gene_snp_causal_eqtl_effects = np.zeros(n_snps)
		gene_snp_causal_eqtl_effects[gene_snp_indices] = np.random.normal(loc=0, scale=np.sqrt(ge_h2_1/snps_per_gene), size=snps_per_gene)

		gene_causal_eqtl_effects.append(gene_snp_causal_eqtl_effects)
		gene_snp_indices_arr.append(gene_snp_indices)
	gene_causal_eqtl_effects = np.asarray(gene_causal_eqtl_effects)

	# Simulate eQTL genotype data
	eqtl_geno_1 = np.random.normal(loc=0,scale=1.0,size=(eqtl_ss_1,n_snps))
	for snp_iter in range(n_snps):
		eqtl_geno_1[:,snp_iter] = (eqtl_geno_1[:,snp_iter] - np.mean(eqtl_geno_1[:,snp_iter]))/np.std(eqtl_geno_1[:,snp_iter])


	# Simulate GWAS genotype data
	gwas_geno = np.random.normal(loc=0,scale=1.0,size=(gwas_ss,n_snps))
	for snp_iter in range(n_snps):
		gwas_geno[:,snp_iter] = (gwas_geno[:,snp_iter] - np.mean(gwas_geno[:,snp_iter]))/np.std(gwas_geno[:,snp_iter])

	# Simulate Gene trait value
	Es = []
	for gene_iter in range(n_genes):
		genetic_gene_trait_1 = np.dot(eqtl_geno_1, gene_causal_eqtl_effects[gene_iter, :])
		E_1 = np.random.normal(loc=genetic_gene_trait_1, scale=np.sqrt(1.0-np.var(genetic_gene_trait_1)))
		E_1 = (E_1 - np.mean(E_1))/np.std(E_1)
		Es.append(E_1)

	# Simulate gwas trait value
	sim_alphas = np.random.normal(loc=0, scale=np.sqrt(sim_alpha_1_sq/n_genes), size=n_genes)
	genetic_trait = np.zeros(gwas_ss)
	for gene_iter in range(n_genes):
		genetic_gene = np.dot(gwas_geno, gene_causal_eqtl_effects[gene_iter, :])
		genetic_trait =genetic_trait + genetic_gene*sim_alphas[gene_iter]
	Y = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0-np.var(genetic_trait)))
	Y = (Y - np.mean(Y))/np.std(Y)

	return sim_alpha_1_sq, Y, Es, gwas_geno, eqtl_geno_1, gene_causal_eqtl_effects, gene_snp_indices_arr


def get_marginal_summary_statistics_no_intercept(Y, G):
	n_snps = G.shape[1]
	beta = []
	beta_se = []
	for snp_iter in range(n_snps):
		olser = sm.OLS(Y, G[:,snp_iter]).fit()
		beta.append(olser.params[0])
		beta_se.append(olser.bse[0])
	return np.asarray(beta), np.asarray(beta_se)

def get_marginal_summary_statistics(Y, G):
	n_snps = G.shape[1]
	beta = []
	beta_se = []
	for snp_iter in range(n_snps):
		olser = sm.OLS(Y, sm.add_constant(G[:,snp_iter])).fit()
		beta.append(olser.params[1])
		beta_se.append(olser.bse[1])
	return np.asarray(beta), np.asarray(beta_se)

def get_marginal_summary_statistics_across_genes(Es, G, gene_snp_indices_arr):
	n_genes = len(Es)
	n_snps = E_geno.shape[1]

	beta_arr = []
	beta_se_arr = []
	# Loop through genes
	for gene_iter in range(n_genes):
		beta = np.zeros(n_snps)
		beta_se = np.ones(n_snps)*1e-7
		gene_snp_indices = gene_snp_indices_arr[gene_iter]

		for gene_snp_index in gene_snp_indices:
			olser = sm.OLS(Es[gene_iter], sm.add_constant(G[:,gene_snp_index])).fit()
			beta[gene_snp_index] = olser.params[1]
			beta_se[gene_snp_index] = olser.bse[1]

		# Add to global array
		beta_arr.append(beta)
		beta_se_arr.append(beta_se)

	beta_arr = np.asarray(beta_arr)
	beta_se_arr = np.asarray(beta_se_arr)

	return beta_arr, beta_se_arr







#########################
# Command line args
#########################
global_simulation_number = sys.argv[1]
n_sims = int(sys.argv[2])
gwas_ss = int(sys.argv[3])
n_snps = int(sys.argv[4])
eqtl_ss_1 = int(sys.argv[5])
med_h2_1 = float(sys.argv[6])
ge_h2_1 = float(sys.argv[7])
output_root = sys.argv[8]
global_simulation_number = int(sys.argv[9])
n_genes = int(sys.argv[10])
snps_per_gene = int(sys.argv[11])


# Set seed
np.random.seed(global_simulation_number+1)


# Open and print header to output file
output_file = output_root + '_effect_est_res_summary.txt'
t = open(output_file,'w')
t.write('method\tsim\tsim_alpha_1\test_alpha_1_no_noise\test_alpha_1\tloss\n')


# Loop through sims
for sim_iter in range(n_sims):
	print(sim_iter)

	# First simulate data
	sim_alpha_1_sq, Y, Es, Y_geno, E_geno, gene_causal_eqtl_effects, gene_snp_indices_arr = simulate_data(gwas_ss, n_snps, eqtl_ss_1, med_h2_1, ge_h2_1, n_genes, snps_per_gene)


	# Next get summary statistics
	gwas_beta, gwas_beta_se = get_marginal_summary_statistics(Y, Y_geno)
	eqtl_beta, eqtl_beta_se = get_marginal_summary_statistics_across_genes(Es, E_geno, gene_snp_indices_arr)


	pdb.set_trace()

	# Print to output
	#t.write(method + '\t' + str(sim_iter) + '\t' + str(sim_alpha_1_sq) + '\t' + str(est_squared_alphas_custom_tf_no_noise) + '\t' + str(est_squared_alphas_custom_tf) + '\n')

t.close()


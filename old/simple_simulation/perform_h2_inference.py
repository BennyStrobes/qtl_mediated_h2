import numpy as np 
import os
import sys
import pdb
import rss_vi_variant_only
import rss_gibbs_variant_only
from sklearn.linear_model import LinearRegression










def load_in_gwas_z_scores(gwas_sum_stats_file):
	f = open(gwas_sum_stats_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(float(data[3]))
	f.close()
	arr = np.asarray(arr)
	return arr

def load_in_gwas_se(gwas_sum_stats_file):
	f = open(gwas_sum_stats_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(float(data[2]))
	f.close()
	arr = np.asarray(arr)
	return arr

def load_in_trait_file(trait_file):
	f = open(trait_file)
	arr = []
	for line in f:
		line = line.rstrip()
		arr.append(float(line))
	f.close()
	return np.asarray(arr)

def estimate_heritability_with_ldsc(LD, marginal_z, N_gwas):
	n_snps = LD.shape[1]

	chi_sq = np.square(marginal_z)

	squared_ld = np.square(LD)

	adj_squared_ld = squared_ld - ((1.0 - squared_ld)/(N_gwas - 2))

	ld_scores = np.sum(adj_squared_ld,axis=0)

	model = LinearRegression(fit_intercept=True)

	ldsc_fit = model.fit(ld_scores.reshape((len(ld_scores),1)), chi_sq)

	ldsc_gene_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas

	return ldsc_gene_h2

def estimate_heritability_with_rough_ldsc(LD, marginal_z, N_gwas):
	n_snps = LD.shape[1]

	chi_sq = np.square(marginal_z)

	squared_ld = np.square(LD)

	adj_squared_ld = squared_ld - ((1.0 - squared_ld)/(N_gwas - 2))

	avg_ld_scores = np.mean(np.sum(adj_squared_ld,axis=0))

	rough_ldsc_gene_h2 = (np.mean(chi_sq) - 1.0)*n_snps/(N_gwas*avg_ld_scores)

	return rough_ldsc_gene_h2




#########################
# Command line args
#########################
simulation_name_string = sys.argv[1]
simulated_gwas_data_dir = sys.argv[2]
simulated_gene_models_dir = sys.argv[3]
mediated_h2_results_dir = sys.argv[4]
processed_genotype_data_dir = sys.argv[5]


# File containing gwas summary statistics
gwas_sum_stats_file = simulated_gwas_data_dir + simulation_name_string + 'simulated_gwas_summary_stats.txt'
gwas_z_scores = load_in_gwas_z_scores(gwas_sum_stats_file)
gwas_se = load_in_gwas_se(gwas_sum_stats_file)
gwas_beta = gwas_z_scores*gwas_se

# Simulated total genetic var
simulated_genetic_var_file = simulated_gwas_data_dir + simulation_name_string + 'simulated_genetic_var.npy'
sim_genetic_var = np.load(simulated_genetic_var_file) + 0.0
print(sim_genetic_var)


# Simulated Trait
trait_file = simulated_gwas_data_dir + simulation_name_string + 'simulated_trait.txt'
trait_vec = load_in_trait_file(trait_file)

# Load in Genotype file
# LD file
ld_file = processed_genotype_data_dir + 'gwas_genotype_LD_1.npy'
LD = np.load(ld_file)

# Genotype file
geno_file = processed_genotype_data_dir + 'gwas_genotype_1.npy'
genotype_mat = np.load(geno_file)

# Gwas sample size
N_gwas = len(trait_vec)

# Perform variant heritability analysis with LDSC
ldsc_h2 = estimate_heritability_with_ldsc(LD, gwas_z_scores, N_gwas)

# Perform variant heritability analysis with rough LDSC
rough_ldsc_h2 = estimate_heritability_with_rough_ldsc(LD, gwas_z_scores, N_gwas)


# Perform variant heritability analysis by optimizing RSS likelihood with Gibbs sampling
#mod = rss_gibbs_variant_only.RSS_GIBBS_VARIANT_ONLY(LD=LD, marginal_beta=gwas_beta, N=N_gwas, X=genotype_mat, Y=trait_vec, max_iter=150, burn_in_iter=100)
mod = rss_gibbs_variant_only.RSS_GIBBS_VARIANT_ONLY(LD=LD, marginal_beta=gwas_beta, N=N_gwas, X=genotype_mat, Y=trait_vec, max_iter=3, burn_in_iter=0)
mod.fit()


# Print results to output
output_file = mediated_h2_results_dir + simulation_name_string + 'h2_variant_only_summary.txt'
t = open(output_file,'w')
t.write('simulated_heritability\testimated_heritability_rss_gibbs\testimated_heritability_rss_gibbs_se\testimated_heritability_ldsc\testimated_heritability_rough_ldsc\n')
t.write(str(sim_genetic_var) + '\t' + str(np.mean(mod.sampled_local_h2s)) + '\t' + str(np.std(mod.sampled_local_h2s)) + '\t' + str(ldsc_h2) + '\t' + str(rough_ldsc_h2) + '\n')
t.close()

print(output_file)








import numpy as np 
import os
import sys
import pdb
import rss_vi_variant_only
import rss_gibbs_variant_only
import rss_gibbs_variant_pred_gene
import rss_gibbs_variant_modeled_gene
import trait_likelihood_gibbs_variant_only
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

def extract_previously_estimsted_causal_eqtl_effects(gene_summary_file):
	causal_eqtl_effects = []
	causal_eqtl_indices = []
	head_count = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
		estimated_causal_eqtl_file = data[10]
		estimated_causal_eqtl_effects = np.load(estimated_causal_eqtl_file)

		# Quick error check
		if len(estimated_causal_eqtl_effects) != len(gene_snp_indices):
			continue

		causal_eqtl_effects.append(estimated_causal_eqtl_effects)
		causal_eqtl_indices.append(gene_snp_indices)

	f.close()
	return causal_eqtl_effects, causal_eqtl_indices


#########################
# Command line args
#########################
simulation_name_string = sys.argv[1]
simulated_gwas_data_dir = sys.argv[2]
simulated_gene_models_dir = sys.argv[3]
eqtl_ss = sys.argv[4]
mediated_h2_results_dir = sys.argv[5]
processed_genotype_data_dir = sys.argv[6]
simulated_expression_data_dir = sys.argv[7]




# File summarizing inferred gene models
gene_summary_file = simulated_gene_models_dir + simulation_name_string + 'model_summaries_' + eqtl_ss + '.txt'

# Extract previously estimated causal eqtl effects
causal_eqtl_effects, causal_eqtl_indices = extract_previously_estimsted_causal_eqtl_effects(gene_summary_file)

# Extract gene expression
gene_expression_file = simulated_expression_data_dir + simulation_name_string + 'eqtl_ss_' + eqtl_ss + '.npy'
gene_expression_mat = np.load(gene_expression_file)
n_genes = gene_expression_mat.shape[0]
gene_expression = []
for gg in range(n_genes):
	gene_expression.append(gene_expression_mat[gg,:])


# File containing gwas summary statistics
gwas_sum_stats_file = simulated_gwas_data_dir + simulation_name_string + 'simulated_gwas_summary_stats.txt'
gwas_z_scores = load_in_gwas_z_scores(gwas_sum_stats_file)
gwas_se = load_in_gwas_se(gwas_sum_stats_file)
gwas_beta = gwas_z_scores*gwas_se

# Simulated total genetic vars
simulated_genetic_var_file = simulated_gwas_data_dir + simulation_name_string + 'simulated_genetic_var.npy'
sim_genetic_var = np.load(simulated_genetic_var_file) + 0.0
simulated_nm_genetic_var_file = simulated_gwas_data_dir + simulation_name_string + 'simulated_nm_genetic_var.npy'
sim_nm_genetic_var = np.load(simulated_nm_genetic_var_file) + 0.0
simulated_mediated_genetic_var_file = simulated_gwas_data_dir + simulation_name_string + 'simulated_mediated_genetic_var.npy'
sim_med_genetic_var = np.load(simulated_mediated_genetic_var_file) + 0.0

# Print what simulated total genetic vars are
print(sim_genetic_var)
print(sim_nm_genetic_var)
print(sim_med_genetic_var)


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

# eqtl genotype file
eqtl_geno_file = processed_genotype_data_dir + 'eqtl_' + str(eqtl_ss) + '_genotype_1.npy'
eqtl_genotype_mat = np.load(eqtl_geno_file)

# Parse eqtl genotypes by gene
gene_genotypes = []
for gg in range(len(causal_eqtl_indices)):
	gene_genotypes.append(eqtl_genotype_mat[:, causal_eqtl_indices[gg]])





# Gwas sample size
N_gwas = len(trait_vec)


print('RSS GIBBS variant-predGene')
# Run RSS-GIBBS with variants and predicted genetic gene expression to assess heritabilities
rss_gibbs_var_pred_gene = rss_gibbs_variant_pred_gene.RSS_GIBBS(LD=LD, marginal_beta=gwas_beta, N=N_gwas, delta=causal_eqtl_effects, delta_indices=causal_eqtl_indices, X=genotype_mat, Y=trait_vec)
rss_gibbs_var_pred_gene.fit(max_iter=175, burn_in_iter=140)
#rss_gibbs_var_pred_gene.fit(max_iter=3, burn_in_iter=0)

print('RSS GIBBS variant-modGene')
# Run RSS-GIBBS with variants and modeled genetic gene expression to assess heritabilities
rss_gibbs_var_modeled_gene = rss_gibbs_variant_modeled_gene.RSS_GIBBS(LD=LD, marginal_beta=gwas_beta, N=N_gwas, gene_expression=gene_expression, delta_indices=causal_eqtl_indices, X=genotype_mat, Y=trait_vec, gene_genotype=gene_genotypes)
rss_gibbs_var_modeled_gene.fit(max_iter=175, burn_in_iter=140)
#rss_gibbs_var_modeled_gene.fit(max_iter=3, burn_in_iter=0)



# Print results to output
output_file = mediated_h2_results_dir + simulation_name_string + 'eqtl_ss_' + eqtl_ss + '_med_h2_summary.txt'
t = open(output_file,'w')
t.write('sim_h2\tsim_nm_h2\tsim_med_h2\t')
t.write('est_nm_h2_rss_gibbs_var_predgene\test_nm_h2_se_rss_gibbs_var_predgene\t')
t.write('est_med_h2_rss_gibbs_var_predgene\test_med_h2_se_rss_gibbs_var_predgene\t')
t.write('est_nm_h2_rss_gibbs_var_modgene\test_nm_h2_se_rss_gibbs_var_modgene\t')
t.write('est_med_h2_rss_gibbs_var_modgene\test_med_h2_se_rss_gibbs_var_modgene\n')

t.write(str(sim_genetic_var) + '\t' + str(sim_nm_genetic_var) + '\t' + str(sim_med_genetic_var) + '\t')
t.write(str(np.mean(rss_gibbs_var_pred_gene.sampled_local_nm_h2s)) + '\t' + str(np.std(rss_gibbs_var_pred_gene.sampled_local_nm_h2s)) + '\t')
t.write(str(np.mean(rss_gibbs_var_pred_gene.sampled_local_med_h2s)) + '\t' + str(np.std(rss_gibbs_var_pred_gene.sampled_local_med_h2s)) + '\t')
t.write(str(np.mean(rss_gibbs_var_modeled_gene.sampled_local_nm_h2s)) + '\t' + str(np.std(rss_gibbs_var_modeled_gene.sampled_local_nm_h2s)) + '\t')
t.write(str(np.mean(rss_gibbs_var_modeled_gene.sampled_local_med_h2s)) + '\t' + str(np.std(rss_gibbs_var_modeled_gene.sampled_local_med_h2s)) + '\n')

t.close()

print(output_file)







import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import bayesian_lmm_ss_h2_med




def load_in_gwas_data(gwas_summary_file):
	rsids = []
	betas = []
	beta_ses = []

	head_count = 0

	f = open(gwas_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count = head_count + 1
			continue

		rsids.append(data[0])
		betas.append(float(data[1]))
		beta_ses.append(float(data[2]))

	f.close()

	return np.asarray(rsids), np.asarray(betas), np.asarray(beta_ses)


def load_in_variant_ld_scores(var_ld_score_file):
	rsids = []
	ldscores = []

	head_count = 0

	f = open(var_ld_score_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count = head_count + 1
			continue
		rsids.append(data[0])
		ldscores.append(float(data[1]))

	f.close()
	return np.asarray(rsids), np.asarray(ldscores)



def create_mapping_from_rsid_to_position(rsids):
	rsid_to_position = {}

	for ii, rsid in enumerate(rsids):
		if rsid in rsid_to_position:
			print('assumption eroror')
			pdb.set_trace()
		rsid_to_position[rsid] = ii
	return rsid_to_position



def load_in_eqtl_data(eqtl_sumstat_file, eqtl_ldscore_file, rsid_to_position):
	# Initialize data structures
	genes = {}
	gene_to_betas = {}
	gene_to_beta_ses = {}
	gene_to_indices = {}
	gene_to_indices2 = {}
	gene_to_ld_scores = {}
	gene_to_cis_snps = {}

	# Add eqtl sumstats
	f = open(eqtl_sumstat_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Load in relevent fields
		gene_name = data[0]
		rsid = data[1]
		eqtl_beta = float(data[2])
		eqtl_beta_se = float(data[3])

		# Add things to global dictionary
		genes[gene_name] = 1
		# beta
		if gene_name not in gene_to_betas:
			gene_to_betas[gene_name] = []
		gene_to_betas[gene_name].append(eqtl_beta)
		# beta se
		if gene_name not in gene_to_beta_ses:
			gene_to_beta_ses[gene_name] = []
		gene_to_beta_ses[gene_name].append(eqtl_beta_se)
		# Position/index
		position = rsid_to_position[rsid]
		if gene_name not in gene_to_indices:
			gene_to_indices[gene_name] = []
		gene_to_indices[gene_name].append(position)
	f.close()

	# Add eqtl ld scores
	f = open(eqtl_ldscore_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Load in relevent fields
		gene_name = data[0]
		rsid = data[1]
		ldscore = float(data[2])
		cis_window_bool = int(data[3])

		# Add to dictionaries
		# ldscore
		if gene_name not in gene_to_ld_scores:
			gene_to_ld_scores[gene_name] = []
		gene_to_ld_scores[gene_name].append(ldscore)
		# Position/index
		position = rsid_to_position[rsid]
		if gene_name not in gene_to_indices2:
			gene_to_indices2[gene_name] = []
		gene_to_indices2[gene_name].append(position)
		# cis window
		if gene_name not in gene_to_cis_snps:
			gene_to_cis_snps[gene_name] = []
		gene_to_cis_snps[gene_name].append(cis_window_bool)
	f.close()


	# Organize results
	beta_arr = []
	beta_se_arr = []
	index_arr = []
	ldscore_arr = []
	cis_snp_arr = []

	gene_arr = np.sort([*genes])

	for gene in gene_arr:
		gene_beta = np.asarray(gene_to_betas[gene])
		gene_beta_ses = np.asarray(gene_to_beta_ses[gene])
		gene_ldscore = np.asarray(gene_to_ld_scores[gene])
		gene_indices1 = np.asarray(gene_to_indices[gene])
		gene_indices2 = np.asarray(gene_to_indices2[gene])
		gene_cis_snps = np.asarray(gene_to_cis_snps[gene])

		# Quick error check
		if np.array_equal(gene_indices1, gene_indices2) == False:
			print('assumption eroror')

		positive_ld_scores_indices = gene_ldscore > .01

		beta_arr.append(gene_beta[positive_ld_scores_indices])
		beta_se_arr.append(gene_beta_ses[positive_ld_scores_indices])
		index_arr.append(gene_indices1[positive_ld_scores_indices])
		ldscore_arr.append(gene_ldscore[positive_ld_scores_indices])
		cis_snp_arr.append(gene_cis_snps[positive_ld_scores_indices])


	return gene_arr, beta_arr, beta_se_arr, ldscore_arr, index_arr, cis_snp_arr









#####################
# Command line args
#####################
simulation_number = int(sys.argv[1])
simulation_name_string = sys.argv[2]
simulated_trait_dir = sys.argv[3]
simulated_gwas_dir = sys.argv[4]
simulation_genotype_dir = sys.argv[5]
simulated_learned_gene_models_dir = sys.argv[6]
gwas_ss = int(sys.argv[7])
eqtl_ss = int(sys.argv[8])
trait_med_h2_inference_dir = sys.argv[9]


####################
# Load in data
####################
# Load in true simulated data parameters
genetic_trait_expr_med_file = simulated_trait_dir + simulation_name_string +'_expression_mediated_trait_values.txt'
genetic_trait_nm_file = simulated_trait_dir +simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'
sim_med_h2 = np.var(np.loadtxt(genetic_trait_expr_med_file))
sim_nm_h2 = np.var(np.loadtxt(genetic_trait_nm_file))
sim_h2 = np.var(np.loadtxt(genetic_trait_nm_file) + np.loadtxt(genetic_trait_expr_med_file))

# Load in GWAS summary statistics
gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results.txt'
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(gwas_summary_file)

# Load in Variant LD-scores
var_ld_score_file = simulation_genotype_dir + 'variant_ref_geno_eqtl_1000_window_3_mb_ld_scores_bias_corrected.txt'
ldscore_rsids, ldscores = load_in_variant_ld_scores(var_ld_score_file)

# Quick error checking
if np.array_equal(gwas_rsids, ldscore_rsids) == False:
	print('rs match assumption eroror')
	pdb.set_trace()

# Create mapping from rsid to position
rsid_to_position = create_mapping_from_rsid_to_position(gwas_rsids)

# load in eqtl data
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_ss) + '_eqtl_sumstats.txt'
eqtl_ldscore_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_ss) + '_eqtl_reference_bias_corrected_ld_scores.txt'
genes, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, cis_snp_arr = load_in_eqtl_data(eqtl_sumstat_file, eqtl_ldscore_file, rsid_to_position)


####################
# Compute h2 using bayesian lmm
####################
# Equal weights
print('start')
mod = bayesian_lmm_ss_h2_med.Bayesian_LMM_SS_h2_med_inference(gwas_beta, gwas_beta_se, ldscores, eqtl_beta, eqtl_beta_se, eqtl_ldscore, eqtl_position, cis_snp_arr)
mod.fit(burn_in_iterations=2, total_iterations=5000, gamma_var_update_version='equal_weights')
#mod.fit(burn_in_iterations=2, total_iterations=5000, gamma_var_update_version='ld_score_weighting')

pdb.set_trace()

#(mod.GG)*mod.sampled_alpha_vars*np.mean(mod.sampled_eqtl_h2s,axis=1)

import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import bayesian_lmm_rss_h2_med



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


def compute_h2_with_ldsc(gwas_beta, gwas_beta_se, ldscores, n_gwas_individuals):
	z_scores = gwas_beta/gwas_beta_se

	chi_squared = np.square(z_scores)

	olser = sm.OLS(chi_squared-1, ldscores).fit()

	h2_est = len(chi_squared)*(olser.params[0]/n_gwas_individuals)

	return h2_est

def load_in_window_info(ld_window_summary_file, genes, gene_info):
	window_names = []
	window_info = {}
	used_snps = {}
	used_genes = {}
	gene_counts = {}
	f = open(ld_window_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Load in relevent info
		window_name = data[0]
		window_start = int(data[1])
		window_end = int(data[2])
		ld_file = data[3]
		window_position_file = data[4]
		window_rsid_file = data[5]
		window_middle_pos_file = data[6]

		# Load in np info
		window_pos = np.load(window_position_file)
		window_middle_indices = np.load(window_middle_pos_file)
		window_rsids = np.load(window_rsid_file,allow_pickle=True)
		middle_rsids = window_rsids[window_middle_indices==1]
		for rsid in middle_rsids:
			if rsid in used_snps:
				print('assumption eroror')
				pdb.set_trace()
			used_snps[rsid] = 1

		# Quick error check
		if window_name in window_info:
			print('assumption erooror')
			pdb.set_trace()
		if len(window_middle_indices) == 0.0:
			print('assumption eroror')

		# Add window info to global data structures
		window_names.append(window_name)

		window_info[window_name] = {}
		window_info[window_name]['ld_file'] = ld_file
		window_info[window_name]['positions'] = window_pos
		window_info[window_name]['middle_indices'] = window_middle_indices == 1
		window_info[window_name]['start_pos'] = window_start
		window_info[window_name]['end_pos'] = window_end
		window_info[window_name]['rsids'] = window_rsids

		

		window_info[window_name]['all_genes'] = []
		window_info[window_name]['middle_genes'] = []
		window_info[window_name]['gene_to_window_indices'] = {}
		window_pos_dicti = {}
		window_middle_pos_dicti = {}
		for ii,ele in enumerate(window_pos):
			window_pos_dicti[ele] = ii
		for ele in window_pos[window_middle_indices == 1]:
			window_middle_pos_dicti[ele] = 1

		for gene in genes:
			middle_gene_snp = int(np.median(gene_info[gene]['indices']))
			if middle_gene_snp in window_pos_dicti:

				gene_window_indices = []
				booler = True
				for gene_index in gene_info[gene]['indices']:
					if gene_index in window_pos_dicti:
						gene_window_indices.append(window_pos_dicti[gene_index])
					else:
						booler = False

				if booler == False:
					continue
				if gene not in gene_counts:
					gene_counts[gene] = 0
				gene_counts[gene] = gene_counts[gene] + 1

				if gene in window_info[window_name]['gene_to_window_indices']:
					print('assumptione ororor')
					pdb.set_trace()
				window_info[window_name]['gene_to_window_indices'][gene] = np.asarray(gene_window_indices)

				window_info[window_name]['all_genes'].append(gene)

				if middle_gene_snp in window_middle_pos_dicti:
					window_info[window_name]['middle_genes'].append(gene)
					if gene in used_genes:
						print('assumption eroror')
						pdb.set_trace()
					used_genes[gene] =1

		window_gene_mapping = {}
		for gene_iter, gene_name in enumerate(window_info[window_name]['all_genes']):
			window_gene_mapping[gene_name] = gene_iter
		window_info[window_name]['gene_name_to_gene_index'] = window_gene_mapping

		tmp_elements = []
		for kk in range(len(window_pos)):
			tmp_elements.append(str(kk))
		for gene in window_info[window_name]['all_genes']:
			tmp_elements.append(gene)

		window_info[window_name]['elements'] = np.asarray(tmp_elements)

	f.close()

	if len(used_genes) != len(genes):
		print('assumption error')
		pdb.set_trace()

	counts = []
	for gene in genes:
		counts.append(gene_counts[gene])


	return np.asarray(window_names), window_info


def load_in_eqtl_data(eqtl_sumstat_file, rsid_to_position, eqtl_ss):
	# Initialize data structures
	genes = {}
	gene_to_betas = {}
	gene_to_beta_ses = {}
	gene_to_indices = {}

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

	# Organize results
	beta_arr = []
	beta_se_arr = []
	index_arr = []

	gene_info = {}
	gene_arr = np.sort([*genes])

	for gene in gene_arr:
		gene_beta = np.asarray(gene_to_betas[gene])
		gene_beta_ses = np.asarray(gene_to_beta_ses[gene])
		gene_indices1 = np.asarray(gene_to_indices[gene])

		if gene in gene_info:
			print('assumtion eroror')
			pdb.set_trace()

		gene_info[gene] = {}
		gene_info[gene]['eqtl_beta'] = gene_beta
		gene_info[gene]['eqtl_beta_se'] = gene_beta_ses
		gene_info[gene]['indices'] = gene_indices1
		gene_info[gene]['N'] = eqtl_ss


	return gene_arr, gene_info

def create_mapping_from_rsid_to_position(rsids):
	rsid_to_position = {}

	for ii, rsid in enumerate(rsids):
		if rsid in rsid_to_position:
			print('assumption eroror')
			pdb.set_trace()
		rsid_to_position[rsid] = ii
	return rsid_to_position



#####################
# Command line args
#####################
simulation_number = int(sys.argv[1])
simulation_name_string = sys.argv[2]
simulated_trait_dir = sys.argv[3]
simulated_gwas_dir = sys.argv[4]
simulation_genotype_dir = sys.argv[5]
n_gwas_individuals = int(sys.argv[6])
n_eqtl_individuals = int(sys.argv[7])
simulated_learned_gene_models_dir = sys.argv[8]
trait_med_h2_inference_dir = sys.argv[9]


# Set random seed
np.random.seed(simulation_number)


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


# Create mapping from rsid to position
rsid_to_position = create_mapping_from_rsid_to_position(gwas_rsids)

# load in eqtl data
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(n_eqtl_individuals) + '_eqtl_sumstats.txt'
genes, gene_info = load_in_eqtl_data(eqtl_sumstat_file, rsid_to_position, n_eqtl_individuals)


# Load in window and LD info
ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_eqtl_1000_window_1_mb_ld_summary.txt'
window_names, window_info = load_in_window_info(ld_window_summary_file, genes, gene_info)

window_subset = np.random.choice(window_names,size=50, replace=False)
if len(np.unique(window_subset)) != 50:
	print('assumption eroror')

temp_output_file = trait_med_h2_inference_dir + simulation_name_string + '_rss_h2_med_inference_tmp_results.txt'
mod = bayesian_lmm_rss_h2_med.Bayesian_LMM_RSS_h2_med_inference(gwas_beta, n_gwas_individuals, window_subset, window_info, gene_info,temp_output_file)
mod.fit(burn_in_iterations=1, total_iterations=10000)




'''
# Quick error checking
if np.array_equal(gwas_rsids, ldscore_rsids) == False:
	print('rs match assumption eroror')
	pdb.set_trace()


####################
# Compute h2 using ldscore regression
####################
ldsc_h2 = compute_h2_with_ldsc(gwas_beta, gwas_beta_se, ldscores, n_gwas_individuals)

'''





'''
# Load in Variant LD-scores
var_ld_score_file = simulation_genotype_dir + 'variant_ref_geno_eqtl_1000_window_3_mb_ld_scores_bias_corrected.txt'
#var_ld_score_file = simulated_gwas_dir + simulation_name_string + '_variant_ref_geno_window_2_mb_ld_scores_bias_corrected.txt'
#var_ld_score_file = '/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/single_tissue_simulation/simulated_gwas/simulation_1_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_variant_ref_geno_window_3_mb_ld_scores_bias_corrected.txt'
ldscore_rsids, ldscores = load_in_variant_ld_scores(var_ld_score_file)

####################
# Compute h2 using bayesian lmm
####################
# Equal weights
mod = bayesian_lmm_ss_h2.Bayesian_LMM_SS_h2_inference(gwas_beta, gwas_beta_se, ldscores)
mod.fit(burn_in_iterations=15000, total_iterations=20000, gamma_var_update_version='equal_weights')
bayesian_lmm_h2_mean = np.mean(mod.sampled_h2)
bayesian_lmm_h2_lb = np.sort(mod.sampled_h2)[int(np.round(len(mod.sampled_h2)*.025))]
bayesian_lmm_h2_ub = np.sort(mod.sampled_h2)[int(np.round(len(mod.sampled_h2)*.975))]

# Variant Ld scoring weights
mod2 = bayesian_lmm_ss_h2.Bayesian_LMM_SS_h2_inference(gwas_beta, gwas_beta_se, ldscores)
mod2.fit(burn_in_iterations=15000, total_iterations=20000, gamma_var_update_version='ld_score_weighting')
bayesian_lmm_h2_mean2 = np.mean(mod2.sampled_h2)
bayesian_lmm_h2_lb2 = np.sort(mod2.sampled_h2)[int(np.round(len(mod2.sampled_h2)*.025))]
bayesian_lmm_h2_ub2 = np.sort(mod2.sampled_h2)[int(np.round(len(mod2.sampled_h2)*.975))]


####################
# Print results to output file
####################
output_file = trait_h2_inference_dir + simulation_name_string + '_h2_inference_summary.txt'
t = open(output_file,'w')
# header
t.write('method\tsim_h2\test_h2\test_h2_95_lb\test_h2_95_ub\n')
# Results
t.write('LDSC\t' + str(sim_h2) + '\t' + str(ldsc_h2) + '\t' + str('NA') + '\t' + str('NA') + '\n')
t.write('bayesian_lmm_equal_weights\t' + str(sim_h2) + '\t' + str(bayesian_lmm_h2_mean) + '\t' + str(bayesian_lmm_h2_lb) + '\t' + str(bayesian_lmm_h2_ub) + '\n')
t.write('bayesian_lmm_ldscore_weights\t' + str(sim_h2) + '\t' + str(bayesian_lmm_h2_mean2) + '\t' + str(bayesian_lmm_h2_lb2) + '\t' + str(bayesian_lmm_h2_ub2) + '\n')

t.close()
'''



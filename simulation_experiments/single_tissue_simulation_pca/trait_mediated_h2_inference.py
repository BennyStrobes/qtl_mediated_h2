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



def convert_to_pc_sum_stat_space(gwas_beta, gwas_rsids, quasi_ld_window_summary_file):
	gwas_beta_pc = []
	var_ld_score_pc = []
	pc_snp_names = []
	f = open(quasi_ld_window_summary_file)
	window_to_q_w_files = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		window_name = data[0]
		window_snp_indices = np.load(data[4])
		window_rsids = np.load(data[5])
		# QUick error checking
		if np.array_equal(window_rsids, gwas_rsids[window_snp_indices]) == False:
			print('assumption eroror')
			pdb.set_trace()
		Q_mat = np.load(data[7])
		w_premult_mat = np.load(data[8])

		window_pc_gwas_beta = np.dot(w_premult_mat, gwas_beta[window_snp_indices])
		window_pc_ld_scores = np.diag(np.dot(Q_mat, np.transpose(Q_mat)))

		gwas_beta_pc.append(window_pc_gwas_beta)
		var_ld_score_pc.append(window_pc_ld_scores)


		if window_name in window_to_q_w_files:
			print('assumption eororor')
			pdb.set_trace()
		window_to_q_w_files[window_name] = (data[7], data[8])

		for ii in range(len(window_pc_ld_scores)):
			pc_snp_names.append(window_name + ':snp_' + str(ii))



	gwas_beta_pc = np.hstack(gwas_beta_pc)
	var_ld_score_pc = np.hstack(var_ld_score_pc)


	return gwas_beta_pc, var_ld_score_pc, len(gwas_beta), np.asarray(pc_snp_names), window_to_q_w_files



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


def extract_eqtl_sumstats_for_specific_gene(sumstats_file, gene_name):
	f = open(sumstats_file)
	sumstats = []
	cis_snps = []
	window_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != gene_name:
			continue
		sumstats.append(data[4])
		cis_snps.append(data[3])
		window_name = data[2]
		window_names.append(window_name)
	f.close()

	unique_window_names = np.unique(window_names)
	if len(unique_window_names) != 1:
		print('assumption eroror')
		pdb.set_trace()

	gene_window_name = unique_window_names[0]

	return np.asarray(sumstats).astype(float), np.asarray(cis_snps).astype(float), gene_window_name

def load_in_eqtl_data(eqtl_sumstat_file, pc_snp_name_to_position, window_name_to_q_mat_and_w_mat_files):
	# First get list of gene names
	gene_names = []
	f = open(eqtl_sumstat_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_names.append(data[0])
	f.close()
	gene_names = np.unique(gene_names)


	gene_arr = []
	beta_arr = []
	ldscore_arr = []
	index_arr = []
	n_cis_snp_arr = []
	# Now loop through genes
	for gene_name in gene_names:
		# Extract summary stats for specific gene
		eqtl_gene_beta, gene_cis_snp_indices, gene_window_name = extract_eqtl_sumstats_for_specific_gene(eqtl_sumstat_file, gene_name)

		# Load in q_mat and w_mat for this window
		q_mat_file, w_mat_file = window_name_to_q_mat_and_w_mat_files[gene_window_name]
		q_mat = np.load(q_mat_file)
		w_mat = np.load(w_mat_file)

		# Get PC space ss
		gene_pc_beta = np.dot(w_mat, eqtl_gene_beta)
		gene_pc_ldscores = np.diag(np.dot(q_mat[:, gene_cis_snp_indices==1], np.transpose(q_mat[:, gene_cis_snp_indices==1])))


		# Get names of snps in gene
		gene_snps = []
		gene_snp_positions = []
		for ii in range(len(gene_pc_ldscores)):
			snp_name = gene_window_name + ':snp_' + str(ii)
			gene_snps.append(snp_name)
			gene_snp_positions.append(pc_snp_name_to_position[snp_name])
		gene_snps = np.asarray(gene_snps)
		gene_snp_positions = np.asarray(gene_snp_positions)

		# Add to global array
		gene_arr.append(gene_name)
		beta_arr.append(gene_pc_beta)
		ldscore_arr.append(gene_pc_ldscores)
		index_arr.append(gene_snp_positions)
		n_cis_snp_arr.append(np.sum(gene_cis_snp_indices))

	return np.asarray(gene_arr), beta_arr, ldscore_arr, index_arr, np.asarray(n_cis_snp_arr)









#####################
# Command line args
#####################
simulation_number= int(sys.argv[1])
simulation_name_string= sys.argv[2]
simulated_trait_dir=sys.argv[3]
simulated_gwas_dir=sys.argv[4]
simulation_genotype_dir=sys.argv[5]
simulated_learned_gene_models_dir=sys.argv[6]
n_gwas_individuals=int(sys.argv[7])
eqtl_sample_size=int(sys.argv[8])
trait_med_h2_inference_dir=sys.argv[9]
cc_hyperparam_str = sys.argv[10]

cc_hyperparam = float(cc_hyperparam_str)




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
N_gwas = n_gwas_individuals

# Convert gwas_beta to gwas_beta_pc and get var_pc_ld_scores
quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_summary.txt' 
gwas_beta_pc, var_ld_score_pc, n_gwas_snps, pc_snp_names, window_name_to_q_mat_and_w_mat_files = convert_to_pc_sum_stat_space(gwas_beta, gwas_rsids, quasi_ld_window_summary_file)


# Create mapping from rsid to position
pc_snp_name_to_position = create_mapping_from_rsid_to_position(pc_snp_names)



# load in eqtl data
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_eqtl_sumstats.txt'
genes, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_n_cis_snps = load_in_eqtl_data(eqtl_sumstat_file, pc_snp_name_to_position, window_name_to_q_mat_and_w_mat_files)

'''
import pickle
f = open('genes.pkl', 'wb')
pickle.dump(genes, f)
f.close()
f = open('eqtl_beta.pkl', 'wb')
pickle.dump(eqtl_beta, f)
f.close()
f = open('eqtl_ldscore.pkl', 'wb')
pickle.dump(eqtl_ldscore, f)
f.close()
f = open('eqtl_position.pkl', 'wb')
pickle.dump(eqtl_position, f)
f.close()
f = open('eqtl_n_cis_snps.pkl', 'wb')
pickle.dump(eqtl_n_cis_snps, f)
f.close()
f = open('genes.pkl', 'rb')
genes = pickle.load(f) 
f.close()
f = open('eqtl_beta.pkl', 'rb')
eqtl_beta = pickle.load(f) 
f.close()
f = open('eqtl_ldscore.pkl', 'rb')
eqtl_ldscore = pickle.load(f) 
f.close()
f = open('eqtl_position.pkl', 'rb')
eqtl_position = pickle.load(f) 
f.close()
f = open('eqtl_n_cis_snps.pkl', 'rb')
eqtl_n_cis_snps = pickle.load(f) 
f.close()
'''


####################
# Compute h2 using bayesian lmm
####################
# Equal weights
print('start')
mod = bayesian_lmm_ss_h2_med.Bayesian_LMM_SS_h2_med_inference(gwas_beta_pc, var_ld_score_pc, N_gwas, n_gwas_snps, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_sample_size, eqtl_n_cis_snps)
mod.fit(burn_in_iterations=6000, total_iterations=12000, v0=0.0, s_sq=0.0, cc=cc_hyperparam)

output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_' + str(eqtl_sample_size) + '_' + cc_hyperparam_str + '_med_h2_results.txt'
t = open(output_file,'w')
t.write('iter\tobs_med_h2\tobs_nm_h2\tmed_h2\talt_med_h2\tnm_h2\teqtl_h2\n')

for sample_iter in range(len(mod.sampled_med_h2)):
	t.write(str(sample_iter) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(mod.sampled_med_h2[sample_iter]) + '\t' + str(mod.sampled_alt_med_h2[sample_iter]) + '\t' + str(mod.sampled_nm_h2[sample_iter]) + '\t' + str(mod.sampled_eqtl_h2s[sample_iter]) + '\n')
t.close()

print(output_file)





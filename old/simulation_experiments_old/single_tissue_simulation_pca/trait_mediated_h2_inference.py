import sys
import numpy as np 
import pandas as pd
import os
import pdb
from scipy.stats import invgamma
import statsmodels.api as sm
import bayesian_lmm_ss_h2_med
import bayesian_lmm_VI_ss_h2_med
import bayesian_lmm_ss_h2_med_no_ldscores
import bayesian_multivariate_lmm_ss_h2_med
import bayesian_lmm_ss_h2_med_shared_gene_variance
import bayesian_lmm_ss_h2_med_cross_gene_gamma_scale
import bayesian_lmm_ss_h2_med_seperated_ld_scores
import bayesian_lmm_ss_two_step_h2_med
import bayesian_lmm_ss_double_dip_two_step_h2_med
import bayesian_lmm_ss_h2_med_vi_variance

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
		ld_file = data[3]

		window_pc_gwas_beta = np.dot(w_premult_mat, gwas_beta[window_snp_indices])
		window_pc_ld_scores = np.diag(np.dot(Q_mat, np.transpose(Q_mat)))


		gwas_beta_pc.append(window_pc_gwas_beta)
		var_ld_score_pc.append(window_pc_ld_scores)


		if window_name in window_to_q_w_files:
			print('assumption eororor')
			pdb.set_trace()
		window_to_q_w_files[window_name] = (data[7], data[8], ld_file)

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



def load_in_eqtl_data_and_ratio_between_two_data_sets(eqtl_sumstat_file, pc_snp_name_to_position, window_name_to_q_mat_and_w_mat_files, eqtl_sample_size):
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
	ldscore_ratio_arr = []
	index_arr = []
	n_cis_snp_arr = []
	# Now loop through genes
	for gene_name in gene_names:
		# Extract summary stats for specific gene
		eqtl_gene_beta, gene_cis_snp_indices, gene_window_name = extract_eqtl_sumstats_for_specific_gene(eqtl_sumstat_file, gene_name)

		# Load in q_mat and w_mat for this window
		q_mat_file, w_mat_file, ld_file = window_name_to_q_mat_and_w_mat_files[gene_window_name]
		q_mat = np.load(q_mat_file)
		w_mat = np.load(w_mat_file)

		gwas_gene_pc_ldscores = np.diag(np.dot(q_mat[:, gene_cis_snp_indices==1], np.transpose(q_mat[:, gene_cis_snp_indices==1])))
		eqtl_ld_file = ld_file.split('geno_gwas')[0] + 'geno_eqtl_' + str(eqtl_sample_size) + ld_file.split('geno_gwas')[1]
		eqtl_ld_mat = np.load(eqtl_ld_file)
		q_mat = np.dot(w_mat, eqtl_ld_mat)

		# Get PC space ss
		gene_pc_beta = np.dot(w_mat, eqtl_gene_beta)
		eqtl_gene_pc_ldscores = np.diag(np.dot(q_mat[:, gene_cis_snp_indices==1], np.transpose(q_mat[:, gene_cis_snp_indices==1])))


		# Get names of snps in gene
		gene_snps = []
		gene_snp_positions = []
		for ii in range(len(eqtl_gene_pc_ldscores)):
			snp_name = gene_window_name + ':snp_' + str(ii)
			gene_snps.append(snp_name)
			gene_snp_positions.append(pc_snp_name_to_position[snp_name])
		gene_snps = np.asarray(gene_snps)
		gene_snp_positions = np.asarray(gene_snp_positions)

		# Add to global array
		gene_arr.append(gene_name)
		beta_arr.append(gene_pc_beta)
		ldscore_arr.append(gwas_gene_pc_ldscores)
		ldscore_ratio_arr.append(eqtl_gene_pc_ldscores/gwas_gene_pc_ldscores)
		index_arr.append(gene_snp_positions)
		n_cis_snp_arr.append(np.sum(gene_cis_snp_indices))

	return np.asarray(gene_arr), beta_arr, ldscore_arr, ldscore_ratio_arr, index_arr, np.asarray(n_cis_snp_arr)


def load_in_eqtl_data(eqtl_sumstat_file, pc_snp_name_to_position, window_name_to_q_mat_and_w_mat_files, eqtl_sample_size, eqtl_ld='out_of_sample'):
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

	print(len(gene_names))


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
		q_mat_file, w_mat_file, ld_file = window_name_to_q_mat_and_w_mat_files[gene_window_name]
		q_mat = np.load(q_mat_file)
		w_mat = np.load(w_mat_file)

		if eqtl_ld == 'in_sample':
			gene_pc_ldscores_out_of_sample = np.diag(np.dot(q_mat[:, gene_cis_snp_indices==1], np.transpose(q_mat[:, gene_cis_snp_indices==1])))
			eqtl_ld_file = ld_file.split('geno_gwas')[0] + 'geno_eqtl_' + str(eqtl_sample_size) + ld_file.split('geno_gwas')[1]
			eqtl_ld_mat = np.load(eqtl_ld_file)
			q_mat = np.dot(w_mat, eqtl_ld_mat)

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

print(sim_h2)
print(sim_med_h2)
print(sim_nm_h2)

# Load in GWAS summary statistics
gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results.txt'
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(gwas_summary_file)
N_gwas = n_gwas_individuals

# Convert gwas_beta to gwas_beta_pc and get var_pc_ld_scores
quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_summary.txt'
#quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_summary.txt'
gwas_beta_pc, var_ld_score_pc, n_gwas_snps, pc_snp_names, window_name_to_q_mat_and_w_mat_files = convert_to_pc_sum_stat_space(gwas_beta, gwas_rsids, quasi_ld_window_summary_file)


# Create mapping from rsid to position
pc_snp_name_to_position = create_mapping_from_rsid_to_position(pc_snp_names)

# load in eqtl data
# Out of sample eqtl ld
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
#eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_big_window_eqtl_sumstats.txt'
genes, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_n_cis_snps = load_in_eqtl_data(eqtl_sumstat_file, pc_snp_name_to_position, window_name_to_q_mat_and_w_mat_files,eqtl_sample_size, eqtl_ld='out_of_sample')



####################
# Compute mediated h2 using bayesian lmm
####################
'''
# Equal weights
print('start')
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_resid_var_True_cc_' + str(cc_hyperparam) + '_alt_3.txt'
mod = bayesian_lmm_ss_h2_med.Bayesian_LMM_SS_h2_med_inference(gwas_beta_pc, var_ld_score_pc, N_gwas, n_gwas_snps, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_sample_size, eqtl_n_cis_snps)
mod.fit(burn_in_iterations=100, total_iterations=1000000, update_gwas_resid_var=True, update_eqtl_resid_var=True, v0=0.0, s_sq=0.0, cc=cc_hyperparam)
t = open(output_file,'w')
t.write('sim_h2\tsim_nm_h2\tsim_med_h2\titer\tnm_h2_v1\tnm_h2_v2\tmed_h2_v1\tmed_h2_v2\ttotal_h2\teqtl_h2_v1\teqtl_h2_v2\n')
for ii, itera in enumerate(mod.sampled_iters):
	t.write(str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(itera) + '\t' + str(mod.sampled_nm_h2_v1[ii]) + '\t' + str(mod.sampled_nm_h2_v2[ii]) + '\t' + str(mod.sampled_med_h2_v1[ii]) + '\t' + str(mod.sampled_med_h2_v2[ii]) + '\t' + str(mod.sampled_total_h2[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v1[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v2[ii]) + '\n')
t.close()
'''


####################
# Compute mediated h2 using bayesian lmm (but with wrong eqtl data)
####################
# load in eqtl data
# Out of sample eqtl ld
info_str = simulation_name_string.split('_chrom')
simulation_number = info_str[0].split('_')[1]
new_simulation_name_string = 'simulation_' + str(int(simulation_number) + 1) + '_chrom' + info_str[1]
eqtl_sumstat_file = simulated_learned_gene_models_dir + new_simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'

genes, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_n_cis_snps = load_in_eqtl_data(eqtl_sumstat_file, pc_snp_name_to_position, window_name_to_q_mat_and_w_mat_files,eqtl_sample_size, eqtl_ld='out_of_sample')

# Equal weights
print('start')
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_resid_var_True_cc_' + str(cc_hyperparam) + '_wrong_eqtls.txt'
mod = bayesian_lmm_ss_h2_med.Bayesian_LMM_SS_h2_med_inference(gwas_beta_pc, var_ld_score_pc, N_gwas, n_gwas_snps, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_sample_size, eqtl_n_cis_snps)
mod.fit(burn_in_iterations=10000, total_iterations=600000, update_gwas_resid_var=True, update_eqtl_resid_var=True, v0=0.0, s_sq=0.0, cc=cc_hyperparam)
t = open(output_file,'w')
t.write('sim_h2\tsim_nm_h2\tsim_med_h2\titer\tnm_h2_v1\tnm_h2_v2\tmed_h2_v1\tmed_h2_v2\ttotal_h2\teqtl_h2_v1\teqtl_h2_v2\n')
for ii, itera in enumerate(mod.sampled_iters):
	t.write(str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(itera) + '\t' + str(mod.sampled_nm_h2_v1[ii]) + '\t' + str(mod.sampled_nm_h2_v2[ii]) + '\t' + str(mod.sampled_med_h2_v1[ii]) + '\t' + str(mod.sampled_med_h2_v2[ii]) + '\t' + str(mod.sampled_total_h2[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v1[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v2[ii]) + '\n')
t.close()


'''
####################
# Compute h2 using bayesian lmm
# with VI variance updates
####################
# Equal weights
print('start')
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_resid_var_True_cc_' + str(cc_hyperparam) + '_vi_variance_alt.txt'
mod = bayesian_lmm_ss_h2_med_vi_variance.Bayesian_LMM_SS_h2_med_inference(gwas_beta_pc, var_ld_score_pc, N_gwas, n_gwas_snps, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_sample_size, eqtl_n_cis_snps)
mod.fit(burn_in_iterations=200, total_iterations=600000, update_gwas_resid_var=True, update_eqtl_resid_var=True, v0=0.0, s_sq=0.0, cc=cc_hyperparam)
t = open(output_file,'w')
t.write('sim_h2\tsim_nm_h2\tsim_med_h2\titer\tnm_h2_v1\tnm_h2_v2\tmed_h2_v1\tmed_h2_v2\ttotal_h2\teqtl_h2_v1\teqtl_h2_v2\n')
for ii, itera in enumerate(mod.sampled_iters):
	t.write(str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(itera) + '\t' + str(mod.sampled_nm_h2_v1[ii]) + '\t' + str(mod.sampled_nm_h2_v2[ii]) + '\t' + str(mod.sampled_med_h2_v1[ii]) + '\t' + str(mod.sampled_med_h2_v2[ii]) + '\t' + str(mod.sampled_total_h2[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v1[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v2[ii]) + '\n')
t.close()
'''


'''
####################
# Compute h2 using bayesian lmm
# with VI variance updates
####################
# Equal weights
print('start')
gene_gamma_shape_param=1e-5
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_resid_var_True_cc_' + str(cc_hyperparam) + '_cross_gene_gamma_scale_' + str(gene_gamma_shape_param) + '.txt'
mod = bayesian_lmm_ss_h2_med_cross_gene_gamma_scale.Bayesian_LMM_SS_h2_med_inference(gwas_beta_pc, var_ld_score_pc, N_gwas, n_gwas_snps, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_sample_size, eqtl_n_cis_snps)
mod.fit(burn_in_iterations=10000, total_iterations=600000, update_resid_var=True, v0=0.0, s_sq=0.0, cc=cc_hyperparam, gene_gamma_shape=gene_gamma_shape_param)
t = open(output_file,'w')
t.write('sim_h2\tsim_nm_h2\tsim_med_h2\titer\tnm_h2_v1\tnm_h2_v2\tmed_h2_v1\tmed_h2_v2\ttotal_h2\teqtl_h2_v1\teqtl_h2_v2\n')
for ii, itera in enumerate(mod.sampled_iters):
	t.write(str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(itera) + '\t' + str(mod.sampled_nm_h2_v1[ii]) + '\t' + str(mod.sampled_nm_h2_v2[ii]) + '\t' + str(mod.sampled_med_h2_v1[ii]) + '\t' + str(mod.sampled_med_h2_v2[ii]) + '\t' + str(mod.sampled_total_h2[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v1[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v2[ii]) + '\n')
t.close()

gene_gamma_shape_param=1e-3
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_resid_var_True_cc_' + str(cc_hyperparam) + '_cross_gene_gamma_scale_' + str(gene_gamma_shape_param) + '.txt'
mod = bayesian_lmm_ss_h2_med_cross_gene_gamma_scale.Bayesian_LMM_SS_h2_med_inference(gwas_beta_pc, var_ld_score_pc, N_gwas, n_gwas_snps, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_sample_size, eqtl_n_cis_snps)
mod.fit(burn_in_iterations=10000, total_iterations=600000, update_resid_var=True, v0=0.0, s_sq=0.0, cc=cc_hyperparam, gene_gamma_shape=gene_gamma_shape_param)
t = open(output_file,'w')
t.write('sim_h2\tsim_nm_h2\tsim_med_h2\titer\tnm_h2_v1\tnm_h2_v2\tmed_h2_v1\tmed_h2_v2\ttotal_h2\teqtl_h2_v1\teqtl_h2_v2\n')
for ii, itera in enumerate(mod.sampled_iters):
	t.write(str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(itera) + '\t' + str(mod.sampled_nm_h2_v1[ii]) + '\t' + str(mod.sampled_nm_h2_v2[ii]) + '\t' + str(mod.sampled_med_h2_v1[ii]) + '\t' + str(mod.sampled_med_h2_v2[ii]) + '\t' + str(mod.sampled_total_h2[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v1[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v2[ii]) + '\n')
t.close()

gene_gamma_shape_param=1e-1
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_resid_var_True_cc_' + str(cc_hyperparam) + '_cross_gene_gamma_scale_' + str(gene_gamma_shape_param) + '.txt'
mod = bayesian_lmm_ss_h2_med_cross_gene_gamma_scale.Bayesian_LMM_SS_h2_med_inference(gwas_beta_pc, var_ld_score_pc, N_gwas, n_gwas_snps, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_sample_size, eqtl_n_cis_snps)
mod.fit(burn_in_iterations=10000, total_iterations=600000, update_resid_var=True, v0=0.0, s_sq=0.0, cc=cc_hyperparam, gene_gamma_shape=gene_gamma_shape_param)
t = open(output_file,'w')
t.write('sim_h2\tsim_nm_h2\tsim_med_h2\titer\tnm_h2_v1\tnm_h2_v2\tmed_h2_v1\tmed_h2_v2\ttotal_h2\teqtl_h2_v1\teqtl_h2_v2\n')
for ii, itera in enumerate(mod.sampled_iters):
	t.write(str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(itera) + '\t' + str(mod.sampled_nm_h2_v1[ii]) + '\t' + str(mod.sampled_nm_h2_v2[ii]) + '\t' + str(mod.sampled_med_h2_v1[ii]) + '\t' + str(mod.sampled_med_h2_v2[ii]) + '\t' + str(mod.sampled_total_h2[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v1[ii]) + '\t' + str(mod.sampled_eqtl_h2s_v2[ii]) + '\n')
t.close()

'''




















































































# OLD. didn't seem to work

'''
####################
# LDSC h2 stuff
####################
h2_v1 = np.sum(np.square(gwas_beta_pc) - (1.0/N_gwas))
model = sm.OLS(np.square(gwas_beta_pc) - (1.0/N_gwas), var_ld_score_pc).fit()
h2_v2 = model.params[0]*np.sum(var_ld_score_pc)

mod = bayesian_lmm_ss_two_step_h2_med.Bayesian_LMM_SS_h2_med_inference(gwas_beta_pc, var_ld_score_pc, N_gwas, n_gwas_snps, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_sample_size, eqtl_n_cis_snps)
mod.fit(burn_in_iterations=100000, total_iterations=110000, update_eqtl_resid_var=True, v0=0.0, s_sq=0.0, cc=cc_hyperparam)
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_two_step_resid_var_True_cc_' + str(cc_hyperparam) + '_alt_2.txt'
t = open(output_file,'w')
t.write('sim_h2\tsim_nm_h2\tsim_med_h2\tnm_h2\tmed_h2\ttotal_h2_v1\ttotal_h2_v2\teqtl_h2\n')
t.write(str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(mod.est_nm_h2) + '\t' + str(mod.est_med_h2) + '\t' + str(h2_v1) + '\t' + str(h2_v2) + '\t' + str(mod.eqtl_h2) + '\n')
t.close()
print(output_file)




####################
# Compute h2 using bayesian lmm two step
####################
# Equal weights
print('start')
mod = bayesian_lmm_ss_double_dip_two_step_h2_med.Bayesian_LMM_SS_h2_med_inference(gwas_beta_pc, var_ld_score_pc, N_gwas, n_gwas_snps, eqtl_beta, eqtl_ldscore, eqtl_position, eqtl_sample_size, eqtl_n_cis_snps)
mod.fit(burn_in_iterations=590000, total_iterations=600000, update_gwas_resid_var=True, update_eqtl_resid_var=True, v0=0.0, s_sq=0.0, cc=cc_hyperparam)
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_double_dip_two_step_resid_var_True_cc_' + str(cc_hyperparam) + '_alt_2.txt'

t = open(output_file,'w')
t.write('sim_h2\tsim_nm_h2\tsim_med_h2\tnm_h2\tmed_h2\tnm_h2_bayes\tmed_h2_bayes\n')
t.write(str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(mod.two_step_nm_h2_est) + '\t' + str(mod.two_step_med_h2_est) + '\t' + str(np.mean(mod.sampled_med_h2_v2)) + '\t' + str(np.mean(mod.sampled_nm_h2_v2)) + '\n')
t.close()
print(output_file)
'''


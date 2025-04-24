import numpy as np 
import os
import sys
import pdb
from sklearn.linear_model import LinearRegression


def extract_causal_eqtl_effects(gene_summary_file):
	causal_eqtl_indices = []
	true_causal_eqtl_effects = []
	head_count = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
		true_causal_eqtl_effect = np.asarray(data[2].split(',')).astype(float)

		causal_eqtl_indices.append(gene_snp_indices)
		true_causal_eqtl_effects.append(true_causal_eqtl_effect)
	f.close()
	return causal_eqtl_indices, true_causal_eqtl_effects


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

def extract_true_sim_gene_ld_scores(LD, causal_eqtl_indices, causal_eqtl_effects):
	# Extract general parameters
	n_snps = LD.shape[0]
	n_genes = len(causal_eqtl_indices)

	##########################
	# Create gene-variant LD mat
	##########################
	n_snps = LD.shape[0]
	gene_variant_ld = np.zeros((n_genes, n_snps))
	for gg in range(n_genes):
		gene_variant_ld[gg,:] = np.dot(LD[:, causal_eqtl_indices[gg]], causal_eqtl_effects[gg])
	# Square LD
	gene_variant_ld_sq = np.square(gene_variant_ld)

	##########################
	# Get gene LD scores
	##########################
	true_gene_ld_scores = np.hstack((np.sum(gene_variant_ld_sq,axis=0)))

	return true_gene_ld_scores


def run_sldsc_to_get_mediated_heritabilities(gwas_chi_sq, var_ld_scores, gene_ld_scores, N_gwas, n_snps, n_genes, average_ge_h2, fit_intercept=True, noise_vec=None):
	if fit_intercept:
		# Get n_eqtl_datasets
		n_eqtl_datasets = gene_ld_scores.shape[0]

		if noise_vec is None:
			noise_vec = np.zeros(n_eqtl_datasets + 1)  # Plus one comes from variant annotation

		# Prepare LDSC data
		covs = np.transpose(np.vstack((np.ones(n_snps), var_ld_scores,gene_ld_scores)))
		noise_vec = np.hstack((np.zeros(1),noise_vec))
		S = np.diag(noise_vec)

		# Run SLDSC
		ldsc_coefs = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(covs), covs) - S*n_snps), np.transpose(covs)), gwas_chi_sq-1)
		
		# Convert SLDSC coefficients to heritability parameters
		nm_h2 = ldsc_coefs[1]*n_snps/N_gwas
		med_h2 = ldsc_coefs[2:]*n_genes*average_ge_h2/N_gwas
	else:
		print('The following scenario has not been tested: this type of SLDSC without an intercept and with noise correction.. PLEASE PROCEED WITH CAUTION')
		pdb.set_trace()
		# Prepare LDSC data
		covs = np.transpose(np.vstack((var_ld_scores, gene_ld_scores)))
		S = np.diag(noise_vec)

		# Run SLDSC
		ldsc_coefs = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(covs), covs) - S*n_snps), np.transpose(covs)), gwas_chi_sq-1)
		
		# Convert SLDSC coefficients to heritability parameters
		nm_h2 = ldsc_coefs[0]*n_snps/N_gwas
		med_h2 = ldsc_coefs[1]*n_genes*average_ge_h2/N_gwas


	return med_h2






#########################
# Command line args
#########################
simulation_name_string = sys.argv[1]
simulated_gwas_data_dir = sys.argv[2]
mediated_h2_results_dir = sys.argv[3]
processed_genotype_data_dir = sys.argv[4]
simulated_eqtl_data_dir = sys.argv[5]
selection = sys.argv[6]
eqtl_ss = sys.argv[7]
simulated_learned_gene_models_dir = sys.argv[8]
n_eqtl_datasets = int(sys.argv[9])





# File containing gwas summary statistics
gwas_sum_stats_file = simulated_gwas_data_dir + simulation_name_string + selection + '_selection_' + 'simulated_gwas_summary_stats.txt'
gwas_z_scores = load_in_gwas_z_scores(gwas_sum_stats_file)
gwas_se = load_in_gwas_se(gwas_sum_stats_file)
gwas_beta = gwas_z_scores*gwas_se
gwas_chi_sq = np.square(gwas_z_scores)

# Simulated total genetic vars
simulated_genetic_var_file = simulated_gwas_data_dir + simulation_name_string+ selection + '_selection_' + 'simulated_genetic_var.npy'
sim_genetic_var = np.load(simulated_genetic_var_file) + 0.0
simulated_nm_genetic_var_file = simulated_gwas_data_dir + simulation_name_string + selection + '_selection_'+ 'simulated_nm_genetic_var.npy'
sim_nm_genetic_var = np.load(simulated_nm_genetic_var_file) + 0.0
simulated_mediated_genetic_var_file = simulated_gwas_data_dir + simulation_name_string+ selection + '_selection_' + 'simulated_mediated_genetic_var.npy'
sim_med_genetic_var = np.load(simulated_mediated_genetic_var_file) + 0.0

# Print what simulated total genetic vars are
print(sim_genetic_var)
print(sim_nm_genetic_var)
print(sim_med_genetic_var)



# Load in Genotype file
ld_file = processed_genotype_data_dir + 'gwas_genotype_LD_1.npy'
LD = np.load(ld_file)

# Simulated Trait
trait_file = simulated_gwas_data_dir + simulation_name_string + selection + '_selection_' + 'simulated_trait.txt'
trait_vec = load_in_trait_file(trait_file)

# Gwas sample size
N_gwas = len(trait_vec)


# Number of snps
n_snps = LD.shape[0]


# Extract variant LD scores
var_ld_scores = np.hstack((np.sum(np.square(LD),axis=0)))



# Initialize global vectors across eqtl datasets keeping track of eqtldata set parameter
global_n_genes = []
global_sim_gene_ld_scores = []
global_noisy_gene_ld_scores = []
global_noise_gene_ld_score_noise_est = []
global_noise_gene_ld_score_noise_true = []
global_posterior_known_param_gene_ld_scores = []
global_posterior_blup_ld_scores = []
global_prediction_blup_ld_scores = []


# Loop through eqtl datasets
for eqtl_dataset_iter in range(n_eqtl_datasets):

	# File summarizing simulated gene models
	gene_summary_file = simulated_eqtl_data_dir + simulation_name_string + 'causal_eqtl_effect_summary_eqtl_dataset_' + str(eqtl_dataset_iter + 1) + '.txt'
	# Extract simulated causal eqtl effects
	causal_eqtl_indices, true_causal_eqtl_effects = extract_causal_eqtl_effects(gene_summary_file)

	# Number of genes
	n_genes = len(causal_eqtl_indices)

	###########################
	# Extract various gene LD Scores for this dataset
	###########################
	# Extract true gene LD scores
	sim_gene_ld_scores = extract_true_sim_gene_ld_scores(LD, causal_eqtl_indices, true_causal_eqtl_effects)
	# Extract estimated gene LD scores according to square of marginal association analysis
	noisy_gene_ld_score_file = simulated_learned_gene_models_dir + simulation_name_string + 'estimated_noisy_gene_related_snp_annotations_' + eqtl_ss+ '_eqtl_dataset_' + str(eqtl_dataset_iter+1) + '.npy'
	noisy_gene_ld_scores = np.load(noisy_gene_ld_score_file)
	# Extract estimated noise on above estimated gene ld scores
	noisy_gene_ld_score_noise_file = simulated_learned_gene_models_dir + simulation_name_string + 'estimated_noisy_gene_related_snp_annotation_noises_' + eqtl_ss+ '_eqtl_dataset_' + str(eqtl_dataset_iter+1) + '.npy'
	noise_gene_ld_score_noise_ests = np.load(noisy_gene_ld_score_noise_file)
	noise_gene_ld_score_noise_est = np.mean(noise_gene_ld_score_noise_ests)  # Take average across snps
	# Extract posterior estimated gene ld scores (known gene h2, known residual variance)
	posterior_t_h2_t_residvar_gene_ld_score_file = simulated_learned_gene_models_dir + simulation_name_string + 'posterior_known_h2_known_resid_var_gene_related_snp_annotations_' + eqtl_ss+ '_eqtl_dataset_' + str(eqtl_dataset_iter+1) + '.npy'
	posterior_t_h2_t_residvar_gene_ld_scores = np.load(posterior_t_h2_t_residvar_gene_ld_score_file)
	# Extract posterior estimated gene ld scores (BLUP)
	posterior_blup_gene_ld_score_file = simulated_learned_gene_models_dir + simulation_name_string + 'posterior_blup_gene_related_snp_annotations_' + eqtl_ss+ '_eqtl_dataset_' + str(eqtl_dataset_iter+1) + '.npy'
	posterior_blup_gene_ld_scores = np.load(posterior_blup_gene_ld_score_file)
	# Extract predicted estimated gene ld scores (BLUP)
	prediction_blup_gene_ld_score_file = simulated_learned_gene_models_dir + simulation_name_string + 'prediction_blup_gene_related_snp_annotations_' + eqtl_ss+ '_eqtl_dataset_' + str(eqtl_dataset_iter+1) + '.npy'
	prediction_blup_gene_ld_scores = np.load(prediction_blup_gene_ld_score_file)

	# Get simulated noise
	sim_noise = np.var(noisy_gene_ld_scores,ddof=1) - np.var(sim_gene_ld_scores,ddof=1)

	# Save to global data structures
	global_n_genes.append(n_genes)
	global_sim_gene_ld_scores.append(sim_gene_ld_scores)
	global_noisy_gene_ld_scores.append(noisy_gene_ld_scores)
	global_noise_gene_ld_score_noise_est.append(noise_gene_ld_score_noise_est)
	global_noise_gene_ld_score_noise_true.append(sim_noise)
	global_posterior_known_param_gene_ld_scores.append(posterior_t_h2_t_residvar_gene_ld_scores)
	global_posterior_blup_ld_scores.append(posterior_blup_gene_ld_scores)
	global_prediction_blup_ld_scores.append(prediction_blup_gene_ld_scores)

# Convert global data structures to numpy arrays
global_n_genes = np.asarray(global_n_genes)
global_sim_gene_ld_scores = np.asarray(global_sim_gene_ld_scores)
global_noisy_gene_ld_scores = np.asarray(global_noisy_gene_ld_scores)
global_noise_gene_ld_score_noise_est = np.asarray(global_noise_gene_ld_score_noise_est)
global_noise_gene_ld_score_noise_true = np.asarray(global_noise_gene_ld_score_noise_true)
global_posterior_known_param_gene_ld_scores = np.asarray(global_posterior_known_param_gene_ld_scores)
global_posterior_blup_ld_scores = np.asarray(global_posterior_blup_ld_scores)
global_prediction_blup_ld_scores = np.asarray(global_prediction_blup_ld_scores)


###########################
# Initialize output file handle
###########################
# Print results to output
output_file = mediated_h2_results_dir + simulation_name_string + selection + '_selection_eqtl_ss_' + str(eqtl_ss)+ 'med_h2_ldsc_style_summary.txt'
t = open(output_file,'w')
# Header
t.write('sim_h2\tsim_nm_h2\tsim_med_h2')
t.write('\tmed_h2_sim_gene_ld_scores')
t.write('\tmed_h2_est_gene_ld_scores')
t.write('\tmed_h2_est_gene_ld_scores_dissattenuated_true_noise')
t.write('\tmed_h2_est_gene_ld_scores_dissattenuated_est_noise')
t.write('\tmed_h2_est_posterior_gene_ld_scores_known_gene_h2_known_residvar')
t.write('\tmed_h2_est_posterior_gene_ld_scores_blup')
t.write('\tmed_h2_est_prediction_gene_ld_scores_blup')
t.write('\n')
# First couple of simulated elements
t.write(str(sim_genetic_var) + '\t' + str(sim_nm_genetic_var) + '\t' + str(sim_med_genetic_var))


###########################
# Run S-LDSC to get mediated heritabilities
###########################
sim_average_ge_h2 = np.asarray([.05]*n_eqtl_datasets)


# true, simulated gene ld scores
med_h2_sim_gene_ld_scores = run_sldsc_to_get_mediated_heritabilities(gwas_chi_sq, var_ld_scores, global_sim_gene_ld_scores, N_gwas, n_snps, global_n_genes, sim_average_ge_h2)
t.write('\t' + ','.join(med_h2_sim_gene_ld_scores.astype(str)))


# estimated, noisy gene ld scores
med_h2_est_gene_ld_scores = run_sldsc_to_get_mediated_heritabilities(gwas_chi_sq, var_ld_scores, global_noisy_gene_ld_scores, N_gwas, n_snps, global_n_genes, sim_average_ge_h2)
t.write('\t' + ','.join(med_h2_est_gene_ld_scores.astype(str)))

# estimated, noisy gene ld scores with dissattenuation correction based on true, simulated variance
sim_noise_vec = np.hstack((np.zeros(1), global_noise_gene_ld_score_noise_true))
med_h2_est_gene_ld_scores_dissatten_true_noise = run_sldsc_to_get_mediated_heritabilities(gwas_chi_sq, var_ld_scores, global_noisy_gene_ld_scores, N_gwas, n_snps, global_n_genes, sim_average_ge_h2, noise_vec=sim_noise_vec)
t.write('\t' + ','.join(med_h2_est_gene_ld_scores_dissatten_true_noise.astype(str)))


# estimated, noisy gene ld scores with dissattenuation correction based on estimated variance
est_noise_vec = np.hstack((np.zeros(1), global_noise_gene_ld_score_noise_est))
med_h2_est_gene_ld_scores_dissatten_est_noise = run_sldsc_to_get_mediated_heritabilities(gwas_chi_sq, var_ld_scores, global_noisy_gene_ld_scores, N_gwas, n_snps, global_n_genes, sim_average_ge_h2, noise_vec=est_noise_vec)
t.write('\t' + ','.join(med_h2_est_gene_ld_scores_dissatten_est_noise.astype(str)))


# posterior gene ld scores (known gene h2, known resid var)
med_h2_posterior_t_h2_t_residvar_gene_ld_scores = run_sldsc_to_get_mediated_heritabilities(gwas_chi_sq, var_ld_scores, global_posterior_known_param_gene_ld_scores, N_gwas, n_snps, global_n_genes, sim_average_ge_h2)
t.write('\t' + ','.join(med_h2_posterior_t_h2_t_residvar_gene_ld_scores.astype(str)))

# posterior gene ld scores (BLUP)
med_h2_posterior_blup_gene_ld_scores = run_sldsc_to_get_mediated_heritabilities(gwas_chi_sq, var_ld_scores, global_posterior_blup_ld_scores, N_gwas, n_snps, global_n_genes, sim_average_ge_h2)
t.write('\t' + ','.join(med_h2_posterior_blup_gene_ld_scores.astype(str)))

# prediction gene ld scores (BLUP)
med_h2_prediction_blup_gene_ld_scores = run_sldsc_to_get_mediated_heritabilities(gwas_chi_sq, var_ld_scores, global_prediction_blup_ld_scores, N_gwas, n_snps, global_n_genes, sim_average_ge_h2)
t.write('\t' + ','.join(med_h2_prediction_blup_gene_ld_scores.astype(str)))

# Close file handles
t.write('\n')
t.close()


print(output_file)







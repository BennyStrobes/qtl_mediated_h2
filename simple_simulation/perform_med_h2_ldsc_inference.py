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
	true_causal_eqtl_effects = []
	est_noise_ratios = []
	true_noise_ratios = []
	causal_eqtl_effect_covs = []
	head_count = 0
	f = open(gene_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
		gene_snp_causal_effects = np.asarray(data[2].split(',')).astype(float)
		estimated_causal_eqtl_file = data[10]
		estimated_causal_eqtl_effects = np.load(estimated_causal_eqtl_file)
		est_noise_ratio = float(data[11])
		true_noise_ratio = float(data[12])
		causal_eqtl_effect_cov = np.load(data[13])

		# Quick error check
		if len(estimated_causal_eqtl_effects) != len(gene_snp_indices):
			continue

		causal_eqtl_effects.append(estimated_causal_eqtl_effects)
		causal_eqtl_indices.append(gene_snp_indices)
		true_causal_eqtl_effects.append(gene_snp_causal_effects)
		est_noise_ratios.append(est_noise_ratio)
		true_noise_ratios.append(true_noise_ratio)
		causal_eqtl_effect_covs.append(causal_eqtl_effect_cov)

	f.close()
	return causal_eqtl_effects, causal_eqtl_indices, true_causal_eqtl_effects, est_noise_ratios, true_noise_ratios, causal_eqtl_effect_covs


def sldsc_with_variants_and_distribution_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_effect_covs, causal_eqtl_indices):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	LD_sq = np.square(LD)

	#####################################################
	# Standardize genetic distribution of gene expression
	#####################################################
	std_causal_eqtl_effects = []
	std_causal_eqtl_covs = []
	for gg in range(n_genes):
		# LD specifically for gene
		R_gene = LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]
		# Calculate gene variance
		causal_eqtl_effects_mean = causal_eqtl_effects[gg]
		#gene_meanT_mean = np.dot(causal_eqtl_effects_mean.reshape(len(causal_eqtl_effects_mean),1), causal_eqtl_effects_mean.reshape(1,len(causal_eqtl_effects_mean)))
		#gene_mean_variance = np.sum(R_gene*gene_meanT_mean)
		gene_mean_variance = np.dot(np.dot(causal_eqtl_effects[gg], R_gene), causal_eqtl_effects[gg])
		gene_cov_variance = np.sum(R_gene*causal_eqtl_effect_covs[gg])
		tot_gene_variance = gene_mean_variance + gene_cov_variance

		# Standardize
		std_causal_eqtl_effects_mean = causal_eqtl_effects_mean/np.sqrt(tot_gene_variance)
		std_causal_eqtl_effect_cov = causal_eqtl_effect_covs[gg]/tot_gene_variance


		# Check that it sums to one
		#std_gene_mean_variance = np.dot(np.dot(std_causal_eqtl_effects_mean, R_gene), std_causal_eqtl_effects_mean)
		#std_gene_cov_variance = np.sum(R_gene*std_causal_eqtl_effect_cov)

		std_causal_eqtl_effects.append(std_causal_eqtl_effects_mean)
		std_causal_eqtl_covs.append(std_causal_eqtl_effect_cov)

	#####################################################
	# Create gene LD scores
	#####################################################
	gene_ld_scores = np.zeros(n_snps)
	for gg in range(n_genes):
		causal_eqtl_effects_mean = std_causal_eqtl_effects[gg]
		gene_meanT_mean = np.dot(causal_eqtl_effects_mean.reshape(len(causal_eqtl_effects_mean),1), causal_eqtl_effects_mean.reshape(1,len(causal_eqtl_effects_mean)))
		
		sub_ld = LD[:, causal_eqtl_indices[gg]]
		for kk in range(n_snps):
			gene_ld_scores[kk] = gene_ld_scores[kk] + np.dot(np.dot(sub_ld[kk,:], gene_meanT_mean + std_causal_eqtl_covs[gg]), sub_ld[kk,:])

		#gene_score = np.diag(np.dot(np.dot(LD[:, causal_eqtl_indices[gg]], (gene_meanT_mean + std_causal_eqtl_covs[gg])), LD[causal_eqtl_indices[gg], :]))

	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores)))
	var_ld_scores = np.hstack((np.sum(LD_sq,axis=0)))
	joint_ld_scores = np.transpose(np.vstack((var_ld_scores,gene_ld_scores)))

	############################
	# Run SLDSC
	############################
	model = LinearRegression(fit_intercept=False)
	ldsc_fit = model.fit(joint_ld_scores, chi_sq-1)

	nm_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas
	med_h2 = ldsc_fit.coef_[1]*n_genes/N_gwas

	return nm_h2, med_h2



def sldsc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, noise_ratios):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	LD_sq = np.square(LD)

	#####################################################
	# Standardize genetically predicted gene expression
	#####################################################
	std_causal_eqtl_effects = []
	for gg in range(n_genes):
		# Calculate gene variance
		gene_variance = np.dot(np.dot(causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), causal_eqtl_effects[gg])

		# Standardize gene effects
		std_causal_eqtl_effect = causal_eqtl_effects[gg]/np.sqrt(gene_variance)

		# Append to array
		std_causal_eqtl_effects.append(std_causal_eqtl_effect)


	##########################
	# Create gene-variant LD mat
	##########################
	gene_variant_ld = np.zeros((n_genes, n_snps))
	for gg in range(n_genes):
		gene_variant_ld[gg,:] = np.dot(LD[:, causal_eqtl_indices[gg]], std_causal_eqtl_effects[gg])
	# Square LD
	gene_variant_ld_sq = np.square(gene_variant_ld)

	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores)))
	var_ld_scores = np.hstack((np.sum(LD_sq,axis=0)))
	gene_ld_scores = np.hstack((np.sum(gene_variant_ld_sq,axis=0)))
	joint_ld_scores = np.transpose(np.vstack((var_ld_scores,gene_ld_scores)))
	
	############################
	# Run SLDSC
	############################
	model = LinearRegression(fit_intercept=False)
	ldsc_fit = model.fit(joint_ld_scores, chi_sq-1)

	nm_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas
	med_h2 = ldsc_fit.coef_[1]*n_genes/N_gwas

	return nm_h2, med_h2





def sldsc_with_variants_and_noisy_genetic_gene_expression(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, noise_ratios):
	# Get number of genes
	n_genes = len(causal_eqtl_indices)
	# Get number of snps
	n_snps = len(gwas_z_scores)

	# Square LD
	LD_sq = np.square(LD)

	#####################################################
	# Standardize genetically predicted gene expression
	#####################################################
	std_causal_eqtl_effects = []
	for gg in range(n_genes):
		# Calculate gene variance
		gene_variance = np.dot(np.dot(causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), causal_eqtl_effects[gg])

		# Standardize gene effects
		std_causal_eqtl_effect = causal_eqtl_effects[gg]/np.sqrt(gene_variance)

		# Append to array
		std_causal_eqtl_effects.append(std_causal_eqtl_effect)

	##########################
	# Compute TWAS Z scores
	##########################
	twas_z_scores = []
	for gg in range(n_genes):
		tmp_numerator = np.dot(std_causal_eqtl_effects[gg], gwas_z_scores[causal_eqtl_indices[gg]])
		# NOTE: Denomenator should always be equal to one cause we standardized. but for completeness..
		tmp_denomenator = np.sqrt(np.dot(np.dot(std_causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg], :][:, causal_eqtl_indices[gg]]), std_causal_eqtl_effects[gg])) 
		twas_z = tmp_numerator/tmp_denomenator
		# Add to global array
		twas_z_scores.append(twas_z)
	twas_z_scores = np.asarray(twas_z_scores)


	##########################
	# Create gene-gene LD mat
	##########################
	gene_gene_ld = np.zeros((n_genes, n_genes))
	for gg in range(n_genes):
		for kk in range(n_genes):
			corry = np.dot(np.dot(std_causal_eqtl_effects[gg], LD[causal_eqtl_indices[gg],:][:,causal_eqtl_indices[kk]]), std_causal_eqtl_effects[kk])
			gene_gene_ld[gg,kk] = corry
			gene_gene_ld[kk,gg] = corry


	##########################
	# Correct gene-gene LD mat for gene noise
	##########################
	noise_ratios = np.asarray(noise_ratios)
	# Set noise ratios larger than 1 to 1
	noise_ratios[noise_ratios > 1.0] = 1.0
	# Correction
	gene_gene_ld = gene_gene_ld - np.diag(noise_ratios)
	# Squared LD
	gene_gene_ld_sq = np.square(gene_gene_ld)

	##########################
	# Create gene-variant LD mat
	##########################
	gene_variant_ld = np.zeros((n_genes, n_snps))
	for gg in range(n_genes):
		gene_variant_ld[gg,:] = np.dot(LD[:, causal_eqtl_indices[gg]], std_causal_eqtl_effects[gg])
	# Square LD
	gene_variant_ld_sq = np.square(gene_variant_ld)

	##########################
	# Prepare LDSC regression data
	##########################
	chi_sq = np.hstack((np.square(gwas_z_scores), np.square(twas_z_scores)))
	var_ld_scores = np.hstack((np.sum(LD_sq,axis=0), np.sum(gene_variant_ld_sq,axis=1)))
	gene_ld_scores = np.hstack((np.sum(gene_variant_ld_sq,axis=0), np.sum(gene_gene_ld_sq,axis=0)))
	joint_ld_scores = np.transpose(np.vstack((var_ld_scores,gene_ld_scores)))
	
	############################
	# Run SLDSC
	############################
	model = LinearRegression(fit_intercept=False)
	ldsc_fit = model.fit(joint_ld_scores, chi_sq-1)

	nm_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas
	med_h2 = ldsc_fit.coef_[1]*n_genes/N_gwas

	return nm_h2, med_h2




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
causal_eqtl_effects, causal_eqtl_indices, true_causal_eqtl_effects, est_noise_ratios, true_noise_ratios, causal_eqtl_effect_covs = extract_previously_estimsted_causal_eqtl_effects(gene_summary_file)


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



# Load in Genotype file
# LD file
ld_file = processed_genotype_data_dir + 'gwas_genotype_LD_1.npy'
LD = np.load(ld_file)

# Genotype file
geno_file = processed_genotype_data_dir + 'gwas_genotype_1.npy'
#genotype_mat = np.load(geno_file)


# Simulated Trait
trait_file = simulated_gwas_data_dir + simulation_name_string + 'simulated_trait.txt'
trait_vec = load_in_trait_file(trait_file)

# Gwas sample size
N_gwas = len(trait_vec)

# Run SLDSC with variants and true genetic gene expression (no true noise)
nm_h2_true_gene_no_noise, med_h2_true_gene_no_noise = sldsc_with_variants_and_noisy_genetic_gene_expression(gwas_z_scores, N_gwas, LD, true_causal_eqtl_effects, causal_eqtl_indices, np.zeros(len(est_noise_ratios)))

# Run SLDSC with variants and genetically predicted gene expression (accounting for zero noise)
nm_h2_pred_gene_no_noise, med_h2_pred_gene_no_noise = sldsc_with_variants_and_noisy_genetic_gene_expression(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, np.zeros(len(est_noise_ratios)))

# Run SLDSC with variants and genetically predicted gene expression (accounting for estimated noise)
nm_h2_pred_gene_pred_noise, med_h2_pred_gene_pred_noise = sldsc_with_variants_and_noisy_genetic_gene_expression(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, est_noise_ratios)

# Run SLDSC with variants and genetically predicted gene expression (accounting for true noise)
nm_h2_pred_gene_true_noise, med_h2_pred_gene_true_noise = sldsc_with_variants_and_noisy_genetic_gene_expression(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, true_noise_ratios)

# Run SLDSC with variants and genetically predicted gene expression (accounting for zero noise) and without modeling twas chi-squared stats
nm_h2_pred_gene_no_noise_no_twas_chi_sq, med_h2_pred_gene_no_noise_no_twas_chi_sq = sldsc_with_variants_and_noisy_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_indices, np.zeros(len(est_noise_ratios)))

# Run SLDSC with variants and genetic distribution of gene expression and without modeling twas chi-squared
nm_h2_distr_gene_no_twas_chi_sq, med_h2_distr_gene_no_twas_chi_sq = sldsc_with_variants_and_distribution_genetic_gene_expression_no_twas_chi_sq(gwas_z_scores, N_gwas, LD, causal_eqtl_effects, causal_eqtl_effect_covs, causal_eqtl_indices)


# Print results to output
output_file = mediated_h2_results_dir + simulation_name_string + 'eqtl_ss_' + eqtl_ss + '_med_h2_ldsc_style_summary.txt'
t = open(output_file,'w')
# Header
t.write('sim_h2\tsim_nm_h2\tsim_med_h2\t')
t.write('nm_h2_true_gene_no_noise\tmed_h2_true_gene_no_noise\t')
t.write('nm_h2_pred_gene_no_noise\tmed_h2_pred_gene_no_noise\t')
t.write('nm_h2_pred_gene_pred_noise\tmed_h2_pred_gene_pred_noise\t')
t.write('nm_h2_pred_gene_true_noise\tmed_h2_pred_gene_true_noise\t')
t.write('nm_h2_pred_gene_no_noise_no_twas_chi_sq\tmed_h2_pred_gene_no_noise_no_twas_chi_sq\t')
t.write('nm_h2_distr_gene_no_twas_chi_sq\tmed_h2_distr_gene_no_twas_chi_sq\n')

# Line
t.write(str(sim_genetic_var) + '\t' + str(sim_nm_genetic_var) + '\t' + str(sim_med_genetic_var) + '\t')
t.write(str(nm_h2_true_gene_no_noise) + '\t' + str(med_h2_true_gene_no_noise) + '\t')
t.write(str(nm_h2_pred_gene_no_noise) + '\t' + str(med_h2_pred_gene_no_noise) + '\t')
t.write(str(nm_h2_pred_gene_pred_noise) + '\t' + str(med_h2_pred_gene_pred_noise) + '\t')
t.write(str(nm_h2_pred_gene_true_noise) + '\t' + str(med_h2_pred_gene_true_noise) + '\t')
t.write(str(nm_h2_pred_gene_no_noise_no_twas_chi_sq) + '\t' + str(med_h2_pred_gene_no_noise_no_twas_chi_sq) + '\t')
t.write(str(nm_h2_distr_gene_no_twas_chi_sq) + '\t' + str(med_h2_distr_gene_no_twas_chi_sq) + '\n')


t.close()

print(output_file)







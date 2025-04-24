import sys
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from bgen import BgenReader

def load_in_alt_allele_genotype_dosage_mat(bfile, window_indices, ref_alt_alleles):
	dosages = []

	for window_index in window_indices:
		var = bfile[window_index]
		dosage = var.minor_allele_dosage
		ma = var.minor_allele

		index_ref_alt_allele = ref_alt_alleles[window_index]
		# Flip dosage if alt-allele is not equal to minor allele
		if index_ref_alt_allele[1] != ma:
			# Quick error check
			if ma != index_ref_alt_allele[0]:
				print('assumptino eroror')
				pdb.set_trace()
			# Flip dosage
			dosage = 2.0 - dosage

		# Append snp dosage to global array
		dosages.append(dosage)

	# Convert to 2d matrix
	dosages = np.asarray(dosages)

	return np.transpose(dosages)


def standardize_genotype_dosage_matrix(genotype_dosage):
	# Quick error checking to make sure there do not exist missing entries
	n_missing = np.sum(np.isnan(genotype_dosage))
	if n_missing != 0:
		print('assumption eroror')
		pdb.set_trace()

	# Now standardize genotype of each snp
	n_snps = genotype_dosage.shape[1]
	# Initialize standardize genotype dosage matrix
	std_genotype_dosage = np.copy(genotype_dosage)
	for snp_iter in range(n_snps):
		# Standardize
		std_genotype_dosage[:, snp_iter] = (genotype_dosage[:,snp_iter] - np.mean(genotype_dosage[:,snp_iter]))/np.std(genotype_dosage[:,snp_iter])

	return std_genotype_dosage

def load_in_ref_alt_allele_arr(pvar_file):
	ref_alt_alleles = []
	f = open(pvar_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 5:
			print('assumptino erooror')
			pdb.set_trace()
		ref_allele = data[3]
		alt_allele = data[4]
		if ref_allele == alt_allele:
			print('assumptino eororor')
			pdb.set_trace()
		ref_alt_alleles.append((ref_allele, alt_allele))
	f.close()
	return ref_alt_alleles
def get_noise_scaler(genetic_gene, random_noise):
	scaler_values = np.arange(0,4,.001)
	diffs = []
	for scaler_value in scaler_values:
		emperical_var = genetic_gene + (random_noise*scaler_value)
		diffs.append(np.abs(np.var(emperical_var) - 1.0))
	diffs = np.asarray(diffs)
	best_diff = np.argmin(diffs)
	return scaler_values[best_diff]


def get_genotype_association_sumstats(trait, genotype_mat):
	n_snps = genotype_mat.shape[1]
	betas = []
	beta_vars = []

	for snp_iter in range(n_snps):
		olser = sm.OLS(trait, genotype_mat[:,snp_iter]).fit()
		beta = olser.params[0]
		beta_se = olser.bse[0]
		betas.append(beta)
		beta_vars.append(np.square(beta_se))

	return np.asarray(betas), np.asarray(beta_vars)

def simulate_gwas_data_and_get_gwas_sum_stats(stand_gwas_geno, causal_eqtl_effect_sizes, nm_h2, med_h2, architecture):
	##################################
	# Simulate causal nm gwas effects
	n_snps = stand_gwas_geno.shape[1]
	n_inidi = stand_gwas_geno.shape[0]
	if architecture == 'polygenic':
		causal_gwas_nm_betas = np.random.normal(loc=0, scale=np.sqrt(nm_h2/n_snps), size=n_snps)
	nm_genetic_trait = np.dot(stand_gwas_geno, causal_gwas_nm_betas)

	##################################
	# Simulate causal mediated gene trait effects
	genetic_gene = np.dot(stand_gwas_geno, causal_eqtl_effect_sizes)
	standardized_genetic_gene = genetic_gene/np.std(genetic_gene)
	#gene_trait_effect = np.random.normal(loc=0, scale=np.sqrt(med_h2))
	gene_trait_effect = np.sqrt(med_h2)
	med_genetic_trait = standardized_genetic_gene*gene_trait_effect

	# Simulate noise
	genetic_trait = med_genetic_trait + nm_genetic_trait
	# Simulate expression
	noise = np.random.normal(loc=0, scale=np.sqrt(1.0-np.var(genetic_trait)), size=n_inidi)
	noise_scaler = get_noise_scaler(genetic_trait, noise)
	trait = genetic_trait + noise*noise_scaler

	trait = (trait - np.mean(trait))/np.std(trait)

	##############################
	# Get gwas effect sizes
	gwas_beta, gwas_beta_vars = get_genotype_association_sumstats(trait, stand_gwas_geno)


	return causal_gwas_nm_betas, gene_trait_effect, gwas_beta, gwas_beta_vars, np.var(nm_genetic_trait), np.var(med_genetic_trait)


def simulate_eqtl_data_and_get_eqtl_sum_stats(stand_eqtl_geno, gene_h2, architecture):
	##################################
	# Simulate causal eqtl effects
	n_snps = stand_eqtl_geno.shape[1]
	n_inidi = stand_eqtl_geno.shape[0]
	if architecture == 'polygenic':
		causal_eqtl_betas = np.zeros(n_snps)
		causal_eqtl_betas[:(int(n_snps/2))] = np.random.normal(loc=0, scale=np.sqrt(gene_h2/(n_snps/2)), size=int(n_snps/2))

	##############################
	# Simulate Gene expression
	# Get genetic component of gene expression
	genetic_gene = np.dot(stand_eqtl_geno, causal_eqtl_betas)
	# Simulate expression
	noise = np.random.normal(loc=0, scale=np.sqrt(1.0-np.var(genetic_gene)), size=n_inidi)
	noise_scaler = get_noise_scaler(genetic_gene, noise)
	expression = genetic_gene + noise*noise_scaler
	expression = (expression - np.mean(expression))/np.std(expression)

	##############################
	# Get eqtl effect sizes
	eqtl_beta, eqtl_beta_vars = get_genotype_association_sumstats(expression, stand_eqtl_geno)

	return causal_eqtl_betas, eqtl_beta, eqtl_beta_vars, np.var(genetic_gene)

def run_single_simulation(simulation_genotype_dir, output_root_dir, window_name, nm_h2, med_h2, gene_h2, architecture, eqtl_sample_size, sim_iter):
	# Load in some pre-computed genotype data
	gwas_ld_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_window_' + window_name + '_ld.npy'
	gwas_ld = np.load(gwas_ld_file)
	rs_id_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_window_' + window_name + '_rsids.npy'
	rsids = np.load(rs_id_file)
	gwas_q_mat_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_window_' + window_name + '_Q_mat.npy'
	gwas_q_mat = np.load(gwas_q_mat_file)
	gwas_w_mat_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_window_' + window_name + '_w_premult_mat.npy'
	gwas_w_mat = np.load(gwas_w_mat_file)
	window_position_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_window_' + window_name + '_window_positions.npy'
	window_positions = np.load(window_position_file)

	Q_ld = np.dot(gwas_q_mat, np.transpose(gwas_q_mat))
	
	# Load in eQTL genotype
	eqtl_plink_stem = simulation_genotype_dir + 'simulated_eqtl_' + str(eqtl_sample_size) + '_data_' + str(1)  # Genotype directory
	# Load in genotype object
	eqtl_genotype_obj = BgenReader(eqtl_plink_stem + '.bgen')
	# Load in ref-alt alleles
	eqtl_ref_alt_alleles = load_in_ref_alt_allele_arr(eqtl_plink_stem + '.pvar')
	# load in actual genotype data
	eqtl_geno = load_in_alt_allele_genotype_dosage_mat(eqtl_genotype_obj, window_positions, eqtl_ref_alt_alleles)
	stand_eqtl_geno = standardize_genotype_dosage_matrix(eqtl_geno)

	# Simulate eQTL data and get eqtl sum-stats
	causal_eqtl_effect_sizes, eqtl_beta, eqtl_beta_var, sim_genetic_gene_variance = simulate_eqtl_data_and_get_eqtl_sum_stats(stand_eqtl_geno, gene_h2, architecture)

	# Load in GWAS genotype
	gwas_plink_stem = simulation_genotype_dir + 'simulated_gwas_data_' + str(1)  # Genotype directory
	# Load in genotype object
	genotype_obj = BgenReader(gwas_plink_stem + '.bgen')
	# Load in ref-alt alleles
	ref_alt_alleles = load_in_ref_alt_allele_arr(gwas_plink_stem + '.pvar')
	# load in actual genotype data
	gwas_geno = load_in_alt_allele_genotype_dosage_mat(genotype_obj, window_positions, ref_alt_alleles)
	stand_gwas_geno = standardize_genotype_dosage_matrix(gwas_geno)

	# Simulate gwas data and get eqtl sum-stats
	causal_gwas_variant_trait_effect_sizes, causal_gwas_gene_trait_effect_sizes, gwas_beta, gwas_beta_var, sim_nm_h2, sim_med_h2 = simulate_gwas_data_and_get_gwas_sum_stats(stand_gwas_geno, causal_eqtl_effect_sizes, nm_h2, med_h2, architecture)

	'''
	np.save('gwas_beta.npy', gwas_beta)
	np.save('eqtl_beta.npy', eqtl_beta)
	np.save('Q_ld.npy', Q_ld)

	# temp loading
	gwas_beta = np.load('gwas_beta.npy')
	eqtl_beta = np.load('eqtl_beta.npy')
	Q_ld = np.load('Q_ld.npy')
	'''


	# Convert gwas and eqtl sum stats to genotype pc space
	pc_gwas_beta = np.dot(gwas_w_mat, gwas_beta)
	pc_eqtl_beta = np.dot(gwas_w_mat, eqtl_beta)

	# Regress pc_gwas_beta on pc_eqtl_beta
	olser = sm.OLS(pc_gwas_beta, pc_eqtl_beta).fit()
	beta_effect = olser.params[0]
	resid_pc_gwas_beta = pc_gwas_beta - (pc_eqtl_beta*beta_effect)

	
	n_snps = len(rsids)
	# Run standard ld score regression
	standard_ldsc_fit = sm.OLS(np.square(pc_gwas_beta) - (1/100000.0), np.diag(Q_ld)).fit()
	ldsc_h2_est = n_snps*standard_ldsc_fit.params[0]

	# do conditional ld score regression
	ld_scores = np.diag(Q_ld)
	alt_ld_scores = np.square(ld_scores)/((gene_h2*ld_scores/n_snps) + (1.0/eqtl_sample_size))

	joint_ld_scores = np.transpose(np.vstack((ld_scores, alt_ld_scores)))

	conditional_ldsc_fit = sm.OLS(np.square(resid_pc_gwas_beta) - (1/100000.0), joint_ld_scores).fit()
	c_ldsc_h2_est = n_snps*conditional_ldsc_fit.params[0]
	c_ldsc_med_h2_est = (-(conditional_ldsc_fit.params[1])*np.square(n_snps))/gene_h2

	print(sim_iter)
	print(sim_med_h2)
	print(c_ldsc_med_h2_est)

	return


def run_single_simulation2(simulation_genotype_dir, output_root_dir, window_name, nm_h2, med_h2, gene_h2, architecture, eqtl_sample_size, sim_iter):
	# Load in some pre-computed genotype data
	gwas_ld_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_window_' + window_name + '_ld.npy'
	gwas_ld = np.load(gwas_ld_file)
	rs_id_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_window_' + window_name + '_rsids.npy'
	rsids = np.load(rs_id_file)
	gwas_q_mat_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_window_' + window_name + '_Q_mat.npy'
	gwas_q_mat = np.load(gwas_q_mat_file)
	gwas_w_mat_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_window_' + window_name + '_w_premult_mat.npy'
	gwas_w_mat = np.load(gwas_w_mat_file)
	window_position_file = simulation_genotype_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld_window_' + window_name + '_window_positions.npy'
	window_positions = np.load(window_position_file)

	Q_ld = np.dot(gwas_q_mat, np.transpose(gwas_q_mat))
	
	# Load in eQTL genotype
	eqtl_plink_stem = simulation_genotype_dir + 'simulated_eqtl_' + str(eqtl_sample_size) + '_data_' + str(1)  # Genotype directory
	# Load in genotype object
	eqtl_genotype_obj = BgenReader(eqtl_plink_stem + '.bgen')
	# Load in ref-alt alleles
	eqtl_ref_alt_alleles = load_in_ref_alt_allele_arr(eqtl_plink_stem + '.pvar')
	# load in actual genotype data
	eqtl_geno = load_in_alt_allele_genotype_dosage_mat(eqtl_genotype_obj, window_positions, eqtl_ref_alt_alleles)
	stand_eqtl_geno = standardize_genotype_dosage_matrix(eqtl_geno)

	# Simulate eQTL data and get eqtl sum-stats
	causal_eqtl_effect_sizes, eqtl_beta, eqtl_beta_var, sim_genetic_gene_variance = simulate_eqtl_data_and_get_eqtl_sum_stats(stand_eqtl_geno, gene_h2, architecture)

	# Load in GWAS genotype
	gwas_plink_stem = simulation_genotype_dir + 'simulated_gwas_data_' + str(1)  # Genotype directory
	# Load in genotype object
	genotype_obj = BgenReader(gwas_plink_stem + '.bgen')
	# Load in ref-alt alleles
	ref_alt_alleles = load_in_ref_alt_allele_arr(gwas_plink_stem + '.pvar')
	# load in actual genotype data
	gwas_geno = load_in_alt_allele_genotype_dosage_mat(genotype_obj, window_positions, ref_alt_alleles)
	stand_gwas_geno = standardize_genotype_dosage_matrix(gwas_geno)

	# Simulate gwas data and get eqtl sum-stats
	causal_gwas_variant_trait_effect_sizes, causal_gwas_gene_trait_effect_sizes, gwas_beta, gwas_beta_var, sim_nm_h2, sim_med_h2 = simulate_gwas_data_and_get_gwas_sum_stats(stand_gwas_geno, causal_eqtl_effect_sizes, nm_h2, med_h2, architecture)

	'''
	np.save('gwas_beta.npy', gwas_beta)
	np.save('eqtl_beta.npy', eqtl_beta)
	np.save('Q_ld.npy', Q_ld)

	# temp loading
	gwas_beta = np.load('gwas_beta.npy')
	eqtl_beta = np.load('eqtl_beta.npy')
	Q_ld = np.load('Q_ld.npy')
	'''


	# Convert gwas and eqtl sum stats to genotype pc space
	pc_gwas_beta = np.dot(gwas_w_mat, gwas_beta)
	pc_eqtl_beta = np.dot(gwas_w_mat, eqtl_beta)

	# Regress pc_gwas_beta on pc_eqtl_beta
	olser = sm.OLS(pc_gwas_beta, pc_eqtl_beta).fit()
	beta_effect = olser.params[0]
	resid_pc_gwas_beta = pc_gwas_beta - (pc_eqtl_beta*beta_effect)

	
	n_snps = len(rsids)
	# Run standard ld score regression
	standard_ldsc_fit = sm.OLS(np.square(pc_gwas_beta) - (1/100000.0), np.diag(Q_ld)).fit()
	ldsc_h2_est = n_snps*standard_ldsc_fit.params[0]


	# do conditional ld score regression
	ld_scores = np.diag(Q_ld)


	term1_numerator = np.square(ld_scores)*np.square(pc_eqtl_beta)/np.square(n_snps)
	term1_denominator = ((gene_h2/np.square(n_snps))*np.square(ld_scores)) + ((2.0*ld_scores)/(n_snps*eqtl_sample_size)) + (1.0/(np.square(eqtl_sample_size)*gene_h2))

	term2_numerator = np.square(ld_scores)/np.square(n_snps)
	term2_denominator = (ld_scores/n_snps) + (1.0/(eqtl_sample_size*gene_h2))


	alt_ld_scores = (term1_numerator/term1_denominator) - (term2_numerator/term2_denominator)


	joint_ld_scores = np.transpose(np.vstack((ld_scores, alt_ld_scores)))

	conditional_ldsc_fit = sm.OLS(np.square(pc_gwas_beta) - (1/100000.0), joint_ld_scores).fit()
	c_ldsc_h2_est = n_snps*conditional_ldsc_fit.params[0]
	c_ldsc_med_h2_est = conditional_ldsc_fit.params[1]

	print(sim_iter)
	print(sim_med_h2)
	print(c_ldsc_med_h2_est)

	return c_ldsc_med_h2_est





#######################
# Command line args
#######################
simulation_genotype_dir = sys.argv[1]
output_root_dir = sys.argv[2]
simulation_number = sys.argv[3]

np.random.seed(int(simulation_number)+2)

#################
# simulation parameters
window_name = '10583_49894177'
nm_h2 = .3
med_h2 = .2
gene_h2 = .1
architecture = 'polygenic'


eqtl_sample_size = 1000

ests = []
for itera in range(200):
	#run_single_simulation(simulation_genotype_dir, output_root_dir, window_name, nm_h2, med_h2, gene_h2, architecture, eqtl_sample_size, int(simulation_number))
	c_ldsc_med_h2_est = run_single_simulation2(simulation_genotype_dir, output_root_dir, window_name, nm_h2, med_h2, gene_h2, architecture, eqtl_sample_size, int(simulation_number))

	ests.append(c_ldsc_med_h2_est)

ests = np.asarray(ests)

meany = np.mean(ests)
error = np.std(ests)/np.sqrt(len(ests))

print('Mean')
print(np.mean(ests))
print('CI')
print('[ ' + str(meany - 1.96*error) + ', ' + str(meany + 1.96*error) + ']')


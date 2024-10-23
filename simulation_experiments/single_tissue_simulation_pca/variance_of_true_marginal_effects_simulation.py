import numpy as np
import os
import sys
import pdb
from bgen import BgenReader
import statsmodels.api as sm




def load_in_genotype_data(genotype_obj):
	dosages = []

	for window_index in range(len(genotype_obj)):
		var = genotype_obj[window_index]
		dosage = var.minor_allele_dosage
		dosages.append(dosage)

	dosages = np.asarray(dosages)

	return dosages


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


def run_eiv_decomp_ld_simulation(simulation_iter, stand_eqtl_genotype, stand_gwas_genotype, eqtl_ld, gwas_ld, w_premult, q_mat,per_snp_h2=1e-5):
	# Extract relevent parameters
	n_snps, n_eqtl_samp = stand_eqtl_genotype.shape
	n_snps2, n_gwas_samp = stand_gwas_genotype.shape

	# Simulate causal betas
	betas = np.random.normal(loc=0.0, scale=np.sqrt(per_snp_h2), size=n_snps)

	# Get marginal gwas and eqtl betas
	true_marginal_gwas_betas = np.dot(gwas_ld, betas)
	true_marginal_eqtl_betas = np.dot(eqtl_ld, betas)

	# Get genetic components
	gwas_genetic_component = np.dot(np.transpose(stand_gwas_genotype), betas)
	eqtl_genetic_component = np.dot(np.transpose(stand_eqtl_genotype), betas)


	true_pc_marginal_gwas_betas = np.dot(w_premult, true_marginal_gwas_betas)
	true_pc_marginal_eqtl_betas = np.dot(w_premult, true_marginal_eqtl_betas)

	eqtl_q_mat = np.dot(w_premult,eqtl_ld)

	Q_variances = (1.0 - np.square(eqtl_q_mat))/(n_eqtl_samp-2)

	expected_var = np.sum(Q_variances,axis=1)*per_snp_h2

	squared_diffs = np.square(true_pc_marginal_gwas_betas - true_pc_marginal_eqtl_betas)


	print(np.var(true_pc_marginal_gwas_betas - true_pc_marginal_eqtl_betas))
	print(np.mean(expected_var))

	return


def run_eiv_decomp_ld_simulation2(simulation_iter, stand_eqtl_genotype, stand_gwas_genotype, eqtl_ld, gwas_ld, w_premult, q_mat,per_snp_h2=1e-5):
	# Extract relevent parameters
	n_snps, n_eqtl_samp = stand_eqtl_genotype.shape
	n_snps2, n_gwas_samp = stand_gwas_genotype.shape

	# Simulate causal betas
	betas = np.zeros(n_snps)
	betas[:200] = np.random.normal(loc=0.0, scale=np.sqrt(per_snp_h2), size=200)

	# Get marginal gwas and eqtl betas
	true_marginal_gwas_betas = np.dot(gwas_ld, betas)
	true_marginal_eqtl_betas = np.dot(eqtl_ld, betas)

	# Get genetic components
	gwas_genetic_component = np.dot(np.transpose(stand_gwas_genotype), betas)
	eqtl_genetic_component = np.dot(np.transpose(stand_eqtl_genotype), betas)

	true_pc_marginal_gwas_betas = np.dot(w_premult, true_marginal_gwas_betas)
	true_pc_marginal_eqtl_betas = np.dot(w_premult, true_marginal_eqtl_betas)


	eqtl_q_mat = np.dot(w_premult,eqtl_ld)

	eqtl_q_mat2 = eqtl_q_mat[:,:200]

	Q_variances = (1.0 - np.square(eqtl_q_mat2))/(n_eqtl_samp-2)

	expected_var = np.sum(Q_variances,axis=1)*per_snp_h2

	squared_diffs = np.square(true_pc_marginal_gwas_betas - true_pc_marginal_eqtl_betas)


	print(np.var(true_pc_marginal_gwas_betas - true_pc_marginal_eqtl_betas))
	print(np.mean(expected_var))

	return




def run_simulation(simulation_iter, stand_eqtl_genotype, stand_gwas_genotype, eqtl_ld, gwas_ld, per_snp_h2=1e-5):
	# Extract relevent parameters
	n_snps, n_eqtl_samp = stand_eqtl_genotype.shape
	n_snps2, n_gwas_samp = stand_gwas_genotype.shape

	# Simulate causal betas
	betas = np.random.normal(loc=0.0, scale=np.sqrt(per_snp_h2), size=n_snps)

	# Get marginal gwas and eqtl betas
	true_marginal_gwas_betas = np.dot(gwas_ld, betas)
	true_marginal_eqtl_betas = np.dot(eqtl_ld, betas)

	# Get genetic components
	gwas_genetic_component = np.dot(np.transpose(stand_gwas_genotype), betas)
	eqtl_genetic_component = np.dot(np.transpose(stand_eqtl_genotype), betas)

	# Double check that marginal gwas and eqtl betas are correct
	true_marginal_gwas_betas2 = np.dot(stand_gwas_genotype, gwas_genetic_component)/n_gwas_samp
	true_marginal_eqtl_betas2 = np.dot(stand_eqtl_genotype, eqtl_genetic_component)/n_eqtl_samp
	if np.max(np.abs(true_marginal_gwas_betas - true_marginal_gwas_betas2)) > 1e-5:
		print('assumption erroror')
		pdb.set_trace()
	if np.max(np.abs(true_marginal_eqtl_betas - true_marginal_eqtl_betas2)) > 1e-5:
		print('assumption erroror')
		pdb.set_trace()

	delta_marginal_betas = true_marginal_gwas_betas - true_marginal_eqtl_betas

	ld_variances = (1.0 - np.square(eqtl_ld))/(n_eqtl_samp-2)
	

	#print(np.var(delta_marginal_betas,ddof=0))
	#print(np.sort(np.sum(ld_variances,axis=0)*per_snp_h2))

	emperical_var = np.var(delta_marginal_betas,ddof=0)
	expected_var = np.mean(np.sum(ld_variances,axis=0)*per_snp_h2)

	return emperical_var, expected_var

def compute_lambda_thresh(lambdas, rho_thresh):
	totaler = np.sum(lambdas)
	cur_total = 0
	lambda_thresh = -1
	for lambda_val in -np.sort(-lambdas):
		cur_total = cur_total + lambda_val
		if cur_total/totaler > rho_thresh:
			if lambda_thresh == -1:
				lambda_thresh = lambda_val


	if lambda_thresh == -1:
		print('assumption eroror')
		pdb.set_trace()

	return lambda_thresh

def eigenvalue_decomp_of_ld(ld_mat):
	# EIG value decomp
	lambdas_full, U_full = np.linalg.eig(ld_mat)
	non_negative_components = lambdas_full > 0.0
	lambdas = lambdas_full[non_negative_components]
	U = U_full[:, non_negative_components]
	real_components = np.iscomplex(lambdas) == False
	lambdas = lambdas[real_components]
	U = U[:, real_components]
	if np.sum(np.iscomplex(lambdas)) > 0:
		print('assumption eroror')
		pdb.set_trace()
	lambdas = lambdas.astype(float)
	U = U.astype(float)

	rho_thresh = 0.99
	lambda_thresh = compute_lambda_thresh(lambdas, rho_thresh)
	thresh_components = lambdas >= lambda_thresh
	lambdas = lambdas[thresh_components]
	U = U[:, thresh_components]


	# Note that reconstruction of ld_mat is achieved with np.dot(np.dot(U, np.diag(lambdas)), np.transpose(U))

	# Compute some relevent quantities
	Q_mat = np.dot(np.diag(lambdas**(.5)), np.transpose(U))
	w_premult = np.dot(np.diag(lambdas**(-.5)), np.transpose(U))

	return w_premult, Q_mat









#####################
# Command line args
####################
true_marginal_effect_simulation_dir = sys.argv[1]  # Output dir
chrom_num = sys.argv[2]
processed_genotype_data_dir = sys.argv[3]


# Load in genotype data
'''
gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype directory
genotype_obj = BgenReader(gwas_plink_stem + '.bgen')
genotype_dosage = load_in_genotype_data(genotype_obj)
# Split into eqtl and genotype
gwas_dosage = genotype_dosage[:,100:]
eqtl_dosage = genotype_dosage[:,:100]
# Standardize
stand_eqtl_genotype = np.transpose(standardize_genotype_dosage_matrix(np.transpose(eqtl_dosage)))
stand_gwas_genotype = np.transpose(standardize_genotype_dosage_matrix(np.transpose(gwas_dosage)))

np.save('eqtl_geno.npy', stand_eqtl_genotype)
np.save('gwas_geno.npy', stand_gwas_genotype)
'''
stand_eqtl_genotype = np.load('eqtl_geno.npy')
stand_gwas_genotype = np.load('gwas_geno.npy')
'''
# Compute ld
eqtl_ld = np.corrcoef(stand_eqtl_genotype)
gwas_ld = np.corrcoef(stand_gwas_genotype)
np.save('eqtl_ld.npy', eqtl_ld)
np.save('gwas_ld.npy', gwas_ld)
'''

eqtl_ld = np.load('eqtl_ld.npy')
gwas_ld = np.load('gwas_ld.npy')

# Compute eigenvalue decomposition of LD
'''
w_premult, q_mat = eigenvalue_decomp_of_ld(gwas_ld)
np.save('eiv_decomp_w_premult.npy', w_premult)
np.save('eiv_decomp_q_mat.npy', q_mat)
'''
w_premult = np.load('eiv_decomp_w_premult.npy')
q_mat = np.load('eiv_decomp_q_mat.npy')


####################### 
# EIV-decomp LD SIMULATION
for simulation_iter in range(2000):
	print(simulation_iter)
	run_eiv_decomp_ld_simulation2(simulation_iter, stand_eqtl_genotype, stand_gwas_genotype, eqtl_ld, gwas_ld, w_premult, q_mat)

####################### 
# LD SIMULATION
'''
emperical_vars = []
expected_vars = []
for simulation_iter in range(2000):
	print(simulation_iter)
	emperical_var, expected_var = run_simulation(simulation_iter, stand_eqtl_genotype, stand_gwas_genotype, eqtl_ld, gwas_ld)
	emperical_vars.append(emperical_var)
	expected_vars.append(expected_var)
emperical_vars = np.asarray(emperical_vars)
expected_vars = np.asarray(expected_vars)
print(np.mean(emperical_vars))
print(np.mean(expected_vars))
'''


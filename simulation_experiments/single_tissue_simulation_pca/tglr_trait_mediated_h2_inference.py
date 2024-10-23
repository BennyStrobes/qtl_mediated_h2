import numpy as np
import os
import sys
import pdb
import statsmodels.api as sm
from bgen import BgenReader




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

def extract_gwas_variant_ldscores(gwas_rsids, quasi_ld_window_summary_file, LD_mat):
	var_ld_score = []
	f = open(quasi_ld_window_summary_file)
	window_to_ld_files = {}
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
		ld_file = data[3]
		ld = np.load(ld_file)

		window_scores = np.diag(np.dot(ld, np.transpose(ld)))

		var_ld_score.append(window_scores)


		if window_name in window_to_ld_files:
			print('assumption eororor')
			pdb.set_trace()
		window_to_ld_files[window_name] = (ld_file, window_snp_indices)

	f.close()

	var_ld_score_old = np.hstack(var_ld_score)


	var_ld_score = np.sum(np.square(LD_mat),axis=0)

	return var_ld_score, window_to_ld_files

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


def load_in_eqtl_data2(eqtl_sumstat_file, eqtl_gene_model_file, window_name_to_ld_mat_files, eqtl_ld_scores, gwas_rsids, LD_mat):
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


	gene_name_to_window = {}
	# Now loop through genes
	for gene_name in gene_names:
		# Extract summary stats for specific gene
		eqtl_gene_beta, gene_cis_snp_indices, gene_window_name = extract_eqtl_sumstats_for_specific_gene(eqtl_sumstat_file, gene_name)
		gene_name_to_window[gene_name] = (gene_window_name, gene_cis_snp_indices)

	n_genes = 0
	f = open(eqtl_gene_model_file)
	gene_vars = []
	final_genes = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_booler = data[1]
		if gene_booler != 'True':
			continue
		n_genes = n_genes + 1

		gene_name = data[0]
		pmces = np.load(data[2])
		cis_snp_names = np.load(data[3])

		window_name, gene_cis_snp_indices_raw = gene_name_to_window[gene_name]

		gene_cis_snp_indices = gene_cis_snp_indices_raw ==1.0

		window_ld_file, window_indices = window_name_to_ld_mat_files[window_name]
		window_ld = np.load(window_ld_file)

		new_cis_snps = gwas_rsids[window_indices][gene_cis_snp_indices]

		# Quick error check
		if np.array_equal(new_cis_snps, cis_snp_names) == False:
			print('assumption oeroror')
			pdb.set_trace()

		gene_ld = window_ld[:, gene_cis_snp_indices][gene_cis_snp_indices,:]

		gene_var_old = np.dot(np.dot(pmces, gene_ld), pmces)

		gene_ld2 = (LD_mat[:, window_indices][:, gene_cis_snp_indices][window_indices,:][gene_cis_snp_indices,:])
		gene_var = np.dot(np.dot(pmces,gene_ld2), pmces)

		gene_vars.append(gene_var)
		final_genes.append(gene_name)


		squared_corrs = np.square(np.dot((LD_mat[:, window_indices][:, gene_cis_snp_indices]), pmces/np.sqrt(gene_var)))

		eqtl_ld_scores = eqtl_ld_scores + squared_corrs

	f.close()

	return eqtl_ld_scores, np.asarray(final_genes), np.asarray(gene_vars)

def load_in_eqtl_data(eqtl_sumstat_file, eqtl_gene_model_file, window_name_to_ld_mat_files, eqtl_ld_scores, gwas_rsids, LD_mat):
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


	gene_name_to_window = {}
	# Now loop through genes
	for gene_name in gene_names:
		# Extract summary stats for specific gene
		eqtl_gene_beta, gene_cis_snp_indices, gene_window_name = extract_eqtl_sumstats_for_specific_gene(eqtl_sumstat_file, gene_name)
		gene_name_to_window[gene_name] = (gene_window_name, gene_cis_snp_indices)

	n_genes = 0
	f = open(eqtl_gene_model_file)
	gene_vars = []
	final_genes = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_booler = data[1]
		if gene_booler != 'True':
			continue
		n_genes = n_genes + 1

		gene_name = data[0]
		pmces = np.load(data[2])
		cis_snp_names = np.load(data[3])

		window_name, gene_cis_snp_indices_raw = gene_name_to_window[gene_name]

		gene_cis_snp_indices = gene_cis_snp_indices_raw ==1.0

		window_ld_file, window_indices = window_name_to_ld_mat_files[window_name]
		window_ld = np.load(window_ld_file)

		new_cis_snps = gwas_rsids[window_indices][gene_cis_snp_indices]

		# Quick error check
		if np.array_equal(new_cis_snps, cis_snp_names) == False:
			print('assumption oeroror')
			pdb.set_trace()

		gene_ld = window_ld[:, gene_cis_snp_indices][gene_cis_snp_indices,:]

		gene_var_old = np.dot(np.dot(pmces, gene_ld), pmces)

		gene_ld2 = (LD_mat[:, window_indices][:, gene_cis_snp_indices][window_indices,:][gene_cis_snp_indices,:])
		gene_var = np.dot(np.dot(pmces,gene_ld2), pmces)

		gene_vars.append(gene_var)
		final_genes.append(gene_name)

		gene_ld2 = window_ld[:, gene_cis_snp_indices]

		squared_corrs = np.square(np.dot((LD_mat[:, window_indices][:, gene_cis_snp_indices]), pmces))

		eqtl_ld_scores = eqtl_ld_scores + squared_corrs

	f.close()

	return eqtl_ld_scores, np.asarray(final_genes), np.asarray(gene_vars)

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


###############################
# Command line args
###############################
simulation_number = sys.argv[1]
simulation_name_string = sys.argv[2]
simulated_trait_dir = sys.argv[3]
simulated_gwas_dir = sys.argv[4]
simulation_genotype_dir = sys.argv[5]
simulated_learned_gene_models_dir = sys.argv[6]
n_gwas_individuals = sys.argv[7]
eqtl_sample_size = sys.argv[8]
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

#print(sim_h2)
#print(sim_med_h2)
#print(sim_nm_h2)

# Load in Genotype data
# Load in genotype object
genotype_stem = simulation_genotype_dir + 'simulated_gwas_data_' + '1'
genotype_obj = BgenReader(genotype_stem + '.bgen')
G_obj_pos = np.asarray(genotype_obj.positions())
G_obj_rsids = np.asarray(genotype_obj.rsids())
# Load in ref-alt alleles
ref_alt_alleles = load_in_ref_alt_allele_arr(genotype_stem + '.pvar')
genotype_dosage = load_in_alt_allele_genotype_dosage_mat(genotype_obj, np.arange(len(G_obj_rsids)), ref_alt_alleles)
G_obj_geno_stand = standardize_genotype_dosage_matrix(genotype_dosage)
LD_mat = np.corrcoef(np.transpose(G_obj_geno_stand))

# Load in GWAS summary statistics
gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results.txt'
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(gwas_summary_file)
N_gwas = n_gwas_individuals
gwas_chi_sq = np.square(gwas_beta/gwas_beta_se)

# Convert gwas_beta to gwas_beta_pc and get var_pc_ld_scores
quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_summary.txt'
var_ld_scores, window_name_to_ld_mat_files = extract_gwas_variant_ldscores(gwas_rsids, quasi_ld_window_summary_file, LD_mat)

# Compute total h2 through ldsc
model = sm.OLS(gwas_chi_sq -1, (var_ld_scores))
results = model.fit()
ldsc_h2_est = results.params[0]*len(gwas_rsids)/float(n_gwas_individuals)





# load in eqtl data
# Out of sample eqtl ld
eqtl_ld_scores = np.copy(var_ld_scores)*0.0
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
eqtl_gene_model_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_multivariate_gene_model_output.txt'
eqtl_ld_scores, gene_names, gene_vars = load_in_eqtl_data(eqtl_sumstat_file, eqtl_gene_model_file, window_name_to_ld_mat_files,eqtl_ld_scores, gwas_rsids, LD_mat)

joint_ld_scores = np.transpose(np.vstack((var_ld_scores, eqtl_ld_scores)))
model = sm.OLS(gwas_chi_sq -1, (joint_ld_scores))
results = model.fit()

ldsc_nm_h2_est = results.params[0]*len(gwas_rsids)/float(n_gwas_individuals)
ldsc_med_h2_est = results.params[1]*len(gene_names)*np.mean(gene_vars)/float(n_gwas_individuals)



eqtl_ld_scores2, gene_names2, gene_vars2 = load_in_eqtl_data2(eqtl_sumstat_file, eqtl_gene_model_file, window_name_to_ld_mat_files,eqtl_ld_scores, gwas_rsids, LD_mat)

joint_ld_scores2 = np.transpose(np.vstack((var_ld_scores, eqtl_ld_scores2)))
model2 = sm.OLS(gwas_chi_sq -1, (joint_ld_scores2))
results2 = model2.fit()

ldsc_nm_h2_est2 = results2.params[0]*len(gwas_rsids)/float(n_gwas_individuals)
ldsc_med_h2_est2 = results2.params[1]*len(gene_names2)/float(n_gwas_individuals)








output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_tglr_estimates.txt'
t = open(output_file,'w')
t.write('method\tsim_tot_h2\tsim_nm_h2\tsim_med_h2\test_tot_h2\test_nm_h2\test_med_h2\n')

t.write('tglr_unscaled\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(ldsc_h2_est) + '\t' + str(ldsc_nm_h2_est) + '\t' + str(ldsc_med_h2_est) + '\n')
t.write('tglr_scaled\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(ldsc_h2_est) + '\t' + str(ldsc_nm_h2_est2) + '\t' + str(ldsc_med_h2_est2) + '\n')

t.close()

print(output_file)





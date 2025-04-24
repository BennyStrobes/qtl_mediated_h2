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


def get_gwas_variant_ld_scores(gwas_beta, gwas_rsids, quasi_ld_window_summary_file):
	var_ld_score = []
	f = open(quasi_ld_window_summary_file)
	window_to_ld_files = {}
	head_count = 0
	tmp_rsids = []
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
		ld_mat = np.load(ld_file)
		tmp_rsids.append(window_rsids)

		squared_ld_mat = np.square(ld_mat)
		squared_adj_ld_mat = squared_ld_mat - ((1.0-squared_ld_mat)/(100000-2.0))

		window_ld_scores = np.sum(squared_adj_ld_mat,axis=0)

		var_ld_score.append(window_ld_scores)

		if window_name in window_to_ld_files:
			print('assumption eororor')
			pdb.set_trace()
		window_to_ld_files[window_name] = (ld_file)


	if np.array_equal(np.hstack(tmp_rsids), gwas_rsids) == False:
		print('assumption erorror')
		pdb.set_trace()



	return np.hstack(var_ld_score), window_to_ld_files

def extract_eqtl_sumstats_for_specific_gene(sumstats_file, gene_name):
	f = open(sumstats_file)
	sumstats = []
	sumstat_ses = []
	cis_snps = []
	window_names = []
	rsids = []
	class_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != gene_name:
			continue
		sumstats.append(data[5])
		cis_snps.append(data[4])
		sumstat_ses.append(data[6])
		window_name = data[3]
		window_names.append(window_name)
		rsids.append(data[2])
		class_names.append(data[1])
	f.close()

	unique_window_names = np.unique(window_names)
	if len(unique_window_names) != 1:
		print('assumption eroror')
		pdb.set_trace()
	class_names = np.unique(class_names)
	if len(class_names) != 1:
		print('assumption eroror')
		pdb.set_trace()

	gene_window_name = unique_window_names[0]

	return np.asarray(sumstats).astype(float), np.asarray(sumstat_ses).astype(float), np.asarray(cis_snps).astype(float), gene_window_name, np.asarray(rsids), class_names[0]



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

def create_mapping_from_rsid_to_position(rsids):
	rsid_to_position = {}

	for ii, rsid in enumerate(rsids):
		if rsid in rsid_to_position:
			print('assumption eroror')
			pdb.set_trace()
		rsid_to_position[rsid] = ii
	return rsid_to_position


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


def create_mapping_form_gene_name_to_multivariate_gene_model_info(multivariate_gene_model_summary_file):
	mapping = {}
	f = open(multivariate_gene_model_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[0]
		if data[1] == 'False':
			continue
		ldsc_h2 = float(data[2])
		eqtl_pred_h2 = float(data[5])
		gene_model_file = data[3]
		cis_snp_names_file = data[4]
		if gene_name in mapping:
			print('assumption erororo')
			pdb.set_trace()

		mapping[gene_name] = (gene_model_file, cis_snp_names_file, ldsc_h2, eqtl_pred_h2)

	f.close()

	return mapping


def load_in_eqtl_data3(eqtl_sumstat_file, snp_name_to_position, window_name_to_ld_files, eqtl_sample_size, num_snps, alt_simulated_learned_gene_models_dir, simulation_name_string, gene_model_type='susie', eqtl_ld='out_of_sample', ldsc_variance=False, single_tissue=False):
	# First get list of gene names
	gene_names = []
	gene_classes = []
	f = open(eqtl_sumstat_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if single_tissue and data[1] != 'tissue0':
			continue
		gene_names.append(data[0])
		gene_classes.append(data[0].split(':')[1])
	f.close()
	gene_names = np.unique(gene_names)
	gene_classes = np.asarray(gene_classes)
	unique_gene_classes = np.sort(np.unique(gene_classes))
	class_name_to_class_index = {}
	index_to_gene_h2 = {}
	for ii, ele in enumerate(unique_gene_classes):
		class_name_to_class_index[ele] = ii
		index_to_gene_h2[ii] = []


	# Initialize matrix of eQTL ld scores
	eqtl_ld_scores = np.zeros((num_snps, len(unique_gene_classes)))
	n_genes = np.zeros(len(unique_gene_classes))


	# Create mapping from gene name to multivariate gene model info
	multivariate_gene_model_summary_file = alt_simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_multivariate_gene_model_output.txt'
	gene_name_to_multivariate_gene_model = create_mapping_form_gene_name_to_multivariate_gene_model_info(multivariate_gene_model_summary_file)


	# Now loop through genes
	for gg_iter, gene_name in enumerate(gene_names):
		# Extract summary stats for specific gene
		eqtl_gene_beta, eqtl_gene_beta_se, gene_cis_snp_indices, gene_window_name, gene_rsids, gene_class_name = extract_eqtl_sumstats_for_specific_gene(eqtl_sumstat_file, gene_name)
		gene_class_index = class_name_to_class_index[gene_class_name]

		# Load in q_mat and w_mat for this window
		ld_file = window_name_to_ld_files[gene_window_name]
		if eqtl_ld == 'out_of_sample':
			ld_mat = np.load(ld_file)
		elif eqtl_ld == 'in_sample_adjusted':
			new_ld_file = ld_file.split('ref_geno_gwas_')[0] + 'ref_geno_eqtl_' + str(eqtl_sample_size) + '_' + ld_file.split('ref_geno_gwas_')[1]
			ld_mat = np.load(new_ld_file)

		# SKip genes with no gene model
		if gene_model_type == 'susie' and gene_name not in gene_name_to_multivariate_gene_model:
			continue

			

		# Get names of snps in gene
		gene_snp_positions = []
		for gene_rsid in gene_rsids:
			gene_snp_positions.append(snp_name_to_position[gene_rsid])
		gene_snp_positions = np.asarray(gene_snp_positions)


		if gene_model_type == 'susie':
			#extract gene model info
			gene_model_file, gene_snp_name_file, gene_ldsc_h2, gene_pred_h2 = gene_name_to_multivariate_gene_model[gene_name]

			# load in gene model
			gene_model = np.load(gene_model_file)

			# get variance of genetically predicted gene
			small_ld_mat = ld_mat[gene_cis_snp_indices==1,:][:, gene_cis_snp_indices==1]
			gene_var = np.dot(np.dot(gene_model, small_ld_mat), gene_model)

			# Computed squared correlations
			snp_gene_squared_correlations = np.square(np.dot(ld_mat[:, gene_cis_snp_indices==1], gene_model))
		elif gene_model_type == 'summary_stats':
			snp_gene_squared_correlations = np.square((eqtl_gene_beta/eqtl_gene_beta_se)/np.sqrt(float(eqtl_sample_size))) - (1.0/float(eqtl_sample_size))
			gene_var = 0.05

		# Update eqtl ld scores matrix
		eqtl_ld_scores[gene_snp_positions, gene_class_index] = eqtl_ld_scores[gene_snp_positions, gene_class_index] + snp_gene_squared_correlations

		# Update n_genes matrix
		n_genes[gene_class_index] = n_genes[gene_class_index] + 1


		index_to_gene_h2[gene_class_index].append(gene_var)



	gene_cis_h2 = np.zeros(len(unique_gene_classes))
	for ii in range(len(unique_gene_classes)):
		gene_cis_h2[ii] = np.mean(index_to_gene_h2[ii])




	return eqtl_ld_scores, n_genes, gene_cis_h2





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
alt_simulated_learned_gene_models_dir = sys.argv[10]



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
# Load in GWAS summary statistics
gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results.txt'
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(gwas_summary_file)

# Get variant ld scores
quasi_ld_window_summary_file = simulation_genotype_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld_summary.txt'
gwas_variant_ld_scores, window_to_ld_files = get_gwas_variant_ld_scores(gwas_beta, gwas_rsids, quasi_ld_window_summary_file)
# Other relevent quantities
N_gwas = n_gwas_individuals
gwas_chi_sq = np.square(gwas_beta/gwas_beta_se)


# Create mapping from rsid to position
snp_name_to_position = create_mapping_from_rsid_to_position(gwas_rsids)

# Compute total h2 through ldsc
model = sm.OLS(gwas_chi_sq -1, (gwas_variant_ld_scores))
results = model.fit()
ldsc_h2_est = results.params[0]*len(gwas_rsids)/float(n_gwas_individuals)


###################################
# 5 tissues (Unpermuted eqtls)
###################################
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_tglr_estimates.txt'
t = open(output_file,'w')
t.write('method\tsim_tot_h2\tsim_nm_h2\tsim_med_h2\test_tot_h2\test_nm_h2\test_total_med_h2\test_per_tissue_med_h2\n')

# Write LDSC to output
t.write('ldsc\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(ldsc_h2_est) + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')


#####################
# load in eqtl data
# Susie PMCES
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
eqtl_ld_scores, n_genes, gene_cis_h2 = load_in_eqtl_data3(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, len(gwas_rsids), alt_simulated_learned_gene_models_dir, simulation_name_string, gene_model_type='susie', eqtl_ld='out_of_sample', ldsc_variance=False)
# Run TGLR
joint_ld_scores = np.hstack((np.transpose(gwas_variant_ld_scores.reshape(1,-1)), eqtl_ld_scores))
model = sm.OLS(gwas_chi_sq -1, joint_ld_scores)
results = model.fit()
# Extract estimated parameters
nm_h2 = results.params[0]*len(gwas_rsids)/float(n_gwas_individuals)
per_tissue_h2 = results.params[1:]*n_genes*gene_cis_h2/float(n_gwas_individuals)
t.write('tglr\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(nm_h2 + np.sum(per_tissue_h2)) + '\t' + str(nm_h2) + '\t' + str(np.sum(per_tissue_h2)) + '\t' + ','.join(per_tissue_h2.astype(str)) + '\n')


#####################
# load in eqtl data
# sumstats
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
eqtl_ld_scores, n_genes, gene_cis_h2 = load_in_eqtl_data3(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, len(gwas_rsids), alt_simulated_learned_gene_models_dir, simulation_name_string, gene_model_type='summary_stats', eqtl_ld='out_of_sample', ldsc_variance=False)
# Run TGLR
joint_ld_scores = np.hstack((np.transpose(gwas_variant_ld_scores.reshape(1,-1)), eqtl_ld_scores))
model = sm.OLS(gwas_chi_sq -1, joint_ld_scores)
results = model.fit()
# Extract estimated parameters
nm_h2 = results.params[0]*len(gwas_rsids)/float(n_gwas_individuals)
per_tissue_h2 = results.params[1:]*n_genes*gene_cis_h2/float(n_gwas_individuals)
t.write('tglr_sumstats\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(nm_h2 + np.sum(per_tissue_h2)) + '\t' + str(nm_h2) + '\t' + str(np.sum(per_tissue_h2)) + '\t' + ','.join(per_tissue_h2.astype(str)) + '\n')



t.close()
print(output_file)


###################################
# 5 tissues (Permute eqtls) 
###################################
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_tglr_estimates_permuted_eqtls.txt'
t = open(output_file,'w')
t.write('method\tsim_tot_h2\tsim_nm_h2\tsim_med_h2\test_tot_h2\test_nm_h2\test_total_med_h2\test_per_tissue_med_h2\n')

# Write LDSC to output
t.write('ldsc\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(ldsc_h2_est) + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')



#####################
# load in eqtl data
# Susie PMCES
perm_simulation_name_string = 'simulation_' + str(int(simulation_number)+1) + '_chrom' + simulation_name_string.split('_chrom')[1]
eqtl_sumstat_file = simulated_learned_gene_models_dir + perm_simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
eqtl_ld_scores, n_genes, gene_cis_h2 = load_in_eqtl_data3(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, len(gwas_rsids), alt_simulated_learned_gene_models_dir, perm_simulation_name_string,gene_model_type='susie', eqtl_ld='out_of_sample', ldsc_variance=False)
# Run TGLR
joint_ld_scores = np.hstack((np.transpose(gwas_variant_ld_scores.reshape(1,-1)), eqtl_ld_scores))
model = sm.OLS(gwas_chi_sq -1, joint_ld_scores)
results = model.fit()
# Extract estimated parameters
nm_h2 = results.params[0]*len(gwas_rsids)/float(n_gwas_individuals)
per_tissue_h2 = results.params[1:]*n_genes*gene_cis_h2/float(n_gwas_individuals)
t.write('tglr\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(nm_h2 + np.sum(per_tissue_h2)) + '\t' + str(nm_h2) + '\t' + str(np.sum(per_tissue_h2)) + '\t' + ','.join(per_tissue_h2.astype(str)) + '\n')


#####################
# load in eqtl data
# Summary stats
perm_simulation_name_string = 'simulation_' + str(int(simulation_number)+1) + '_chrom' + simulation_name_string.split('_chrom')[1]
eqtl_sumstat_file = simulated_learned_gene_models_dir + perm_simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
eqtl_ld_scores, n_genes, gene_cis_h2 = load_in_eqtl_data3(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, len(gwas_rsids), alt_simulated_learned_gene_models_dir, perm_simulation_name_string,gene_model_type='summary_stats', eqtl_ld='out_of_sample', ldsc_variance=False)
# Run TGLR
joint_ld_scores = np.hstack((np.transpose(gwas_variant_ld_scores.reshape(1,-1)), eqtl_ld_scores))
model = sm.OLS(gwas_chi_sq -1, joint_ld_scores)
results = model.fit()
# Extract estimated parameters
nm_h2 = results.params[0]*len(gwas_rsids)/float(n_gwas_individuals)
per_tissue_h2 = results.params[1:]*n_genes*gene_cis_h2/float(n_gwas_individuals)
t.write('tglr_sumstats\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(nm_h2 + np.sum(per_tissue_h2)) + '\t' + str(nm_h2) + '\t' + str(np.sum(per_tissue_h2)) + '\t' + ','.join(per_tissue_h2.astype(str)) + '\n')




t.close()
print(output_file)




'''
###################################
# 1 tissues (Unpermuted eqtls)
###################################
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_tglr_estimates_single_tissue.txt'
t = open(output_file,'w')
t.write('method\tsim_tot_h2\tsim_nm_h2\tsim_med_h2\test_tot_h2\test_nm_h2\test_total_med_h2\test_per_tissue_med_h2\n')

# Write LDSC to output
t.write('ldsc\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(ldsc_h2_est) + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')

# load in eqtl data
# Out of sample eqtl ld
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
eqtl_ld_scores, n_genes, gene_cis_h2 = load_in_eqtl_data3(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, len(gwas_rsids), alt_simulated_learned_gene_models_dir, simulation_name_string, eqtl_ld='out_of_sample', ldsc_variance=False, single_tissue=True)

# Run TGLR
joint_ld_scores = np.hstack((np.transpose(gwas_variant_ld_scores.reshape(1,-1)), eqtl_ld_scores))
model = sm.OLS(gwas_chi_sq -1, joint_ld_scores)
results = model.fit()
# Extract estimated parameters
nm_h2 = results.params[0]*len(gwas_rsids)/float(n_gwas_individuals)
per_tissue_h2 = results.params[1:]*n_genes*gene_cis_h2/float(n_gwas_individuals)
t.write('tglr\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(nm_h2 + np.sum(per_tissue_h2)) + '\t' + str(nm_h2) + '\t' + str(np.sum(per_tissue_h2)) + '\t' + ','.join(per_tissue_h2.astype(str)) + '\n')

# load in eqtl data
# Out of sample eqtl ld
eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
eqtl_ld_scores, n_genes, gene_cis_h2 = load_in_eqtl_data3(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, len(gwas_rsids), alt_simulated_learned_gene_models_dir, simulation_name_string, eqtl_ld='out_of_sample', ldsc_variance=True,single_tissue=True)

# Run TGLR
joint_ld_scores = np.hstack((np.transpose(gwas_variant_ld_scores.reshape(1,-1)), eqtl_ld_scores))
model = sm.OLS(gwas_chi_sq -1, joint_ld_scores)
results = model.fit()
# Extract estimated parameters
nm_h2 = results.params[0]*len(gwas_rsids)/float(n_gwas_individuals)
per_tissue_h2 = results.params[1:]*n_genes*gene_cis_h2/float(n_gwas_individuals)
t.write('tglr_ldsc_h2_reweight\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(nm_h2 + np.sum(per_tissue_h2)) + '\t' + str(nm_h2)+ '\t' + str(np.sum(per_tissue_h2)) + '\t' + ','.join(per_tissue_h2.astype(str)) + '\n')

t.close()
print(output_file)
'''

###################################
# 1 tissues (Permuted eqtls)
###################################
'''
output_file = trait_med_h2_inference_dir + simulation_name_string + '_eqtl_SS_' + str(eqtl_sample_size) + '_tglr_estimates_single_tissue_permuted_eqtls.txt'
t = open(output_file,'w')
t.write('method\tsim_tot_h2\tsim_nm_h2\tsim_med_h2\test_tot_h2\test_nm_h2\test_total_med_h2\test_per_tissue_med_h2\n')

# Write LDSC to output
t.write('ldsc\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(ldsc_h2_est) + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')

# load in eqtl data
# Out of sample eqtl ld
perm_simulation_name_string = 'simulation_' + str(int(simulation_number)+1) + '_chrom' + simulation_name_string.split('_chrom')[1]
eqtl_sumstat_file = simulated_learned_gene_models_dir + perm_simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
eqtl_ld_scores, n_genes, gene_cis_h2 = load_in_eqtl_data3(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, len(gwas_rsids), alt_simulated_learned_gene_models_dir, perm_simulation_name_string, eqtl_ld='out_of_sample', ldsc_variance=True, single_tissue=True)


# Run TGLR
joint_ld_scores = np.hstack((np.transpose(gwas_variant_ld_scores.reshape(1,-1)), eqtl_ld_scores))
model = sm.OLS(gwas_chi_sq -1, joint_ld_scores)
results = model.fit()
# Extract estimated parameters
nm_h2 = results.params[0]*len(gwas_rsids)/float(n_gwas_individuals)
per_tissue_h2 = results.params[1:]*n_genes*gene_cis_h2/float(n_gwas_individuals)
t.write('tglr\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(nm_h2 + np.sum(per_tissue_h2)) + '\t' + str(nm_h2) + '\t' + str(np.sum(per_tissue_h2)) + '\t' + ','.join(per_tissue_h2.astype(str)) + '\n')

# load in eqtl data
# Out of sample eqtl ld
perm_simulation_name_string = 'simulation_' + str(int(simulation_number)+1) + '_chrom' + simulation_name_string.split('_chrom')[1]
eqtl_sumstat_file = simulated_learned_gene_models_dir + perm_simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
eqtl_ld_scores, n_genes, gene_cis_h2 = load_in_eqtl_data3(eqtl_sumstat_file, snp_name_to_position, window_to_ld_files ,eqtl_sample_size, len(gwas_rsids), alt_simulated_learned_gene_models_dir, perm_simulation_name_string, eqtl_ld='out_of_sample', ldsc_variance=True, single_tissue=True)

# Run TGLR
joint_ld_scores = np.hstack((np.transpose(gwas_variant_ld_scores.reshape(1,-1)), eqtl_ld_scores))
model = sm.OLS(gwas_chi_sq -1, joint_ld_scores)
results = model.fit()
# Extract estimated parameters
nm_h2 = results.params[0]*len(gwas_rsids)/float(n_gwas_individuals)
per_tissue_h2 = results.params[1:]*n_genes*gene_cis_h2/float(n_gwas_individuals)
t.write('tglr_ldsc_h2_reweight\t' + str(sim_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t' + str(nm_h2 + np.sum(per_tissue_h2)) + '\t' + str(nm_h2)+ '\t' + str(np.sum(per_tissue_h2)) + '\t' + ','.join(per_tissue_h2.astype(str)) + '\n')

t.close()
print(output_file)
'''




import numpy as np
import os
import sys
import pdb
import time
import pickle
import gzip
#import joint_ldsc
#import joint_ldsc_gibbs
import argparse
import gibbs_ash
import statsmodels.api as sm
from linearmodels import IV2SLS, IVGMM, IVGMMCUE, IVLIML

def linear_regression(XX, YY, intercept=False):
	if intercept:
		print('not yet implemented. add column of 1s to X')
		pdb.set_trace()

	return np.dot(np.dot(np.linalg.pinv(np.dot(np.transpose(XX), XX)), np.transpose(XX)), YY)

def load_in_gwas_data(gwas_summary_file, chromosome_dictionary):
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
		# Filter out variants not in chromosome list
		line_chrom_num = data[1]
		if line_chrom_num not in chromosome_dictionary:
			print('filtered snp')
			continue

		rsids.append(data[0])
		betas.append(float(data[5]))
		beta_ses.append(float(data[6]))

	f.close()

	return np.asarray(rsids), np.asarray(betas), np.asarray(beta_ses)

def load_in_regression_snps_for_gene(regression_snp_file):
	f = open(regression_snp_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[1])
	f.close()
	return np.asarray(arr)

def load_in_eqtl_data(rsid_to_position, gene_ldscore_filestem, gene_ldscore_filesuffix, eqtl_dataset_names, chrom_arr):
	genes = []
	gene_info = {}

	for chrom_num in chrom_arr:
		gene_summary_file = gene_ldscore_filestem + str(chrom_num) + gene_ldscore_filesuffix
		f = open(gene_summary_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue

			# Each line is a gene. Extract relevent fields
			ensamble_id = data[1]
			genes.append(ensamble_id)
			squared_ld_mat_file = data[7]
			n_low_dim_snp_file = data[8]
			regression_snp_file = data[6]
			if ensamble_id in gene_info:
				print('assumption eroror')
				pdb.set_trace()
			# Load in regression snps
			gene_regression_snps = load_in_regression_snps_for_gene(regression_snp_file)
			# Get corresponding indices of regression snps
			regression_snp_indices = []
			for rsid in gene_regression_snps:
				regression_snp_indices.append(rsid_to_position[rsid])
			regression_snp_indices = np.asarray(regression_snp_indices)
			# Create mapping from rsid to gene position
			rsid_to_gene_position = {}
			for ii, val in enumerate(gene_regression_snps):
				rsid_to_gene_position[val] = ii


			gene_info[ensamble_id] = {}
			gene_info[ensamble_id]['squared_ld_file'] = squared_ld_mat_file
			gene_info[ensamble_id]['n_low_dimensional_snps'] = np.load(n_low_dim_snp_file)
			gene_info[ensamble_id]['regression_snp_indices'] = regression_snp_indices
			gene_info[ensamble_id]['rsid_to_gene_position'] = rsid_to_gene_position
			gene_info[ensamble_id]['n_gene_regression_snps'] = len(gene_regression_snps)
			gene_info[ensamble_id]['squared_sumstats'] = np.zeros((len(gene_regression_snps), len(eqtl_dataset_names)))
			gene_info[ensamble_id]['beta_sumstats'] = np.zeros((len(gene_regression_snps), len(eqtl_dataset_names)))
			gene_info[ensamble_id]['beta_se_sumstats'] = np.zeros((len(gene_regression_snps), len(eqtl_dataset_names)))
			gene_info[ensamble_id]['eqtl_category_names'] = np.copy(eqtl_dataset_names)
		f.close()


	return np.asarray(genes), gene_info


def get_gwas_variant_ld_scores(variant_ldscore_filestem, variant_ldscore_filesuffix, chrom_arr, non_mediated_annotation_version):

	rsids = []
	ldscores = []

	for chrom_num in chrom_arr:
		filer = variant_ldscore_filestem + chrom_num + variant_ldscore_filesuffix
		f = open(filer)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				anno_names = np.asarray(data[6:])
				continue
			rsid = data[1]
			ldscore = np.asarray(data[6:]).astype(float)
			rsids.append(rsid)
			ldscores.append(ldscore)
		f.close()
	ldscores = np.asarray(ldscores)

	if non_mediated_annotation_version == 'genotype_intercept':
		ldscores = ldscores[:,:1]
		anno_names = anno_names[:1]

	return ldscores, np.asarray(rsids), anno_names


def create_mapping_from_rsid_to_position(rsids):
	rsid_to_position = {}

	for ii, rsid in enumerate(rsids):
		if rsid in rsid_to_position:
			print('assumption eroror')
			pdb.set_trace()
		rsid_to_position[rsid] = ii
	return rsid_to_position

def fill_in_eqtl_sumstats(gene_info, eqtl_dataset_files, eqtl_dataset_names, genes, chrom_dicti, rsid_to_variant_stdev):
	for ii, eqtl_dataset_name in enumerate(eqtl_dataset_names):
		eqtl_sumstat_file = eqtl_dataset_files[ii]
		f = open(eqtl_sumstat_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			line_chrom_num = data[2]
			if line_chrom_num not in chrom_dicti:
				print('skipped variant-gene because wrong chromosome')
				pdb.set_trace()
			ens_id = data[0]
			rsid = data[1]
			beta = float(data[6])
			beta_se = float(data[7])
			stdev = rsid_to_variant_stdev[rsid]
			# Update betas and standard errors according to this
			beta = beta*stdev
			beta_se = beta_se*stdev
			gene_row_index = gene_info[ens_id]['rsid_to_gene_position'][rsid]
			gene_info[ens_id]['squared_sumstats'][gene_row_index, ii] = np.square(beta) - np.square(beta_se)
			gene_info[ens_id]['beta_sumstats'][gene_row_index, ii] = beta
			gene_info[ens_id]['beta_se_sumstats'][gene_row_index, ii] = beta_se
		f.close()

	# Filter out missing eqtl categories for each gene
	for ens_id in genes:
		valid_categories = []
		for ii, category_name in enumerate(gene_info[ens_id]['eqtl_category_names']):
			gene_info[ens_id]['squared_sumstats'][:, ii]
			if np.array_equal(gene_info[ens_id]['squared_sumstats'][:, ii], np.zeros(len(gene_info[ens_id]['squared_sumstats'][:, ii]))):
				valid_categories.append(False)
			else:
				valid_categories.append(True)
		valid_categories = np.asarray(valid_categories)
		gene_info[ens_id]['squared_sumstats'] = gene_info[ens_id]['squared_sumstats'][:,valid_categories]
		gene_info[ens_id]['beta_sumstats'] = gene_info[ens_id]['beta_sumstats'][:,valid_categories]
		gene_info[ens_id]['beta_se_sumstats'] = gene_info[ens_id]['beta_se_sumstats'][:,valid_categories]
		gene_info[ens_id]['eqtl_category_names'] = gene_info[ens_id]['eqtl_category_names'][valid_categories]
	return gene_info

def extract_number_of_reference_snps(variant_M_filestem, chrom_arr, non_mediated_annotation_version):
	n_snps = []
	for chrom_num in chrom_arr:
		filer = variant_M_filestem + str(chrom_num) + '_M.txt'
		f = open(filer)
		head_count = 0
		tmp_arr = []
		for line in f:
			line = line.rstrip()
			if head_count == 0:
				head_count = head_count + 1
				continue
			data = line.split('\t')
			tmp_arr.append(float(data[1]))
		f.close()
		n_snps.append(np.asarray(tmp_arr))
	n_snps = np.asarray(n_snps)
	final_n_snps = np.sum(n_snps,axis=0)

	# Filter columns if using genotype intercept
	if non_mediated_annotation_version == 'genotype_intercept':
		final_n_snps = final_n_snps[:1]

	return final_n_snps

def extract_chromosome_names(chromosome_file):
	if chromosome_file is None:
		tmp_chrom_arr = np.arange(1,23)
		arr = []
		dicti = {}
		for ele in tmp_chrom_arr:
			arr.append(str(ele))
			dicti[str(ele)] = 1
	else:
		arr = []
		dicti = {}
		f = open(chromosome_file)
		for line in f:
			line = line.rstrip()
			arr.append(line)
			dicti[line] = 1
		f.close()
	return np.asarray(arr), dicti

def load_in_validation_eqtl_dataset_summary_file(eqtl_summary_file):
	names = []
	sample_sizes = []
	filers = []
	f = open(eqtl_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		names.append(data[0])
		sample_sizes.append(float(data[1]))
		filers.append(data[2].split('.txt')[0] + '_repeat.txt')
	f.close()
	return np.asarray(names), np.asarray(sample_sizes), np.asarray(filers)


def load_in_eqtl_dataset_summary_file(eqtl_summary_file):
	names = []
	sample_sizes1 = []
	sample_sizes2 = []
	sample_sizesfull = []
	filers1 = []
	filers2 = []
	filersfull = []
	f = open(eqtl_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		names.append(data[0])
		sample_sizes1.append(float(data[1]))
		sample_sizes2.append(float(data[2]))
		sample_sizesfull.append(float(data[3]))
		filers1.append(data[4])
		filers2.append(data[5])
		filersfull.append(data[6])
	f.close()
	return np.asarray(names), np.asarray(sample_sizes1), np.asarray(sample_sizes2), np.asarray(filers1), np.asarray(filers2), np.asarray(sample_sizesfull), np.asarray(filersfull)


def update_gene_info_to_include_eqtl_dataset_names(gene_info, squared_eqtl_effect_threshold):
	# Get list of genes
	genes = np.asarray([*gene_info])
	# loop through genes
	for gene in genes:
		gene_info[gene]['eqtl_dataset_names'] = np.copy(gene_info[gene]['eqtl_category_names'])

	##############################
	# Estimate gene cis h2
	##############################

	# loop through genes
	for gene in genes:
		info = gene_info[gene]
		# Load in squared ld 
		squared_ld = np.load(info['squared_ld_file'])

		# Initialize gene cis_h2 vec
		gene_cis_h2_vec = np.zeros(len(info['eqtl_category_names']))

		# Load in eqtl sumstats
		eqtl_sq_sumstats = info['squared_sumstats']

		# Load in gene eqtl categories
		gene_eqtl_categories = info['eqtl_category_names']

		# Loop through eQTL categories
		# Initialize seperately for each category
		for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):


			temp_Y = eqtl_sq_sumstats[:, eqtl_category_iter]
			valid_indices = (np.abs(temp_Y) < squared_eqtl_effect_threshold) & (np.isnan(temp_Y) == False)
			
			weights = linear_regression(squared_ld[valid_indices,:], temp_Y[valid_indices])


			gene_est_cis_h2 = np.sum(weights*info['n_low_dimensional_snps'])
			#weights = linear_regression(np.sum(squared_ld,axis=1).reshape(-1,1), eqtl_sq_sumstats[:, eqtl_category_iter])
			#gene_est_cis_h2 = weights[0]*np.sum(info['n_low_dimensional_snps'])
			
			gene_cis_h2_vec[eqtl_category_iter] = gene_est_cis_h2
		# Update dictionary
		gene_info[gene]['gene_cis_h2_vec'] = gene_cis_h2_vec



	return gene_info

def bin_assignment(gene_cis_h2, dataset_bins, bin_names):
	n_bins = len(bin_names)
	bin_index_assignment = None
	for bin_iter in range(len(dataset_bins)-1):
		if gene_cis_h2 >= dataset_bins[bin_iter] and gene_cis_h2 < dataset_bins[(bin_iter+1)]:
			bin_index_assignment = bin_iter
	if bin_index_assignment is None:
		print('assumption error: gene not assigned to bin')
		pdb.set_trace()

	return bin_names[bin_index_assignment]

def update_gene_info_to_bin_genes_by_random(gene_info, eqtl_dataset_names, n_expr_cis_h2_bins, squared_eqtl_effect_threshold):
	##############################
	# First estimate gene cis h2
	##############################
	# Get list of genes
	genes = np.asarray([*gene_info])

	# Keep track of estimated heritability
	gene_cis_h2 = {}

	# loop through genes
	for gene in genes:
		info = gene_info[gene]
		# Load in squared ld 
		squared_ld = np.load(info['squared_ld_file'])

		# Initialize gene cis_h2 vec
		gene_cis_h2_vec = np.zeros(len(info['eqtl_category_names']))

		# Load in eqtl sumstats
		eqtl_sq_sumstats = info['squared_sumstats']

		# Load in gene eqtl categories
		gene_eqtl_categories = info['eqtl_category_names']

		# Loop through eQTL categories
		# Initialize seperately for each category
		for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):


			temp_Y = eqtl_sq_sumstats[:, eqtl_category_iter]
			valid_indices = (np.abs(temp_Y) < squared_eqtl_effect_threshold) & (np.isnan(temp_Y) == False)
			
			weights = linear_regression(squared_ld[valid_indices,:], temp_Y[valid_indices])

			gene_est_cis_h2 = np.sum(weights*info['n_low_dimensional_snps'])
			#weights = linear_regression(np.sum(squared_ld,axis=1).reshape(-1,1), eqtl_sq_sumstats[:, eqtl_category_iter])
			#gene_est_cis_h2 = weights[0]*np.sum(info['n_low_dimensional_snps'])
			
			gene_cis_h2_vec[eqtl_category_iter] = gene_est_cis_h2
		# Update dictionary
		gene_info[gene]['gene_cis_h2_vec'] = gene_cis_h2_vec



	bin_names = []
	for bin_iter in range(n_expr_cis_h2_bins):
		bin_names.append('bin' + str(bin_iter))
	bin_names = np.asarray(bin_names)

	# create list of new eqtl categories
	new_eqtl_categories = []
	new_eqtl_datasets = []
	for eqtl_dataset_name in eqtl_dataset_names:
		for bin_name in bin_names:
			new_eqtl_categories.append(eqtl_dataset_name + ':' + bin_name)
			new_eqtl_datasets.append(eqtl_dataset_name)
	new_eqtl_categories = np.asarray(new_eqtl_categories)
	new_eqtl_datasets = np.asarray(new_eqtl_datasets)


	##############################
	# Update gene_info file to now contain bin_info
	##############################
	# loop through genes
	for gene in genes:
		info = gene_info[gene]

		# Load in gene eqtl categories
		gene_eqtl_categories = info['eqtl_category_names']


		# Initialize vector to keep track of update gene_eqtl_categories
		updated_gene_eqtl_categories = []
		# Loop through eQTL categories
		for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):

			# Bin assignment
			bin_name = np.random.choice(bin_names)

			# Update gnee eqtl category
			updated_gene_eqtl_categories.append(eqtl_category_name + ':' + bin_name)
		# Put into nice array
		updated_gene_eqtl_categories = np.asarray(updated_gene_eqtl_categories)

		# Update gene info file
		gene_info[gene]['eqtl_dataset_names'] = np.copy(gene_info[gene]['eqtl_category_names'])
		gene_info[gene]['eqtl_category_names'] = updated_gene_eqtl_categories

	return gene_info, new_eqtl_categories, new_eqtl_datasets



def update_gene_info_to_bin_genes_by_est_cis_h2(gene_info, eqtl_dataset_names, n_expr_cis_h2_bins, squared_eqtl_effect_threshold):
	##############################
	# First estimate gene cis h2
	##############################
	# Get list of genes
	genes = np.asarray([*gene_info])
	# Mapping from dataset name to index
	dataset_name_to_index = {}
	dataset_name_to_cis_h2s = {}
	for ii, val in enumerate(eqtl_dataset_names):
		dataset_name_to_index[val] = ii
		dataset_name_to_cis_h2s[val] = []

	# Keep track of estimated heritability
	gene_cis_h2 = {}

	# loop through genes
	for gene in genes:
		info = gene_info[gene]
		# Load in squared ld 
		squared_ld = np.load(info['squared_ld_file'])

		# Initialize gene cis_h2 vec
		gene_cis_h2_vec = np.zeros(len(info['eqtl_category_names']))

		# Load in eqtl sumstats
		eqtl_sq_sumstats = info['squared_sumstats']

		# Load in gene eqtl categories
		gene_eqtl_categories = info['eqtl_category_names']

		# Loop through eQTL categories
		# Initialize seperately for each category
		for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):


			temp_Y = eqtl_sq_sumstats[:, eqtl_category_iter]
			valid_indices = (np.abs(temp_Y) < squared_eqtl_effect_threshold) & (np.isnan(temp_Y) == False)
			
			weights = linear_regression(squared_ld[valid_indices,:], temp_Y[valid_indices])

			gene_est_cis_h2 = np.sum(weights*info['n_low_dimensional_snps'])
			#weights = linear_regression(np.sum(squared_ld,axis=1).reshape(-1,1), eqtl_sq_sumstats[:, eqtl_category_iter])
			#gene_est_cis_h2 = weights[0]*np.sum(info['n_low_dimensional_snps'])
			
			gene_cis_h2_vec[eqtl_category_iter] = gene_est_cis_h2
			dataset_name_to_cis_h2s[eqtl_category_name].append(gene_est_cis_h2)
		# Update dictionary
		gene_cis_h2[gene] = gene_cis_h2_vec

	##############################
	# Bin cis_h2s in each eqtl dataset
	##############################
	dataset_to_bins = {}
	quantiles = np.linspace(0, 1, n_expr_cis_h2_bins + 1)
	for eqtl_dataset_name in eqtl_dataset_names:
		bins = np.quantile(dataset_name_to_cis_h2s[eqtl_dataset_name], quantiles)
		bins[-1] = bins[-1] + 1  # Done to make upper bound <= instead of <
		bins[0] = bins[0] - 1  # Done to make lower bound >= instead of >
		dataset_to_bins[eqtl_dataset_name] = bins
	# Get names of bins
	bin_names = []
	for bin_iter in range(n_expr_cis_h2_bins):
		bin_names.append('bin' + str(bin_iter))
	bin_names = np.asarray(bin_names)

	# create list of new eqtl categories
	new_eqtl_categories = []
	new_gene_eqtl_data_set_names = []
	for eqtl_dataset_name in eqtl_dataset_names:
		for bin_name in bin_names:
			new_eqtl_categories.append(eqtl_dataset_name + ':' + bin_name)
			new_gene_eqtl_data_set_names.append(eqtl_dataset_name)
	new_eqtl_categories = np.asarray(new_eqtl_categories)
	new_gene_eqtl_data_set_names = np.asarray(new_gene_eqtl_data_set_names)

	##############################
	# Update gene_info file to now contain bin_info
	##############################
	# loop through genes
	for gene in genes:
		info = gene_info[gene]

		# Get gene cis_h2 vec
		gene_cis_h2_vec = gene_cis_h2[gene]

		# Load in gene eqtl categories
		gene_eqtl_categories = info['eqtl_category_names']

		# quick error check
		if len(gene_cis_h2_vec) != len(gene_eqtl_categories):
			print('assumptionerororo')
			pdb.set_trace()

		# Initialize vector to keep track of update gene_eqtl_categories
		updated_gene_eqtl_categories = []
		# Loop through eQTL categories
		for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):

			# Bin assignment
			bin_name = bin_assignment(gene_cis_h2_vec[eqtl_category_iter], dataset_to_bins[eqtl_category_name], bin_names)
			
			# Update gnee eqtl category
			updated_gene_eqtl_categories.append(eqtl_category_name + ':' + bin_name)
		# Put into nice array
		updated_gene_eqtl_categories = np.asarray(updated_gene_eqtl_categories)

		# Update gene info file
		gene_info[gene]['eqtl_dataset_names'] = np.copy(gene_info[gene]['eqtl_category_names'])
		gene_info[gene]['eqtl_category_names'] = updated_gene_eqtl_categories
		gene_info[gene]['gene_cis_h2_vec'] = gene_cis_h2_vec

	return gene_info, new_eqtl_categories, new_gene_eqtl_data_set_names


def extract_gene_ld_scores(genes, gene_info, eqtl_category_names, gene_ldscore_type, n_regression_snps, squared_eqtl_effect_threshold):
	# Initialize output data types
	gene_ldscores = np.zeros((n_regression_snps, len(eqtl_category_names)))
	n_genes = np.zeros(len(eqtl_category_names))
	cis_snp_h2_dicti = {}
	eqtl_category_name_to_index = {}
	for ii, eqtl_category_name in enumerate(eqtl_category_names):
		cis_snp_h2_dicti[eqtl_category_name] = []
		eqtl_category_name_to_index[eqtl_category_name] = ii

	# loop through genes
	for gene in genes:
		info = gene_info[gene]
		# Load in squared ld 
		squared_ld = np.load(info['squared_ld_file'])

		# Load in eqtl sumstats
		eqtl_sq_sumstats = info['squared_sumstats']

		# Load in gene eqtl categories
		gene_eqtl_categories = info['eqtl_category_names']

		# Loop through eQTL categories
		# Initialize seperately for each category
		for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):
			temp_Y = eqtl_sq_sumstats[:, eqtl_category_iter]
			valid_indices = (np.abs(temp_Y) < squared_eqtl_effect_threshold) & (np.isnan(temp_Y) == False)
			
			weights = linear_regression(squared_ld[valid_indices,:], temp_Y[valid_indices])
			gene_est_cis_h2 = np.sum(weights*info['n_low_dimensional_snps'])

			# Keep track of cis snps h2s

			cis_snp_h2_dicti[eqtl_category_name].append(gene_est_cis_h2)

			global_eqtl_category_index = eqtl_category_name_to_index[eqtl_category_name]

			if gene_ldscore_type == 'squared_marginal_sumstats' or gene_ldscore_type == 'ashr_style_pred':
				gene_ldscore = np.copy(eqtl_sq_sumstats[:, eqtl_category_iter])
			elif gene_ldscore_type == 'ldsc_style_pred':
				gene_ldscore = np.dot(squared_ld, weights)
			else:
				print('not yet implemented')
				pdb.set_trace()

			gene_ldscores[info['regression_snp_indices'], global_eqtl_category_index] = gene_ldscores[info['regression_snp_indices'], global_eqtl_category_index] + gene_ldscore

			n_genes[global_eqtl_category_index] = n_genes[global_eqtl_category_index] + 1
	
	# Put gene cis-snp h2s into array
	gene_cis_snp_h2s = []
	for eqtl_category_name in eqtl_category_names:
		gene_cis_snp_h2s.append(np.mean(cis_snp_h2_dicti[eqtl_category_name]))
	gene_cis_snp_h2s = np.asarray(gene_cis_snp_h2s)


	return gene_ldscores, gene_cis_snp_h2s, n_genes



def extract_gene_ld_scores_for_stdExpr_analysis(genes, gene_info, gene_info_validation, eqtl_category_names, gene_ldscore_type, n_regression_snps, squared_eqtl_effect_threshold):
	# Initialize output data types
	gene_ldscores_training = np.zeros((n_regression_snps, len(eqtl_category_names)))
	gene_ldscores_validation = np.zeros((n_regression_snps, len(eqtl_category_names)))
	n_genes = np.zeros(len(eqtl_category_names))
	cis_snp_h2_dicti = {}
	eqtl_category_name_to_index = {}
	for ii, eqtl_category_name in enumerate(eqtl_category_names):
		cis_snp_h2_dicti[eqtl_category_name] = []
		eqtl_category_name_to_index[eqtl_category_name] = ii

	if gene_ldscore_type == 'ashr_style_pred':
		gene_info = run_ashr_style_posterior_on_eqtl_data(gene_info, genes, eqtl_category_names)

	# loop through genes
	for gene in genes:
		info = gene_info[gene]
		info_validation = gene_info_validation[gene]
		# Load in eqtl sumstats
		if gene_ldscore_type == 'ashr_style_pred':
			eqtl_sq_sumstats = info['ashr_squared_sumstats']
		else:
			eqtl_sq_sumstats = info['squared_sumstats']
		eqtl_sq_sumstats_validation = info_validation['squared_sumstats']


		# Load in gene eqtl categories
		gene_eqtl_categories = info['eqtl_category_names']

		# Loop through eQTL categories
		# Initialize seperately for each category
		for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):


			gene_est_cis_h2 = info['gene_cis_h2_vec'][eqtl_category_iter]
			gene_est_cis_h2_validation = info_validation['gene_cis_h2_vec'][eqtl_category_iter]


			# Keep track of cis snps h2s
			cis_snp_h2_dicti[eqtl_category_name].append(gene_est_cis_h2_validation)

			global_eqtl_category_index = eqtl_category_name_to_index[eqtl_category_name]

			if gene_ldscore_type == 'squared_marginal_sumstats' or gene_ldscore_type == 'ashr_style_pred':
				gene_ldscore = np.copy(eqtl_sq_sumstats[:, eqtl_category_iter])
				gene_ldscore_validation = np.copy(eqtl_sq_sumstats_validation[:, eqtl_category_iter])
			elif gene_ldscore_type == 'ldsc_style_pred':
				print('not yet implemented')
				pdb.set_trace()
				gene_ldscore = np.dot(squared_ld, weights)
			else:
				print('not yet implemented')
				pdb.set_trace()

			if np.array_equal(info['regression_snp_indices'], info_validation['regression_snp_indices']) == False:
				print('assumption eroorroror')
				pdb.set_trace()

			gene_ldscores_training[info['regression_snp_indices'], global_eqtl_category_index] = gene_ldscores_training[info['regression_snp_indices'], global_eqtl_category_index] + gene_ldscore
			gene_ldscores_validation[info_validation['regression_snp_indices'], global_eqtl_category_index] = gene_ldscores_validation[info_validation['regression_snp_indices'], global_eqtl_category_index] + gene_ldscore_validation
			n_genes[global_eqtl_category_index] = n_genes[global_eqtl_category_index] + 1
	
	# Put gene cis-snp h2s into array
	gene_cis_snp_h2s = []
	for eqtl_category_name in eqtl_category_names:
		gene_cis_snp_h2s.append(np.mean(cis_snp_h2_dicti[eqtl_category_name]))
	gene_cis_snp_h2s = np.asarray(gene_cis_snp_h2s)


	return gene_ldscores_training, gene_ldscores_validation, gene_cis_snp_h2s, n_genes

def block_bootstrap_the_regression_to_get_nonparametric_fstat(YY, XX, ww, step1_regression_method, n_bootstraps=100, n_blocks=100):
	# Dimension of space
	NN, KK = XX.shape
	# Split indices into blocks
	blocks = np.array_split(np.arange(NN), n_blocks)
	
	# Storage for bootstrap coefficients
	sampled_coefs = []

	for i in range(n_bootstraps):
		# Sample blocks with replacement
		sampled_block_idxs = np.random.choice(n_blocks, size=n_blocks, replace=True)
		# Gather indices for this bootstrap sample
		sample_indices = np.concatenate([blocks[idx] for idx in sampled_block_idxs])
		
		# Get bootstrap sample
		Y_boot = YY[sample_indices]
		X_boot = XX[sample_indices, :]
		ww_boot = ww[sample_indices]
				
		# Fit OLS
		if step1_regression_method == 'all_snps':
			model = sm.WLS(Y_boot, X_boot, weights=ww_boot).fit()
		elif step1_regression_method == 'non_zero_snps':
			valid_indices = Y_boot != 0.0
			model = sm.WLS(Y_boot[valid_indices], X_boot[valid_indices, :], weights=ww_boot[valid_indices]).fit()
		else:
			print('assumtpioner roror: Not yet implemented')
			pdb.set_trace()
		sampled_coefs.append(model.params)

	# Put into nice matrix
	sampled_coefs = np.asarray(sampled_coefs)

	mean_beta = np.mean(sampled_coefs,axis=0)
	cov_est = np.cov(np.transpose(sampled_coefs))

	mean_beta_minus_intercept = mean_beta[1:]
	cov_est_minus_intercept = cov_est[1:,1:]

	F_stat = np.dot(np.dot(mean_beta_minus_intercept, np.linalg.inv(cov_est_minus_intercept)), mean_beta_minus_intercept)/len(mean_beta_minus_intercept)

	return F_stat


def single_run_of_uncalibrated_mesc(gwas_variant_ld_scores, gene_ldscores_training1, gene_ldscores_validation1,gene_ldscores_training2, gene_ldscores_validation2, gwas_E_beta_sq, regression_weights, squared_eqtl_effect_threshold, n_reference_snps, n_genetic_genes1, n_genetic_genes2, full_eqtl_dataset_names, eqtl_dataset_names):
	##############
	# Data set 1
	# Filter regression variables
	gwas_variant_ld_scores1, gene_ldscores_training1, gene_ldscores_validation1, gwas_E_beta_sq1, regression_weights1 = filter_regression_variables(gwas_variant_ld_scores, gene_ldscores_training1, gene_ldscores_validation1, gwas_E_beta_sq, regression_weights, squared_eqtl_effect_threshold)


	##############
	# Data set 2
	# Filter regression variables
	gwas_variant_ld_scores2, gene_ldscores_training2, gene_ldscores_validation2, gwas_E_beta_sq2, regression_weights2 = filter_regression_variables(gwas_variant_ld_scores, gene_ldscores_training2, gene_ldscores_validation2, gwas_E_beta_sq, regression_weights, squared_eqtl_effect_threshold)


	# Run Mesc in data set 1
	est_nm_h2s1, est_med_h2s1 = run_mesc(gwas_E_beta_sq1, gwas_variant_ld_scores1, n_reference_snps, gene_ldscores_training1, n_genetic_genes1, regression_weights1)
	est_dataset_med_h2s1 = get_per_dataset_med_h2(est_med_h2s1, full_eqtl_dataset_names, eqtl_dataset_names)
	est_nm_h2s2, est_med_h2s2 = run_mesc(gwas_E_beta_sq2, gwas_variant_ld_scores2, n_reference_snps, gene_ldscores_training2, n_genetic_genes2, regression_weights2)
	est_dataset_med_h2s2 = get_per_dataset_med_h2(est_med_h2s2, full_eqtl_dataset_names, eqtl_dataset_names)

	# Average two
	est_nm_h2s = (est_nm_h2s1 + est_nm_h2s2)/2.0
	est_med_h2s = (est_med_h2s1 + est_med_h2s2)/2.0
	est_dataset_med_h2s = (est_dataset_med_h2s1 + est_dataset_med_h2s2)/2.0

	return est_nm_h2s, est_med_h2s, est_dataset_med_h2s


def single_run_of_calibrated_mesc(gwas_variant_ld_scores, gene_ldscores_training1, gene_ldscores_validation1,gene_ldscores_training2, gene_ldscores_validation2, gwas_E_beta_sq, regression_weights, squared_eqtl_effect_threshold, n_reference_snps, n_genetic_genes1, n_genetic_genes2, full_eqtl_dataset_names, eqtl_dataset_names, step1_regression_method, inference_approach):
	##############
	# Data set 1
	# Filter regression variables
	gwas_variant_ld_scores1, gene_ldscores_training1, gene_ldscores_validation1, gwas_E_beta_sq1, regression_weights1 = filter_regression_variables(gwas_variant_ld_scores, gene_ldscores_training1, gene_ldscores_validation1, gwas_E_beta_sq, regression_weights, squared_eqtl_effect_threshold)
	# Use validation data to extract calibrated gene ldscores
	if inference_approach == '2SLS':
		calibrated_gene_ldscores1 = use_validation_data_to_get_calibrated_gene_ldscores(gwas_variant_ld_scores1, gene_ldscores_training1, gene_ldscores_validation1, regression_weights1, step1_regression_method)
	elif inference_approach == 'JIVE':
		calibrated_gene_ldscores1 = use_validation_data_to_get_calibrated_gene_ldscores_JIVE(gwas_variant_ld_scores1, gene_ldscores_training1, gene_ldscores_validation1, regression_weights1, step1_regression_method)


	##############
	# Data set 2
	# Filter regression variables
	gwas_variant_ld_scores2, gene_ldscores_training2, gene_ldscores_validation2, gwas_E_beta_sq2, regression_weights2 = filter_regression_variables(gwas_variant_ld_scores, gene_ldscores_training2, gene_ldscores_validation2, gwas_E_beta_sq, regression_weights, squared_eqtl_effect_threshold)
	##############################
	# Use validation data to extract calibrated gene ldscores
	if inference_approach == '2SLS':
		calibrated_gene_ldscores2 = use_validation_data_to_get_calibrated_gene_ldscores(gwas_variant_ld_scores2, gene_ldscores_training2, gene_ldscores_validation2, regression_weights2, step1_regression_method)
	elif inference_approach == 'JIVE':
		calibrated_gene_ldscores2 = use_validation_data_to_get_calibrated_gene_ldscores_JIVE(gwas_variant_ld_scores2, gene_ldscores_training2, gene_ldscores_validation2, regression_weights2, step1_regression_method)

	# Run Mesc in data set 1
	est_nm_h2s1, est_med_h2s1 = run_mesc(gwas_E_beta_sq1, gwas_variant_ld_scores1, n_reference_snps, calibrated_gene_ldscores1, n_genetic_genes1, regression_weights1)
	est_dataset_med_h2s1 = get_per_dataset_med_h2(est_med_h2s1, full_eqtl_dataset_names, eqtl_dataset_names)
	est_nm_h2s2, est_med_h2s2 = run_mesc(gwas_E_beta_sq2, gwas_variant_ld_scores2, n_reference_snps, calibrated_gene_ldscores2, n_genetic_genes2, regression_weights2)
	est_dataset_med_h2s2 = get_per_dataset_med_h2(est_med_h2s2, full_eqtl_dataset_names, eqtl_dataset_names)



	#ivolsmod = IV2SLS(gwas_E_beta_sq1, instruments=np.hstack((np.ones(gwas_variant_ld_scores1.shape),gene_ldscores_training1)), endog=gene_ldscores_validation1, exog=gwas_variant_ld_scores1, weights=regression_weights1).fit()
	#ivolsmod = IVLIML(gwas_E_beta_sq2, instruments=np.hstack((np.ones(gwas_variant_ld_scores2.shape),gene_ldscores_training2)), endog=gene_ldscores_validation2, exog=gwas_variant_ld_scores2, weights=regression_weights2).fit()
	#ivolsmod = IVGMM(gwas_E_beta_sq2, instruments=np.hstack((np.ones(gwas_variant_ld_scores2.shape),gene_ldscores_training2)), endog=gene_ldscores_validation2, exog=gwas_variant_ld_scores2).fit(cov_type='kernel')
	#tmp_params = np.asarray(ivolsmod.params)


	# Average two
	est_nm_h2s = (est_nm_h2s1 + est_nm_h2s2)/2.0
	est_med_h2s = (est_med_h2s1 + est_med_h2s2)/2.0
	est_dataset_med_h2s = (est_dataset_med_h2s1 + est_dataset_med_h2s2)/2.0

	return est_nm_h2s, est_med_h2s, est_dataset_med_h2s


def run_uncalibrated_mesc(gwas_variant_ld_scores_full, gene_ldscores_training1_full, gene_ldscores_validation1_full, gene_ldscores_training2_full, gene_ldscores_validation2_full, gwas_E_beta_sq_full, regression_weights_full, squared_eqtl_effect_threshold, full_eqtl_dataset_names, eqtl_dataset_names, n_reference_snps, n_genetic_genes1, n_genetic_genes2, avg_eqtl_h2s_1, avg_eqtl_h2s_2, n_blocks=100):
	# Dimension of space
	NN = len(gwas_variant_ld_scores_full)

	# Split indices into blocks
	genomic_blocks = np.array_split(np.arange(NN), n_blocks)	

	# Keep track of jacknifed samples
	jk_nm_h2s = []
	jk_med_h2s = []
	jk_dataset_med_h2s = []
	jk_total_med_h2s = []

	# Run jacknife iterations
	for jk_iter in range(n_blocks):
		# Get blocks corresponding to this JK iteration
		all_block_indices = np.arange(n_blocks)
		jk_block_indices = np.concatenate((all_block_indices[:jk_iter], all_block_indices[(jk_iter+1):]))

		# Gather variant indices for this jacknifed sample
		sample_indices = np.concatenate([genomic_blocks[idx] for idx in jk_block_indices])

		est_nm_h2s, est_med_h2s, est_dataset_med_h2s = single_run_of_uncalibrated_mesc(gwas_variant_ld_scores_full[sample_indices], gene_ldscores_training1_full[sample_indices, :], gene_ldscores_validation1_full[sample_indices, :],gene_ldscores_training2_full[sample_indices,:], gene_ldscores_validation2_full[sample_indices, :], gwas_E_beta_sq_full[sample_indices], regression_weights_full[sample_indices], squared_eqtl_effect_threshold, n_reference_snps, n_genetic_genes1, n_genetic_genes2, full_eqtl_dataset_names, eqtl_dataset_names)

		jk_nm_h2s.append(np.sum(est_nm_h2s))
		jk_med_h2s.append(est_med_h2s)
		jk_dataset_med_h2s.append(est_dataset_med_h2s)
		jk_total_med_h2s.append(np.sum(est_med_h2s))

	jk_nm_h2s = np.asarray(jk_nm_h2s)
	jk_med_h2s = np.asarray(jk_med_h2s)
	jk_dataset_med_h2s = np.asarray(jk_dataset_med_h2s)
	jk_total_med_h2s = np.asarray(jk_total_med_h2s)
	jk_total_h2s = jk_nm_h2s + jk_total_med_h2s



	uncalibrated_mesc_res = {}
	uncalibrated_mesc_res['nm_h2_jk_mean'] = np.mean(jk_nm_h2s)
	uncalibrated_mesc_res['nm_h2_jk_se'] = np.sqrt(np.var(jk_nm_h2s)*(n_blocks-1))
	uncalibrated_mesc_res['med_h2_jk_mean'] = np.mean(jk_med_h2s,axis=0)
	uncalibrated_mesc_res['med_h2_jk_se'] = np.sqrt(np.var(jk_med_h2s,axis=0)*(n_blocks-1))
	uncalibrated_mesc_res['dataset_med_h2_jk_mean'] = np.mean(jk_dataset_med_h2s,axis=0)
	uncalibrated_mesc_res['dataset_med_h2_jk_se'] = np.sqrt(np.var(jk_dataset_med_h2s,axis=0)*(n_blocks-1))
	uncalibrated_mesc_res['total_med_h2_jk_mean'] = np.mean(jk_total_med_h2s)
	uncalibrated_mesc_res['total_med_h2_jk_se'] = np.sqrt(np.var(jk_total_med_h2s)*(n_blocks-1))
	uncalibrated_mesc_res['total_h2_jk_mean'] = np.mean(jk_total_h2s)
	uncalibrated_mesc_res['total_h2_jk_se'] = np.sqrt(np.var(jk_total_h2s)*(n_blocks-1))


	est_nm_h2s, est_med_h2s, est_dataset_med_h2s = single_run_of_uncalibrated_mesc(gwas_variant_ld_scores_full, gene_ldscores_training1_full, gene_ldscores_validation1_full,gene_ldscores_training2_full, gene_ldscores_validation2_full, gwas_E_beta_sq_full, regression_weights_full, squared_eqtl_effect_threshold, n_reference_snps, n_genetic_genes1, n_genetic_genes2, full_eqtl_dataset_names, eqtl_dataset_names)
	
	uncalibrated_mesc_res['nm_h2'] = np.sum(est_nm_h2s)
	uncalibrated_mesc_res['med_h2'] = est_med_h2s
	uncalibrated_mesc_res['dataset_med_h2'] = est_dataset_med_h2s
	uncalibrated_mesc_res['total_med_h2'] = np.mean(np.sum(est_med_h2s))
	uncalibrated_mesc_res['total_h2'] = np.sum(est_nm_h2s) + np.sum(est_med_h2s)

	agg_mean_eqtl_h2 = np.mean([np.mean(avg_eqtl_h2s_1), np.mean(avg_eqtl_h2s_2)])
	uncalibrated_mesc_res['eqtl_h2'] = agg_mean_eqtl_h2

	return uncalibrated_mesc_res


def run_calibrated_mesc_via_bootstrapping(gwas_variant_ld_scores_full, gene_ldscores_training1_full, gene_ldscores_validation1_full, gene_ldscores_training2_full, gene_ldscores_validation2_full, gwas_E_beta_sq_full, regression_weights_full, squared_eqtl_effect_threshold, full_eqtl_dataset_names, eqtl_dataset_names, n_reference_snps, n_genetic_genes1, n_genetic_genes2, avg_eqtl_h2s_1, avg_eqtl_h2s_2, step1_regression_method, n_bootstraps=50, n_blocks=100):
	# Dimension of space
	NN = len(gwas_variant_ld_scores_full)

	# Split indices into blocks
	blocks = np.array_split(np.arange(NN), n_blocks)	

	# Keep track of bootstrapped samples
	bs_nm_h2s = []
	bs_med_h2s = []
	bs_dataset_med_h2s = []
	bs_total_med_h2s = []

	# Run bootstrap iterations
	for i in range(n_bootstraps):
		# Sample blocks with replacement
		sampled_block_idxs = np.random.choice(n_blocks, size=n_blocks, replace=True)
		# Gather indices for this bootstrap sample
		sample_indices = np.concatenate([blocks[idx] for idx in sampled_block_idxs])

		est_nm_h2s, est_med_h2s, est_dataset_med_h2s = single_run_of_calibrated_mesc(gwas_variant_ld_scores_full[sample_indices], gene_ldscores_training1_full[sample_indices, :], gene_ldscores_validation1_full[sample_indices, :],gene_ldscores_training2_full[sample_indices,:], gene_ldscores_validation2_full[sample_indices, :], gwas_E_beta_sq_full[sample_indices], regression_weights_full[sample_indices], squared_eqtl_effect_threshold, n_reference_snps, n_genetic_genes1, n_genetic_genes2, full_eqtl_dataset_names, eqtl_dataset_names, step1_regression_method)

		bs_nm_h2s.append(np.sum(est_nm_h2s))
		bs_med_h2s.append(est_med_h2s)
		bs_dataset_med_h2s.append(est_dataset_med_h2s)
		bs_total_med_h2s.append(np.sum(est_med_h2s))

	bs_nm_h2s = np.asarray(bs_nm_h2s)
	bs_med_h2s = np.asarray(bs_med_h2s)
	bs_dataset_med_h2s = np.asarray(bs_dataset_med_h2s)
	bs_total_med_h2s = np.asarray(bs_total_med_h2s)

	calibrated_mesc_res = {}

	calibrated_mesc_res['nm_h2_bs_avg'] = np.mean(bs_nm_h2s)
	calibrated_mesc_res['nm_h2_bs_std'] = np.std(bs_nm_h2s,ddof=1)
	calibrated_mesc_res['med_h2_bs_avg'] = np.mean(bs_med_h2s,axis=0)
	calibrated_mesc_res['med_h2_bs_std'] = np.std(bs_med_h2s,axis=0, ddof=1)
	calibrated_mesc_res['dataset_med_h2_bs_avg'] = np.mean(bs_dataset_med_h2s,axis=0)
	calibrated_mesc_res['dataset_med_h2_bs_std'] = np.std(bs_dataset_med_h2s,axis=0, ddof=1)
	calibrated_mesc_res['total_med_h2_bs_avg'] = np.mean(bs_total_med_h2s)
	calibrated_mesc_res['total_med_h2_bs_std'] = np.std(bs_total_med_h2s,ddof=1)


	est_nm_h2s, est_med_h2s, est_dataset_med_h2s = single_run_of_calibrated_mesc(gwas_variant_ld_scores_full, gene_ldscores_training1_full, gene_ldscores_validation1_full,gene_ldscores_training2_full, gene_ldscores_validation2_full, gwas_E_beta_sq_full, regression_weights_full, squared_eqtl_effect_threshold, n_reference_snps, n_genetic_genes1, n_genetic_genes2, full_eqtl_dataset_names, eqtl_dataset_names, step1_regression_method)
	
	calibrated_mesc_res['nm_h2'] = np.sum(est_nm_h2s)
	calibrated_mesc_res['med_h2'] = est_med_h2s
	calibrated_mesc_res['dataset_med_h2'] = est_dataset_med_h2s
	calibrated_mesc_res['total_med_h2'] = np.mean(np.sum(est_med_h2s))



	agg_mean_eqtl_h2 = np.mean([np.mean(avg_eqtl_h2s_1), np.mean(avg_eqtl_h2s_2)])
	calibrated_mesc_res['eqtl_h2'] = agg_mean_eqtl_h2


	return calibrated_mesc_res



def run_calibrated_mesc(gwas_variant_ld_scores_full, gene_ldscores_training1_full, gene_ldscores_validation1_full, gene_ldscores_training2_full, gene_ldscores_validation2_full, gwas_E_beta_sq_full, regression_weights_full, squared_eqtl_effect_threshold, full_eqtl_dataset_names, eqtl_dataset_names, n_reference_snps, n_genetic_genes1, n_genetic_genes2, avg_eqtl_h2s_1, avg_eqtl_h2s_2, step1_regression_method, inference_approach, n_blocks=100):
	# Dimension of space
	NN = len(gwas_variant_ld_scores_full)

	# Split indices into blocks
	genomic_blocks = np.array_split(np.arange(NN), n_blocks)	

	# Keep track of Jacknifed samples
	jk_nm_h2s = []
	jk_med_h2s = []
	jk_dataset_med_h2s = []
	jk_total_med_h2s = []
	'''
	# Run jacknife iterations
	for jk_iter in range(n_blocks):
		# Get blocks corresponding to this JK iteration
		all_block_indices = np.arange(n_blocks)
		jk_block_indices = np.concatenate((all_block_indices[:jk_iter], all_block_indices[(jk_iter+1):]))

		# Gather variant indices for this jacknifed sample
		sample_indices = np.concatenate([genomic_blocks[idx] for idx in jk_block_indices])

		est_nm_h2s, est_med_h2s, est_dataset_med_h2s = single_run_of_calibrated_mesc(gwas_variant_ld_scores_full[sample_indices], gene_ldscores_training1_full[sample_indices, :], gene_ldscores_validation1_full[sample_indices, :],gene_ldscores_training2_full[sample_indices,:], gene_ldscores_validation2_full[sample_indices, :], gwas_E_beta_sq_full[sample_indices], regression_weights_full[sample_indices], squared_eqtl_effect_threshold, n_reference_snps, n_genetic_genes1, n_genetic_genes2, full_eqtl_dataset_names, eqtl_dataset_names, step1_regression_method, inference_approach)

		jk_nm_h2s.append(np.sum(est_nm_h2s))
		jk_med_h2s.append(est_med_h2s)
		jk_dataset_med_h2s.append(est_dataset_med_h2s)
		jk_total_med_h2s.append(np.sum(est_med_h2s))

	jk_nm_h2s = np.asarray(jk_nm_h2s)
	jk_med_h2s = np.asarray(jk_med_h2s)
	jk_dataset_med_h2s = np.asarray(jk_dataset_med_h2s)
	jk_total_med_h2s = np.asarray(jk_total_med_h2s)
	jk_total_h2s = jk_nm_h2s + jk_total_med_h2s

	calibrated_mesc_res = {}
	calibrated_mesc_res['nm_h2_jk_mean'] = np.mean(jk_nm_h2s)
	calibrated_mesc_res['nm_h2_jk_se'] = np.sqrt(np.var(jk_nm_h2s)*(n_blocks-1))
	calibrated_mesc_res['med_h2_jk_mean'] = np.mean(jk_med_h2s,axis=0)
	calibrated_mesc_res['med_h2_jk_se'] = np.sqrt(np.var(jk_med_h2s,axis=0)*(n_blocks-1))
	calibrated_mesc_res['dataset_med_h2_jk_mean'] = np.mean(jk_dataset_med_h2s,axis=0)
	calibrated_mesc_res['dataset_med_h2_jk_se'] = np.sqrt(np.var(jk_dataset_med_h2s,axis=0)*(n_blocks-1))
	calibrated_mesc_res['total_med_h2_jk_mean'] = np.mean(jk_total_med_h2s)
	calibrated_mesc_res['total_med_h2_jk_se'] = np.sqrt(np.var(jk_total_med_h2s)*(n_blocks-1))
	calibrated_mesc_res['total_h2_jk_mean'] = np.mean(jk_total_h2s)
	calibrated_mesc_res['total_h2_jk_se'] = np.sqrt(np.var(jk_total_h2s)*(n_blocks-1))
	'''

	est_nm_h2s, est_med_h2s, est_dataset_med_h2s = single_run_of_calibrated_mesc(gwas_variant_ld_scores_full, gene_ldscores_training1_full, gene_ldscores_validation1_full,gene_ldscores_training2_full, gene_ldscores_validation2_full, gwas_E_beta_sq_full, regression_weights_full, squared_eqtl_effect_threshold, n_reference_snps, n_genetic_genes1, n_genetic_genes2, full_eqtl_dataset_names, eqtl_dataset_names, step1_regression_method, inference_approach)
	
	calibrated_mesc_res['nm_h2'] = np.sum(est_nm_h2s)
	calibrated_mesc_res['med_h2'] = est_med_h2s
	calibrated_mesc_res['dataset_med_h2'] = est_dataset_med_h2s
	calibrated_mesc_res['total_med_h2'] = np.sum(est_med_h2s)
	calibrated_mesc_res['total_h2'] = np.sum(est_nm_h2s) + np.sum(est_med_h2s)


	agg_mean_eqtl_h2 = np.mean([np.mean(avg_eqtl_h2s_1), np.mean(avg_eqtl_h2s_2)])
	calibrated_mesc_res['eqtl_h2'] = agg_mean_eqtl_h2


	return calibrated_mesc_res



def use_validation_data_to_get_fstatistics(gwas_variant_ld_scores, gene_ldscores, val_gene_ldscores, regression_weights, step1_regression_method):
	intercept = np.ones((gene_ldscores.shape[0], 1))
	X_mat = np.hstack((intercept, gwas_variant_ld_scores, gene_ldscores))

	calibrated_gene_ldscores = np.zeros(val_gene_ldscores.shape)

	parametric_fstats = []
	bootstrapped_fstats = []
	for jj in range(5):

		if step1_regression_method == 'all_snps':
			model = sm.WLS(val_gene_ldscores[:,jj], X_mat, weights=regression_weights).fit()
			coef = np.copy(model.params)
			pred_gene_ldscores = np.dot(X_mat, coef)
		elif step1_regression_method == 'non_zero_snps':
			non_zero_indices = val_gene_ldscores[:,jj] != 0.0
			model = sm.WLS(val_gene_ldscores[non_zero_indices,jj], X_mat[non_zero_indices,:], weights=regression_weights[non_zero_indices]).fit()
			coef = np.copy(model.params)
			small_pred_gene_ldscores = np.dot(X_mat[non_zero_indices, :], coef)
			pred_gene_ldscores = np.zeros(len(val_gene_ldscores[:,jj]))
			pred_gene_ldscores[non_zero_indices] = small_pred_gene_ldscores
		else:
			print('note yet implemented')
			pdb.set_trace()
	
		calibrated_gene_ldscores[:, jj] = np.copy(pred_gene_ldscores)

		parametric_fstat = model.fvalue
		bootstrapped_fstat = block_bootstrap_the_regression_to_get_nonparametric_fstat(val_gene_ldscores[:,jj], X_mat, regression_weights, step1_regression_method)
		print(bootstrapped_fstat)
		parametric_fstats.append(parametric_fstat)
		bootstrapped_fstats.append(bootstrapped_fstat)

	return bootstrapped_fstats, parametric_fstats


def use_validation_data_to_get_calibrated_gene_ldscores(gwas_variant_ld_scores, gene_ldscores, val_gene_ldscores, regression_weights, step1_regression_method):
	intercept = np.ones((gene_ldscores.shape[0], 1))
	X_mat = np.hstack((intercept, gwas_variant_ld_scores, gene_ldscores))

	calibrated_gene_ldscores = np.zeros(val_gene_ldscores.shape)

	parametric_fstats = []
	bootstrapped_fstats = []
	for jj in range(val_gene_ldscores.shape[1]):
		if step1_regression_method == 'all_snps':
			model = sm.WLS(val_gene_ldscores[:,jj], X_mat, weights=regression_weights).fit()
			coef = np.copy(model.params)
			pred_gene_ldscores = np.dot(X_mat, coef)
			calibrated_gene_ldscores[:, jj] = np.copy(pred_gene_ldscores)
		elif step1_regression_method == 'non_zero_snps':
			non_zero_indices = val_gene_ldscores[:,jj] != 0.0
			model = sm.WLS(val_gene_ldscores[non_zero_indices,jj], X_mat[non_zero_indices,:], weights=regression_weights[non_zero_indices]).fit()
			coef = np.copy(model.params)
			small_pred_gene_ldscores = np.dot(X_mat[non_zero_indices, :], coef)
			pred_gene_ldscores = np.zeros(len(val_gene_ldscores[:,jj]))
			pred_gene_ldscores[non_zero_indices] = small_pred_gene_ldscores
			calibrated_gene_ldscores[:, jj] = np.copy(pred_gene_ldscores)
		else:
			print('not yet implemented')
			pdb.set_trace()


	return calibrated_gene_ldscores

def use_validation_data_to_get_calibrated_gene_ldscores_JIVE(gwas_variant_ld_scores, gene_ldscores, val_gene_ldscores, regression_weights, step1_regression_method, n_blocks=99):
	if step1_regression_method != 'all_snps':
		print('assumption eroror: not yet implemented')
		pdb.set_trace()


	intercept = np.ones((gene_ldscores.shape[0], 1))
	X_mat = np.hstack((intercept, gwas_variant_ld_scores, gene_ldscores))
	calibrated_gene_ldscores = np.zeros(val_gene_ldscores.shape)

	NN = val_gene_ldscores.shape[0]

	# Split indices into blocks
	genomic_blocks = np.array_split(np.arange(NN), n_blocks)	


	# Lol: chatgpt helped speed this computation up.
	# But verified it works
	N, K = val_gene_ldscores.shape
	_, J = X_mat.shape
	n_blocks = len(genomic_blocks)

	# Otherwise, just regress without intercept:
	X_aug = X_mat.copy()
	J_aug = J

	# 2) Form weighted design and weighted outcomes:
	sqrt_w = np.sqrt(regression_weights)       # (N,)
	Xw     = X_aug * sqrt_w[:, None]           # (N, J_aug)
	Yw     = val_gene_ldscores * sqrt_w[:, None]  # (N, K)

	# 3) Full cross‐products:
	XtX_full = Xw.T @ Xw       # (J_aug, J_aug)
	XtY_full = Xw.T @ Yw       # (J_aug, K)

	# 4) Precompute each block’s cross‐products:
	XtX_blocks = []
	XtY_blocks = []
	for blk in genomic_blocks:
		Xw_blk = Xw[blk, :]       # (|block|, J_aug)
		Yw_blk = Yw[blk, :]       # (|block|, K)
		XtX_blocks.append(Xw_blk.T @ Xw_blk)  # (J_aug, J_aug)
		XtY_blocks.append(Xw_blk.T @ Yw_blk)  # (J_aug, K)

	# 5) Allocate storage for all jackknife‐betas:
	betas = np.zeros((n_blocks, J_aug, K))

	# 6) Loop over blocks and solve:
	for b in range(n_blocks):
		XtX_jack = XtX_full - XtX_blocks[b]    # (J_aug, J_aug)
		XtY_jack = XtY_full - XtY_blocks[b]    # (J_aug, K)

		# Solve (J_aug × J_aug) · (J_aug × K) = (J_aug × K):
		#   beta[:, k] = (X^T W X)^{-1} (X^T W y_k),
		# for all k = 0,…,K-1, in one go.
		beta_b = np.linalg.solve(XtX_jack, XtY_jack)  # shape (J_aug, K)
		betas[b, :, :] = beta_b

	# 7) Build a blank prediction array of shape (N, K):
	calibrated_gene_ldscores = np.zeros((N, K))

	# 8) For each block, take its rows from the unweighted design X_aug,
	#    multiply by the block‐specific beta, and store into predicted[][]:
	for b, blk in enumerate(genomic_blocks):
		# blk is an array of row‐indices to predict (the “left‐out” samples)
		X_blk = X_aug[blk, :]        # shape = (|blk|, J_aug)
		beta_b = betas[b, :, :]      # shape = (J_aug, K)

		# Matrix‐multiply to get (|blk| × K) predictions for block b:
		pred_blk = X_blk @ beta_b    # shape = (|blk|, K)

		# Insert into the right rows of the full N×K prediction matrix:
		calibrated_gene_ldscores[blk, :] = pred_blk


	return calibrated_gene_ldscores


def run_mesc(gwas_E_beta_sq, gwas_variant_ld_scores, n_reference_snps, gene_ldscores, n_genetic_genes, regression_weights):

	# Get number of non-mediated variant annotation
	n_non_mediated_variant_anno = gwas_variant_ld_scores.shape[1]

	# Input X for mesc regression
	XX = np.hstack((gwas_variant_ld_scores, gene_ldscores))

	# Run regression
	'''
	coef = linear_regression(XX, gwas_E_beta_sq)
	'''
	model = sm.WLS(gwas_E_beta_sq, XX, weights=regression_weights).fit()
	coef = np.copy(model.params)

	# Convert to heritabilities
	nm_h2s = coef[:n_non_mediated_variant_anno]*n_reference_snps
	med_h2s = coef[n_non_mediated_variant_anno:]*n_genetic_genes

	return nm_h2s, med_h2s

def run2SLS_package(gwas_E_beta_sq, gwas_variant_ld_scores, gene_ldscores, val_gene_ldscores, n_reference_snps, n_genetic_genes):
	ivolsmod = IV2SLS(gwas_E_beta_sq, instruments=np.hstack((np.ones(gwas_variant_ld_scores.shape),gene_ldscores)), endog=val_gene_ldscores, exog=gwas_variant_ld_scores).fit()
	tmp_params = np.asarray(ivolsmod.params)

	#ivolsmod = IVGMM(gwas_E_beta_sq, instruments=np.hstack((np.ones(gwas_variant_ld_scores.shape),gene_ldscores)), endog=val_gene_ldscores, exog=gwas_variant_ld_scores).fit()

	# Get number of non-mediated variant annotation
	n_non_mediated_variant_anno = gwas_variant_ld_scores.shape[1]	
	# Convert to heritabilities
	nm_h2s = tmp_params[:n_non_mediated_variant_anno]*n_reference_snps
	med_h2s = tmp_params[n_non_mediated_variant_anno:]*n_genetic_genes

	print(ivolsmod)
	ivolsmod = IVLIML(gwas_E_beta_sq, instruments=np.hstack((np.ones(gwas_variant_ld_scores.shape),gene_ldscores)), endog=val_gene_ldscores, exog=gwas_variant_ld_scores).fit()
	print(ivolsmod)

	return nm_h2s, med_h2s

def runIVGMM_package(gwas_E_beta_sq, gwas_variant_ld_scores, gene_ldscores, val_gene_ldscores, n_reference_snps, n_genetic_genes):
	ivolsmod = IVGMM(gwas_E_beta_sq, instruments=np.hstack((np.ones(gwas_variant_ld_scores.shape),gene_ldscores)), endog=val_gene_ldscores, exog=gwas_variant_ld_scores).fit()
	tmp_params = np.asarray(ivolsmod.params)

	# Get number of non-mediated variant annotation
	n_non_mediated_variant_anno = gwas_variant_ld_scores.shape[1]	
	# Convert to heritabilities
	nm_h2s = tmp_params[:n_non_mediated_variant_anno]*n_reference_snps
	med_h2s = tmp_params[n_non_mediated_variant_anno:]*n_genetic_genes


	return nm_h2s, med_h2s




def run_ashr_style_posterior_on_eqtl_data(gene_info, genes, eqtl_category_names):

	######################
	# Convert gene info betas and beta_ses into a long vector
	######################
	for eqtl_category_name in eqtl_category_names:
		# Create arrays to keep track of snp gene info
		snp_gene_info_gene_arr = []
		snp_gene_info_gene_index_arr = []
		snp_gene_info_gene_column_index_arr = []

		betas = []
		beta_ses = []

		for gene in genes:

			gene_info[gene]['ashr_squared_sumstats'] = np.copy(gene_info[gene]['squared_sumstats'])

			for kk, gene_eqtl_category_name in enumerate(gene_info[gene]['eqtl_category_names']):
				if gene_eqtl_category_name != eqtl_category_name:
					continue
				gene_eqtl_betas = np.copy(gene_info[gene]['beta_sumstats'][:,kk])
				gene_eqtl_beta_ses = np.copy(gene_info[gene]['beta_se_sumstats'][:,kk])

				for ii in range(len(gene_eqtl_betas)):
					if np.isnan(gene_eqtl_betas[ii]) or np.isnan(gene_eqtl_beta_ses[ii]):
						continue
					betas.append(gene_eqtl_betas[ii])
					beta_ses.append(gene_eqtl_beta_ses[ii])
					snp_gene_info_gene_arr.append(gene)
					snp_gene_info_gene_index_arr.append(ii)
					snp_gene_info_gene_column_index_arr.append(kk)
		# Convert to nice arrays
		betas = np.copy(np.asarray(betas))
		betas_ses = np.copy(np.asarray(beta_ses))
		snp_gene_info_gene_arr = np.asarray(snp_gene_info_gene_arr)
		snp_gene_info_gene_index_arr = np.asarray(snp_gene_info_gene_index_arr)
		snp_gene_info_gene_column_index_arr = np.asarray(snp_gene_info_gene_column_index_arr)

		######################
		# Get posteriors using betas and betas_ses
		######################
		#posterior_pred_beta_sq = np.square(betas) - np.square(betas_ses)
		# Ash style model
		obj = gibbs_ash.model(burn_in_iter=10, max_iter=20, alpha_0=1e-10, variance_grid_lb_divisor=10)
		obj.fit(betas, betas_ses)
		# Updata posteriors
		posterior_pred_beta_sq = np.square(obj.sampled_beta_mean) + obj.sampled_beta_var



		######################
		# Re-fill in gene info with updated posteriors
		######################
		for ii, gene_name in enumerate(snp_gene_info_gene_arr):
			gene_index = snp_gene_info_gene_index_arr[ii]
			column_index = snp_gene_info_gene_column_index_arr[ii]

			# error checking
			if gene_info[gene_name]['squared_sumstats'][gene_index, column_index] != np.square(betas[ii]) - np.square(beta_ses[ii]):
				print('assumption erororr')
				pdb.set_trace()

			gene_info[gene_name]['ashr_squared_sumstats'][gene_index, column_index] = posterior_pred_beta_sq[ii]



	return gene_info

def create_mapping_from_rsid_to_variant_standard_deviation(variant_stdev_filestem, chrom_arr):
	dicti = {}

	for chrom_num in chrom_arr:
		filer = variant_stdev_filestem + chrom_num + '.txt'
		f = open(filer)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			rsid = data[0]
			stdev = float(data[5])
			if rsid in dicti:
				print('assumption eororro')
				pdb.set_trace()
			dicti[rsid] = stdev
		f.close()
	return dicti



def filter_regression_variables(gwas_variant_ld_scores, gene_ldscores_training, gene_ldscores_validation, gwas_E_beta_sq, regression_weights_full, squared_eqtl_effect_threshold):
	filter1 = np.sum(np.isnan(gene_ldscores_training),axis=1) ==0
	filter2 = np.sum(np.abs(gene_ldscores_training) > squared_eqtl_effect_threshold,axis=1) ==0
	filter3 = np.sum(np.isnan(gene_ldscores_validation),axis=1) ==0
	filter4 = np.sum(np.abs(gene_ldscores_validation) > squared_eqtl_effect_threshold,axis=1) ==0
	big_filter = (filter1) & (filter2) & (filter3) & (filter4)

	n_filtered = np.sum(big_filter==False)
	#print('Filtered out ' + str(n_filtered) + ' regression snps')

	return gwas_variant_ld_scores[big_filter], gene_ldscores_training[big_filter, :], gene_ldscores_validation[big_filter, :], gwas_E_beta_sq[big_filter], regression_weights_full[big_filter]


def get_per_dataset_med_h2(med_h2, full_eqtl_dataset_names, eqtl_dataset_names):
	dataset_med_h2s = np.zeros(len(eqtl_dataset_names))

	mapping = {}
	for ii, dataset_name in enumerate(eqtl_dataset_names):
		mapping[dataset_name] = ii


	for ii, category_med_h2 in enumerate(med_h2):
		dataset_name = full_eqtl_dataset_names[ii]

		dataset_med_h2s[mapping[dataset_name]] = dataset_med_h2s[mapping[dataset_name]] + category_med_h2

	return dataset_med_h2s



def get_regression_snp_ldscores(regression_snp_ldscore_filestem, regression_snp_ldscore_filesuffix, chrom_arr):
	regression_snp_ldscores = []
	regrression_snp_rsids = []
	for chrom_num in chrom_arr:
		filer = regression_snp_ldscore_filestem + chrom_num + regression_snp_ldscore_filesuffix
		f = open(filer)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			rsid = data[1]
			ldscore = float(data[6])
			regrression_snp_rsids.append(rsid)
			regression_snp_ldscores.append(ldscore)

		f.close()
	return np.asarray(regression_snp_ldscores), np.asarray(regrression_snp_rsids)


def get_bootstrapped_f_stats_from_1st_stage_least_squares(gwas_variant_ld_scores_full, gene_ldscores_training_full, gene_ldscores_validation_full, gwas_E_beta_sq_full, regression_weights_full, squared_eqtl_effect_threshold, step1_regression_method):
	gwas_variant_ld_scores1, gene_ldscores_training1, gene_ldscores_validation1, gwas_E_beta_sq1, regression_weights1 = filter_regression_variables(gwas_variant_ld_scores_full, gene_ldscores_training_full, gene_ldscores_validation_full, gwas_E_beta_sq_full, regression_weights_full, squared_eqtl_effect_threshold)

	bootstrapped_fstats1, parametric_fstats1 = use_validation_data_to_get_fstatistics(gwas_variant_ld_scores1, gene_ldscores_training1, gene_ldscores_validation1, regression_weights1, step1_regression_method)

	return bootstrapped_fstats1, parametric_fstats1

def print_mesc_regression_results(mesc_obj, output_file):
	t = open(output_file,'w')
	t.write('parameter_name\tparameter_estimate\tJK_mean\tJK_SE\n')

	t.write('total_h2\t' + str(mesc_obj['total_h2']) + '\t' + str(mesc_obj['total_h2_jk_mean']) + '\t' + str(mesc_obj['total_h2_jk_se']) + '\n')

	t.write('nm_h2\t' + str(mesc_obj['nm_h2']) + '\t' + str(mesc_obj['nm_h2_jk_mean']) + '\t' + str(mesc_obj['nm_h2_jk_se']) + '\n')
	t.write('total_med_h2\t' + str(mesc_obj['total_med_h2']) + '\t' + str(mesc_obj['total_med_h2_jk_mean']) + '\t' + str(mesc_obj['total_med_h2_jk_se']) + '\n')

	for ii, eqtl_cat_med_h2 in enumerate(mesc_obj['med_h2']):
		t.write('category_med_h2_' + str(ii) + '\t' + str(eqtl_cat_med_h2) + '\t' + str(mesc_obj['med_h2_jk_mean'][ii]) + '\t' + str(mesc_obj['med_h2_jk_se'][ii]) + '\n')

	for ii, eqtl_dataset_med_h2 in enumerate(mesc_obj['dataset_med_h2']):
		t.write('dataset_med_h2_' + str(ii) + '\t' + str(eqtl_dataset_med_h2) + '\t' + str(mesc_obj['dataset_med_h2_jk_mean'][ii]) + '\t' + str(mesc_obj['dataset_med_h2_jk_se'][ii]) + '\n')

	t.write('eqtl_h2\t' + str(mesc_obj['eqtl_h2']) + '\t' + str(mesc_obj['eqtl_h2']) + '\t' + str(0.0) + '\n')


	t.close()
	return

def print_fstatistics(bootstrapped_fstats1, parametric_fstats1, bootstrapped_fstats2, parametric_fstats2, output_file):
	t = open(output_file,'w')
	t.write('eqtl_category\tfold\tfstat_method\tvalue\n')
	for ii, val in enumerate(bootstrapped_fstats1):
		t.write('category' + str(ii) + '\t' + '1' + '\t' + 'bootstrapped\t' + str(val) + '\n')
	for ii, val in enumerate(bootstrapped_fstats2):
		t.write('category' + str(ii) + '\t' + '2' + '\t' + 'bootstrapped\t' + str(val) + '\n')
	for ii, val in enumerate(parametric_fstats1):
		t.write('category' + str(ii) + '\t' + '1' + '\t' + 'parametric\t' + str(val) + '\n')
	for ii, val in enumerate(parametric_fstats2):
		t.write('category' + str(ii) + '\t' + '2' + '\t' + 'parametric\t' + str(val) + '\n')
	t.close()
	return


def extract_mesc_gene_scores(rsid_to_position, eqtl_summary_file, mesc_expression_score_dir, replication_string, chrom_arr):
	f = open(eqtl_summary_file)
	head_count = 0
	global_gene_ldscores = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_ldscore_vec = np.zeros(len(rsid_to_position))
		file_stem = mesc_expression_score_dir + data[4].split('/')[-1].split('_replicate')[0]
		for chrom_num_str in chrom_arr:
			file_name = file_stem + '_' + replication_string + '_' + chrom_num_str + '.' + chrom_num_str + '.expscore.gz'

			g = gzip.open(file_name)
			head_count2 = 0
			for line2 in g:
				if head_count2 == 0:
					head_count2 = head_count2 + 1
					continue
				line2 = line2.decode('utf-8').rstrip()
				data2 = line2.split()
				rsid = data2[1]
				gene_score = np.sum(np.asarray(data2[3:]).astype(float))
				tmp_pos = rsid_to_position[rsid]
				gene_ldscore_vec[tmp_pos] = gene_ldscore_vec[tmp_pos] + gene_score
			g.close()
		global_gene_ldscores.append(gene_ldscore_vec)
	f.close()
	global_gene_ldscores = np.transpose(np.asarray(global_gene_ldscores))
	return global_gene_ldscores


def get_temp_gene_ldscores(eqtl_dataset_file, rsid_to_position, rsid_to_variant_stdev):
	gene_dicti = {}
	tmp_gene_ldscores = np.zeros(len(rsid_to_position))
	f = open(eqtl_dataset_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensid = data[0]
		rsid = data[1]
		variant_pos = rsid_to_position[rsid]
		variant_sdev = rsid_to_variant_stdev[rsid]
		beta = float(data[6])
		beta_se = float(data[7])
		# Update betas and standard errors according to this
		beta = beta*variant_sdev
		beta_se = beta_se*variant_sdev

		tmp_gene_ldscores[variant_pos] = tmp_gene_ldscores[variant_pos] + np.square(beta) - np.square(beta_se)

		gene_dicti[ensid] = 1
	f.close()

	return tmp_gene_ldscores, len(gene_dicti)


def load_in_eqtl_ldscores_for_single_replicate(eqtl_training_data_summary_file, training_data_eqtl_ldscores_type, standardize_gene_training_scores, rsid_to_position, rsid_to_variant_stdev, eqtl_cis_snp_h2_summary_file, replicate_index):
	scale_variant_effects_by_ref_sdev = False
	if training_data_eqtl_ldscores_type == 'MarginalSS':
		scale_variant_effects_by_ref_sdev = True

	#####################################################
	# First generate mapping from study:gene to cis-snp h2
	#####################################################
	study_gene_to_cis_snp_h2 = {}
	f = open(eqtl_cis_snp_h2_summary_file)
	study_names = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		study_name = data[0]
		study_names.append(study_name)
		study_cis_h2_summary_file = data[(3+replicate_index)]
		head_count2 = 0
		gg = open(study_cis_h2_summary_file)
		for line2 in gg:
			line2 = line2.rstrip()
			data2 = line2.split('\t')
			if head_count2 == 0:
				head_count2 = head_count2 + 1
				continue
			gene_name = data2[0]
			cis_snp_h2 = data2[1]
			study_gene = study_name + ':' + gene_name
			if study_gene in study_gene_to_cis_snp_h2:
				print('assumptione rororor')
				pdb.set_trace()
			study_gene_to_cis_snp_h2[study_gene] = float(cis_snp_h2)
		gg.close()
	f.close()
	study_names = np.asarray(study_names)

	#####################################################
	# Generate gene LD scores
	#####################################################
	gene_ldscores = np.zeros((len(rsid_to_position), len(study_names)))
	n_genes = np.zeros(len(study_names))
	avg_cis_h2s = np.zeros(len(study_names))
	f = open(eqtl_training_data_summary_file)
	head_count = 0
	study_counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		study_name = data[0]
		if study_name != study_names[study_counter]:
			print('assumption error')
			pdb.set_trace()
		sumstat_file = data[(4+replicate_index)]
		gg = open(sumstat_file)
		head_count2 = 0
		study_genes = {}
		study_cis_h2s = []
		for line2 in gg:
			line2 = line2.rstrip()
			data2 = line2.split('\t')
			if head_count2 == 0:
				head_count2 = head_count2 + 1
				continue
			ens_id = data2[0]
			rsid = data2[1]
			beta = float(data2[6])
			beta_se = float(data2[7])
			if scale_variant_effects_by_ref_sdev:
				variant_sdev = rsid_to_variant_stdev[rsid]
				# Update betas and standard errors according to this
				beta = beta*variant_sdev
				beta_se = beta_se*variant_sdev
			variant_pos = rsid_to_position[rsid]
			gene_ldscores[variant_pos, study_counter] = gene_ldscores[variant_pos, study_counter] + np.square(beta) - np.square(beta_se)
			if ens_id in study_genes:
				continue
			study_genes[ens_id] = 1
			study_cis_h2s.append(study_gene_to_cis_snp_h2[study_name + ':' + ens_id])
		gg.close()
		if len(study_cis_h2s) != len(study_genes):
			print('assumption eroror')
			pdb.set_trace()
		n_genes[study_counter] = len(study_genes)
		avg_cis_h2s[study_counter] = np.mean(np.asarray(study_cis_h2s))
		study_counter = study_counter + 1
	f.close()

	return gene_ldscores, n_genes, avg_cis_h2s, study_names







def load_in_eqtl_ldscores(eqtl_training_data_summary_file, training_data_eqtl_ldscores_type, standardize_gene_training_scores, rsid_to_position, rsid_to_variant_stdev, eqtl_cis_snp_h2_summary_file, training=True):
	# Do first for replicate 1
	eqtl_ldscores_rep1, n_genes_rep1, avg_cis_h2s_rep1, dataset_names = load_in_eqtl_ldscores_for_single_replicate(eqtl_training_data_summary_file, training_data_eqtl_ldscores_type, standardize_gene_training_scores, rsid_to_position, rsid_to_variant_stdev, eqtl_cis_snp_h2_summary_file, 0)

	# Do for second replicate
	eqtl_ldscores_rep2, n_genes_rep2, avg_cis_h2s_rep2, dataset_names = load_in_eqtl_ldscores_for_single_replicate(eqtl_training_data_summary_file, training_data_eqtl_ldscores_type, standardize_gene_training_scores, rsid_to_position, rsid_to_variant_stdev, eqtl_cis_snp_h2_summary_file, 1)

	if training:
		return eqtl_ldscores_rep1, eqtl_ldscores_rep2, avg_cis_h2s_rep1, n_genes_rep1, avg_cis_h2s_rep2, n_genes_rep2, dataset_names
	else:
		return eqtl_ldscores_rep2, eqtl_ldscores_rep1, avg_cis_h2s_rep2, n_genes_rep2, avg_cis_h2s_rep1, n_genes_rep1, dataset_names



#####################
# Parse command line arguments
#####################
parser = argparse.ArgumentParser()
parser.add_argument('--gwas-sumstat-file', default='None', type=str,
					help='Full path to GWAS summary stat file')
parser.add_argument('--gwas-N', default='None', type=float,
					help='GWAS sample size')
parser.add_argument('--variant-ldscore-filestem', default='None', type=str,
					help='Per chromosome variant-ldscore-filestem')
parser.add_argument('--variant-ldscore-filesuffix', default='.txt', type=str,
					help='Per chromosome variant-ldscore-file suffix')
parser.add_argument('--variant-M-filestem', default='None', type=str,
					help='Per chromosome variant-M file filestems')
parser.add_argument('--regression-snp-ldscore-filestem', default=None, type=str,
					help='Per chromosome ldscore filestems (corresponding to only regression snps)')
parser.add_argument('--regression-snp-ldscore-filesuffix', default='.txt', type=str,
					help='Per chromosome ldscore filesuffix (corresponding to only regression snps)')
parser.add_argument('--variant-stdev-filestem', default='None', type=str,
					help='Per chromosome file containing standard deviation of each snp')
parser.add_argument('--chromosome-file', default=None, type=str,
					help='File containing chromosomes to run analysis on. If None, set to autosomal chromosomes')
parser.add_argument('--eqtl-training-data-summary-file', default=None, type=str,
					help='File containing one line for each eqtl data set')
parser.add_argument('--eqtl-validation-data-summary-file', default=None, type=str,
					help='File containing one line for each eqtl data set')
parser.add_argument('--standardize-gene-training-scores', default=False, action='store_true',
					help='Boolean on whether or not to Standardize expression scores used for training')
parser.add_argument('--standardize-gene-validation-scores', default=False, action='store_true',
					help='Boolean on whether or not to Standardize expression scores used for validation')
parser.add_argument('--non-mediated-annotation-version', default="full_anno", type=str,
					help='Takes as values "genotype_intercept" or "full_anno"')
parser.add_argument('--gene-trait-architecture', default="stdExpr", type=str,
					help='Takes as values "linear" or "stdExpr"')
parser.add_argument('--step1-regression-method', default="all_snps", type=str,   ## SHOULD GET RID OF?
					help='Takes as values "all_snps" or "non_zero_snps"')
parser.add_argument('--n-expr-cis-h2-bins', default=4, type=int,
					help='Only used if args.gene_trait_architecture=="stdExpr"')
parser.add_argument('--inference-approach', default='2SLS', type=str,
					help='How to perform optimization')
parser.add_argument('--jacknife', default=False, action='store_true',
					help='Boolean on whether or not to Jacknife the estimtes')
parser.add_argument('--squared-eqtl-effect-threshold', default=1.0, type=float,
					help='Filter out any squared eQTL effect with absolute value greater than this value')
parser.add_argument('--mesc-expression-score-dir', default=None, type=str,
					help='Directory containing mesc expression scores')
parser.add_argument('--training-data-eqtl-ldscores-type', default="MarginalSS", type=str,
					help='type of expression scores used for training data')
parser.add_argument('--validation-data-eqtl-ldscores-type', default="MarginalSS", type=str,
					help='type of expression scores used for validation data')
parser.add_argument('--eqtl-cis-snp-h2-summary-file', default="None", type=str,
					help='File containing estimated per-gene cis-snp h2s')
parser.add_argument('--output-stem', default=None, type=str,
					help='Output file stem to save data to')
args = parser.parse_args()
np.random.seed(1)


#####################
# Load in relevent data
#####################
# Names of chromosomes to run analysis on
chrom_arr, chrom_dicti = extract_chromosome_names(args.chromosome_file)

# Extract number of reference snps
n_reference_snps = extract_number_of_reference_snps(args.variant_M_filestem, chrom_arr, args.non_mediated_annotation_version)

# Load in GWAS summary statistics
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(args.gwas_sumstat_file, chrom_dicti)
gwas_E_beta_sq = np.square(gwas_beta) - np.square(gwas_beta_se)

# Create mapping from rsid to variant standard deviation
rsid_to_variant_stdev = create_mapping_from_rsid_to_variant_standard_deviation(args.variant_stdev_filestem, chrom_arr)

# Load in GWAS variant ld scores
gwas_variant_ld_scores_full, gwas_rsids_tmp, non_mediated_annotation_names = get_gwas_variant_ld_scores(args.variant_ldscore_filestem, args.variant_ldscore_filesuffix, chrom_arr, args.non_mediated_annotation_version)

# Load in regression snp ldscores
if args.regression_snp_ldscore_filestem is None:
	regression_weights_full = np.ones(len(gwas_variant_ld_scores_full))
	regression_snp_rsids = np.copy(gwas_rsids_tmp)
else:
	regression_snp_ldscores_full, regression_snp_rsids = get_regression_snp_ldscores(args.regression_snp_ldscore_filestem, args.regression_snp_ldscore_filesuffix, chrom_arr)
	regression_weights_full = 1.0/regression_snp_ldscores_full

# QUick error checking
if np.array_equal(gwas_rsids, gwas_rsids_tmp) == False or np.array_equal(gwas_rsids, regression_snp_rsids) == False:
	print('assumption eroror')
	pdb.set_trace()

# Create mapping from rsid to position
rsid_to_position = create_mapping_from_rsid_to_position(gwas_rsids)


# Quick concern
if args.validation_data_eqtl_ldscores_type != 'MarginalSS':
	print('assumption eroror: need to edit cis-h2 estimates and gene counts I think')
	pdb.set_trace()

# Load in eQTL LD scores for training data
gene_ldscores_training1, gene_ldscores_training2, avg_eqtl_h2s_train_1, n_genes_train_1, avg_eqtl_h2s_train_2, n_genes_train_2, eqtl_dataset_names = load_in_eqtl_ldscores(args.eqtl_training_data_summary_file, args.training_data_eqtl_ldscores_type, args.standardize_gene_training_scores, rsid_to_position, rsid_to_variant_stdev, args.eqtl_cis_snp_h2_summary_file, training=True)

# Load in eQTL LD scores for validation data
gene_ldscores_validation1, gene_ldscores_validation2, avg_eqtl_h2s_1, n_genes_1, avg_eqtl_h2s_2, n_genes_2, eqtl_dataset_names = load_in_eqtl_ldscores(args.eqtl_validation_data_summary_file, args.validation_data_eqtl_ldscores_type, args.standardize_gene_validation_scores, rsid_to_position, rsid_to_variant_stdev, args.eqtl_cis_snp_h2_summary_file, training=False)
full_eqtl_dataset_names = np.copy(eqtl_dataset_names)
tmp = np.copy(avg_eqtl_h2s_1)
avg_eqtl_h2s_1 = np.copy(avg_eqtl_h2s_2)
avg_eqtl_h2s_2 = np.copy(avg_eqtl_h2s_1)

###########
# TO DO
## 2. Get working for raw mesc ld scores
## 3. Get working for processed mesc ld scores
## 4. Get working for standardized run.





'''
##############################
# Load in eqtl Training data
eqtl_dataset_names, eqtl_dataset_Ns_training, eqtl_dataset_Ns_validation, eqtl_dataset_files_training, eqtl_dataset_files_validation, eqtl_dataset_Ns_full, eqtl_dataset_files_full = load_in_eqtl_dataset_summary_file(args.eqtl_summary_file)
genes_training, gene_info_training = load_in_eqtl_data(rsid_to_position, args.gene_ldscore_filestem, args.gene_ldscore_filesuffix, eqtl_dataset_names, chrom_arr)
gene_info_training = fill_in_eqtl_sumstats(gene_info_training, eqtl_dataset_files_training, eqtl_dataset_names, genes_training, chrom_dicti, rsid_to_variant_stdev)


# Create cis h2 bins
if args.gene_trait_architecture == 'linear':
	eqtl_category_names = np.copy(eqtl_dataset_names)
	full_eqtl_dataset_names = np.copy(eqtl_dataset_names)
	gene_info_training = update_gene_info_to_include_eqtl_dataset_names(gene_info_training, args.squared_eqtl_effect_threshold)
elif args.gene_trait_architecture == 'stdExpr':
	gene_info_training, eqtl_category_names, full_eqtl_dataset_names = update_gene_info_to_bin_genes_by_est_cis_h2(gene_info_training, eqtl_dataset_names, args.n_expr_cis_h2_bins, args.squared_eqtl_effect_threshold)
elif args.gene_trait_architecture == 'random_bins':
	gene_info_training, eqtl_category_names, full_eqtl_dataset_names = update_gene_info_to_bin_genes_by_random(gene_info_training, eqtl_dataset_names, args.n_expr_cis_h2_bins, args.squared_eqtl_effect_threshold)
else:
	print('assumption error: ' + str(args.gene_trait_architecture) + ' not an implemented option')
	pdb.set_trace()


##############################
# Load in eqtl Validation data
genes_validation, gene_info_validation = load_in_eqtl_data(rsid_to_position, args.gene_ldscore_filestem, args.gene_ldscore_filesuffix, eqtl_dataset_names, chrom_arr)
gene_info_validation = fill_in_eqtl_sumstats(gene_info_validation, eqtl_dataset_files_validation, eqtl_dataset_names, genes_validation, chrom_dicti,rsid_to_variant_stdev)

if args.gene_trait_architecture == 'linear':
	gene_info_validation = update_gene_info_to_include_eqtl_dataset_names(gene_info_validation, args.squared_eqtl_effect_threshold)
elif args.gene_trait_architecture == 'stdExpr':
	gene_info_validation, eqtl_category_names_validation, full_eqtl_dataset_names = update_gene_info_to_bin_genes_by_est_cis_h2(gene_info_validation, eqtl_dataset_names, args.n_expr_cis_h2_bins, args.squared_eqtl_effect_threshold)
elif args.gene_trait_architecture == 'random_bins':
	gene_info_validation, eqtl_category_names_validation, full_eqtl_dataset_names = update_gene_info_to_bin_genes_by_random(gene_info_validation, eqtl_dataset_names, args.n_expr_cis_h2_bins, args.squared_eqtl_effect_threshold)
else:
	print('erroror: not yet implemented')
	pdb.set_trace()




# Data set 1
gene_ldscores_training1, gene_ldscores_validation1, avg_eqtl_h2s_1, n_genes_1 = extract_gene_ld_scores_for_stdExpr_analysis(genes_training, gene_info_training, gene_info_validation, eqtl_category_names, args.gene_ldscore_type, len(gwas_variant_ld_scores_full), args.squared_eqtl_effect_threshold)
# Data set 2
gene_ldscores_training2, gene_ldscores_validation2, avg_eqtl_h2s_2, n_genes_2 = extract_gene_ld_scores_for_stdExpr_analysis(genes_training, gene_info_validation, gene_info_training, eqtl_category_names, args.gene_ldscore_type, len(gwas_variant_ld_scores_full), args.squared_eqtl_effect_threshold)
if args.gene_trait_architecture == 'linear':
	tmp1 = np.copy(avg_eqtl_h2s_1)
	tmp2 = np.copy(avg_eqtl_h2s_2)
	avg_eqtl_h2s_1 = np.copy(tmp2)
	avg_eqtl_h2s_2 = np.copy(tmp1)
'''


# NOTE: Need eqtl_dataset_names and full_eqtl_dataset_names
# eqtl_dataset_names is ordered list of datasets
# full_eqtl_dataset_names is ordered list of of number of columns in gene_ldscores_validation (with each element being one of the eqtl_dataset_names)


# NOTE: eQTL estimation not in bootstrap. probs ok for now
un_calibrated_mesc_obj = run_uncalibrated_mesc(gwas_variant_ld_scores_full, gene_ldscores_training1, gene_ldscores_validation1, gene_ldscores_training2, gene_ldscores_validation2, gwas_E_beta_sq, regression_weights_full, args.squared_eqtl_effect_threshold, full_eqtl_dataset_names, eqtl_dataset_names, n_reference_snps, n_genes_1*avg_eqtl_h2s_1, n_genes_2*avg_eqtl_h2s_2, avg_eqtl_h2s_1, avg_eqtl_h2s_2)
print_mesc_regression_results(un_calibrated_mesc_obj, args.output_stem + '_uncalibrated_mesc_jk_results.txt')
print(args.output_stem + '_uncalibrated_mesc_jk_results.txt')


if args.training_data_eqtl_ldscores_type != "MarginalSS":
	mesc_gene_scores_training1 = extract_mesc_gene_scores(rsid_to_position, args.eqtl_summary_file, args.mesc_expression_score_dir, 'replicate1', chrom_arr)
	mesc_gene_scores_training2 = extract_mesc_gene_scores(rsid_to_position, args.eqtl_summary_file, args.mesc_expression_score_dir, 'replicate2', chrom_arr)
	if args.step1_gene_ldscores == 'mescLassoPlusMarginalSS':
		gene_ldscores_training1 = np.hstack((gene_ldscores_training1,mesc_gene_scores_training1))
		gene_ldscores_training2 = np.hstack((gene_ldscores_training2,mesc_gene_scores_training2))
	elif args.step1_gene_ldscores == 'mescLasso':
		gene_ldscores_training1 = np.copy(mesc_gene_scores_training1)
		gene_ldscores_training2 = np.copy(mesc_gene_scores_training2)
	else:
		print('not yet implemented')
		pdb.set_trace()



# Get f-statistics
bootstrapped_fstats1, parametric_fstats1 = get_bootstrapped_f_stats_from_1st_stage_least_squares(gwas_variant_ld_scores_full, gene_ldscores_training1, gene_ldscores_validation1, gwas_E_beta_sq, regression_weights_full, args.squared_eqtl_effect_threshold, args.step1_regression_method)
bootstrapped_fstats2, parametric_fstats2 = get_bootstrapped_f_stats_from_1st_stage_least_squares(gwas_variant_ld_scores_full, gene_ldscores_training2, gene_ldscores_validation2, gwas_E_beta_sq, regression_weights_full, args.squared_eqtl_effect_threshold, args.step1_regression_method)
print_fstatistics(bootstrapped_fstats1, parametric_fstats1, bootstrapped_fstats2, parametric_fstats2, args.output_stem + '_step_1_f_stats.txt')
print(args.output_stem + '_step_1_f_stats.txt')



calibrated_mesc_obj = run_calibrated_mesc(gwas_variant_ld_scores_full, gene_ldscores_training1, gene_ldscores_validation1, gene_ldscores_training2, gene_ldscores_validation2, gwas_E_beta_sq, regression_weights_full, args.squared_eqtl_effect_threshold, full_eqtl_dataset_names, eqtl_dataset_names, n_reference_snps, n_genes_1*avg_eqtl_h2s_1, n_genes_2*avg_eqtl_h2s_2, avg_eqtl_h2s_1, avg_eqtl_h2s_2, args.step1_regression_method, args.inference_approach)
print_mesc_regression_results(calibrated_mesc_obj, args.output_stem + '_calibrated_mesc_jk_results.txt')
print(args.output_stem + '_calibrated_mesc_jk_results.txt')


###########
# TO DO:
# 2. JK for F-stats
# 3. Fix Uncalibrated MESC to use the full data
















'''
##############
# Data set 1
# Filter regression variables

gwas_variant_ld_scores1, gene_ldscores_training1, gene_ldscores_validation1, gwas_E_beta_sq1, regression_weights1 = filter_regression_variables(gwas_variant_ld_scores_full, gene_ldscores_training1, gene_ldscores_validation1, gwas_E_beta_sq, regression_weights_full, args.squared_eqtl_effect_threshold)
# Use validation data to extract calibrated gene ldscores
calibrated_gene_ldscores1, bootstrapped_fstats1, parametric_fstats1 = use_validation_data_to_get_calibrated_gene_ldscores(gwas_variant_ld_scores1, gene_ldscores_training1, gene_ldscores_validation1, regression_weights1)

print(parametric_fstats1)
print(bootstrapped_fstats1)

##############
# Data set 2
# Filter regression variables
gwas_variant_ld_scores2, gene_ldscores_training2, gene_ldscores_validation2, gwas_E_beta_sq2, regression_weights2 = filter_regression_variables(gwas_variant_ld_scores_full, gene_ldscores_training2, gene_ldscores_validation2, gwas_E_beta_sq, regression_weights_full, args.squared_eqtl_effect_threshold)
##############################
# Use validation data to extract calibrated gene ldscores
calibrated_gene_ldscores2, bootstrapped_fstats2, parametric_fstats2 = use_validation_data_to_get_calibrated_gene_ldscores(gwas_variant_ld_scores2, gene_ldscores_training2, gene_ldscores_validation2, regression_weights2)




#####################
# Open output file handle
#####################
output_file = args.output_stem + '_calibrated_mesc_results.txt'
t = open(output_file,'w')
t.write('method\test_med_h2\test_nm_h2\test_mean_eqtl_h2\tper_dataset_med_h2\tper_category_med_h2\n')



####################################
# Run non-calibrated approach
# Run MESC
# Data set 1
est_nm_h2s1, est_med_h2s1 = run_mesc(gwas_E_beta_sq1, gwas_variant_ld_scores1, n_reference_snps, gene_ldscores_training1, n_genes_1*avg_eqtl_h2s_1, regression_weights1)
est_dataset_med_h2s1 = get_per_dataset_med_h2(est_med_h2s1, full_eqtl_dataset_names, eqtl_dataset_names)
# Print to output file
t.write('two_step_joint_ldsc_rep1' + '\t' + str(np.sum(est_med_h2s1)) + '\t' + str(np.sum(est_nm_h2s1)) + '\t' + str(np.mean(avg_eqtl_h2s_1)) + '\t' + ','.join(est_dataset_med_h2s1.astype(str)) + '\t' + ','.join(est_med_h2s1.astype(str)) + '\n')
# Data set 2
est_nm_h2s2, est_med_h2s2 = run_mesc(gwas_E_beta_sq2, gwas_variant_ld_scores2, n_reference_snps, gene_ldscores_training2, n_genes_2*avg_eqtl_h2s_2, regression_weights2)
est_dataset_med_h2s2 = get_per_dataset_med_h2(est_med_h2s2, full_eqtl_dataset_names, eqtl_dataset_names)
t.write('two_step_joint_ldsc_rep2' + '\t' + str(np.sum(est_med_h2s2)) + '\t' + str(np.sum(est_nm_h2s2)) + '\t' + str(np.mean(avg_eqtl_h2s_2)) + '\t' + ','.join(est_dataset_med_h2s2.astype(str)) + '\t' + ','.join(est_med_h2s2.astype(str)) + '\n')

# Average two
est_nm_h2s = (est_nm_h2s1 + est_nm_h2s2)/2.0
est_med_h2s = (est_med_h2s1 + est_med_h2s2)/2.0
est_dataset_med_h2s = (est_dataset_med_h2s1 + est_dataset_med_h2s2)/2.0
agg_mean_eqtl_h2 = np.mean([np.mean(avg_eqtl_h2s_1), np.mean(avg_eqtl_h2s_2)])
t.write('two_step_joint_ldsc' + '\t' + str(np.sum(est_med_h2s)) + '\t' + str(np.sum(est_nm_h2s)) + '\t' + str(agg_mean_eqtl_h2) + '\t' + ','.join(est_dataset_med_h2s.astype(str)) + '\t' + ','.join(est_med_h2s.astype(str)) + '\n')



####################################
# Run calibrated approach
# Run MESC
est_nm_h2s1, est_med_h2s1 = run_mesc(gwas_E_beta_sq1, gwas_variant_ld_scores1, n_reference_snps, calibrated_gene_ldscores1, n_genes_1*avg_eqtl_h2s_1, regression_weights1)
est_dataset_med_h2s1 = get_per_dataset_med_h2(est_med_h2s1, full_eqtl_dataset_names, eqtl_dataset_names)
est_nm_h2s2, est_med_h2s2 = run_mesc(gwas_E_beta_sq2, gwas_variant_ld_scores2, n_reference_snps, calibrated_gene_ldscores2, n_genes_2*avg_eqtl_h2s_2, regression_weights2)
est_dataset_med_h2s2 = get_per_dataset_med_h2(est_med_h2s2, full_eqtl_dataset_names, eqtl_dataset_names)


# Print to output file
t.write('calibrated_two_step_joint_ldsc_rep1' + '\t' + str(np.sum(est_med_h2s1)) + '\t' + str(np.sum(est_nm_h2s1)) + '\t' + str(np.mean(avg_eqtl_h2s_1)) + '\t' + ','.join(est_dataset_med_h2s1.astype(str)) + '\t' + ','.join(est_med_h2s1.astype(str)) + '\n')
t.write('calibrated_two_step_joint_ldsc_rep2' + '\t' + str(np.sum(est_med_h2s2)) + '\t' + str(np.sum(est_nm_h2s2)) + '\t' + str(np.mean(avg_eqtl_h2s_2)) + '\t' + ','.join(est_dataset_med_h2s2.astype(str)) + '\t' + ','.join(est_med_h2s2.astype(str)) + '\n')
# Average two
est_nm_h2s = (est_nm_h2s1 + est_nm_h2s2)/2.0
est_med_h2s = (est_med_h2s1 + est_med_h2s2)/2.0
est_dataset_med_h2s = (est_dataset_med_h2s1 + est_dataset_med_h2s2)/2.0
agg_mean_eqtl_h2 = np.mean([np.mean(avg_eqtl_h2s_1), np.mean(avg_eqtl_h2s_2)])
t.write('calibrated_two_step_joint_ldsc' + '\t' + str(np.sum(est_med_h2s)) + '\t' + str(np.sum(est_nm_h2s)) + '\t' + str(agg_mean_eqtl_h2) + '\t' + ','.join(est_dataset_med_h2s.astype(str)) + '\t' + ','.join(est_med_h2s.astype(str)) + '\n')




t.close()
print(output_file)

'''



'''
#################
# Temp saving
np.save('genes.npy', genes)
np.save('gwas_variant_ld_scores.npy', gwas_variant_ld_scores)
np.save('gwas_E_beta_sq.npy', gwas_E_beta_sq)
np.save('n_reference_snps.npy', n_reference_snps)
np.save('gene_ldscores.npy', gene_ldscores)
np.save('val_gene_ldscores.npy', val_gene_ldscores)
np.save('avg_eqtl_h2s.npy', avg_eqtl_h2s)
np.save('val_avg_eqtl_h2s.npy', val_avg_eqtl_h2s)
np.save('n_genes.npy', n_genes)
'''

'''
#################
# Temp loading
genes = np.load('genes.npy')
gwas_variant_ld_scores = np.load('gwas_variant_ld_scores.npy')
gwas_E_beta_sq = np.load('gwas_E_beta_sq.npy')
n_reference_snps = np.load('n_reference_snps.npy')
gene_ldscores = np.load('gene_ldscores.npy')
val_gene_ldscores = np.load('val_gene_ldscores.npy')
avg_eqtl_h2s = np.load('avg_eqtl_h2s.npy')
val_avg_eqtl_h2s = np.load('val_avg_eqtl_h2s.npy')
n_genes = np.load('n_genes.npy')
####################
'''




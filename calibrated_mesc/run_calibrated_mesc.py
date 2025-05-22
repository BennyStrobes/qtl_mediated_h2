import numpy as np
import os
import sys
import pdb
import time
import pickle
#import joint_ldsc
#import joint_ldsc_gibbs
import argparse
import gibbs_ash
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
	filers1 = []
	filers2 = []
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
		filers1.append(data[3])
		filers2.append(data[4])
	f.close()
	return np.asarray(names), np.asarray(sample_sizes1), np.asarray(sample_sizes2), np.asarray(filers1), np.asarray(filers2)


def update_gene_info_to_include_eqtl_dataset_names(gene_info):
	# Get list of genes
	genes = np.asarray([*gene_info])
	# loop through genes
	for gene in genes:
		gene_info[gene]['eqtl_dataset_names'] = np.copy(gene_info[gene]['eqtl_category_names'])

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

def update_gene_info_to_bin_genes_by_random(gene_info, eqtl_dataset_names, n_expr_cis_h2_bins):
	##############################
	# First estimate gene cis h2
	##############################
	# Get list of genes
	genes = np.asarray([*gene_info])

	bin_names = []
	for bin_iter in range(n_expr_cis_h2_bins):
		bin_names.append('bin' + str(bin_iter))
	bin_names = np.asarray(bin_names)

	# create list of new eqtl categories
	new_eqtl_categories = []
	for eqtl_dataset_name in eqtl_dataset_names:
		for bin_name in bin_names:
			new_eqtl_categories.append(eqtl_dataset_name + ':' + bin_name)
	new_eqtl_categories = np.asarray(new_eqtl_categories)


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

	return gene_info, new_eqtl_categories



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
	for eqtl_dataset_name in eqtl_dataset_names:
		for bin_name in bin_names:
			new_eqtl_categories.append(eqtl_dataset_name + ':' + bin_name)
	new_eqtl_categories = np.asarray(new_eqtl_categories)

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

	return gene_info, new_eqtl_categories


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


def use_validation_data_to_get_calibrated_gene_ldscores(gwas_variant_ld_scores, gene_ldscores, val_gene_ldscores):
	intercept = np.ones((gene_ldscores.shape[0], 1))
	X_mat = np.hstack((intercept, gwas_variant_ld_scores, gene_ldscores))

	calibrated_gene_ldscores = np.zeros(gene_ldscores.shape)


	for jj in range(val_gene_ldscores.shape[1]):

		coef = linear_regression(X_mat, val_gene_ldscores[:,jj])

		pred_gene_ldscores = np.dot(X_mat, coef)
	
		calibrated_gene_ldscores[:, jj] = np.copy(pred_gene_ldscores)

	return calibrated_gene_ldscores

def run_mesc(gwas_E_beta_sq, gwas_variant_ld_scores, n_reference_snps, gene_ldscores, n_genetic_genes):

	# Get number of non-mediated variant annotation
	n_non_mediated_variant_anno = gwas_variant_ld_scores.shape[1]

	# Input X for mesc regression
	XX = np.hstack((gwas_variant_ld_scores, gene_ldscores))

	# Run regression
	coef = linear_regression(XX, gwas_E_beta_sq)

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

			for kk, gene_eqtl_category_name in enumerate(gene_info[gene]['eqtl_category_names']):
				if gene_eqtl_category_name != eqtl_category_name:
					continue
				gene_eqtl_betas = np.copy(gene_info[gene]['beta_sumstats'][:,kk])
				gene_eqtl_beta_ses = np.copy(gene_info[gene]['beta_se_sumstats'][:,kk])

				for ii in range(len(gene_eqtl_betas)):
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
		obj = gibbs_ash.model(burn_in_iter=600, max_iter=1000, alpha_0=1e-10, variance_grid_lb_divisor=10)
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

			gene_info[gene_name]['squared_sumstats'][gene_index, column_index] = posterior_pred_beta_sq[ii]



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



def filter_regression_variables(gwas_variant_ld_scores, gene_ldscores_training, gene_ldscores_validation, gwas_E_beta_sq, squared_eqtl_effect_threshold):
	filter1 = np.sum(np.isnan(gene_ldscores_training),axis=1) ==0
	filter2 = np.sum(np.abs(gene_ldscores_training) > squared_eqtl_effect_threshold,axis=1) ==0
	filter3 = np.sum(np.isnan(gene_ldscores_validation),axis=1) ==0
	filter4 = np.sum(np.abs(gene_ldscores_validation) > squared_eqtl_effect_threshold,axis=1) ==0
	big_filter = (filter1) & (filter2) & (filter3) & (filter4)

	n_filtered = np.sum(big_filter==False)
	print('Filtered out ' + str(n_filtered) + ' regression snps')

	return gwas_variant_ld_scores[big_filter], gene_ldscores_training[big_filter, :], gene_ldscores_validation[big_filter, :], gwas_E_beta_sq[big_filter]




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
parser.add_argument('--gene-ldscore-filestem', default='None', type=str,
                    help='Per chromosome gene-ldscore-filestem')
parser.add_argument('--gene-ldscore-filesuffix', default='.txt', type=str,
                    help='Per chromosome gene-ldscore-file suffix')
parser.add_argument('--variant-stdev-filestem', default='None', type=str,
                    help='Per chromosome file containing standard deviation of each snp')
parser.add_argument('--chromosome-file', default=None, type=str,
                    help='File containing chromosomes to run analysis on. If None, set to autosomal chromosomes')
parser.add_argument('--eqtl-summary-file', default=None, type=str,
                    help='File containing one line for each eqtl data set')
parser.add_argument('--non-mediated-annotation-version', default="full_anno", type=str,
                    help='Takes as values "genotype_intercept" or "full_anno"')
parser.add_argument('--gene-trait-architecture', default="stdExpr", type=str,
                    help='Takes as values "linear" or "stdExpr"')
parser.add_argument('--n-expr-cis-h2-bins', default=4, type=int,
                    help='Only used if args.gene_trait_architecture=="stdExpr"')
parser.add_argument('--inference-approach', default='mle', type=str,
                    help='How to perform optimization')
parser.add_argument('--jacknife', default=False, action='store_true',
                    help='Boolean on whether or not to Jacknife the estimtes')
parser.add_argument('--fixed-variance', default=False, action='store_true',
                    help='Boolean on whether the variance parameters should be fixed')
parser.add_argument('--per-dataset-variance', default=False, action='store_true',
                    help='Boolean on whether the variance parameters should be fixed')
parser.add_argument('--gene-ldscore-type', default="sqaured_marginal_sumstats", type=str,
                    help='Takes as values "sqaured_marginal_sumstats" or other more complicated things')
parser.add_argument('--squared-eqtl-effect-threshold', default=1.0, type=float,
                    help='Filter out any squared eQTL effect with absolute value greater than this value')
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
gwas_variant_ld_scores, gwas_rsids_tmp, non_mediated_annotation_names = get_gwas_variant_ld_scores(args.variant_ldscore_filestem, args.variant_ldscore_filesuffix, chrom_arr, args.non_mediated_annotation_version)

# QUick error checking
if np.array_equal(gwas_rsids, gwas_rsids_tmp) == False:
	print('assumption eroror')
	pdb.set_trace()


# Create mapping from rsid to position
rsid_to_position = create_mapping_from_rsid_to_position(gwas_rsids)



##############################
# Load in eqtl Training data
eqtl_dataset_names, eqtl_dataset_Ns_training, eqtl_dataset_Ns_validation, eqtl_dataset_files_training, eqtl_dataset_files_validation = load_in_eqtl_dataset_summary_file(args.eqtl_summary_file)
genes_training, gene_info_training = load_in_eqtl_data(rsid_to_position, args.gene_ldscore_filestem, args.gene_ldscore_filesuffix, eqtl_dataset_names, chrom_arr)
gene_info_training = fill_in_eqtl_sumstats(gene_info_training, eqtl_dataset_files_training, eqtl_dataset_names, genes_training, chrom_dicti, rsid_to_variant_stdev)

# Create cis h2 bins
if args.gene_trait_architecture == 'linear':
	eqtl_category_names = np.copy(eqtl_dataset_names)
	gene_info_training = update_gene_info_to_include_eqtl_dataset_names(gene_info_training)
elif args.gene_trait_architecture == 'stdExpr':
	gene_info_training, eqtl_category_names = update_gene_info_to_bin_genes_by_est_cis_h2(gene_info_training, eqtl_dataset_names, args.n_expr_cis_h2_bins, args.squared_eqtl_effect_threshold)
elif args.gene_trait_architecture == 'random_bins':
	gene_info, eqtl_category_names = update_gene_info_to_bin_genes_by_random(gene_info, eqtl_dataset_names, args.n_expr_cis_h2_bins)
else:
	print('assumption error: ' + str(args.gene_trait_architecture) + ' not an implemented option')
	pdb.set_trace()


if args.gene_ldscore_type == 'ashr_style_pred':
	print('not yet ready')
	gene_ldscores, avg_eqtl_h2s, n_genes = extract_gene_ld_scores(genes, gene_info, eqtl_category_names, args.gene_ldscore_type, len(gwas_variant_ld_scores))
	gene_info = run_ashr_style_posterior_on_eqtl_data(gene_info, genes, eqtl_category_names)
	gene_ldscores, fake_avg_eqtl_h2s, fake_n_genes = extract_gene_ld_scores(genes, gene_info, eqtl_category_names, args.gene_ldscore_type, len(gwas_variant_ld_scores))
else:
	gene_ldscores_training, avg_eqtl_h2s_training, n_genes_training = extract_gene_ld_scores(genes_training, gene_info_training, eqtl_category_names, args.gene_ldscore_type, len(gwas_variant_ld_scores), args.squared_eqtl_effect_threshold)



##############################
# Load in eqtl Validation data
genes_validation, gene_info_validation = load_in_eqtl_data(rsid_to_position, args.gene_ldscore_filestem, args.gene_ldscore_filesuffix, eqtl_dataset_names, chrom_arr)
gene_info_validation = fill_in_eqtl_sumstats(gene_info_validation, eqtl_dataset_files_validation, eqtl_dataset_names, genes_validation, chrom_dicti,rsid_to_variant_stdev)

if args.gene_trait_architecture == 'linear':
	gene_info_validation = update_gene_info_to_include_eqtl_dataset_names(gene_info_validation)
else:
	print('erroror: not yet implemented')
	pdb.set_trace()


gene_ldscores_validation, avg_eqtl_h2s_validation, n_genes_validation = extract_gene_ld_scores(genes_validation, gene_info_validation, eqtl_category_names, "squared_marginal_sumstats", len(gwas_variant_ld_scores), args.squared_eqtl_effect_threshold)


# Quick error check
if np.array_equal(genes_training, genes_validation) == False or np.array_equal(n_genes_validation, n_genes_training) == False:
	print('assumption eroroor')

##############
# Filter regression variables
gwas_variant_ld_scores, gene_ldscores_training, gene_ldscores_validation, gwas_E_beta_sq = filter_regression_variables(gwas_variant_ld_scores, gene_ldscores_training, gene_ldscores_validation, gwas_E_beta_sq, args.squared_eqtl_effect_threshold)


##############################
# Use validation data to extract calibrated gene ldscores
calibrated_gene_ldscores1 = use_validation_data_to_get_calibrated_gene_ldscores(gwas_variant_ld_scores, gene_ldscores_training, gene_ldscores_validation)
calibrated_gene_ldscores2 = use_validation_data_to_get_calibrated_gene_ldscores(gwas_variant_ld_scores, gene_ldscores_validation, gene_ldscores_training)



#####################
# Open output file handle
#####################
output_file = args.output_stem + '_calibrated_mesc_results.txt'
t = open(output_file,'w')
t.write('method\test_med_h2\test_nm_h2\test_mean_eqtl_h2\tper_dataset_med_h2\tper_category_med_h2\n')


if args.gene_trait_architecture != 'linear':
	print('assumptiontoeoror: Need to update print statements')
	pdb.set_trace()
####################################
# Run non-calibrated approach
# Run MESC
# training data
est_nm_h2s1, est_med_h2s1 = run_mesc(gwas_E_beta_sq, gwas_variant_ld_scores, n_reference_snps, gene_ldscores_training, n_genes_training*avg_eqtl_h2s_training)
# Print to output file
t.write('two_step_joint_ldsc_rep1' + '\t' + str(np.sum(est_med_h2s1)) + '\t' + str(np.sum(est_nm_h2s1)) + '\t' + str(np.mean(avg_eqtl_h2s_training)) + '\t' + ','.join(est_med_h2s1.astype(str)) + '\t' + ','.join(est_med_h2s1.astype(str)) + '\n')
# Validation data
est_nm_h2s2, est_med_h2s2 = run_mesc(gwas_E_beta_sq, gwas_variant_ld_scores, n_reference_snps, gene_ldscores_validation, n_genes_validation*avg_eqtl_h2s_validation)
t.write('two_step_joint_ldsc_rep2' + '\t' + str(np.sum(est_med_h2s2)) + '\t' + str(np.sum(est_nm_h2s2)) + '\t' + str(np.mean(avg_eqtl_h2s_validation)) + '\t' + ','.join(est_med_h2s2.astype(str)) + '\t' + ','.join(est_med_h2s2.astype(str)) + '\n')

# Average two
est_nm_h2s = (est_nm_h2s1 + est_nm_h2s2)/2.0
est_med_h2s = (est_med_h2s1 + est_med_h2s2)/2.0
agg_mean_eqtl_h2 = np.mean([np.mean(avg_eqtl_h2s_validation), np.mean(avg_eqtl_h2s_training)])
t.write('two_step_joint_ldsc' + '\t' + str(np.sum(est_med_h2s)) + '\t' + str(np.sum(est_nm_h2s)) + '\t' + str(agg_mean_eqtl_h2) + '\t' + ','.join(est_med_h2s.astype(str)) + '\t' + ','.join(est_med_h2s.astype(str)) + '\n')



####################################
# Run calibrated approach
# Run MESC
est_nm_h2s1, est_med_h2s1 = run_mesc(gwas_E_beta_sq, gwas_variant_ld_scores, n_reference_snps, calibrated_gene_ldscores1, n_genes_training*avg_eqtl_h2s_training)
est_nm_h2s2, est_med_h2s2 = run_mesc(gwas_E_beta_sq, gwas_variant_ld_scores, n_reference_snps, calibrated_gene_ldscores2, n_genes_validation*avg_eqtl_h2s_validation)

'''
# Alt way to get identical results
alt_est_nm_h2s, alt_est_med_h2s = run2SLS_package(gwas_E_beta_sq, gwas_variant_ld_scores, gene_ldscores, val_gene_ldscores, n_reference_snps, n_genes*avg_eqtl_h2s)
if np.max(np.abs(est_med_h2s - alt_est_med_h2s)) > 1e-12:
	print('assumption eroorr')
	pdb.set_trace()
'''


#XX = np.hstack((gwas_variant_ld_scores, calibrated_gene_ldscores))
#coef = linear_regression(XX, gwas_E_beta_sq)
#dicti = {'Y': gwas_E_beta_sq, 'X2': val_gene_ldscores[:,0], 'X3':val_gene_ldscores[:,1], 'X4':val_gene_ldscores[:,2], 'X5':val_gene_ldscores[:,3], 'X6':val_gene_ldscores[:,4], 'X1':gwas_variant_ld_scores[:,0], 'Z2': gene_ldscores[:,0], 'Z3':gene_ldscores[:,1], 'Z4':gene_ldscores[:,2], 'Z5':gene_ldscores[:,3], 'Z6':gene_ldscores[:,4]}
#df = pd.DataFrame(dicti)
#mod = IV2SLS.from_formula("Y ~ X1 + [X2 ~ X1 + Z2 + Z3 + Z4 + Z5 + Z6]+ [X3 ~ X1 + Z2 + Z3 + Z4 + Z5 + Z6] + [X4 ~ X1 + Z2 + Z3 + Z4 + Z5 + Z6] + [X5 ~ X1 + Z2 + Z3 + Z4 + Z5 + Z6]+ [X6 ~ X1 + Z2 + Z3 + Z4 + Z5 + Z6]",data=df)
#ivolsmod = IV2SLS(gwas_E_beta_sq, instruments=np.hstack((np.ones(gwas_variant_ld_scores.shape),gwas_variant_ld_scores,gene_ldscores)), endog=val_gene_ldscores, exog=gwas_variant_ld_scores).fit()
#ivolsmod = IVLIML(gwas_E_beta_sq, instruments=np.hstack((np.ones(gwas_variant_ld_scores.shape),gwas_variant_ld_scores,gene_ldscores)), exog=val_gene_ldscores, endog=gwas_variant_ld_scores).fit()
#iv_mod = IV2SLS(dependent=gwas_E_beta_sq,exog=calibrated_gene_ldscores, endog=gwas_variant_ld_scores,instruments=gene_ldscores).fit()



# Print to output file
t.write('calibrated_two_step_joint_ldsc_rep1' + '\t' + str(np.sum(est_med_h2s1)) + '\t' + str(np.sum(est_nm_h2s1)) + '\t' + str(np.mean(avg_eqtl_h2s_training)) + '\t' + ','.join(est_med_h2s1.astype(str)) + '\t' + ','.join(est_med_h2s1.astype(str)) + '\n')
t.write('calibrated_two_step_joint_ldsc_rep2' + '\t' + str(np.sum(est_med_h2s2)) + '\t' + str(np.sum(est_nm_h2s2)) + '\t' + str(np.mean(avg_eqtl_h2s_validation)) + '\t' + ','.join(est_med_h2s2.astype(str)) + '\t' + ','.join(est_med_h2s2.astype(str)) + '\n')
# Average two
est_nm_h2s = (est_nm_h2s1 + est_nm_h2s2)/2.0
est_med_h2s = (est_med_h2s1 + est_med_h2s2)/2.0
agg_mean_eqtl_h2 = np.mean([np.mean(avg_eqtl_h2s_validation), np.mean(avg_eqtl_h2s_training)])
t.write('calibrated_two_step_joint_ldsc' + '\t' + str(np.sum(est_med_h2s)) + '\t' + str(np.sum(est_nm_h2s)) + '\t' + str(agg_mean_eqtl_h2) + '\t' + ','.join(est_med_h2s.astype(str)) + '\t' + ','.join(est_med_h2s.astype(str)) + '\n')





'''
####################################
# Run calibrated GMM approach
# Run MESC
est_nm_h2s, est_med_h2s = runIVGMM_package(gwas_E_beta_sq, gwas_variant_ld_scores, gene_ldscores, val_gene_ldscores, n_reference_snps, n_genes*avg_eqtl_h2s)
# Print to output file
t.write('gmm_calibrated_two_step_joint_ldsc' + '\t' + str(np.sum(est_med_h2s)) + '\t' + str(np.sum(est_nm_h2s)) + '\t' + str(np.mean(avg_eqtl_h2s)) + '\t' + ','.join(est_med_h2s.astype(str)) + '\t' + ','.join(est_med_h2s.astype(str)) + '\n')
t.write('gmm_calibrated_two_step_joint_ldsc_unscaled' + '\t' + str(np.sum(est_med_h2s)) + '\t' + str(np.sum(est_nm_h2s)) + '\t' + str(np.mean(avg_eqtl_h2s)) + '\t' + ','.join((est_med_h2s/avg_eqtl_h2s).astype(str)) + '\t' + ','.join((est_med_h2s/avg_eqtl_h2s).astype(str)) + '\n')
'''



'''
#####################
# Run optimization
#####################
# Run two step joint ldsc
obj = joint_ldsc.med_h2(version='two_step')
obj.fit(genes, gene_info, gwas_variant_ld_scores, gwas_E_beta_sq, eqtl_category_names, eqtl_dataset_names, n_reference_snps)
t.write('two_step_joint_ldsc' + '\t' + str(np.sum(obj.category_med_h2)) + '\t' + str(obj.nm_h2) + '\t' + str(obj.avg_eqtl_h2) + '\t' + ','.join(obj.dataset_med_h2.astype(str)) + '\t' + ','.join(obj.category_med_h2.astype(str)) + '\n')
obj = joint_ldsc.med_h2(max_iter=1000,convergence_thresh=1e-10, fixed_variance=args.fixed_variance, per_dataset_variance=args.per_dataset_variance)
if args.inference_approach == 'mle':
	# Run joint ldsc
	obj = joint_ldsc.med_h2(max_iter=1000,convergence_thresh=1e-10, fixed_variance=args.fixed_variance, per_dataset_variance=args.per_dataset_variance)
	obj.fit(genes, gene_info, gwas_variant_ld_scores,gwas_E_beta_sq, eqtl_category_names, eqtl_dataset_names, n_reference_snps)
	t.write('joint_ldsc' + '\t' + str(np.sum(obj.category_med_h2)) + '\t' + str(obj.nm_h2) + '\t' + str(obj.avg_eqtl_h2) + '\t' + ','.join(obj.dataset_med_h2.astype(str)) + '\t' + ','.join(obj.category_med_h2.astype(str)) + '\n')
elif args.inference_approach == 'gibbs':
	# Run joint ldsc
	obj = joint_ldsc_gibbs.med_h2(burn_in_iter=500, max_iter=800,fixed_variance=args.fixed_variance, per_dataset_variance=args.per_dataset_variance)
	obj.fit(genes, gene_info, gwas_variant_ld_scores,gwas_E_beta_sq, eqtl_category_names, eqtl_dataset_names, n_reference_snps)
	t.write('joint_ldsc_gibbs' + '\t' + str(np.sum(obj.category_med_h2)) + '\t' + str(obj.nm_h2) + '\t' + str(obj.avg_eqtl_h2) + '\t' + ','.join(obj.dataset_med_h2.astype(str)) + '\t' + ','.join(obj.category_med_h2.astype(str)) + '\n')
'''


t.close()
print(output_file)





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




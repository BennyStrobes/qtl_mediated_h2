import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import os
import pdb
import numpy as np
from pandas_plink import read_plink1_bin
import pickle
import pandas as pd
import pyreadr
import gzip
import time




def get_pseudotissue_names(gtex_pseudotissue_file):
	f = open(gtex_pseudotissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)



def create_gene_model_df(mesc_lasso_file, plink_eqtl_genotype_bim_file, plink_reference_genotype_bim_file):
	# First create mapping from rsid to a1 and a2 for eqtl data
	rsid_to_eqtl_alleles = {}
	f = open(plink_eqtl_genotype_bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		rsid = data[1]
		a1 = data[4]
		a2 = data[5]
		if rsid in rsid_to_eqtl_alleles:
			print('assumption eroror')
			pdb.set_trace()
		rsid_to_eqtl_alleles[rsid] = a1 + '_' + a2
	f.close()

	# First create mapping from rsid to a1 and a2 for reference data
	rsid_to_ref_alleles = {}
	f = open(plink_reference_genotype_bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		rsid = data[1]
		a1 = data[4]
		a2 = data[5]
		if rsid in rsid_to_ref_alleles:
			print('assumption eroror')
			pdb.set_trace()
		rsid_to_ref_alleles[rsid] = a1 + '_' + a2
	f.close()

	# Now create data structure containing lasso predicted causal eqtl effect sizes
	gene_model_df = {}
	f = open(mesc_lasso_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Quick error check
		if len(data) != 4:
			print('assumption eroror')
			pdb.set_trace()

		gene_id = data[0]
		rsid = data[2]
		eqtl_effect = float(data[3])

		# Correct for allele flipping between data sets
		eqtl_alleles = rsid_to_eqtl_alleles[rsid]
		ref_alleles = rsid_to_ref_alleles[rsid]

		# Correct if alleles are flipped
		if eqtl_alleles != ref_alleles:

			# Flip sign of eqtl effect
			eqtl_effect = eqtl_effect*-1.0

			# Quick error check to make sure its a flip
			eqtl_allele_info = eqtl_alleles.split('_')
			if ref_alleles != eqtl_allele_info[1] + '_' + eqtl_allele_info[0]:
				print('assumption errror')
				pdb.set_trace()
		
		# Now add to data structure
		if gene_id not in gene_model_df:
			gene_model_df[gene_id] = []
		gene_model_df[gene_id].append((rsid, eqtl_effect))
	f.close()
	return gene_model_df



def create_mapping_from_rsid_to_snpid(rsids, snp_ids):
	mapping = {}
	for ii, rsid in enumerate(rsids):
		if rsid in mapping:
			print('assumption eroror')
			pdb.set_trace()
		mapping[rsid] = snp_ids[ii]
	return mapping

def create_mapping_from_snpid_to_rsid(rsids, snp_ids):
	mapping = {}
	for ii, rsid in enumerate(rsids):
		if snp_ids[ii] in mapping:
			print('assumption eororor')
			pdb.set_trace()
		mapping[snp_ids[ii]] = rsid
	return mapping

def create_mapping_rsid_to_reference_index(rsids):
	mapping = {}
	for ii, rsid in enumerate(rsids):
		mapping[rsid] = ii
	return mapping

def get_gene_model_snp_positions(variant_names):
	positions = []
	for variant_name in variant_names:
		positions.append(variant_name.split('_')[1])
	return np.asarray(positions).astype(int)

def project_gene_model_onto_full_window(gene_susie_mu, gene_snpid_names, window_snpid_to_reference_index, n_window_snps, sign_aware_model=True):
	# First get number of components
	n_components = gene_susie_mu.shape[0]

	# Initialize gene model on full window
	window_susie_mu = np.zeros((n_components, n_window_snps))

	# Loop through gene snpid names
	for ii, gene_snpid_name in enumerate(gene_snpid_names):
		# Get window position and sign
		if gene_snpid_name not in window_snpid_to_reference_index:
			print('skip snp ' + gene_snpid_name)
			continue
		window_position, sign = window_snpid_to_reference_index[gene_snpid_name]

		if sign_aware_model:
			window_susie_mu[:, window_position] = (gene_susie_mu[:, ii])*sign
		else:
			window_susie_mu[:, window_position] = (gene_susie_mu[:, ii])

	return window_susie_mu

def get_window_regression_snp_indices(G_window_snpids, regression_snp_id_to_regression_snp_position):
	window_regression_snp_indices = []
	global_regression_snp_positions = []

	for snpid in G_window_snpids:
		if snpid in regression_snp_id_to_regression_snp_position:
			window_regression_snp_indices.append(True)
			global_regression_snp_positions.append(regression_snp_id_to_regression_snp_position[snpid])
		else:
			window_regression_snp_indices.append(False)

	# Put into nice array format
	window_regression_snp_indices = np.asarray(window_regression_snp_indices)
	global_regression_snp_positions = np.asarray(global_regression_snp_positions)

	# Quick error checking
	if np.sum(window_regression_snp_indices) != len(global_regression_snp_positions):
		print('assumption eroror')
		pdb.set_trace()

	return window_regression_snp_indices, global_regression_snp_positions

def compute_gene_variance(susie_mu, susie_mu_sd, susie_alpha, ld):
	gene_var = 0.0

	# Component level eqtl effect sizes for this gene		
	gene_component_effect_sizes = (susie_mu)*susie_alpha

	# eQTL effect sizes for this gene
	gene_eqtl_effect_sizes = np.sum(gene_component_effect_sizes,axis=0)


	num_susie_components = susie_mu.shape[0]
	for k_index in range(num_susie_components):
		gene_var = gene_var + np.sum((np.square(susie_mu[k_index,:]) + np.square(susie_mu_sd[k_index,:]))*np.diag(ld)*susie_alpha[k_index,:])
		eqtl_component_pmces = (susie_mu[k_index,:])*(susie_alpha[k_index,:])
		gene_var = gene_var - np.dot(np.dot(eqtl_component_pmces,ld), eqtl_component_pmces)
	gene_var = gene_var + np.dot(np.dot(gene_eqtl_effect_sizes,ld), gene_eqtl_effect_sizes)
				
	return gene_var


def extract_ld_annotations_for_this_gene_region(ld, susie_mu, susie_mu_sd, susie_alpha, variant_indices, n_ref_panel_samples):
	# Number of variants
	num_var = ld.shape[0]

	# Component level eqtl effect sizes for this gene		
	gene_component_effect_sizes = susie_mu*susie_alpha

	# eQTL effect sizes for this gene
	gene_eqtl_effect_sizes = np.sum(gene_component_effect_sizes,axis=0)

	# Compute squared eqtl effect sizes for this gene
	gene_squared_eqtl_effect_sizes = np.sum((np.square(susie_mu) + np.square(susie_mu_sd))*susie_alpha,axis=0) + gene_eqtl_effect_sizes*gene_eqtl_effect_sizes - np.sum(gene_component_effect_sizes*gene_component_effect_sizes,axis=0)

	# E[beta_k*beta_j]
	cross_terms = np.dot(np.reshape(gene_eqtl_effect_sizes, (num_var,1)), np.reshape(gene_eqtl_effect_sizes, (1,num_var))) - np.dot(np.transpose(gene_component_effect_sizes), gene_component_effect_sizes)

	# Compute gene variance
	gene_variance = compute_gene_variance(susie_mu, susie_mu_sd, susie_alpha, ld)

	# Compute ld scores for diagonal piece
	diagonal_ld_scores = np.sum((np.square(ld[variant_indices,:])*gene_squared_eqtl_effect_sizes),axis=1)
	
	# Comptute ld scores for off diagonal elements
	np.fill_diagonal(cross_terms, np.zeros(num_var))
	non_diagonal_ld_scores = np.sum(np.dot(cross_terms, ld[:, variant_indices])*ld[:,variant_indices],axis=0)

	# Generate complete ld scores
	ld_scores = (diagonal_ld_scores + non_diagonal_ld_scores)/gene_variance

	# Get adjusted ld scores
	adj_ld_scores = ld_scores - ((1.0-ld_scores)/(n_ref_panel_samples-2.0))

	return ld_scores, adj_ld_scores

def extract_ld_annotations_for_this_gene_region_with_eqtl_point_estimate(ld, gene_eqtl_effect_sizes, variant_indices, n_ref_panel_samples, gene_variance, version='B'):
	#gene_variance = np.dot(np.dot(gene_eqtl_effect_sizes, ld), gene_eqtl_effect_sizes)
	if version == 'A':
		# Number of variants
		num_var = ld.shape[0]

		# Compute squared eqtl effect sizes for this gene
		gene_squared_eqtl_effect_sizes = np.square(gene_eqtl_effect_sizes)

		# E[beta_k*beta_j]
		cross_terms = np.dot(np.reshape(gene_eqtl_effect_sizes, (num_var,1)), np.reshape(gene_eqtl_effect_sizes, (1,num_var)))


		diagonal_ld_scores = np.sum((np.square(ld[variant_indices,:])*gene_squared_eqtl_effect_sizes),axis=1)

		# Temp
		np.fill_diagonal(cross_terms, np.zeros(num_var))
		non_diagonal_ld_scores = np.sum(np.dot(cross_terms, ld[:, variant_indices])*ld[:,variant_indices],axis=0)

		# Generate complete ld scores
		ld_scores = (diagonal_ld_scores + non_diagonal_ld_scores)/gene_variance
	elif version == 'B':
		standardized_gene_eqtl_effect_sizes = gene_eqtl_effect_sizes/np.sqrt(gene_variance)
		ld_scores = np.square(np.dot(ld[variant_indices,:], standardized_gene_eqtl_effect_sizes))

	# Get adjusted ld scores
	adj_ld_scores = ld_scores - ((1.0-ld_scores)/(n_ref_panel_samples-2.0))

	return ld_scores, adj_ld_scores

def extract_gene_regression_weights_for_this_gene_region_with_eqtl_point_estimate(ld, gene_eqtl_effect_sizes, variant_indices, n_ref_panel_samples):
	gene_variance = np.dot(np.dot(gene_eqtl_effect_sizes, ld), gene_eqtl_effect_sizes)
	standardized_gene_eqtl_effect_sizes = gene_eqtl_effect_sizes/np.sqrt(gene_variance)

	standardized_gene_eqtl_effect_sizes_shaped = np.reshape(standardized_gene_eqtl_effect_sizes, (len(standardized_gene_eqtl_effect_sizes), 1))

	beta_beta_t = np.dot(standardized_gene_eqtl_effect_sizes_shaped, np.transpose(standardized_gene_eqtl_effect_sizes_shaped))

	eqtl_cov = np.dot(np.dot(ld[variant_indices,:], beta_beta_t), ld[:,variant_indices])

	# Remove diagonal elements
	np.fill_diagonal(eqtl_cov,0.0)

	# Square and sum
	np.sum(np.square(eqtl_cov), axis=0)

	weights = np.sum(np.square(eqtl_cov), axis=1)

	return weights

def extract_ld_annotations_for_this_gene_region_with_sample_correlation(pmces, sample_geno, variant_indices, n_ref_panel_samples):
	stand_sample_geno = np.copy(sample_geno)
	for kk in range(stand_sample_geno.shape[1]):
		stand_sample_geno[:, kk] = (sample_geno[:, kk] - np.mean(sample_geno[:,kk]))/np.std(sample_geno[:,kk])

	pred_expr = np.dot(stand_sample_geno, pmces)

	full = np.hstack((np.reshape(pred_expr, (len(pred_expr), 1)), stand_sample_geno[:, variant_indices]))

	tmper = np.corrcoef(np.transpose(full))

	ld_scores = np.square(tmper[0,1:])
	adj_ld_scores = ld_scores - ((1.0-ld_scores)/(n_ref_panel_samples-2.0))
	return ld_scores, adj_ld_scores




def load_in_regression_snp_ids(variant_level_ld_score_file, rsid_to_snpid, snpid_to_reference_index):
	rsids = []
	snpids = []
	snp_to_regression_snp_position = {}
	f = gzip.open(variant_level_ld_score_file)
	head_count = 0
	snp_counter = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[1]
		# This filter occurs because of weird filter at the beginnong on removing repeat snps (only allow this to happen once per genome)
		snp_id = rsid_to_snpid[rsid]
		snp_id_info = snp_id.split('_')
		snp_id_alt = snp_id_info[0] + '_' + snp_id_info[1] + '_' + snp_id_info[3] + '_' + snp_id_info[2] + '_' + snp_id_info[4]
		# Quick error check
		if snp_id not in snpid_to_reference_index:
			print('assumption eroror')
			pdb.set_trace()
	
		# Add to arrays
		rsids.append(rsid)
		snpids.append(snp_id)
		if snp_id in snp_to_regression_snp_position or snp_id_alt in snp_to_regression_snp_position:
			print('assumption eroror')
			pdb.set_trace()
		snp_to_regression_snp_position[snp_id] = snp_counter
		snp_to_regression_snp_position[snp_id_alt] = snp_counter
		snp_counter = snp_counter + 1
	f.close()

	return np.asarray(rsids), np.asarray(snpids), snp_to_regression_snp_position

def load_in_regression_rsids(variant_level_ld_score_file):
	rsids = []
	rsid_to_regression_snp_position = {}
	f = gzip.open(variant_level_ld_score_file)
	head_count = 0
	snp_counter = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[1]
		# Add to arrays
		rsids.append(rsid)
		if rsid in rsid_to_regression_snp_position:
			print('assumption eroror')
			pdb.set_trace()
		rsid_to_regression_snp_position[rsid] = snp_counter
		snp_counter = snp_counter + 1
	f.close()

	return np.asarray(rsids), rsid_to_regression_snp_position



def get_non_repeat_columns_from_plink_obj(G_obj, variant_level_ld_score_file):
	f = gzip.open(variant_level_ld_score_file)
	head_count = 0
	hm3_rs_ids = {}
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[1]
		if rsid in hm3_rs_ids:
			print('assumption erorro')
			pdb.set_trace()
		hm3_rs_ids[rsid] = 1
	f.close()

	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	G_obj_rs = np.asarray(G_obj.snp)

	n_var = len(G_obj_pos)
	used = {}
	valid_columns = []

	# Go through hm3 snps first
	for ii in range(n_var):
		snp_id1 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a0[ii] + '_' + G_obj_a1[ii]
		snp_id2 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a1[ii] + '_' + G_obj_a0[ii]
		rsid = G_obj_rs[ii]
		if rsid not in hm3_rs_ids:
			continue
		if snp_id1 in used or snp_id2 in used:
			print('assumption error')
			pdb.set_trace()
		used[snp_id1] = 1
		used[snp_id2] = 1
		valid_columns.append(ii)
	for ii in range(n_var):
		snp_id1 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a0[ii] + '_' + G_obj_a1[ii]
		snp_id2 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a1[ii] + '_' + G_obj_a0[ii]
		rsid = G_obj_rs[ii]
		if rsid in hm3_rs_ids:
			continue
		if snp_id1 in used or snp_id2 in used:
			continue
		used[snp_id1] = 1
		used[snp_id2] = 1
		valid_columns.append(ii)
	valid_columns = np.sort(np.asarray(valid_columns))
	return valid_columns

def extract_protein_chrom_num(fit_file):
	f = open(fit_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		p_chrom_num = data[1]
		break
	f.close()
	return p_chrom_num


def extract_protein_to_protein_model_file_mapping(protein_models_dir, chrom_num):
	arr = []
	mapping = {}
	for file_name in os.listdir(protein_models_dir):

		fit_file = protein_models_dir + file_name + '/' + 'susie_fit.tsv'
		if os.path.isfile(fit_file) == False:
			continue

		protein_chrom_num = extract_protein_chrom_num(fit_file)
		if protein_chrom_num != chrom_num:
			continue

		arr.append(file_name)
		if file_name in mapping:
			print('assumpiotn oeroro')
			pdb.set_trace()
		mapping[file_name] = protein_models_dir + file_name + '/' + 'susie_fit.tsv'
	return np.unique(np.asarray(arr)), mapping

def get_gene_names(protein_models_dir, cis_snp_dir):
	gene_names = []
	miss_count = 0
	for file_name in os.listdir(protein_models_dir):
		if file_name.endswith('_sbrcout.txt') == False:
			continue
		gene_name = file_name.split('_sbrcout.t')[0]

		if os.path.isfile(cis_snp_dir + gene_name + '.txt'):
			gene_names.append(gene_name)
	return np.asarray(gene_names)

def extract_dictionary_of_cis_snps_for_this_gene(cis_snp_file):
	f = open(cis_snp_file)
	dicti = {}
	for line in f:
		line = line.rstrip()
		dicti[line] = 1
	f.close()

	return dicti

def load_in_pqtl_data(protein_models_dir, gene_names, geno_rsid_to_a1_a2, eqtl_class, cis_snp_dir):
	# Example gene name
	gene_name = gene_names[0]
	snp_list = []
	flip_values = []
	f = open(protein_models_dir + gene_name + '_sbrcout.txt')
	head_count =0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[0]
		a1 = data[1]
		if rsid not in geno_rsid_to_a1_a2:
			continue
		ref_a0, ref_a1 = geno_rsid_to_a1_a2[rsid]

		if ref_a1 == a1:
			flip = 1.0
		elif ref_a0 == a1:
			flip = -1.0
		else:
			print('allele mismatch')
			continue
		snp_list.append(rsid)
		flip_values.append(flip)
	f.close()
	snp_list = np.asarray(snp_list)
	flip_values = np.asarray(flip_values)

	if len(np.unique(snp_list)) != len(snp_list):
		print('assumption eroror')
		pdb.set_trace()

	# Create dictionary mapping from snp to array of length number of genes where element of array has snp-gene effect size
	# Initialize dictionaryu
	dicti = {}
	flips_dicti = {}
	for ii,snp in enumerate(snp_list):
		dicti[snp] = []
		flips_dicti[snp] = flip_values[ii]

	# Fill in dictionary for each gene
	for gene_name in gene_names:
		print(gene_name)
		
		# First extract dictionary of cis snps for this gene
		cis_snp_file = cis_snp_dir + gene_name + '.txt'
		cis_snps = extract_dictionary_of_cis_snps_for_this_gene(cis_snp_file)


		f = open(protein_models_dir + gene_name + '_sbrcout.txt')
		head_count =0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			rsid = data[0]
			if rsid not in dicti:
				continue

			cis_annotation_scaler = 1.0
			if eqtl_class == 'cis_pqtl' and rsid not in cis_snps:
				cis_annotation_scaler = 0.0
			elif eqtl_class == 'trans_pqtl' and rsid in cis_snps:
				cis_annotation_scaler = 0.0

			coef = cis_annotation_scaler*float(data[2])*flips_dicti[rsid]
			dicti[rsid].append(coef)
		f.close()

	for rsid in snp_list:
		dicti[rsid] = np.asarray(dicti[rsid])
		if len(dicti[rsid]) != len(gene_names):
			print('assumption eroror')
			pdb.set_trace()

	return dicti, snp_list









####################################################
# Command line args
####################################################
eqtl_class = sys.argv[1]
reference_genotype_stem = sys.argv[2]
protein_models_dir = sys.argv[3]
ldsc_baseline_ld_annotation_stem = sys.argv[4]
tglr_expression_score_dir = sys.argv[5]
method = sys.argv[6]
cis_snp_dir = sys.argv[7]

# Get names of genes/ proteins
gene_names = get_gene_names(protein_models_dir, cis_snp_dir)

# Get snp names
chrom_nums = np.arange(1,23)
geno_rsid_to_a1_a2 = {}
# Loop through chromosomes
for chrom_num in chrom_nums:
	genotype_stem = reference_genotype_stem + str(chrom_num)
	variant_level_ld_score_file = ldsc_baseline_ld_annotation_stem + str(chrom_num) + '.l2.ldscore.gz'

	# Load in Reference Genotype data
	G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)

	#G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	# For our purposes, a0 is the effect allele
	# For case of plink package, a0 is the first column in the plink bim file
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	# RSids
	G_obj_rsids = np.asarray(G_obj.snp)
	# Centimorgan distances
	G_obj_cm = np.asarray(G_obj.cm)
	n_samp = len(np.asarray(G_obj.iid))

	for ii, rsid in enumerate(G_obj_rsids):
		if rsid in geno_rsid_to_a1_a2:
			print('assumption erororor')
			pdb.set_trace()
		geno_rsid_to_a1_a2[rsid] = (G_obj_a0[ii], G_obj_a1[ii])

	
# Load in pqtl data (ie. load in causal pqtl effect sizes for variants)
print('load in PQTL data')
pqtl_effect_sizes, pqtl_variant_ids = load_in_pqtl_data(protein_models_dir, gene_names, geno_rsid_to_a1_a2, eqtl_class, cis_snp_dir)



###################
###################
# Nice, now get genetically predicted protein levels for each protein
print('Get genetically predicted protein levels')
genetic_protein = np.zeros((len(gene_names), n_samp))


for chrom_num in chrom_nums:
	print(chrom_num)
	genotype_stem = reference_genotype_stem + str(chrom_num)
	variant_level_ld_score_file = ldsc_baseline_ld_annotation_stem + str(chrom_num) + '.l2.ldscore.gz'

	# Load in Reference Genotype data
	G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)

	G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	# For our purposes, a0 is the effect allele
	# For case of plink package, a0 is the first column in the plink bim file
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	# RSids
	G_obj_rsids = np.asarray(G_obj.snp)
	# Centimorgan distances
	G_obj_cm = np.asarray(G_obj.cm)
	n_samp = len(np.asarray(G_obj.iid))

	# Put geno into organized dictionary
	genotype_obj = {'G': G_obj_geno, 'rsid': G_obj_rsids, 'position': G_obj_pos, 'cm': G_obj_cm, 'a0': G_obj_a0, 'a1': G_obj_a1}

	for ii, rsid in enumerate(G_obj_rsids):
		if rsid not in pqtl_effect_sizes:
			continue

		genotype = G_obj_geno[:,ii]
		variant_gene_effects = pqtl_effect_sizes[rsid]

		pred_genotype = np.dot(variant_gene_effects.reshape((len(variant_gene_effects),1)), genotype.reshape((1, len(genotype))))
		
		genetic_protein = genetic_protein + pred_genotype

# Standardize genetic protein
for gg in range(genetic_protein.shape[0]):
	genetic_protein[gg,:] = (genetic_protein[gg,:] - np.mean(genetic_protein[gg,:]))/np.std(genetic_protein[gg,:])



###################
###################
# Nice, now compute gene ld scores on each chromosome
print('Get genee ld scores')


for chrom_num in chrom_nums:
	print(chrom_num)
	genotype_stem = reference_genotype_stem + str(chrom_num)
	variant_level_ld_score_file = ldsc_baseline_ld_annotation_stem + str(chrom_num) + '.l2.ldscore.gz'

	# Load in Reference Genotype data
	G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)

	G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	# For our purposes, a0 is the effect allele
	# For case of plink package, a0 is the first column in the plink bim file
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	# RSids
	G_obj_rsids = np.asarray(G_obj.snp)
	# Centimorgan distances
	G_obj_cm = np.asarray(G_obj.cm)
	n_samp = len(np.asarray(G_obj.iid))

	# Put geno into organized dictionary
	genotype_obj = {'G': G_obj_geno, 'rsid': G_obj_rsids, 'position': G_obj_pos, 'cm': G_obj_cm, 'a0': G_obj_a0, 'a1': G_obj_a1}


	# Create mapping from rsid to reference index
	rsid_to_reference_index = create_mapping_rsid_to_reference_index(genotype_obj['rsid'])
	
	# Get regression rsids
	regression_rsids, regression_rsid_to_regression_snp_position = load_in_regression_rsids(variant_level_ld_score_file)
	n_regression_snps = len(regression_rsids)

	# Initialize vector to keep track of gene level ld scores for regression snps on this chromosome
	gene_level_ld_scores = np.zeros(n_regression_snps)


	for regression_index, regression_rsid in enumerate(regression_rsids):
		# Quick error check
		if regression_index != regression_rsid_to_regression_snp_position[regression_rsid]:
			print('assumption eroror')
			pdb.set_trace()


		ref_index = rsid_to_reference_index[regression_rsid]
		genotype = G_obj_geno[:,ref_index]
		std_genotype = (genotype-np.mean(genotype))/np.std(genotype)
		'''
		t1 = time.time()
		for gene_index in range(genetic_protein.shape[0]):
			corry = np.corrcoef(genotype, genetic_protein[gene_index,:])[0,1]
			squared_corry = np.square(corry)
			squared_corry_adj = squared_corry - ((1.0-squared_corry)/(n_samp-2.0))

			gene_level_ld_scores[regression_index] = gene_level_ld_scores[regression_index] + squared_corry_adj
		'''
		squared_corrz = np.square(np.dot(genetic_protein,std_genotype)/len(std_genotype))
		squared_corrz_adj = squared_corrz - ((1.0-squared_corrz)/(n_samp-2.0))
		tmp_gene_ldscore = np.sum(squared_corrz_adj)
		gene_level_ld_scores[regression_index] = gene_level_ld_scores[regression_index] + tmp_gene_ldscore

	# Save gene level ld scores to output to output file
	output_file = tglr_expression_score_dir + eqtl_class + '_' + method + '_' + str(chrom_num) + '_' + 'tglr_gene_ld_scores.txt'
	np.savetxt(output_file, gene_level_ld_scores, fmt="%s")

	# Save number of genes to output also
	n_genes = genetic_protein.shape[0]
	output_file2 = tglr_expression_score_dir + eqtl_class + '_' + method + '_' + str(chrom_num) + '_' + 'tglr_n_genes.txt'
	t = open(output_file2,'w')
	t.write(str(n_genes) + '\n')
	t.close()


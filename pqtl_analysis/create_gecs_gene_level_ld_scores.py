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


def extract_causal_eqtl_effects_for_single_gene(gene_info_file, ref_genotype_obj, rsid_to_reference_index, regression_rsid_to_regression_snp_position, method):
	gene_rsids = []
	gene_eqtl_effect_sizes = []
	gene_ref_positions = []
	f = open(gene_info_file)
	head_count = 0
	misses = 0
	total = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header = np.copy(data)
			continue
		if len(data) != 51 and len(data) != 52:
			print('assumptinoer orro')
			pdb.set_trace()

		rsid = data[0]
		snp_a0 = data[3]
		snp_a1 = data[4]
		if method == 'susie_inf_pmces':
			effect_size = float(data[46])
		ref_index = rsid_to_reference_index[rsid]
		total = total + 1
		sign_flip = 1.0
		if ref_genotype_obj['a0'][ref_index] == snp_a0 and ref_genotype_obj['a1'][ref_index] == snp_a1:
			sign_flip = 1.0
		elif ref_genotype_obj['a1'][ref_index] == snp_a0 and ref_genotype_obj['a0'][ref_index] == snp_a1:
			sign_flip = -1.0
		else:
			misses = misses + 1
			continue
		gene_rsids.append(rsid)
		gene_eqtl_effect_sizes.append(effect_size*sign_flip)
		gene_ref_positions.append(ref_index)

	print(str(misses) + ' : ' + str(total))
	gene_rsids = np.asarray(gene_rsids)
	gene_eqtl_effect_sizes = np.asarray(gene_eqtl_effect_sizes)
	gene_ref_positions = np.asarray(gene_ref_positions)

	if np.array_equal(gene_ref_positions, np.sort(gene_ref_positions)) == False:
		print('assumption erororor')
		pdb.set_trace()


	n_snps = len(ref_genotype_obj['rsid'])
	gene_snp_boolean_arr = np.asarray([False]*n_snps)
	gene_snp_boolean_arr[gene_ref_positions] = True


	# Get pred expression
	# First need to standardize expression
	eqtl_geno = ref_genotype_obj['G'][:, gene_ref_positions]
	for snp_iter in range(eqtl_geno.shape[1]):
		eqtl_geno[:, snp_iter] = (eqtl_geno[:, snp_iter] - np.mean(eqtl_geno[:, snp_iter]))/np.std(eqtl_geno[:, snp_iter])
	# Now compute genetically predicted expression
	pred_expr = np.dot(eqtl_geno, gene_eqtl_effect_sizes)

	standardized_gene_eqtl_effect_sizes = gene_eqtl_effect_sizes/np.sqrt(np.var(pred_expr))

	return standardized_gene_eqtl_effect_sizes, gene_snp_boolean_arr


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

		#fit_file = protein_models_dir + file_name + '/' + 'susie_fit.tsv'
		fit_file = protein_models_dir + file_name
		if fit_file.endswith('.susie_fit.tsv') == False:
			continue

		protein_chrom_num = extract_protein_chrom_num(fit_file)
		if protein_chrom_num != chrom_num:
			continue

		arr.append(file_name)
		if file_name in mapping:
			print('assumpiotn oeroro')
			pdb.set_trace()
		#mapping[file_name] = protein_models_dir + file_name + '/' + 'susie_fit.tsv'
		mapping[file_name] = fit_file
	return np.unique(np.asarray(arr)), mapping

def create_mapping_from_rsid_to_annovec(variant_level_annotation_file):
	f = open(variant_level_annotation_file)
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header_vec = np.asarray(data[5:])
			continue
		if len(data) != 97:
			print('assumption erororo')
			pdb.set_trace()
		rsid = data[2]
		anno_vec = np.asarray(data[5:]).astype(float)
		if rsid in mapping:
			print('repeat variant assumption eroror')
			pdb.set_trace()
		mapping[rsid] = anno_vec
	return mapping, header_vec

def create_snp_annotation_matrix(rsid_to_anno_vec, G_obj_rsids):
	anno_mat = []
	for rsid in G_obj_rsids:
		anno_mat.append(rsid_to_anno_vec[rsid])
	return np.asarray(anno_mat)

def check_if_variant_pair_in_single_cis_gene(snp_pos1, snp_pos2, gene_windows):
	per_gene_boolers = (snp_pos1 >= gene_windows[:,0]) & (snp_pos1 <= gene_windows[:,1]) & (snp_pos2 >= gene_windows[:,0]) & (snp_pos2 <= gene_windows[:,1])
	total_bool = np.sum(per_gene_boolers)
	final_bool = 'False'
	if total_bool > 0:
		final_bool = 'True'
	return final_bool


####################################################
# Command line args
####################################################
chrom_num = sys.argv[1]
tissue_name = sys.argv[2]
reference_genotype_stem = sys.argv[3]
protein_models_dir = sys.argv[4]
ldsc_baseline_ld_annotation_stem = sys.argv[5]
gecs_expression_score_dir = sys.argv[6]
method = sys.argv[7]



print('CHROM' + str(chrom_num))

genotype_stem = reference_genotype_stem + str(chrom_num)
variant_level_ld_score_file = ldsc_baseline_ld_annotation_stem + str(chrom_num) + '.l2.ldscore.gz'
variant_level_annotation_file = ldsc_baseline_ld_annotation_stem + str(chrom_num) + '.annot'


# Load in Reference Genotype data
G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)
# For some reason, this file has some repeats. First get columns that don't correspond to repeats
#valid_variant_columns = get_non_repeat_columns_from_plink_obj(G_obj, variant_level_ld_score_file)
#G_obj = G_obj[:, valid_variant_columns]

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

# Put geno into organized dictionary
genotype_obj = {'G': G_obj_geno, 'rsid': G_obj_rsids, 'position': G_obj_pos, 'cm': G_obj_cm, 'a0': G_obj_a0, 'a1': G_obj_a1}

# Create rsid to ref alleles mapping
rsid_to_alleles = {}
rsid_to_pos = {}
for ii, rsid in enumerate(G_obj_rsids):
	if rsid in rsid_to_alleles:
		print('assumption erororo')
		pdb.set_trace()
	rsid_to_alleles[rsid] = (G_obj_a0[ii], G_obj_a1[ii])
	rsid_to_pos[rsid] = G_obj_pos[ii]

# Create mapping from rsid to annotion vector
rsid_to_anno_vec, anno_header = create_mapping_from_rsid_to_annovec(variant_level_annotation_file)

# Create snp matrix of annotations
anno_mat = create_snp_annotation_matrix(rsid_to_anno_vec, G_obj_rsids)

# Add intercept to anno_vec
anno_mat = np.hstack((np.ones((anno_mat.shape[0],1)), anno_mat))
anno_header = np.hstack((['genotype_intercept'], anno_header))

# Create mapping from rsid to reference index
rsid_to_reference_index = create_mapping_rsid_to_reference_index(genotype_obj['rsid'])


# Get regression rsids
regression_rsids, regression_rsid_to_regression_snp_position = load_in_regression_rsids(variant_level_ld_score_file)

# Fill in regression snps
regression_snp_boolean = []
for rsid in G_obj_rsids:
	if rsid in regression_rsid_to_regression_snp_position:
		regression_snp_boolean.append(True)
	else:
		regression_snp_boolean.append(False)
regression_snp_boolean = np.asarray(regression_snp_boolean)
genotype_obj['regression_snp_boolean'] = regression_snp_boolean


# Get snp index corresponding to first regression snp
window_start = np.where(genotype_obj['position'] == np.min(genotype_obj['position'][regression_snp_boolean]))[0]
if len(window_start) != 1:
	print('assumption errror')
active_window_start_position = genotype_obj['position'][window_start[0]] -1 


# Create dictionary of genes (proteins)
# Also filter to proteins on this chromosome
protein_arr, protein_to_model_file = extract_protein_to_protein_model_file_mapping(protein_models_dir, chrom_num)

# Loop through genes
n_genes = len(protein_arr)
gene_arr = []
gene_windows = []
for g_counter, gene_name in enumerate(protein_arr):
	# Load in causal eqtls for this gene
	gene_info_file = protein_to_model_file[gene_name]
	# Extract causal eQTL effect sizes for this gene
	gene_pmces, gene_snp_indices = extract_causal_eqtl_effects_for_single_gene(gene_info_file, genotype_obj, rsid_to_reference_index, regression_rsid_to_regression_snp_position, method)
	gene_start = np.min(G_obj_pos[gene_snp_indices])
	gene_end = np.max(G_obj_pos[gene_snp_indices])
	gene_arr.append((gene_name, gene_pmces, gene_snp_indices))
	gene_windows.append([gene_start, gene_end])
gene_windows = np.asarray(gene_windows)


# Initialize list of data to save
snp_pair_names = []
pairwise_ld_scores = []
pairwise_anno_ld_scores = []
pairwise_eqtl_ld_scores = []
pairwise_ld = []



# Size of active windows
active_window_size = 1000000
interaction_window_size = 2000000

counter = 0
# Loop through genomic windows
while active_window_start_position < np.max(genotype_obj['position']):

	print(active_window_start_position)
	# Extract active window snps
	active_window_end_position = active_window_start_position + active_window_size
	active_window_snp_indices = (G_obj_pos >= active_window_start_position) & (G_obj_pos < active_window_end_position)

	# Extract interaction window snps
	interaction_window_start_position = active_window_end_position + 0
	interaction_window_end_position = interaction_window_start_position + interaction_window_size
	interaction_window_snp_indices = (G_obj_pos >= interaction_window_start_position) & (G_obj_pos < interaction_window_end_position)

	# Joint active plus interaction window
	active_plus_interaction_window_snp_indices = active_window_snp_indices + interaction_window_snp_indices
	
	# Skip windows without active window snp indices
	if len(active_plus_interaction_window_snp_indices) == 0 or len(active_window_snp_indices) == 0:
		active_window_start_position = active_window_start_position + active_window_size
		continue
	if np.sum(active_plus_interaction_window_snp_indices) == 0:
		active_window_start_position = active_window_start_position + active_window_size
		continue	


	# Extract left flanking window snps
	left_flanking_window_end_cm_value = np.min(G_obj_cm[active_plus_interaction_window_snp_indices])
	left_flanking_window_start_cm_value = left_flanking_window_end_cm_value - 1.0
	left_flanking_window_snp_indices = (G_obj_cm > left_flanking_window_start_cm_value) & (G_obj_cm < left_flanking_window_end_cm_value)
	
	# Extract right flanking window snps
	right_flanking_window_start_cm_value = np.max(G_obj_cm[active_plus_interaction_window_snp_indices])
	right_flanking_window_end_cm_value = right_flanking_window_start_cm_value + 1.0
	right_flanking_window_snp_indices = (G_obj_cm > right_flanking_window_start_cm_value) & (G_obj_cm < right_flanking_window_end_cm_value)

	# Quick error checking
	if np.sum(1.0*active_window_snp_indices + 1.0*interaction_window_snp_indices + 1.0*left_flanking_window_snp_indices+ 1.0*right_flanking_window_snp_indices > 1.0) != 0.0:
		print('assumption error')
		pdb.set_trace()

	counter = counter + 1

	# Get window snp indices
	window_snp_indices = left_flanking_window_snp_indices + active_window_snp_indices + interaction_window_snp_indices + right_flanking_window_snp_indices

	# Extract genotype of window
	window_geno = G_obj_geno[:, window_snp_indices]
	window_rsid = G_obj_rsids[window_snp_indices]
	window_snp_cm = G_obj_cm[window_snp_indices]
	window_snp_pos = G_obj_pos[window_snp_indices]
	window_snp_a0 = G_obj_a0[window_snp_indices]
	window_snp_a1 = G_obj_a1[window_snp_indices]
	window_regression_snps = regression_snp_boolean[window_snp_indices]
	window_anno_mat = anno_mat[window_snp_indices,:]

	# get window subsets
	window_active_indices = active_window_snp_indices[window_snp_indices]
	window_interaction_indices = interaction_window_snp_indices[window_snp_indices]
	window_active_regression_indices = window_active_indices*window_regression_snps
	window_active_plus_interaction_regression_indices = (window_active_indices + window_interaction_indices)*window_regression_snps


	# Skip window cause no interaction
	if np.sum(window_active_regression_indices) == 0 or np.sum(window_active_plus_interaction_regression_indices) <= 1:
		active_window_start_position = active_window_start_position + active_window_size
		continue

	# Get window LD
	window_LD = np.corrcoef(np.transpose(window_geno))

	# Create matrix of dimension n_snpsXn_snps where each element contains the CM difference between snp pairs
	window_cm_diff_mat = np.abs(window_snp_cm[:, None] - window_snp_cm)
	# Use matrix to set LD between snp-pairs with CM > 1 to 0.0
	window_LD[window_cm_diff_mat > 1.0] = 0.0
	# Delete window_cm_diff_mat for memory purposes

	# Compute ld_t_ld for pairwise combo of window_active_regression_indices and window_active_plus_interaction_regression_indices
	window_ld_t_ld = np.dot(window_LD[window_active_regression_indices, :], window_LD[:, window_active_plus_interaction_regression_indices])
	# Compute ld for subset
	window_LD_subset = window_LD[window_active_regression_indices, :][:, window_active_plus_interaction_regression_indices]
	# Compute pairwise distance between snps
	window_snp_pair_pos_diff_subset = np.abs(window_snp_pos[:, None] - window_snp_pos)[window_active_regression_indices,:][:, window_active_plus_interaction_regression_indices]
	window_snp_pair_cm_diff_subset = window_cm_diff_mat[window_active_regression_indices,:][:, window_active_plus_interaction_regression_indices]


	window_eqtl_eqtl_LD_subset = np.copy(window_LD_subset)*0.0
	for gene_tuple in gene_arr:
		gene_pmces = gene_tuple[1]
		gene_snp_indices = gene_tuple[2]
		# Skip genes with no overlap with the window
		if np.sum(np.sum(gene_snp_indices[window_snp_indices])) == 0:
			continue
		n_snps = len(G_obj_rsids)
		causal_eqtl_effects = np.zeros(n_snps)
		causal_eqtl_effects[gene_snp_indices] = gene_pmces
		window_causal_eqtl_effects = causal_eqtl_effects[window_snp_indices]
		marginal_term_a = np.dot(window_LD[window_active_regression_indices, :], window_causal_eqtl_effects)
		marginal_term_b = np.dot(window_causal_eqtl_effects, window_LD[:, window_active_plus_interaction_regression_indices])
		window_eqtl_eqtl_LD_subset = window_eqtl_eqtl_LD_subset + np.dot(marginal_term_a.reshape(-1,1), marginal_term_b.reshape(1,-1))


	# Compute matrices into vectors
	# Get upper triangle indices 
	upper_triangle_indices = np.triu_indices(window_ld_t_ld.shape[0], k=1, m=window_ld_t_ld.shape[1])
	# Use to conver to vectors
	window_ld_t_ld_vec = window_ld_t_ld[upper_triangle_indices]
	window_ld_vec = window_LD_subset[upper_triangle_indices]
	window_snp_pair_pos_diff_vec = window_snp_pair_pos_diff_subset[upper_triangle_indices]
	window_snp_pair_cm_diff_vec = window_snp_pair_cm_diff_subset[upper_triangle_indices]
	window_eqtl_ld_vec = window_eqtl_eqtl_LD_subset[upper_triangle_indices]
	window_snp_pairs = window_rsid[window_active_regression_indices][upper_triangle_indices[0]] + ':' + window_rsid[window_active_plus_interaction_regression_indices][upper_triangle_indices[1]]
	
	# Filter
	filtered_snp_pairs = (np.abs(window_snp_pair_pos_diff_vec) < interaction_window_size) & (np.abs(window_snp_pair_cm_diff_vec) < 1.0) & (np.abs(window_ld_t_ld_vec) < 1.0) & (np.abs(window_ld_vec) < .1)

	snp_pair_names.append(window_snp_pairs[filtered_snp_pairs])
	pairwise_ld_scores.append(window_ld_t_ld_vec[filtered_snp_pairs])
	pairwise_ld.append(window_ld_vec[filtered_snp_pairs])
	pairwise_eqtl_ld_scores.append(window_eqtl_ld_vec[filtered_snp_pairs])

	# Run analysis for annotation weighted ld scores
	tmp_pairwise_anno_ld_scores = []
	for anno_iter in range(len(anno_header)):
		anno_ld_ld_t = np.dot(np.transpose(np.transpose((window_LD[window_active_regression_indices, :]))*((anno_mat[window_snp_indices,anno_iter])[:, np.newaxis])), window_LD[:, window_active_plus_interaction_regression_indices])
		tmp_pairwise_anno_ld_scores.append(anno_ld_ld_t[upper_triangle_indices][filtered_snp_pairs])
	tmp_pairwise_anno_ld_scores = np.transpose(np.asarray(tmp_pairwise_anno_ld_scores))
	pairwise_anno_ld_scores.append(tmp_pairwise_anno_ld_scores)


	# Start next iter
	active_window_start_position = active_window_start_position + active_window_size

	'''
	if counter > 5:
		break
	'''

snp_pair_names = np.hstack(snp_pair_names)
pairwise_ld_scores = np.hstack(pairwise_ld_scores)
pairwise_ld = np.hstack(pairwise_ld)
pairwise_eqtl_ld_scores = np.hstack(pairwise_eqtl_ld_scores)
pairwise_anno_ld_scores = np.vstack(pairwise_anno_ld_scores)

# Save gene level ld scores to output to output file
output_file = gecs_expression_score_dir + tissue_name + '_' + method + '_' + chrom_num + '_' + 'gecs_scores.txt'
t = open(output_file,'w')
t.write('snp_pair\tsnp1\tsnp2\tsnp1_A1\tsnp1_A2\tsnp2_A1\tsnp2_A2\tsnp_pair_in_cis_gene\tLD\teqtl_LD_score\t' + '\t'.join(anno_header) + '\n')


for ii, snp_pair_name in enumerate(snp_pair_names):
	snp_name_info = snp_pair_name.split(':')
	if len(snp_name_info) != 2:
		print('assumption eororor')
		pdb.set_trace()
	snp1 = snp_name_info[0]
	snp2 = snp_name_info[1]
	allele_string = '\t'.join(rsid_to_alleles[snp1]) + '\t' + '\t'.join(rsid_to_alleles[snp2])

	variant_pair_in_single_cis_gene = check_if_variant_pair_in_single_cis_gene(rsid_to_pos[snp1], rsid_to_pos[snp2], gene_windows)

	t.write(snp_pair_name + '\t' + snp1 + '\t' + snp2 + '\t' + allele_string + '\t' + variant_pair_in_single_cis_gene + '\t' + str(pairwise_ld[ii]) + '\t' + str(pairwise_eqtl_ld_scores[ii]) + '\t' + '\t'.join(pairwise_anno_ld_scores[ii,:].astype(str)) + '\n')
t.close()
print(output_file)

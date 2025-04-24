import sys
import os
import pdb
import numpy as np
from pandas_plink import read_plink1_bin
import pickle
import pandas as pd
import gzip
import time

def get_all_tissue_names(causal_eqtl_summary_file):
	dicti = {}
	f = open(causal_eqtl_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		dicti[data[0]] = 1
	f.close()
	return np.asarray([*dicti])


def create_distance_stratefication_gene_model_arr(causal_eqtl_summary_file, chrom_num, plink_reference_bim_file, gene_annotation_file, all_distances, n_bins=10, max_min_ratio_thresh=100):
	# Create mapping from gene name to gene class
	gene_name_to_gene_info = create_mapping_from_gene_name_to_gene_info(gene_annotation_file)

	splits = np.array_split(all_distances, n_bins)
	category_names = []
	category_ranges = []
	for ii, grouper in enumerate(splits):
		category_names.append('eqtl_distance_group_' + str(ii))
		category_ranges.append((np.min(grouper), np.max(grouper)))
	category_names = np.asarray(category_names)
	category_counts = np.zeros(len(category_names))

	# Create mapping from rsid to reference and alternate alleles
	f = open(plink_reference_bim_file)
	rsid_mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != chrom_num:
			print('assumption error')
			pdb.set_trace()
		rsid1 = data[1] + ':' + data[4] + '_' + data[5]
		rsid2 = data[1] + ':' + data[5] + '_' + data[4]
		if rsid1 in rsid_mapping or rsid2 in rsid_mapping:
			print('assumption eororro')
			pdb.set_trace()
		rsid_mapping[rsid1] = 1.0
		rsid_mapping[rsid2] = -1.0
	f.close()

	# Create gene array
	gene_arr = []
	f = open(causal_eqtl_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Skip genes not on the chromosome of interest
		line_chrom_num = data[2]
		if line_chrom_num != chrom_num:
			continue
		# Skip genes that dont pass boolean
		gene_boolean = data[4]
		if gene_boolean != 'Pass':
			continue

		# Get gene TSS
		ensamble_id = data[1]
		gene_info = gene_name_to_gene_info[ensamble_id]
		gene_tss = int(gene_info[3])

		# Extract relevent fields
		gene_id = data[1] + ':' + data[0]
		variant_info_file = data[5]
		tissue_name = data[0]

		susie_alpha_file = data[6]
		susie_mu_file = data[7]

		susie_alpha = np.load(susie_alpha_file)
		susie_mu = np.load(susie_mu_file)

		# Loop through components
		for kk in range(susie_alpha.shape[0]):

			if np.min(susie_alpha[kk,:]) == 0.0:
				max_min_ratio = np.max(susie_alpha[kk,:])/(1e-20)
			else:
				max_min_ratio = np.max(susie_alpha[kk,:])/np.min(susie_alpha[kk,:])
			if max_min_ratio < max_min_ratio_thresh:
				continue

			component_pmces = susie_alpha[kk,:]*susie_mu[kk,:]
			component_gene_id = gene_id + ':' + str(kk)
			g = open(variant_info_file)
			gene_rsids = []
			snp_gene_dist = []
			var_counter = 0
			for line2 in g:
				line2 = line2.rstrip()
				data2 = line2.split('\t')
				snp_pos = float(data2[3])
				snp_gene_dist.append(np.abs(snp_pos- gene_tss))
				full_rsid = data2[1] + ':' + data2[4] + '_' + data2[5]
				gene_rsids.append(data2[1])
				component_pmces[var_counter] = component_pmces[var_counter]*rsid_mapping[full_rsid]
				var_counter = var_counter + 1
			g.close()

			if var_counter != len(component_pmces):
				print('assumption eroror')
				pdb.set_trace()

			# Compute expected component distance
			snp_gene_dist = np.asarray(snp_gene_dist)
			component_expected_snp_gene_distance = np.sum(snp_gene_dist*susie_alpha[kk,:])
			
			# Get category number
			category_number = None
			for ii, ranger in enumerate(category_ranges):
				if component_expected_snp_gene_distance >= ranger[0] and component_expected_snp_gene_distance <= ranger[1]:
					if category_number is not None:
						print('assumption eroror')
						pdb.set_trace()
					category_number = ii
			if category_number is None:
				print('assumption eroror')
				pdb.set_trace()

			# Add to global gene array
			category_counts[category_number] = category_counts[category_number] + 1
			gene_arr.append((component_gene_id,category_number, np.asarray(gene_rsids), component_pmces))
	f.close()
	return gene_arr, category_names, category_counts



def create_qtl_rank_stratefication_gene_model_arr(causal_eqtl_summary_file, chrom_num, plink_reference_bim_file, max_min_ratio_thresh=100):
	category_names = []
	category_names.append('eqtl_rank_1')
	category_names.append('eqtl_rank_2')
	category_names.append('eqtl_rank_3')
	category_names = np.asarray(category_names)
	category_counts = np.zeros(len(category_names))

	# Create mapping from rsid to reference and alternate alleles
	f = open(plink_reference_bim_file)
	rsid_mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != chrom_num:
			print('assumption error')
			pdb.set_trace()
		rsid1 = data[1] + ':' + data[4] + '_' + data[5]
		rsid2 = data[1] + ':' + data[5] + '_' + data[4]
		if rsid1 in rsid_mapping or rsid2 in rsid_mapping:
			print('assumption eororro')
			pdb.set_trace()
		rsid_mapping[rsid1] = 1.0
		rsid_mapping[rsid2] = -1.0
	f.close()

	# Create gene array
	gene_arr = []
	f = open(causal_eqtl_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Skip genes not on the chromosome of interest
		line_chrom_num = data[2]
		if line_chrom_num != chrom_num:
			continue
		# Skip genes that dont pass boolean
		gene_boolean = data[4]
		if gene_boolean != 'Pass':
			continue

		# Get gene TSS
		ensamble_id = data[1]

		# Extract relevent fields
		gene_id = data[1] + ':' + data[0]
		variant_info_file = data[5]
		tissue_name = data[0]

		susie_alpha_file = data[6]
		susie_mu_file = data[7]

		susie_alpha = np.load(susie_alpha_file)
		susie_mu = np.load(susie_mu_file)

		# Loop through components
		valid_components = []
		valid_component_effect_magnitudes = []
		for kk in range(susie_alpha.shape[0]):

			if np.min(susie_alpha[kk,:]) == 0.0:
				max_min_ratio = np.max(susie_alpha[kk,:])/(1e-20)
			else:
				max_min_ratio = np.max(susie_alpha[kk,:])/np.min(susie_alpha[kk,:])
			if max_min_ratio < max_min_ratio_thresh:
				continue
			valid_components.append(kk)
			valid_component_effect_magnitudes.append(np.abs(susie_mu[kk, np.argmax(susie_alpha[kk,:])]))
		valid_components = np.asarray(valid_components)
		valid_component_effect_magnitudes = np.asarray(valid_component_effect_magnitudes)

		if len(valid_components) == 0:
			continue

		for ii, kk in enumerate(valid_components[np.argsort(-valid_component_effect_magnitudes)]):
			component_pmces = susie_alpha[kk,:]*susie_mu[kk,:]
			component_gene_id = gene_id + ':' + str(kk)
			g = open(variant_info_file)
			gene_rsids = []
			var_counter = 0
			for line2 in g:
				line2 = line2.rstrip()
				data2 = line2.split('\t')
				snp_pos = float(data2[3])
				full_rsid = data2[1] + ':' + data2[4] + '_' + data2[5]
				gene_rsids.append(data2[1])
				component_pmces[var_counter] = component_pmces[var_counter]*rsid_mapping[full_rsid]
				var_counter = var_counter + 1
			g.close()

			if var_counter != len(component_pmces):
				print('assumption eroror')
				pdb.set_trace()
			
			# Get category number
			category_number = None

			if ii == 0:
				category_number = 0
			elif ii == 1:
				category_number = 1
			elif ii >= 2:
				category_number = 2

			# Add to global gene array
			category_counts[category_number] = category_counts[category_number] + 1
			gene_arr.append((component_gene_id,category_number, np.asarray(gene_rsids), component_pmces))
	f.close()
	return gene_arr, category_names, category_counts





def create_tissue_stratefication_gene_model_arr(causal_eqtl_summary_file, chrom_num, plink_reference_bim_file):
	# Create mapping from tissue names to category numbers
	tissue_names = np.sort(get_all_tissue_names(causal_eqtl_summary_file))
	category_names = []
	tissue_to_category_number_mapping = {}
	for tissue_num, tissue_name in enumerate(tissue_names):
		if tissue_name in tissue_to_category_number_mapping:
			print('assumption eroror')
			pdb.set_trace()
		tissue_to_category_number_mapping[tissue_name] = tissue_num
		category_names.append('eQTL_' + tissue_name)
	category_names = np.asarray(category_names)
	category_counts = np.zeros(len(category_names))

	# Create mapping from rsid to reference and alternate alleles
	f = open(plink_reference_bim_file)
	rsid_mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != chrom_num:
			print('assumption error')
			pdb.set_trace()
		rsid1 = data[1] + ':' + data[4] + '_' + data[5]
		rsid2 = data[1] + ':' + data[5] + '_' + data[4]
		if rsid1 in rsid_mapping or rsid2 in rsid_mapping:
			print('assumption eororro')
			pdb.set_trace()
		rsid_mapping[rsid1] = 1.0
		rsid_mapping[rsid2] = -1.0
	f.close()
	
	# Create gene array
	gene_arr = []
	f = open(causal_eqtl_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Skip genes not on the chromosome of interest
		line_chrom_num = data[2]
		if line_chrom_num != chrom_num:
			continue
		# Skip genes that dont pass boolean
		gene_boolean = data[4]
		if gene_boolean != 'Pass':
			continue

		# Extract relevent fields
		gene_id = data[1] + ':' + data[0]
		variant_info_file = data[5]
		susie_pmces_file = data[9]
		tissue_name = data[0]

		susie_pmces = np.load(susie_pmces_file)
		init_susie_pmces = np.copy(susie_pmces)
		g = open(variant_info_file)
		gene_rsids = []
		var_counter = 0
		for line2 in g:
			line2 = line2.rstrip()
			data2 = line2.split('\t')
			full_rsid = data2[1] + ':' + data2[4] + '_' + data2[5]
			gene_rsids.append(data2[1])
			susie_pmces[var_counter] = susie_pmces[var_counter]*rsid_mapping[full_rsid]
			var_counter = var_counter + 1
		g.close()
		gene_rsids = np.asarray(gene_rsids)
		if var_counter != len(susie_pmces):
			print('assumption eroror')
			pdb.set_trace()
		category_number = tissue_to_category_number_mapping[tissue_name]
		gene_arr.append((gene_id,category_number, gene_rsids, susie_pmces))
		category_counts[category_number] = category_counts[category_number] + 1
	f.close()
	return gene_arr, category_names, category_counts

def create_mapping_from_gene_name_to_gene_info(gene_annotation_file):
	f = open(gene_annotation_file)
	mapping = {}
	for line in f:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		data = line.split('\t')
		if len(data) != 9:
			print('assumption eroror')
			pdb.set_trace()
		if data[2] != 'gene':
			continue
		ensamble_id = 'null'
		gene_type = 'null'
		gene_info = data[8].split(';')
		for info in gene_info:
			if info.startswith('gene_id'):
				ensamble_id = info.split('"')[1]
			elif info.startswith(' gene_type'):
				gene_type = info.split('"')[1]
		if ensamble_id == 'null' or gene_type == 'null':
			print('assumption eroror')
			pdb.set_trace()
		gene_chrom_num = data[0]
		gene_strand = data[6]
		if float(data[3]) > float(data[4]):
			print('assumption erroror')
			pdb.set_trace()
		if gene_strand == '+':
			tss = data[3]
		elif gene_strand == '-':
			tss = data[4]
		else:
			print('assumption error')


		# Add to info
		if ensamble_id not in mapping:
			mapping[ensamble_id] = (gene_type, gene_chrom_num, gene_strand, tss)
		else:
			if mapping[ensamble_id][0] != gene_type:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][1] != gene_chrom_num:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][2] != gene_strand:
				print('assumption eroror')
				pdb.set_trace()
			if mapping[ensamble_id][3] != tss:
				print('assumption eroror')
				pdb.set_trace()
	f.close()
	return mapping

def extract_gene_snp_distances_all_chrom(causal_eqtl_summary_file, gene_annotation_file, max_min_ratio_thresh=100):
	# Create mapping from gene name to gene class
	gene_name_to_gene_info = create_mapping_from_gene_name_to_gene_info(gene_annotation_file)

	
	# Create gene array
	gene_arr = []
	f = open(causal_eqtl_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Skip genes that dont pass boolean
		gene_boolean = data[4]
		if gene_boolean != 'Pass':
			continue

		# Get gene TSS
		ensamble_id = data[1]
		gene_info = gene_name_to_gene_info[ensamble_id]

		gene_tss = int(gene_info[3])

		# Extract relevent fields
		gene_id = data[1] + ':' + data[0]
		variant_info_file = data[5]
		susie_pmces_file = data[9]
		tissue_name = data[0]
		susie_alpha_file = data[6]
		susie_mu_file = data[7]

		susie_alpha = np.load(susie_alpha_file)


		# Loop through components
		for kk in range(susie_alpha.shape[0]):
			if np.min(susie_alpha[kk,:]) == 0.0:
				max_min_ratio = np.max(susie_alpha[kk,:])/(1e-20)
			else:
				max_min_ratio = np.max(susie_alpha[kk,:])/np.min(susie_alpha[kk,:])
			if max_min_ratio < max_min_ratio_thresh:
				continue

			component_gene_id = gene_id + ':' + str(kk)
			g = open(variant_info_file)
			gene_rsids = []
			snp_gene_dist = []
			var_counter = 0
			for line2 in g:
				line2 = line2.rstrip()
				data2 = line2.split('\t')
				snp_pos = float(data2[3])
				snp_gene_dist.append(np.abs(snp_pos- gene_tss))
				var_counter = var_counter + 1
			g.close()

			# Compute expected component distance
			snp_gene_dist = np.asarray(snp_gene_dist)
			component_expected_snp_gene_distance = np.sum(snp_gene_dist*susie_alpha[kk,:])
			# Add to global gene array
			gene_arr.append((component_gene_id,component_expected_snp_gene_distance))
	f.close()
	return gene_arr





def create_tissue_aggregation_gene_model_arr(causal_eqtl_summary_file, chrom_num, plink_reference_bim_file):
	# Create mapping from tissue names to category numbers
	tissue_names = np.sort(get_all_tissue_names(causal_eqtl_summary_file))
	category_names = []
	tissue_to_category_number_mapping = {}
	for tissue_num, tissue_name in enumerate(tissue_names):
		if tissue_name in tissue_to_category_number_mapping:
			print('assumption eroror')
			pdb.set_trace()
		tissue_to_category_number_mapping[tissue_name] = tissue_num
		category_names.append('eQTL_' + tissue_name)
	category_names = np.asarray(['eQTL_aggregate'])
	category_counts = np.zeros(1)

	# Create mapping from rsid to reference and alternate alleles
	f = open(plink_reference_bim_file)
	rsid_mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] != chrom_num:
			print('assumption error')
			pdb.set_trace()
		rsid1 = data[1] + ':' + data[4] + '_' + data[5]
		rsid2 = data[1] + ':' + data[5] + '_' + data[4]
		if rsid1 in rsid_mapping or rsid2 in rsid_mapping:
			print('assumption eororro')
			pdb.set_trace()
		rsid_mapping[rsid1] = 1.0
		rsid_mapping[rsid2] = -1.0
	f.close()
	
	# Create gene array
	gene_arr = []
	f = open(causal_eqtl_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Skip genes not on the chromosome of interest
		line_chrom_num = data[2]
		if line_chrom_num != chrom_num:
			continue
		# Skip genes that dont pass boolean
		gene_boolean = data[4]
		if gene_boolean != 'Pass':
			continue

		# Extract relevent fields
		gene_id = data[1] + ':' + data[0]
		variant_info_file = data[5]
		susie_pmces_file = data[9]
		tissue_name = data[0]

		susie_pmces = np.load(susie_pmces_file)
		init_susie_pmces = np.copy(susie_pmces)
		g = open(variant_info_file)
		gene_rsids = []
		var_counter = 0
		for line2 in g:
			line2 = line2.rstrip()
			data2 = line2.split('\t')
			full_rsid = data2[1] + ':' + data2[4] + '_' + data2[5]
			gene_rsids.append(data2[1])
			susie_pmces[var_counter] = susie_pmces[var_counter]*rsid_mapping[full_rsid]
			var_counter = var_counter + 1
		g.close()
		gene_rsids = np.asarray(gene_rsids)
		if var_counter != len(susie_pmces):
			print('assumption eroror')
			pdb.set_trace()
		category_number = 0
		gene_arr.append((gene_id,category_number, gene_rsids, susie_pmces))
		category_counts[category_number] = category_counts[category_number] + 1
	f.close()
	return gene_arr, category_names, category_counts

def create_mapping_rsid_to_reference_index(rsids):
	mapping = {}
	for ii, rsid in enumerate(rsids):
		if rsid in mapping:
			print('assumption erororor')
			pdb.set_trace()
		mapping[rsid] = ii
	return mapping


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



def compute_gene_level_ld_scores_for_single_gene(gene_rsids, gene_eqtl_effect_sizes, gene_level_ld_scores, gene_category, ref_genotype_obj, rsid_to_reference_index, regression_rsid_to_regression_snp_position):
	gene_ref_positions = []
	for gene_rsid in gene_rsids:
		gene_ref_positions.append(rsid_to_reference_index[gene_rsid])
	gene_ref_positions = np.asarray(gene_ref_positions)

	# Get n_ref_panel_samples
	n_ref_panel_samples = ref_genotype_obj['G'].shape[0]

	# Now get ref positions of gene snp id lb and ub
	ref_pos_gene_snp_lb = np.min(gene_ref_positions)
	ref_pos_gene_snp_ub = np.max(gene_ref_positions)

	# Gene window cm lb and ub
	gene_window_cm_lb = ref_genotype_obj['cm'][ref_pos_gene_snp_lb] - 1.0
	gene_window_cm_ub = ref_genotype_obj['cm'][ref_pos_gene_snp_ub] + 1.0

	# Get indices corresponding to ref genotype for this gene window
	ref_genotype_indices_for_window = (ref_genotype_obj['cm'] >= gene_window_cm_lb) & (ref_genotype_obj['cm'] < gene_window_cm_ub)

	# Get pred expression
	# First need to standardize expression
	eqtl_geno = ref_genotype_obj['G'][:, gene_ref_positions]
	for snp_iter in range(eqtl_geno.shape[1]):
		eqtl_geno[:, snp_iter] = (eqtl_geno[:, snp_iter] - np.mean(eqtl_geno[:, snp_iter]))/np.std(eqtl_geno[:, snp_iter])
	# Now compute genetically predicted expression
	pred_expr = np.dot(eqtl_geno, gene_eqtl_effect_sizes)


	# Get window rsids and window genotype
	window_rsids = ref_genotype_obj['rsid'][ref_genotype_indices_for_window]
	window_geno = ref_genotype_obj['G'][:, ref_genotype_indices_for_window]

	# Loop through window rsids looking for ones that are regression rsids
	for ii, window_rsid in enumerate(window_rsids):
		if window_rsid not in regression_rsid_to_regression_snp_position:
			continue
		# Get squared correlation
		sqared_correlation = np.square(np.corrcoef(pred_expr, window_geno[:,ii])[0,1])
		# Get adjusted squared correlation
		adj_squared_correlation = sqared_correlation - ((1.0-sqared_correlation)/(n_ref_panel_samples-2.0))

		# Update gene level ld score for this regression snp
		regression_snp_pos = regression_rsid_to_regression_snp_position[window_rsid]
		gene_level_ld_scores[regression_snp_pos, gene_category] = gene_level_ld_scores[regression_snp_pos, gene_category] + adj_squared_correlation

	return gene_level_ld_scores


def extract_snp_gene_distances(tmp_gene_model_arr, distance_output_file):
	t = open(distance_output_file,'w')
	dist_arr = []
	t.write('snp_gene_pair\tdistance\n')
	for snp_gene in tmp_gene_model_arr:
		t.write(snp_gene[0] + '\t' + str(snp_gene[1]) + '\n')
		dist_arr.append(snp_gene[1])
	t.close()
	return np.asarray(dist_arr)


#######################
# Command line args
#######################
chrom_num = sys.argv[1]
causal_eqtl_summary_file = sys.argv[2]
gwas_genotype_stem = sys.argv[3]
ldsc_baseline_ld_annotation_stem = sys.argv[4]
gene_annotation_file = sys.argv[5]
tglr_expression_score_dir = sys.argv[6]
version = sys.argv[7]

print('CHROM' + str(chrom_num) + '\t' + version)

variant_level_ld_score_file = ldsc_baseline_ld_annotation_stem + str(chrom_num) + '.l2.ldscore.gz'

# Create array of genes (together with their stratefication category)
plink_reference_bim_file = gwas_genotype_stem + str(chrom_num) + '.bim'
if version == 'tissue_stratification':
	gene_model_arr, gene_category_names, gene_category_counts = create_tissue_stratefication_gene_model_arr(causal_eqtl_summary_file, chrom_num, plink_reference_bim_file)
elif version == 'tissue_aggregation':
	gene_model_arr, gene_category_names, gene_category_counts = create_tissue_aggregation_gene_model_arr(causal_eqtl_summary_file, chrom_num, plink_reference_bim_file)
elif version == 'distance_stratification':
	tmp_gene_model_arr = extract_gene_snp_distances_all_chrom(causal_eqtl_summary_file, gene_annotation_file)
	distance_output_file = tglr_expression_score_dir + 'snp_gene_distance_summary_' + version + '_chrom_' + str(chrom_num) + '.txt'
	ordered_snp_gene_distances = np.sort(extract_snp_gene_distances(tmp_gene_model_arr, distance_output_file))
	#np.save('snp_gene_dist.npy', ordered_snp_gene_distances)
	#ordered_snp_gene_distances = np.load('snp_gene_dist.npy')
	gene_model_arr, gene_category_names, gene_category_counts = create_distance_stratefication_gene_model_arr(causal_eqtl_summary_file, chrom_num, plink_reference_bim_file, gene_annotation_file, ordered_snp_gene_distances)
elif version == 'qtl_rank_stratification':
	gene_model_arr, gene_category_names, gene_category_counts = create_qtl_rank_stratefication_gene_model_arr(causal_eqtl_summary_file, chrom_num, plink_reference_bim_file)
else:
	print('not yet implemented')
	pdb.set_trace()


# Load in Reference Genotype data
G_obj = read_plink1_bin(gwas_genotype_stem + str(chrom_num) + '.bed', gwas_genotype_stem + str(chrom_num) + '.bim', gwas_genotype_stem + str(chrom_num) + '.fam', verbose=False)



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
genotype_obj = {'G': G_obj_geno, 'rsid': G_obj_rsids, 'position': G_obj_pos, 'cm': G_obj_cm}


# Create mapping from rsid to reference index
rsid_to_reference_index = create_mapping_rsid_to_reference_index(genotype_obj['rsid'])


# Get regression rsids
regression_rsids, regression_rsid_to_regression_snp_position = load_in_regression_rsids(variant_level_ld_score_file)
n_regression_snps = len(regression_rsids)

# Initialize matrix to keep track of gene level ld scores for regression snps on this chromosome
gene_level_ld_scores = np.zeros((n_regression_snps, len(gene_category_names)))


# Loop through genes
for gene_info in gene_model_arr:
	# Extract relevent fields
	gene_name = gene_info[0]
	gene_category = gene_info[1]
	gene_rsids = gene_info[2]
	gene_pmces = gene_info[3]

	# Update gene level ld scores with this gene
	gene_level_ld_scores = compute_gene_level_ld_scores_for_single_gene(gene_rsids, gene_pmces, gene_level_ld_scores, gene_category, genotype_obj, rsid_to_reference_index, regression_rsid_to_regression_snp_position)


# Save gene level ld scores to output to output file
output_file = tglr_expression_score_dir + 'tglr_gene_ld_scores_' + version + '_chrom_' + str(chrom_num) + '.txt'
np.savetxt(output_file, gene_level_ld_scores, fmt="%s")

# Save number of genes to output also
output_file2 = tglr_expression_score_dir +'tglr_n_genes_' + version + '_chrom_' + str(chrom_num) + '.txt'
t = open(output_file2,'w')
for category_index, category_name in enumerate(gene_category_names):
	t.write(category_name + '\t' + str(gene_category_counts[category_index]) + '\n')
t.close()

print(output_file)
print(output_file2)

import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
from pandas_plink import read_plink1_bin
import statsmodels.api as sm

def mean_impute_and_standardize_genotype(G_obj_geno):
	# Fill in missing values
	G_obj_geno_stand = np.copy(G_obj_geno)
	ncol = G_obj_geno_stand.shape[1]
	n_missing = []
	for col_iter in range(ncol):
		nan_indices = np.isnan(G_obj_geno[:,col_iter])
		non_nan_mean = np.mean(G_obj_geno[nan_indices==False, col_iter])
		G_obj_geno_stand[nan_indices, col_iter] = non_nan_mean
		n_missing.append(np.sum(nan_indices))
	n_missing = np.asarray(n_missing)

	# Quick error check
	if np.sum(np.std(G_obj_geno_stand,axis=0) == 0) > 0:
		print('no variance genotype assumption error')
		pdb.set_trace()

	# Standardize genotype
	G_obj_geno_stand = (G_obj_geno_stand -np.mean(G_obj_geno_stand,axis=0))/np.std(G_obj_geno_stand,axis=0)

	return G_obj_geno_stand




def construct_variant_ld_scores_based_on_distance(genotype_stem, output_stem, window_size=1):
	# Load in genotype object
	genotype_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)

	snp_pos = np.asarray(genotype_obj['pos'])
	snp_rsids = np.asarray(genotype_obj['snp'])
	snp_variant_ids = np.asarray(genotype_obj['variant'])

	inner_window_start_pos = np.min(snp_pos) - 1

	mb = 1000000.0

	used_snps = {}

	output_file = output_stem + '.txt'
	t = open(output_file,'w')
	output_file2 = output_stem + '_bias_corrected.txt'
	t2 = open(output_file2,'w')	
	t.write('variant_id\tld_score\n')
	t2.write('variant_id\tld_score\n')

	window_counter = 0

	while inner_window_start_pos < np.max(snp_pos):

		print(window_counter)
		window_counter = window_counter + 1

		# Define windows
		inner_window_end_pos = inner_window_start_pos + window_size*mb
		outer_window_start_pos = inner_window_start_pos - (window_size*mb)
		outer_window_end_pos = inner_window_end_pos + (window_size*mb)

		window_position_indices = (snp_pos >= outer_window_start_pos) & (snp_pos < outer_window_end_pos)

		if np.sum(window_position_indices) == 0:
			inner_window_start_pos = inner_window_start_pos + (window_size*mb)
			continue

		# Extract genotype data for this window
		window_variant_genotype = np.asarray(genotype_obj.sel(variant=snp_variant_ids[window_position_indices]))
		# standardize genotype
		window_std_variant_genotype = mean_impute_and_standardize_genotype(window_variant_genotype)

		# Compute LD matrix for window
		ref_ss =window_std_variant_genotype.shape[0]
		squared_ld_mat = np.square(np.corrcoef(np.transpose(window_std_variant_genotype)))
		bias_corrected_squared_ld_mat = squared_ld_mat - ((1.0-squared_ld_mat)/(ref_ss-2.0))

		# Get into window orientation
		window_snp_pos = snp_pos[window_position_indices]
		window_rsids = snp_rsids[window_position_indices]

		# Loop through snps until we find a middle snp
		for tmp_snp_iter, tmp_snp_pos in enumerate(window_snp_pos):
			if tmp_snp_pos >= inner_window_start_pos and tmp_snp_pos < inner_window_end_pos:
				# This is a middle snp
				tmp_snp_rsid = window_rsids[tmp_snp_iter]

				if tmp_snp_rsid in used_snps:
					print('assumption eroror')
					pdb.set_trace()
				used_snps[tmp_snp_rsid] = 1

				tmp_snp_indices = (window_snp_pos >= tmp_snp_pos - (window_size*mb)) & (window_snp_pos < tmp_snp_pos + (window_size*mb))
				
				snp_ld_score = np.sum(squared_ld_mat[tmp_snp_iter, tmp_snp_indices])
				bias_corrected_snp_ld_score = np.sum(bias_corrected_squared_ld_mat[tmp_snp_iter, tmp_snp_indices])

				t.write(tmp_snp_rsid + '\t' + str(snp_ld_score) + '\n')
				t2.write(tmp_snp_rsid + '\t' + str(bias_corrected_snp_ld_score) + '\n')

		inner_window_start_pos = inner_window_start_pos + (window_size*mb)

	t.close()
	t2.close()


	return















##############################
# Command line argumemnts
##############################
simulation_number = int(sys.argv[1])
chrom_num = sys.argv[2]
simulation_name_string = sys.argv[3]
processed_genotype_data_dir = sys.argv[4]
simulated_gwas_dir = sys.argv[5]



ref_genotype_stem = processed_genotype_data_dir +'simulated_eqtl_1000_data_1'
output_stem = simulated_gwas_dir + simulation_name_string + '_variant_ref_geno_window_1_mb_ld_scores'
construct_variant_ld_scores_based_on_distance(ref_genotype_stem, output_stem, window_size=1)


ref_genotype_stem = processed_genotype_data_dir +'simulated_eqtl_1000_data_1'
output_stem = simulated_gwas_dir + simulation_name_string + '_variant_ref_geno_window_2_mb_ld_scores'
construct_variant_ld_scores_based_on_distance(ref_genotype_stem, output_stem, window_size=2)


in_sample_genotype_stem = processed_genotype_data_dir + 'simulated_gwas_data_1'


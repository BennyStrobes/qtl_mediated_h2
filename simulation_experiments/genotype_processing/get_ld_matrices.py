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
	G_obj_geno_stand3 = (G_obj_geno_stand -np.mean(G_obj_geno_stand,axis=0))/np.std(G_obj_geno_stand,axis=0)


	return G_obj_geno_stand3



def construct_variant_ld_mat_based_on_distance(genotype_stem, output_stem, window_size=1):
	# Load in genotype object
	genotype_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)

	snp_pos = np.asarray(genotype_obj['pos'])
	snp_rsids = np.asarray(genotype_obj['snp'])
	snp_variant_ids = np.asarray(genotype_obj['variant'])

	inner_window_start_pos = np.min(snp_pos) - 1

	mb = 1000000.0

	used_snps = {}

	output_file = output_stem + '_summary.txt'
	t = open(output_file,'w')
	t.write('window_name\twindow_start\twindow_end\tld_file\twindow_snp_position_file\trsid_file\tmiddle_boolean_file\n')

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
		ld_mat = np.corrcoef(np.transpose(window_std_variant_genotype))

		# Get into window orientation
		window_snp_pos = snp_pos[window_position_indices]
		window_rsids = snp_rsids[window_position_indices]

		middle_rsids_bool = []
		# Loop through snps until we find a middle snp
		for tmp_snp_iter, tmp_snp_pos in enumerate(window_snp_pos):
			if tmp_snp_pos >= inner_window_start_pos and tmp_snp_pos < inner_window_end_pos:
				# This is a middle snp
				middle_rsids_bool.append(1.0)
			else:
				middle_rsids_bool.append(0.0)
		middle_rsids_bool = np.asarray(middle_rsids_bool)
		print(len(middle_rsids_bool))

		window_global_positions = np.where(window_position_indices==True)[0]


		# Start saving things to output file
		window_name = 'window_' + str(int(outer_window_start_pos)) + '_' + str(int(outer_window_end_pos))

		# Save LD mat to output
		ld_mat_file = output_stem + '_' + window_name + '_ld.npy'
		np.save(ld_mat_file, ld_mat)

		# Save rsids to outtput file
		rsid_file = output_stem + '_' + window_name + '_rsids.npy'
		np.save(rsid_file, window_rsids)

		# Save middle snp boolean to output file
		middle_boolean_file = output_stem + '_' + window_name + '_middle_boolean.npy'
		np.save(middle_boolean_file, middle_rsids_bool)

		# Save window positions to output
		window_position_file = output_stem + '_' + window_name + '_window_positions.npy'
		np.save(window_position_file, window_global_positions)


		t.write(window_name + '\t' + str(int(outer_window_start_pos)) + '\t' + str(int(outer_window_end_pos)) + '\t' + ld_mat_file + '\t' + window_position_file + '\t' + rsid_file + '\t' + middle_boolean_file + '\n')	
		inner_window_start_pos = inner_window_start_pos + (window_size*mb)

	t.close()
	return


##############################
# Command line argumemnts
##############################
processed_genotype_data_dir = sys.argv[1]






window_size = 1

ref_genotype_stem = processed_genotype_data_dir +'simulated_eqtl_1000_data_1'
output_stem = processed_genotype_data_dir + 'variant_ref_geno_eqtl_1000_window_' + str(window_size) + '_mb_ld'
construct_variant_ld_mat_based_on_distance(ref_genotype_stem, output_stem, window_size=window_size)

in_sample_genotype_stem = processed_genotype_data_dir +'simulated_gwas_data_1'
output_stem = processed_genotype_data_dir + 'variant_ref_geno_gwas_window_' + str(window_size) + '_mb_ld'
construct_variant_ld_mat_based_on_distance(in_sample_genotype_stem, output_stem, window_size=window_size)

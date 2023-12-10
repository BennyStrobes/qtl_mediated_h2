import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
from pandas_plink import read_plink1_bin




def get_nvar(window, bim_file):
	nvar = 0
	window_start = int(window.split(':')[0])
	window_end = int(window.split(':')[1])
	f = open(bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		pos = float(data[3])
		if pos >= window_start and pos < window_end:
			nvar = nvar + 1

	f.close()
	return nvar


def get_windows(bim_file):
	start = 10000000000
	end = -1
	f = open(bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		pos = float(data[3])
		if pos < start:
			start = pos
		if pos > end:
			end = pos
	f.close()
	windows = []
	cur_start = start
	while cur_start < end:
		window = str(int(cur_start)) + ':' + str(int(cur_start+7000000))
		cur_start = cur_start+7000000
		windows.append(window)

	windows = np.asarray(windows)

	'''
	for window in windows:
		nvar = get_nvar(window, bim_file)
		print(window + '\t' + str(nvar))
	'''

	return windows[2:5]




def extract_dictionary_list_of_window_rsids(bim_file, window_start, window_end):
	dicti = {}
	f = open(bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		rs_id = data[1]
		snp_pos = int(data[3])
		if snp_pos >= window_start and snp_pos < window_end: # in window
			# Quick error check
			if rs_id in dicti:
				print('assumption eroror')
				pdb.set_trace()
			dicti[rs_id] = 1
	f.close()
	return dicti


def load_in_ordered_array_of_rsids_from_bim_file(genotype_bim):
	f = open(genotype_bim)
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		arr.append(data[1])
	f.close()
	arr = np.asarray(arr)

	# Quick error check 
	if len(np.unique(arr)) != len(arr):
		print('assumption eroror')
		pdb.set_trace()

	return arr


def extract_variants_in_window(gwas_plink_stem, window_rsids):
	# Load in ordered array of rsids
	genotype_bim = gwas_plink_stem + '.bim'
	ordered_rsids = load_in_ordered_array_of_rsids_from_bim_file(genotype_bim)

	# Load in variant names in plink format
	window_variant_names = []
	for snp_iter, snp_rsid in enumerate(ordered_rsids):
		# Skip snps that are not hapmap3 snps
		if snp_rsid not in window_rsids:
			continue
		window_variant_names.append('variant' + str(snp_iter))
	window_variant_names = np.asarray(window_variant_names)

	# Load in genotype object
	genotype_obj = read_plink1_bin(gwas_plink_stem + '.bed', gwas_plink_stem + '.bim', gwas_plink_stem + '.fam', verbose=False)


	# Extract genotype data for this batch of snps
	window_variant_genotype = np.asarray(genotype_obj.sel(variant=window_variant_names))

	# mean impute genotype
	nvary = window_variant_genotype.shape[1]
	for var_iter in range(nvary):
		variant_genotype = window_variant_genotype[:,var_iter]
		std_variant_genotype = np.copy(variant_genotype)
		nan_indices = np.isnan(variant_genotype)
		non_nan_mean = np.mean(variant_genotype[nan_indices==False])
		std_variant_genotype[nan_indices] = non_nan_mean

		std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)
		window_variant_genotype[:,var_iter] = std_variant_genotype

	print(np.sum(np.isnan(window_variant_genotype)))

	return window_variant_genotype



#######################
# Command line args
#######################
chrom_num = sys.argv[1]
processed_genotype_data_dir = sys.argv[2]

# Bim file for this chromosomse
bim_file = processed_genotype_data_dir + 'simulated_gwas_data_' + chrom_num + '.bim'

windows = get_windows(bim_file)


# Loop through windows
window_iter = 1
for window_name in windows:
	print(window_name)
	if window_iter > 1:
		continue

	# Extract relevent info from this window
	window_start = int(window_name.split(':')[0])
	window_end = int(window_name.split(':')[1])

	# Extract dictionary list of window rsids
	window_rsids = extract_dictionary_list_of_window_rsids(bim_file, window_start, window_end)

	# Extract gwas genotype matrix
	gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype files
	window_gwas_genotype = extract_variants_in_window(gwas_plink_stem, window_rsids)
	# Save gwas genotype
	output_file = processed_genotype_data_dir + 'gwas_genotype_' + str(window_iter) + '.npy'
	np.save(output_file, window_gwas_genotype)

	# Save variant LD
	LD = np.corrcoef(np.transpose(window_gwas_genotype))
	# Save LD genotype
	output_file = processed_genotype_data_dir + 'gwas_genotype_LD_' + str(window_iter) + '.npy'
	np.save(output_file, LD)

	# LD LD transpose
	'''
	ld_ld_t = np.dot(LD, np.transpose(LD))
	# Save 
	output_file = processed_genotype_data_dir + 'gwas_genotype_LD_LD_t_' + str(window_iter) + '.npy'
	np.save(output_file, ld_ld_t)
	'''


	# Iterate through eqtl sample sizes
	eqtl_sss = [100, 300, 500, 1000, 5000]
	for eqtl_ss in eqtl_sss:
		print(eqtl_ss)
		eqtl_plink_stem = processed_genotype_data_dir + 'simulated_eqtl_' + str(eqtl_ss) + '_data_' + str(chrom_num)
		window_eqtl_genotype = extract_variants_in_window(eqtl_plink_stem, window_rsids)
		# Save eqtl genotype
		output_file = processed_genotype_data_dir + 'eqtl_' + str(eqtl_ss) + '_genotype_' + str(window_iter) + '.npy'
		np.save(output_file, window_eqtl_genotype)

		# Save eqtl genotype LD
		output_file = processed_genotype_data_dir + 'eqtl_' + str(eqtl_ss) + '_genotype_LD_' + str(window_iter) + '.npy'
		eqtl_LD = np.corrcoef(np.transpose(window_eqtl_genotype))
		np.save(output_file, eqtl_LD)


	window_iter = window_iter + 1





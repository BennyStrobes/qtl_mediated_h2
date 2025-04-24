import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import pandas as pd
from sklearn.decomposition import PCA
import gzip


def extract_downsample_names(downsampled_sample_names_file):
	dicti = {}
	arr = []
	f = open(downsampled_sample_names_file)

	for line in f:
		line = line.rstrip()
		arr.append(line)
		dicti[line] = 1
	f.close()
	return np.asarray(arr), dicti


def generate_expression_pcs(expression_file, pc_file):
	# Load in expression data
	expr_full = np.transpose(np.loadtxt(expression_file, dtype=str,delimiter='\t'))
	expr = expr_full[3:,1:].astype(float)
	samples = expr_full[3:,0]
	genes = expr_full[0,1:]

	# Get number of PCs
	# Based on section section 3.5.1 of GTEx v8 paper
	sample_size = len(samples)
	if sample_size < 150:
		num_expression_pcs = 15
	elif sample_size >= 150 and sample_size < 250:
		num_expression_pcs = 30
	elif sample_size >= 250 and sample_size < 350:
		num_expression_pcs = 45
	elif sample_size >= 350:
		num_expression_pcs = 60
	else:
		print('assumption error')
		pdb.set_trace()

	# This way of computing pcs is faster
	_pca = PCA(n_components=num_expression_pcs, svd_solver='arpack')
	expr_pc_loadings = _pca.fit_transform(expr)

	# print to output file 
	cov_names = []
	for cov_iter in range(num_expression_pcs):
		cov_names.append('PC' + str(cov_iter))
	cov_names = np.asarray(cov_names)
	t = open(pc_file, 'w')
	t.write('FID\tIID\t' + '\t'.join(cov_names) + '\n')
	for sample_num, sample_name in enumerate(samples):
		t.write('0' + '\t' + sample_name + '\t' + '\t'.join(expr_pc_loadings[sample_num,:].astype(str)) + '\n')
	t.close()
	return


def downsample_gene_expression_matrix(downsample_names_arr, downsample_names_dicti, full_data_expression_file, downsampled_expression_file):
	f = open(full_data_expression_file)
	t = open(downsampled_expression_file,'w')

	downsample_names_dicti['GENE'] = 1
	downsample_names_dicti['CHR'] = 1
	downsample_names_dicti['GENE_COORD'] = 1

	# For header
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			valid_indices = []
			for ii, val in enumerate(data):
				if val in downsample_names_dicti:
					valid_indices.append(ii)
			valid_indices = np.asarray(valid_indices)

			subset_data = np.asarray(data)[valid_indices]

			# Quick error check
			if np.array_equal(subset_data[3:], downsample_names_arr) == False:
				print('assumption eroorro')
				pdb.set_trace()

			# Print to output
			t.write('\t'.join(subset_data) + '\n')
			continue

		subset_data = np.asarray(data)[valid_indices]

		# Standardize subsetted expression
		un_norm_expr = subset_data[3:].astype(float)
		norm_expr = (un_norm_expr - np.mean(un_norm_expr))/(np.std(un_norm_expr))
		t.write('\t'.join(subset_data[:3]) + '\t' + '\t'.join(norm_expr.astype(str)) + '\n')

	f.close()
	t.close()
	return




###########################
# Command line args
###########################
full_data_expression_file = sys.argv[1]  # Non-downsampled gene expression mat
downsampled_sample_names_file = sys.argv[2]  # file containing list of samples to downsample to
downsampled_expression_file = sys.argv[3]  # output file
downsampled_expression_pc_file = sys.argv[4]  # output file


# Extract array of sample names to downsample to
downsample_names_arr, downsample_names_dicti = extract_downsample_names(downsampled_sample_names_file)


# Downsample gene expression matrix
downsample_gene_expression_matrix(downsample_names_arr, downsample_names_dicti, full_data_expression_file, downsampled_expression_file)

# Generate expression pcs
generate_expression_pcs(downsampled_expression_file, downsampled_expression_pc_file)


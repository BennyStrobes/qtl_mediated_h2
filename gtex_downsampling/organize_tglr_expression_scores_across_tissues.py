import numpy as np
import sys
import pdb
import os
import gzip


def get_tissue_names(tissue_info_file):
	arr = []

	head_count= 0

	f = open(tissue_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)




def get_ordered_regression_variant_strings(ldsc_genotype_intercept_annotation_stem, chrom_num):
	arr = []
	ldscore_file = ldsc_genotype_intercept_annotation_stem + str(chrom_num) + '.l2.ldscore.gz'
	f = gzip.open(ldscore_file)
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		stringer = data[0] + '\t' + data[1] + '\t' + data[2]
		arr.append(stringer)
	f.close()
	return np.asarray(arr)


def get_multi_tissue_gene_ld_scores(tissue_names, tglr_expression_score_dir, mesc_run_name, chrom_num):
	global_arr = []
	n_gene_arr = []
	for tissue_name in tissue_names:
		tissue_arr = []

		tissue_file = tglr_expression_score_dir + mesc_run_name + '_' + tissue_name + '_' + str(chrom_num) + '_tglr_gene_ld_scores.txt'
		f = open(tissue_file)
		for line in f:
			line = line.rstrip()
			tissue_arr.append(line)
		f.close()

		n_gene_file = tglr_expression_score_dir + mesc_run_name + '_' + tissue_name + '_' + str(chrom_num) + '_tglr_n_genes.txt'
		f = open(n_gene_file)
		counter = 0
		for line in f:
			line = line.rstrip()
			n_gene_arr.append(line)
			if counter > 0:
				print('asumption eroro')
				pdb.set_trace()
			counter = counter + 1
		f.close()

		tissue_arr = np.asarray(tissue_arr)
		global_arr.append(tissue_arr)


	n_gene_arr = np.asarray(n_gene_arr)
	global_arr = np.transpose(np.asarray(global_arr))

	if len(n_gene_arr) != global_arr.shape[1]:
		print('assumption eroro')
		pdb.set_trace()

	return global_arr, n_gene_arr

def make_m_file(m_file, n_genes):
	t = open(m_file,'w')

	t.write('\t'.join(n_genes) + '\n')

	t.close()
	return


#######################
# Command line args
#######################
tissue_info_file = sys.argv[1]
mesc_run_name = sys.argv[2]
ldsc_genotype_intercept_annotation_stem = sys.argv[3]
tglr_expression_score_dir = sys.argv[4]
expression_score_annotation_stem = sys.argv[5]


# First extract names of tissues
tissue_names = get_tissue_names(tissue_info_file)


for chrom_num in range(1,23):
	# Extract ordered regression variants
	ordered_regression_variant_strings = get_ordered_regression_variant_strings(ldsc_genotype_intercept_annotation_stem, chrom_num)


	# Extract multi-tissue gene LD scores
	multitissue_gene_ld_scores, n_genes = get_multi_tissue_gene_ld_scores(tissue_names, tglr_expression_score_dir, mesc_run_name, chrom_num)

	# Error checking
	if len(ordered_regression_variant_strings) != multitissue_gene_ld_scores.shape[0]:
		print('assumption eroror')
		pdb.set_trace()


	# Make Gene expression ld score file
	expression_score_output_file = expression_score_annotation_stem + str(chrom_num) + '.l2.ldscore'
	t = open(expression_score_output_file,'w')

	# Header
	t.write('CHR\tSNP\tBP')
	for tissue_name in tissue_names:
		t.write('\t' + 'gene_score_' + tissue_name)
	t.write('\n')

	for ii, ordered_regression_variant_string in enumerate(ordered_regression_variant_strings):
		t.write(ordered_regression_variant_string + '\t' + '\t'.join(multitissue_gene_ld_scores[ii,:]) + '\n')
	t.close()


	# Make M files
	m_file = expression_score_annotation_stem + str(chrom_num) + '.l2.M'
	make_m_file(m_file, n_genes)

	m_5_50_file = expression_score_annotation_stem + str(chrom_num) + '.l2.M_5_50'
	make_m_file(m_5_50_file, n_genes)



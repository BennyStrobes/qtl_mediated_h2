import numpy as np
import sys
import pdb
import os
import gzip




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

def get_gene_ld_scores(tglr_expression_score_dir, qtl_version, chrom_num):
	# Get gene LD scores
	gene_ld_score_file = tglr_expression_score_dir + 'tglr_gene_ld_scores_' + qtl_version + '_chrom_' + str(chrom_num) + '.txt'
	gene_ld_scores = []
	f = open(gene_ld_score_file)
	for line in f:
		line = line.rstrip()
		data = np.asarray(line.split())
		gene_ld_scores.append(data)
	f.close()
	gene_ld_scores = np.asarray(gene_ld_scores)

	# Get number of genes in each category as well as qtl category names
	n_genes_file = tglr_expression_score_dir + 'tglr_n_genes_' + qtl_version + '_chrom_' + str(chrom_num) + '.txt'
	categories = []
	n_genes = []
	f = open(n_genes_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		categories.append(data[0])
		n_genes.append(data[1])
	f.close()
	categories = np.asarray(categories)
	n_genes = np.asarray(n_genes)

	return gene_ld_scores, categories, n_genes

def make_m_file(m_file, n_genes):
	t = open(m_file,'w')

	t.write('\t'.join(n_genes) + '\n')

	t.close()
	return

#######################
# Command line args
#######################
qtl_version = sys.argv[1]
ldsc_genotype_intercept_annotation_stem = sys.argv[2]
tglr_expression_score_dir = sys.argv[3]
expression_score_annotation_stem = sys.argv[4]


for chrom_num in range(1,23):
	print(chrom_num)
	# Extract ordered regression variants
	ordered_regression_variant_strings = get_ordered_regression_variant_strings(ldsc_genotype_intercept_annotation_stem, chrom_num)

	# Extract gene LD scores as well as the number of genes, and names of qtl categories
	gene_ld_scores, qtl_category_names, n_genes = get_gene_ld_scores(tglr_expression_score_dir, qtl_version, chrom_num)
	

	# Error checking
	if len(ordered_regression_variant_strings) != gene_ld_scores.shape[0]:
		print('assumption eroror')
		pdb.set_trace()

	# Make Gene expression ld score file
	expression_score_output_file = expression_score_annotation_stem + str(chrom_num) + '.l2.ldscore'
	t = open(expression_score_output_file,'w')

	# Header
	t.write('CHR\tSNP\tBP')
	for qtl_category_name in qtl_category_names:
		t.write('\t' + 'gene_score_' + qtl_category_name)
	t.write('\n')

	for ii, ordered_regression_variant_string in enumerate(ordered_regression_variant_strings):
		t.write(ordered_regression_variant_string + '\t' + '\t'.join(gene_ld_scores[ii,:]) + '\n')
	t.close()


	# Make M files
	m_file = expression_score_annotation_stem + str(chrom_num) + '.l2.M'
	make_m_file(m_file, n_genes)

	m_5_50_file = expression_score_annotation_stem + str(chrom_num) + '.l2.M_5_50'
	make_m_file(m_5_50_file, n_genes)


import numpy as np
import os
import sys
import pdb
import gzip

def extract_gene_ld_scores(expression_score_dir, mesc_run_name):
	gene_ld_scores = []
	rsids = []

	for chrom_num in range(1,23):
		expr_score_file = expression_score_dir + mesc_run_name + '_meta_analyzed_scores.' + str(chrom_num) + '.expscore.gz'

		f = gzip.open(expr_score_file)
		head_count = 0
		for line in f:
			line = line.decode('utf-8').rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				anno_names = np.asarray(data[3:])
				continue
			rsid = data[1]
			gene_ld_scores.append(np.asarray(data[3:]).astype(float))
			rsids.append(rsid)
		f.close()

	return np.asarray(gene_ld_scores), np.asarray(rsids), anno_names






def extract_variant_ld_scores(ldsc_baseline_ld_annotation_stem):
	var_ld_scores = []
	rsids = []

	for chrom_num in range(1,23):
		filer = ldsc_baseline_ld_annotation_stem + str(chrom_num) + '.l2.ldscore.gz'
		f = gzip.open(filer)
		head_count = 0
		for line in f:
			line = line.decode('utf-8').rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				anno_names = np.asarray(data[3:])
				continue
			rsid = data[1]
			var_ld_scores.append(np.asarray(data[3:]).astype(float))
			rsids.append(rsid)

	return np.asarray(var_ld_scores), np.asarray(rsids), anno_names




def get_shared_rsids(gene_ld_score_rsids, var_ld_score_rsids):
	gene_var = {}
	for var in gene_ld_score_rsids:
		gene_var[var] = 1
	shared = {}
	for var in var_ld_score_rsids:
		if var in gene_var:
			shared[var] = 1
	return shared

def get_shared_indices(shared_rsids, rsid_list):
	indices = []
	for rsid in rsid_list:
		if rsid in shared_rsids:
			indices.append(True)
		else:
			indices.append(False)
	return np.asarray(indices)

def genomic_bootstrap_correlation(a1, a2, n_genomic_bins=200, n_bootstraps=100):
	chunks = np.array_split(np.arange(len(a1)),n_genomic_bins)

	bs_corrs = []
	for bs_iter in range(n_bootstraps):
		bs_sample = np.random.choice(np.arange(n_genomic_bins),size=n_genomic_bins)

		full_indices = []
		for chunk_index in bs_sample:
			full_indices.append(chunks[chunk_index])
		full_indices = np.hstack(full_indices)

		bs_corr = np.corrcoef(a1[full_indices], a2[full_indices])[0,1]
		bs_corrs.append(bs_corr)

	bs_corrs = np.asarray(bs_corrs)

	return np.std(bs_corrs)



######################
# Command line args
######################
mesc_run_name = sys.argv[1]
expression_score_dir = sys.argv[2]
ldsc_baseline_ld_annotation_stem = sys.argv[3]
gene_ldscore_variant_ld_score_corr_dir = sys.argv[4]


####
# Extract gene ldscores
gene_ld_scores, gene_ld_score_rsids, gene_ld_score_anno_names = extract_gene_ld_scores(expression_score_dir, mesc_run_name)

###
# Extract variant ldscores
var_ld_scores, var_ld_score_rsids, var_ld_score_anno_names = extract_variant_ld_scores(ldsc_baseline_ld_annotation_stem)


# Get shared rsids
shared_rsids = get_shared_rsids(gene_ld_score_rsids, var_ld_score_rsids)
shared_gene_indices = get_shared_indices(shared_rsids, gene_ld_score_rsids)
shared_var_indices = get_shared_indices(shared_rsids, var_ld_score_rsids)



if np.array_equal(gene_ld_score_rsids[shared_gene_indices], var_ld_score_rsids[shared_var_indices]) == False:
	print('assumption eroror')
	pdb.set_trace()

# Filter to shared rsids
gene_ld_score_rsids = gene_ld_score_rsids[shared_gene_indices]
gene_ld_scores = gene_ld_scores[shared_gene_indices]
agg_gene_ld_scores = np.sum(gene_ld_scores,axis=1)

var_ld_score_rsids = var_ld_score_rsids[shared_var_indices]
var_ld_scores = var_ld_scores[shared_var_indices]

output_file = gene_ldscore_variant_ld_score_corr_dir + mesc_run_name + 'gene_ld_score_variant_ld_score_correlations.txt'
t = open(output_file,'w')
t.write('gene_anno\tvariant_anno\tcorrelation\tcorrelation_se\n')

gene_anno_name = 'aggregate'
for var_anno_index, var_anno_name in enumerate(var_ld_score_anno_names):
	corry = np.corrcoef(var_ld_scores[:, var_anno_index], agg_gene_ld_scores)[0,1]
	#corry_se_parametric = np.sqrt((1-np.square(corry))/(len(agg_gene_ld_scores) - 2))
	corry_se = genomic_bootstrap_correlation(var_ld_scores[:, var_anno_index], agg_gene_ld_scores)
	t.write(gene_anno_name + '\t' + var_anno_name + '\t' + str(corry) + '\t' + str(corry_se) + '\n')


for gene_anno_index, gene_anno_name in enumerate(gene_ld_score_anno_names):
	for var_anno_index, var_anno_name in enumerate(var_ld_score_anno_names):
		corry = np.corrcoef(var_ld_scores[:, var_anno_index], gene_ld_scores[:, gene_anno_index])[0,1]
		#corry_se_parametric = np.sqrt((1-np.square(corry))/(len(agg_gene_ld_scores) - 2))
		corry_se = genomic_bootstrap_correlation(var_ld_scores[:, var_anno_index], gene_ld_scores[:, gene_anno_index])
		t.write(gene_anno_name + '\t' + var_anno_name + '\t' + str(corry) + '\t' + str(corry_se) + '\n')

t.close()



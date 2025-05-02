import numpy as np
import os
import sys
import pdb
from bgen import BgenReader
import argparse
import gzip
import time



def create_dictionary_mapping_from_rsid_to_cm_position(sldsc_anno_file):
	dicti = {}
	f = gzip.open(sldsc_anno_file)
	head_count = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid =data[2]
		if rsid in dicti:
			print('assumption error')
			pdb.set_trace()
		cm_position = float(data[3])
		dicti[rsid] = cm_position
	f.close()
	return dicti



def create_mapping_from_hm3_rsid_to_regression_index(variant_ld_score_file):
	dicti = {}
	f = open(variant_ld_score_file)
	counter = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[1]
		if rsid in dicti:
			print('assumptioneroror')
			pdb.set_trace()
		dicti[rsid] = counter
		counter = counter + 1
	f.close()

	return dicti

def extract_snp_cm_pos(ordered_rsids, rsid_to_cm):
	cm_arr = []
	for rsid in ordered_rsids:
		cm_arr.append(rsid_to_cm[rsid])
	return np.asarray(cm_arr)

def make_gene_regression_snp_summary_file(gene_regression_snp_summary_file, chrom_num, tmp_rsids, tmp_poss, tmp_cms, tmp_a0s, tmp_a1s, tmp_cis_snp_booleans, hm3_rsid_to_regression_index):
	t = open(gene_regression_snp_summary_file,'w')
	t.write('chr\trsid\tcm\tpos\tref_allele\talt_allele\tcis_snp_indicator\n')

	for snp_index, tmp_rsid in enumerate(tmp_rsids):
		tmp_cm = tmp_cms[snp_index]
		tmp_pos = tmp_poss[snp_index]
		tmp_a0 = tmp_a0s[snp_index]
		tmp_a1 = tmp_a1s[snp_index]
		tmp_cis_snp_boolean = tmp_cis_snp_booleans[snp_index]
		regression_snp_index = hm3_rsid_to_regression_index[tmp_rsid]

		t.write(str(chrom_num) + '\t' + tmp_rsid + '\t' + str(tmp_cm) + '\t' + str(tmp_pos) + '\t' + tmp_a0 + '\t' + tmp_a1 + '\t' + str(tmp_cis_snp_boolean*1.0) + '\n')
	t.close()
	return

def compute_lambda_thresh(lambdas, rho_thresh):
	totaler = np.sum(lambdas)
	cur_total = 0
	lambda_thresh = -1
	for lambda_val in -np.sort(-lambdas):
		cur_total = cur_total + lambda_val
		if cur_total/totaler > rho_thresh:
			if lambda_thresh == -1:
				lambda_thresh = lambda_val


	if lambda_thresh == -1:
		print('assumption eroror')
		pdb.set_trace()

	return lambda_thresh

def pairwise_correlations(Z, X):
	X_centered = X - X.mean(axis=1, keepdims=True)
	Z_centered = Z - Z.mean(axis=1, keepdims=True)

	# Normalize to unit variance (L2 norm)
	X_norm = np.linalg.norm(X_centered, axis=1, keepdims=True)
	Z_norm = np.linalg.norm(Z_centered, axis=1, keepdims=True)

	# Compute dot product and divide by norms to get Pearson correlation
	# Result: (K x L) matrix of correlations
	correlations = (X_centered @ Z_centered.T) / (X_norm * Z_norm.T)

	return correlations



#####################
# Parse command line arguments
#####################
parser = argparse.ArgumentParser()
parser.add_argument('--chrom', default='None', type=str,
                    help='Chromosome number (e.g. 1)')
parser.add_argument('--gene-position-file', default='None', type=str,
                    help='Absolute path to gene position file')
parser.add_argument('--bgen-file', default='None', type=str,
                    help='Absolute path to reference bgen genotype file')
parser.add_argument('--variant-ld-score-file', default='None', type=str,
                    help='Absolute path to variant ld score file')
parser.add_argument('--gene-snp-representation', default='None', type=str,
                    help='Which type of low dimensional representation to use')
parser.add_argument('--sldsc-annotation-file', default='None', type=str,
                    help='Absolute path to gziped sldsc annotation file')
parser.add_argument('--gene-ld-output-root', default='None', type=str,
                    help='Output root')
parser.add_argument('--cm-window-size', default=1.0, type=float,
                    help='size of CM window around genetic elements to use')
parser.add_argument('--min-cis-snps', default=5, type=int,
                    help='Minimimum number of allowable cis snps for a gene')
parser.add_argument('--min-regr-snps', default=5, type=int,
                    help='Minimimum number of allowable regression snps a gene')
args = parser.parse_args()



# CM window size (impacts memory/speed trast)
cm_window_size = args.cm_window_size

# Create dictionary mapping from rsid to CM position
rsid_to_cm = create_dictionary_mapping_from_rsid_to_cm_position(args.sldsc_annotation_file)

# Create dictionary mapping from HM3 rsid to regression index
hm3_rsid_to_regression_index = create_mapping_from_hm3_rsid_to_regression_index(args.variant_ld_score_file)


# Load in genotype object
genotype_obj = BgenReader(args.bgen_file)
snp_pos = np.asarray(genotype_obj.positions())
ordered_rsids = np.asarray(genotype_obj.rsids())
snp_cm = extract_snp_cm_pos(ordered_rsids, rsid_to_cm)
snp_integers = np.arange(len(ordered_rsids))


# Output summary file
output_summary_file = args.gene_ld_output_root + '_summary_file.txt'
t = open(output_summary_file,'w')


# Now loop through genes
f = open(args.gene_position_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\t' + 'regression_snp_summary_file\tlow_dimensional_squared_ld_file\tn_snps_per_low_dim_file\tregression_snp_ld_file\n')
		continue

	# Extract relevent fields
	gene_name = data[1]
	print(gene_name)
	gene_tss = int(data[2])
	gene_cis_start = int(data[4])
	gene_cis_end = int(data[5])

	# Extract indices of cis snps for gene
	cis_snp_indices = (snp_pos >= gene_cis_start) & (snp_pos <= gene_cis_end)

	# Skip genes that have no cis snps
	if np.sum(cis_snp_indices) < args.min_cis_snps:
		print('skipped gene because too few cis snps')
		continue


	# Get cis snps
	cis_window_pos = snp_pos[cis_snp_indices]
	cis_window_rsids = ordered_rsids[cis_snp_indices]
	cis_window_cm = snp_cm[cis_snp_indices]
	cis_window_snp_integers = snp_integers[cis_snp_indices]  # Used getting dosage

	# Now get ref positions of gene snp id lb and ub
	cis_window_cm_lb = np.min(cis_window_cm)
	cis_window_cm_ub = np.max(cis_window_cm)

	# Gene window cm lb and ub
	big_window_cm_lb = cis_window_cm_lb - cm_window_size
	big_window_cm_ub = cis_window_cm_ub + cm_window_size

	# Extract indices of "big window" for gene
	big_window_snp_indices = (snp_cm >= big_window_cm_lb) & (snp_cm <= big_window_cm_ub)

	# Extract snps in big window
	big_window_pos = snp_pos[big_window_snp_indices]
	big_window_rsids = ordered_rsids[big_window_snp_indices]
	big_window_cm = snp_cm[big_window_snp_indices]
	big_window_snp_integers = snp_integers[big_window_snp_indices]  # Used getting dosage
	# Extract which snps are hm3 snps
	big_window_hm3_boolean = []
	for rsid in big_window_rsids:
		if rsid in hm3_rsid_to_regression_index:
			big_window_hm3_boolean.append(True)
		else:
			big_window_hm3_boolean.append(False)
	big_window_hm3_boolean = np.asarray(big_window_hm3_boolean)
	# Extract cis snps in big window
	big_window_cis_snp_indices = cis_snp_indices[big_window_snp_indices]

	if np.sum(big_window_hm3_boolean) < args.min_regr_snps:
		print('skipped gene because too few regression snps')
		continue

	# Construct LD
	big_window_a0s = []
	big_window_a1s = []
	big_window_geno_mat = []
	for snp_integer in big_window_snp_integers:
		var = genotype_obj[snp_integer]
		#dosage = var.minor_allele_dosage # Note, not standardized but ok for ld
		dosage = var.alt_dosage # Note, not standardized but ok for ld
		dosage = (dosage - np.mean(dosage))/np.std(dosage)
		big_window_geno_mat.append(dosage)
		big_window_a0s.append(var.alleles[0])
		big_window_a1s.append(var.alleles[1]) # THIS IS ALT ALLELE! Note: for some reason this is flipped with the .pvar file, but at least the flip is consistent. So all should be good.

	big_window_geno_mat = np.asarray(big_window_geno_mat)
	geno_ss = big_window_geno_mat.shape[1]
	big_window_a0s = np.asarray(big_window_a0s)
	big_window_a1s = np.asarray(big_window_a1s)
	LD = np.corrcoef(big_window_geno_mat)

	# Need to return a file summarizing regression snps (one line per snp include snp name, but also cis indicator as well as regression index)
	# Make gene regression snp summmary file
	gene_regression_snp_summary_file = args.gene_ld_output_root + '_' + gene_name + '_regression_snp_summary.txt'
	make_gene_regression_snp_summary_file(gene_regression_snp_summary_file, args.chrom, big_window_rsids[big_window_hm3_boolean], big_window_pos[big_window_hm3_boolean], big_window_cm[big_window_hm3_boolean], big_window_a0s[big_window_hm3_boolean], big_window_a1s[big_window_hm3_boolean], big_window_cis_snp_indices[big_window_hm3_boolean], hm3_rsid_to_regression_index)

	# Need to save a npy file with regression snp LD
	gene_regression_snp_LD_file = args.gene_ld_output_root + '_' + gene_name + '_regression_snp_LD.npy'
	np.save(gene_regression_snp_LD_file, LD[big_window_hm3_boolean, :][:, big_window_hm3_boolean])


	# Need to return regression snp X latent snp squared ld matrix
	if args.gene_snp_representation.startswith('pca'):
		rho_thresh = float('.' + args.gene_snp_representation.split('_')[1])
		# Subset LD to cis snps
		gene_cis_snp_ld = LD[:, big_window_cis_snp_indices][big_window_cis_snp_indices,:]
		lambdas_full, U_full = np.linalg.eig(gene_cis_snp_ld)
		non_negative_components = lambdas_full > 0.0
		lambdas = lambdas_full[non_negative_components]
		U = U_full[:, non_negative_components]
		real_components = np.iscomplex(lambdas) == False
		lambdas = lambdas[real_components]
		U = U[:, real_components]
		if np.sum(np.iscomplex(lambdas)) > 0:
			print('assumption eroror')
			pdb.set_trace()
		lambdas = lambdas.astype(float)
		U = U.astype(float)
		lambda_thresh = compute_lambda_thresh(lambdas, rho_thresh)
		thresh_components = lambdas >= lambda_thresh
		lambdas = lambdas[thresh_components]
		U = U[:, thresh_components]

		# PC sample loadings
		Z = np.dot(np.transpose(U), big_window_geno_mat[big_window_cis_snp_indices, :])

		# squared PC correlations
		squared_ld_pc_snps = np.square(pairwise_correlations(Z, big_window_geno_mat[big_window_hm3_boolean,:]))

		# Save low dimensional ld file
		regression_snp_low_dimensional_snp_squared_ld_file = args.gene_ld_output_root + '_' + gene_name + '_regression_snp_by_low_dimensional_squared_ld_mat.npy'
		np.save(regression_snp_low_dimensional_snp_squared_ld_file, squared_ld_pc_snps)

		# number of low dimensions
		n_low_dimensions = squared_ld_pc_snps.shape[1]
		n_snps_per_low_dim = np.ones(n_low_dimensions)
		# Save n_snps per low dimension
		n_snps_per_low_dim_file = args.gene_ld_output_root + '_' + gene_name + '_n_snps_per_low_dimension.npy'
		np.save(n_snps_per_low_dim_file, n_snps_per_low_dim)
	elif args.gene_snp_representation.startswith('bins'):
		# Create LD matrix that is regression snp by cis snp
		small_ld_sq = np.square(LD[big_window_hm3_boolean, :][:, big_window_cis_snp_indices])
		squared_adj_LD = small_ld_sq - ((1.0-small_ld_sq)/(geno_ss-2.0))
		n_bins = int(args.gene_snp_representation.split('_')[1])
		n_cis_snps = np.sum(big_window_cis_snp_indices)
		temp_vec = np.arange(n_cis_snps)
		snp_chunks = np.array_split(temp_vec,n_bins)

		# Construct low dimensional ld score mat
		low_dimensional_ld_score_mat = []
		n_snps_per_low_dim = []
		for snp_chunk in snp_chunks:
			n_snps_per_low_dim.append(len(snp_chunk))
			low_dimensional_ld_score_mat.append(np.sum(squared_adj_LD[:, snp_chunk],axis=1))
		n_snps_per_low_dim = np.asarray(n_snps_per_low_dim)
		low_dimensional_ld_score_mat = np.transpose(np.asarray(low_dimensional_ld_score_mat))

		# Save low dimensional ld file
		regression_snp_low_dimensional_snp_squared_ld_file = args.gene_ld_output_root + '_' + gene_name + '_regression_snp_by_low_dimensional_squared_ld_mat.npy'
		np.save(regression_snp_low_dimensional_snp_squared_ld_file, low_dimensional_ld_score_mat)

		# Save n_snps per low dimension
		n_snps_per_low_dim_file = args.gene_ld_output_root + '_' + gene_name + '_n_snps_per_low_dimension.npy'
		np.save(n_snps_per_low_dim_file, n_snps_per_low_dim)


	t.write(line + '\t' + gene_regression_snp_summary_file + '\t' + regression_snp_low_dimensional_snp_squared_ld_file + '\t' + n_snps_per_low_dim_file + '\t' + gene_regression_snp_LD_file + '\n')




t.close()
f.close()














import numpy as np
import os
import sys
import pdb
import time
import argparse
from bgen import BgenReader
import gzip


def create_dictionary_list_of_hm3_rsids(hm3_rs_id_file):
	dicti = {}
	f = open(hm3_rs_id_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		dicti[line] = 1
	f.close()
	return dicti


def extract_array_of_chromosomes(chromosome_file):
	if chromosome_file == 'None':
		chrom_arr = np.arange(1,23)
	else:
		chrom_arr = []
		f = open(chromosome_file)
		for line in f:
			line = line.rstrip()
			chrom_arr.append(int(line))
		chrom_arr = np.asarray(chrom_arr)
		f.close()
	return chrom_arr

def create_dictionary_mapping_from_rsid_to_cm_position(sldsc_annotation_file):
	cm_dicti = {}
	head_count = 0
	f = gzip.open(sldsc_annotation_file)
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[2]
		cm = float(data[3])
		if rsid in cm_dicti:
			print('assumption erororr')
			pdb.set_trace()
		cm_dicti[rsid] = cm
	f.close()
	return cm_dicti

def create_causal_eqtl_effect_mapping(lasso_summary_file, tmp_chrom, rsid_to_sdev, dont_standardize_snp_effects_bool):
	dicti = {}
	head_count = 0
	genes = []
	f = open(lasso_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		chrom_num = data[1]
		if chrom_num != tmp_chrom:
			continue
		rsid = data[2]
		effect = float(data[6])
		genes.append(gene_id)

		if dont_standardize_snp_effects_bool == False:
			snp_sdev = rsid_to_sdev[rsid]
			effect = effect*snp_sdev

		if gene_id not in dicti:
			dicti[gene_id] = {}
			dicti[gene_id]['rsids'] = []
			dicti[gene_id]['effects'] = []
		dicti[gene_id]['rsids'].append(rsid)
		dicti[gene_id]['effects'].append(effect)

	return dicti, np.unique(genes)

def load_in_ref_alt_allele_arr(pvar_file):
	ref_alt_alleles = []
	f = open(pvar_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 5:
			print('assumptino erooror')
			pdb.set_trace()
		ref_allele = data[3]
		alt_allele = data[4]
		if ref_allele == alt_allele:
			print('assumptino eororor')
			pdb.set_trace()
		ref_alt_alleles.append((ref_allele, alt_allele))
	f.close()
	return ref_alt_alleles


def extract_snp_cm_pos(ordered_rsids, rsid_to_cm):
	cm_arr = []
	for rsid in ordered_rsids:
		cm_arr.append(rsid_to_cm[rsid])
	return np.asarray(cm_arr)

def create_mapping_from_rsid_to_snp_index(G_obj_rsids):
	dicti = {}
	for ii, val in enumerate(G_obj_rsids):
		dicti[val] = ii
	return dicti

def create_mapping_from_gene_name_to_tss(filer):
	f = open(filer)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count== 0:
			head_count = head_count + 1
			continue
		ens_id = data[1]
		tss = int(data[2])
		if ens_id in dicti:
			print('assumption erororo')
			pdb.set_trace()
		dicti[ens_id] = tss
	f.close()
	return dicti

def get_boolean_vector_of_hm3_indices(ordered_rsids, hm3_rsids):
	booler = []
	for rsid in ordered_rsids:
		if rsid in hm3_rsids:
			booler.append(True)
		else:
			booler.append(False)
	return np.asarray(booler)


def load_in_alt_allele_genotype_dosage_mat(bfile, window_indices):
	dosages = []

	for window_index in window_indices:
		var = bfile[window_index]
		dosage = var.alt_dosage
		#ma = var.minor_allele

		# Append snp dosage to global array
		dosages.append(dosage)

	# Convert to 2d matrix
	dosages = np.asarray(dosages)

	return np.transpose(dosages)


def correlation_matrix(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
	"""
	Compute the (K × L) matrix of Pearson correlations between columns of X (shape N×K)
	and columns of Y (shape N×L).

	Parameters
	----------
	X : np.ndarray
		A 2D array of shape (N, K).
	Y : np.ndarray
		A 2D array of shape (N, L).

	Returns
	-------
	corr : np.ndarray
		A 2D array of shape (K, L) where corr[i, j] is the Pearson correlation
		between X[:, i] and Y[:, j].

	Raises
	------
	ValueError
		If X and Y do not have the same number of rows.
	"""
	if X.ndim != 2 or Y.ndim != 2:
		raise ValueError("Both X and Y must be 2D arrays.")
	N, K = X.shape
	Ny, L = Y.shape
	if N != Ny:
		raise ValueError(f"X has {N} rows but Y has {Ny} rows; they must match.")

	# 1. Center each column
	Xc = X - X.mean(axis=0, keepdims=True)   # shape (N, K)
	Yc = Y - Y.mean(axis=0, keepdims=True)   # shape (N, L)

	# 2. Compute cross‐covariance matrix (K × L)
	cov_XY = (Xc.T @ Yc) / (N - 1)           # shape (K, L)

	# 3. Compute column‐wise standard deviations (ddof=1 to match / (N-1))
	std_X = Xc.std(axis=0, ddof=1)           # shape (K,)
	std_Y = Yc.std(axis=0, ddof=1)           # shape (L,)

	# 4. Form outer product of std_X and std_Y, then normalize
	denom = np.outer(std_X, std_Y)           # shape (K, L)
	corr = cov_XY / denom

	return corr

def load_in_alt_allele_sdevs(bfile, window_indices, rsids):
	sdevs = []
	dicti = {}
	for ii,window_index in enumerate(window_indices):
		var = bfile[window_index]
		dosage = var.alt_dosage
		#ma = var.minor_allele


		# Compute standard deviation
		rsid = rsids[ii]
		sdev = np.std(dosage)
		if rsid in dicti:
			print('assumption erororo')
			pdb.set_trace()
		dicti[rsid] = sdev

	return dicti


#####################
# Parse command line arguments
#####################
parser = argparse.ArgumentParser()

parser.add_argument('--study-names-file', default='None', type=str,
					help='File containing list of studies to run analysis for')
parser.add_argument('--bgen', default='None', type=str,
					help='Bgen file')
parser.add_argument('--chromosome-file', default='None', type=str,
					help='Filename containing list of chromosomes to run (no header). If None, then defaults to 22 chromosomes')
parser.add_argument('--hm3-file-stem', default='None', type=str,
					help='Filename containing stem to files containing hm3 snps')
parser.add_argument('--ldsc-annotation-file-stem', default='None', type=str,
					help='Filename containing stem to ldsc annotation files')
parser.add_argument('--gene-summary-file-stem', default='None', type=str,
					help='Filename containing stem to gene-to-tss mapping')
parser.add_argument('--gene-summary-file-suffix', default='None', type=str,
					help='Filename containing suffix to gene-to-tss mapping')
parser.add_argument('--cis-window', default=500000.0, type=float,
					help='BP region around gene to call cis-eqtls')
parser.add_argument('--cm-window', default=1.0, type=float,
					help='CM window used to cutoff LD')
parser.add_argument('--standardize', default=False, action='store_true',
					help='Boolean on whether or not to standardize cis predicted gene expression')
parser.add_argument('--dont-standardize-snp-effects', default=False, action='store_true',
					help='Boolean on whether or not to standardize snp effects or not')
parser.add_argument('--output-dir', default=None, type=str,
					help='Output directory stem to save data to')
args = parser.parse_args()




# Extract array of chromosomes
chromosome_arr = extract_array_of_chromosomes(args.chromosome_file)



# Open file handles for each study
t = {}
study_names = []
study_lasso_files = []
head_count = 0
f = open(args.study_names_file)
for line in f:
	# Skip header
	if head_count == 0:
		head_count = head_count + 1
		continue
	line = line.rstrip()
	data = line.split('\t')
	study_name = data[0]
	study_names.append(study_name)
	study_lasso_files.append(data[1])
	if args.standardize:
		output_file = args.output_dir + study_name + '_standardized_pred_eqtl_sumstats.txt'
	else:
		output_file = args.output_dir + study_name + '_pred_eqtl_sumstats.txt'
	t[study_name] = open(output_file,'w')
	t[study_name].write('gene_id\trsid\tchr\tpos\ta1\ta2\tbeta\tbeta_se\tz\tin_sample_sdev\n')
f.close()
study_names = np.asarray(study_names)
study_lasso_files = np.asarray(study_lasso_files)

# Loop through chromosomes
for chrom_num in chromosome_arr:

	# Create dictionary list of hm3rsids
	hm3_rsid_file = args.hm3_file_stem + str(chrom_num) + '.txt'
	hm3_rsids = create_dictionary_list_of_hm3_rsids(hm3_rsid_file)

	# Create dictionary mapping from rsid to CM position
	ldsc_annotation_file = args.ldsc_annotation_file_stem + str(chrom_num) + '.annot.gz'
	rsid_to_cm = create_dictionary_mapping_from_rsid_to_cm_position(ldsc_annotation_file)

	# Load in genotype object
	ref_alt_alleles = load_in_ref_alt_allele_arr(args.bgen + str(chrom_num) + '.pvar')
	ref_alt_alleles = np.asarray(ref_alt_alleles)
	bgen_file = args.bgen + str(chrom_num) + '.bgen'
	genotype_obj = BgenReader(bgen_file)
	snp_pos = np.asarray(genotype_obj.positions())
	ordered_rsids = np.asarray(genotype_obj.rsids())
	snp_cm = extract_snp_cm_pos(ordered_rsids, rsid_to_cm)
	snp_integers = np.arange(len(ordered_rsids))
	rsid_to_sdev = load_in_alt_allele_sdevs(genotype_obj, snp_integers, ordered_rsids)


	# Create causal eqtl effect mapping
	lasso_effect_dicti = {}
	all_genes = []
	for study_iter, study_name in enumerate(study_names):
		study_lasso_file = study_lasso_files[study_iter]
		
		study_dicti, study_genes = create_causal_eqtl_effect_mapping(study_lasso_file, str(chrom_num), rsid_to_sdev, args.dont_standardize_snp_effects)
		lasso_effect_dicti[study_name] = study_dicti
		all_genes.append(study_genes)
	all_genes = np.unique(np.hstack(all_genes))

	# Create mapping from gene name to tss
	gene_name_to_tss = create_mapping_from_gene_name_to_tss(args.gene_summary_file_stem + str(chrom_num) + args.gene_summary_file_suffix)


	# Get boolean vector of hm3 indices
	hm3_boolean = get_boolean_vector_of_hm3_indices(ordered_rsids, hm3_rsids)

	# Create mapping from rsid to snp index
	rsid_to_snp_index = create_mapping_from_rsid_to_snp_index(ordered_rsids)

	# Now loop through genes
	for ii,gene_id in enumerate(all_genes):
		print(ii)
		print(gene_id)

		# TSS of this gene
		gene_tss = gene_name_to_tss[gene_id]

		# Extract cis regression snps
		cis_start_pos = gene_tss - args.cis_window
		cis_end_pos = gene_tss + args.cis_window
		cis_regr_snp_indices = (snp_pos >= cis_start_pos) & (snp_pos <= cis_end_pos) & (hm3_boolean)
		cis_rsids = ordered_rsids[cis_regr_snp_indices]
		regression_snp_indices_raw = snp_integers[cis_regr_snp_indices]
		regression_snp_genotype_dosage = load_in_alt_allele_genotype_dosage_mat(genotype_obj, regression_snp_indices_raw)
		regression_snp_alleles = ref_alt_alleles[regression_snp_indices_raw,:]
		regression_snp_pos = snp_pos[regression_snp_indices_raw]

		# Now loop through studies
		for study_name in study_names:
			# Ignore studies for which gene has no estimated effect
			if gene_id not in lasso_effect_dicti[study_name]:
				continue

			# Gene has estimated causal effect in this study
			causal_eqtl_rsids = np.asarray(lasso_effect_dicti[study_name][gene_id]['rsids'])
			causal_eqtl_effects = np.asarray(lasso_effect_dicti[study_name][gene_id]['effects'])


			causal_snp_indices = []
			for rsid in causal_eqtl_rsids:
				causal_snp_indices.append(rsid_to_snp_index[rsid])
			causal_snp_indices = np.asarray(causal_snp_indices)


			causal_snp_genotype_dosage = load_in_alt_allele_genotype_dosage_mat(genotype_obj, causal_snp_indices)

			R_sub = correlation_matrix(regression_snp_genotype_dosage, causal_snp_genotype_dosage)

			# Standardize
			tiny_R = np.corrcoef(np.transpose(causal_snp_genotype_dosage))
			gene_var = np.dot(np.dot(causal_eqtl_effects, tiny_R), causal_eqtl_effects)

			if args.standardize:
				std_causal_eqtl_effects = causal_eqtl_effects/np.sqrt(gene_var)
				marginal_effect_est = np.dot(R_sub, std_causal_eqtl_effects)
			else:
				marginal_effect_est = np.dot(R_sub, causal_eqtl_effects)

			for regression_snp_index, regression_snp in enumerate(cis_rsids):
				t[study_name].write(gene_id + '\t' + regression_snp + '\t' + str(chrom_num) + '\t' + str(regression_snp_pos[regression_snp_index]) + '\t' + str(regression_snp_alleles[regression_snp_index,1]) + '\t' + str(regression_snp_alleles[regression_snp_index,0]) + '\t' + str(marginal_effect_est[regression_snp_index]) + '\t' + '0.0' + '\t' + 'NA' + '\t' + 'NA\n')





# Close output file handles
for study_name in study_names:
	t[study_name].close()
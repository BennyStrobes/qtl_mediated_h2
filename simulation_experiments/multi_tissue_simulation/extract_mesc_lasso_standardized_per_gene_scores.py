import numpy as np
import os
import sys
import pdb
from bgen import BgenReader
import gzip
import time



def create_dictionary_list_of_hm3_rsids(hm3_rs_id_file):
	dicti = {}
	f = open(hm3_rs_id_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		dicti[line] = 1
	f.close()
	return dicti

def create_dictionary_mapping_from_rsid_to_cm_position(bim_file):
	dicti = {}
	f = open(bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		rsid =data[1]
		if rsid in dicti:
			print('assumption error')
			pdb.set_trace()
		cm_position = float(data[2])
		dicti[rsid] = cm_position
	f.close()
	return dicti

def create_dictionary_mapping_from_rsid_to_cm_position_and_annotation_vector(sldsc_annotation_file, invalid_annotations):
	cm_dicti = {}
	anno_dicti = {}
	head_count = 0
	f = gzip.open(sldsc_annotation_file)
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			anno_names = np.asarray(data[4:])
			valid_anno_indices = []
			for ii, anno_name in enumerate(anno_names):
				if anno_name not in invalid_annotations:
					valid_anno_indices.append(ii)
			valid_anno_indices = np.asarray(valid_anno_indices)
			anno_names = anno_names[valid_anno_indices]
			continue
		rsid = data[2]
		cm = float(data[3])
		anno_vec = np.asarray(data[4:]).astype(float)

		if rsid in cm_dicti or rsid in anno_vec:
			print('assumption erororr')
			pdb.set_trace()
		cm_dicti[rsid] = cm
		anno_dicti[rsid] = anno_vec[valid_anno_indices]

	f.close()
	return cm_dicti, anno_dicti, anno_names


def extract_snp_cm_pos(ordered_rsids, rsid_to_cm):
	cm_arr = []
	for rsid in ordered_rsids:
		cm_arr.append(rsid_to_cm[rsid])
	return np.asarray(cm_arr)

def create_causal_eqtl_effect_mapping(lasso_summary_file, tmp_chrom):
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
		effect = float(data[3])
		genes.append(gene_id)

		if gene_id not in dicti:
			dicti[gene_id] = {}
			dicti[gene_id]['rsids'] = []
			dicti[gene_id]['effects'] = []
		dicti[gene_id]['rsids'].append(rsid)
		dicti[gene_id]['effects'].append(effect)

	return dicti, np.unique(genes)

def create_mapping_from_rsid_to_snp_index(G_obj_rsids):
	dicti = {}
	for ii, val in enumerate(G_obj_rsids):
		dicti[val] = ii
	return dicti

def extract_regression_snp_indices(regression_snp_summary_file, rsid_to_snp_index):
	f = open(regression_snp_summary_file)
	indices = []
	snp_names = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[1]
		indices.append(rsid_to_snp_index[rsid])
		snp_names.append(rsid)
	f.close()
	return np.asarray(indices), np.asarray(snp_names)



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

####################
# Command line args
####################
simulation_number = sys.argv[1]
chrom_string = sys.argv[2]
simulation_name_string = sys.argv[3]
simulation_genotype_dir = sys.argv[4]
mesc_expression_score_dir = sys.argv[5]

eqtl_ss_arr = ['100', '200', '300', '1000']

# Create causal eqtl effect mapping
t = {}
study_names = []
for tissue_iter in range(5):
	for replicate_name in ['replicate1', 'replicate2']:
		for eqtl_ss in eqtl_ss_arr:
			study_name = 'tissue' + str(tissue_iter) + '_' + replicate_name + '_' + eqtl_ss
			study_names.append(study_name)

			output_file = tissue_sumstat1_output_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_iter) + '_' + str(eqtl_ss) + '_' + replicate_name + '_mesc_lasso_standardized_pred_eqtl_sumstats.txt'
			t[study_name] = open(output_file,'w')
			t[study_name].write('gene_id\trsid\tchr\tpos\ta1\ta2\tbeta\tbeta_se\tz\tin_sample_sdev\n')

study_names = np.asarray(study_names)


cm_window_size = 1.0

for chrom_num in range(1,3):
	hm3_rsid_file = '/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/genotype_processing/' + 'hm3_rsids_chr' + str(chrom_num) + '.txt'

	# Create dictionary list of hm3rsids
	hm3_rsids = create_dictionary_list_of_hm3_rsids(hm3_rsid_file)

	ldsc_annotation_file = '/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/' + 'baselineLD.' + str(chrom_num) + '.annot.gz'
	invalid_annotations = {}
	invalid_annotations['GTEx_eQTL_MaxCPP'] = 1
	invalid_annotations['BLUEPRINT_H3K27acQTL_MaxCPP'] = 1
	invalid_annotations['BLUEPRINT_H3K4me1QTL_MaxCPP'] = 1
	invalid_annotations['BLUEPRINT_DNA_methylation_MaxCPP'] = 1
	rsid_to_cm, rsid_to_anno, anno_names = create_dictionary_mapping_from_rsid_to_cm_position_and_annotation_vector(ldsc_annotation_file, invalid_annotations)

	# Create causal eqtl effect mapping
	lasso_effect_dicti = {}
	all_genes = []
	for tissue_iter in range(5):
		for replicate_name in ['replicate1', 'replicate2']:
			for eqtl_ss in eqtl_ss_arr:
				lasso_summary_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_iter) + '_' + str(eqtl_ss) + '_' + replicate_name + '_' + str(chrom_num) + '.' + str(chrom_num) + '.lasso'
				study_dicti, study_genes = create_causal_eqtl_effect_mapping(lasso_summary_file, str(chrom_num))
				study_name = 'tissue' + str(tissue_iter) + '_' + replicate_name + '_' + eqtl_ss
				lasso_effect_dicti[study_name] = study_dicti
				all_genes.append(study_genes)
	all_genes = np.unique(np.hstack(all_genes))

	# Load in genotype object
	ref_alt_alleles = load_in_ref_alt_allele_arr(simulation_genotype_dir + 'simulated_reference_genotype_data_' + str(chrom_num) + '.pvar')
	bgen_file = simulation_genotype_dir + 'simulated_reference_genotype_data_' + str(chrom_num) + '.bgen'
	genotype_obj = BgenReader(bgen_file)
	snp_pos = np.asarray(genotype_obj.positions())
	ordered_rsids = np.asarray(genotype_obj.rsids())
	snp_cm = extract_snp_cm_pos(ordered_rsids, rsid_to_cm)
	snp_integers = np.arange(len(ordered_rsids))

	# Create mapping from rsid to snp index
	rsid_to_snp_index = create_mapping_from_rsid_to_snp_index(ordered_rsids)


	# Now loop through genes
	for ii,gene_id in enumerate(all_genes):
		print(ii)
		print(gene_id)

		regression_snp_summary_file = simulation_genotype_dir + 'gene_level_ld_chr' + str(chrom_num) + '_bins_10_' + gene_id + '_regression_snp_summary.txt'
		regression_snp_indices_raw, regression_snps = extract_regression_snp_indices(regression_snp_summary_file, rsid_to_snp_index)
		regression_snp_genotype_dosage = load_in_alt_allele_genotype_dosage_mat(genotype_obj, regression_snp_indices_raw)
		regression_snp_alleles = np.asarray(ref_alt_alleles)[regression_snp_indices_raw,:]
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
			std_causal_eqtl_effects = causal_eqtl_effects/np.sqrt(gene_var)


			marginal_effect_est = np.dot(R_sub, std_causal_eqtl_effects)

			for regression_snp_index, regression_snp in enumerate(regression_snps):
				t[study_name].write(gene_id + '\t' + regression_snp + '\t' + str(chrom_num) + '\t' + str(regression_snp_pos[regression_snp_index]) + '\t' + str(regression_snp_alleles[regression_snp_index,1]) + '\t' + str(regression_snp_alleles[regression_snp_index,0]) + '\t' + str(marginal_effect_est[regression_snp_index]) + '\t' + '0.0' + '\t' + 'NA' + '\t' + 'NA\n')


for study_name in study_names:
	t[study_name].close()

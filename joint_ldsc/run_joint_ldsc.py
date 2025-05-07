import numpy as np
import os
import sys
import pdb
import time
import pickle
import joint_ldsc
import argparse



def load_in_gwas_data(gwas_summary_file, chromosome_dictionary):
	rsids = []
	betas = []
	beta_ses = []

	head_count = 0

	f = open(gwas_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count = head_count + 1
			continue
		# Filter out variants not in chromosome list
		line_chrom_num = data[1]
		if line_chrom_num not in chromosome_dictionary:
			print('filtered snp')
			continue

		rsids.append(data[0])
		betas.append(float(data[5]))
		beta_ses.append(float(data[6]))

	f.close()

	return np.asarray(rsids), np.asarray(betas), np.asarray(beta_ses)

def load_in_regression_snps_for_gene(regression_snp_file):
	f = open(regression_snp_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[1])
	f.close()
	return np.asarray(arr)

def load_in_eqtl_data(rsid_to_position, gene_ldscore_filestem, gene_ldscore_filesuffix, eqtl_dataset_names, chrom_arr):
	genes = []
	gene_info = {}

	for chrom_num in chrom_arr:
		gene_summary_file = gene_ldscore_filestem + str(chrom_num) + gene_ldscore_filesuffix
		f = open(gene_summary_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue

			# Each line is a gene. Extract relevent fields
			ensamble_id = data[1]
			genes.append(ensamble_id)
			squared_ld_mat_file = data[7]
			n_low_dim_snp_file = data[8]
			regression_snp_file = data[6]
			if ensamble_id in gene_info:
				print('assumption eroror')
				pdb.set_trace()
			# Load in regression snps
			gene_regression_snps = load_in_regression_snps_for_gene(regression_snp_file)
			# Get corresponding indices of regression snps
			regression_snp_indices = []
			for rsid in gene_regression_snps:
				regression_snp_indices.append(rsid_to_position[rsid])
			regression_snp_indices = np.asarray(regression_snp_indices)
			# Create mapping from rsid to gene position
			rsid_to_gene_position = {}
			for ii, val in enumerate(gene_regression_snps):
				rsid_to_gene_position[val] = ii


			gene_info[ensamble_id] = {}
			gene_info[ensamble_id]['squared_ld_file'] = squared_ld_mat_file
			gene_info[ensamble_id]['n_low_dimensional_snps'] = np.load(n_low_dim_snp_file)
			gene_info[ensamble_id]['regression_snp_indices'] = regression_snp_indices
			gene_info[ensamble_id]['rsid_to_gene_position'] = rsid_to_gene_position
			gene_info[ensamble_id]['n_gene_regression_snps'] = len(gene_regression_snps)
			gene_info[ensamble_id]['squared_sumstats'] = np.zeros((len(gene_regression_snps), len(eqtl_dataset_names)))
			gene_info[ensamble_id]['eqtl_category_names'] = np.copy(eqtl_dataset_names)
		f.close()


	return np.asarray(genes), gene_info


def get_gwas_variant_ld_scores(variant_ldscore_filestem, variant_ldscore_filesuffix, chrom_arr, non_mediated_annotation_version):

	rsids = []
	ldscores = []

	for chrom_num in chrom_arr:
		filer = variant_ldscore_filestem + chrom_num + variant_ldscore_filesuffix
		f = open(filer)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				anno_names = np.asarray(data[6:])
				continue
			rsid = data[1]
			ldscore = np.asarray(data[6:]).astype(float)
			rsids.append(rsid)
			ldscores.append(ldscore)
		f.close()
	ldscores = np.asarray(ldscores)

	if non_mediated_annotation_version == 'genotype_intercept':
		ldscores = ldscores[:,:1]
		anno_names = anno_names[:1]

	return ldscores, np.asarray(rsids), anno_names


def create_mapping_from_rsid_to_position(rsids):
	rsid_to_position = {}

	for ii, rsid in enumerate(rsids):
		if rsid in rsid_to_position:
			print('assumption eroror')
			pdb.set_trace()
		rsid_to_position[rsid] = ii
	return rsid_to_position

def fill_in_eqtl_sumstats(gene_info, eqtl_dataset_files, eqtl_dataset_names, genes, chrom_dicti):
	for ii, eqtl_dataset_name in enumerate(eqtl_dataset_names):
		eqtl_sumstat_file = eqtl_dataset_files[ii]
		f = open(eqtl_sumstat_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			line_chrom_num = data[2]
			if line_chrom_num not in chrom_dicti:
				print('skipped variant-gene because wrong chromosome')
				pdb.set_trace()
			ens_id = data[0]
			rsid = data[1]
			beta = float(data[6])
			beta_se = float(data[7])
			gene_row_index = gene_info[ens_id]['rsid_to_gene_position'][rsid]
			gene_info[ens_id]['squared_sumstats'][gene_row_index, ii] = np.square(beta) - np.square(beta_se)
		f.close()

	# Filter out missing eqtl categories for each gene
	for ens_id in genes:
		valid_categories = []
		for ii, category_name in enumerate(gene_info[ens_id]['eqtl_category_names']):
			gene_info[ens_id]['squared_sumstats'][:, ii]
			if np.array_equal(gene_info[ens_id]['squared_sumstats'][:, ii], np.zeros(len(gene_info[ens_id]['squared_sumstats'][:, ii]))):
				valid_categories.append(False)
			else:
				valid_categories.append(True)
		valid_categories = np.asarray(valid_categories)
		gene_info[ens_id]['squared_sumstats'] = gene_info[ens_id]['squared_sumstats'][:,valid_categories]
		gene_info[ens_id]['eqtl_category_names'] = gene_info[ens_id]['eqtl_category_names'][valid_categories]
	return gene_info

def extract_number_of_reference_snps(variant_M_filestem, chrom_arr, non_mediated_annotation_version):
	n_snps = []
	for chrom_num in chrom_arr:
		filer = variant_M_filestem + str(chrom_num) + '_M.txt'
		f = open(filer)
		head_count = 0
		tmp_arr = []
		for line in f:
			line = line.rstrip()
			if head_count == 0:
				head_count = head_count + 1
				continue
			data = line.split('\t')
			tmp_arr.append(float(data[1]))
		f.close()
		n_snps.append(np.asarray(tmp_arr))
	n_snps = np.asarray(n_snps)
	final_n_snps = np.sum(n_snps,axis=0)

	# Filter columns if using genotype intercept
	if non_mediated_annotation_version == 'genotype_intercept':
		final_n_snps = final_n_snps[:1]

	return final_n_snps

def extract_chromosome_names(chromosome_file):
	if chromosome_file is None:
		tmp_chrom_arr = np.arange(1,23)
		arr = []
		dicti = {}
		for ele in tmp_chrom_arr:
			arr.append(str(ele))
			dicti[str(ele)] = 1
	else:
		arr = []
		dicti = {}
		f = open(chromosome_file)
		for line in f:
			line = line.rstrip()
			arr.append(line)
			dicti[line] = 1
		f.close()
	return np.asarray(arr), dicti

def load_in_eqtl_dataset_summary_file(eqtl_summary_file):
	names = []
	sample_sizes = []
	filers = []
	f = open(eqtl_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		names.append(data[0])
		sample_sizes.append(float(data[1]))
		filers.append(data[2])
	f.close()
	return np.asarray(names), np.asarray(sample_sizes), np.asarray(filers)



#####################
# Parse command line arguments
#####################
parser = argparse.ArgumentParser()
parser.add_argument('--gwas-sumstat-file', default='None', type=str,
                    help='Full path to GWAS summary stat file')
parser.add_argument('--gwas-N', default='None', type=float,
                    help='GWAS sample size')
parser.add_argument('--variant-ldscore-filestem', default='None', type=str,
                    help='Per chromosome variant-ldscore-filestem')
parser.add_argument('--variant-ldscore-filesuffix', default='.txt', type=str,
                    help='Per chromosome variant-ldscore-file suffix')
parser.add_argument('--variant-M-filestem', default='None', type=str,
                    help='Per chromosome variant-M file filestems')
parser.add_argument('--gene-ldscore-filestem', default='None', type=str,
                    help='Per chromosome gene-ldscore-filestem')
parser.add_argument('--gene-ldscore-filesuffix', default='.txt', type=str,
                    help='Per chromosome gene-ldscore-file suffix')
parser.add_argument('--chromosome-file', default=None, type=str,
                    help='File containing chromosomes to run analysis on. If None, set to autosomal chromosomes')
parser.add_argument('--eqtl-summary-file', default=None, type=str,
                    help='File containing one line for each eqtl data set')
parser.add_argument('--non-mediated-annotation-version', default="full_anno", type=str,
                    help='Takes as values "genotype_intercept" or "full_anno"')
parser.add_argument('--gene-trait-architecture', default="stdExpr", type=str,
                    help='Takes as values "linear" or "stdExpr"')
parser.add_argument('--jacknife', default=False, action='store_true',
                    help='Boolean on whether or not to Jacknife the estimtes')
parser.add_argument('--fixed-variance', default=False, action='store_true',
                    help='Boolean on whether the variance parameters should be fixed')
parser.add_argument('--output-stem', default=None, type=str,
                    help='Output file stem to save data to')
args = parser.parse_args()


#####################
# Load in relevent data
#####################
# Names of chromosomes to run analysis on
chrom_arr, chrom_dicti = extract_chromosome_names(args.chromosome_file)

# Extract number of reference snps
n_reference_snps = extract_number_of_reference_snps(args.variant_M_filestem, chrom_arr, args.non_mediated_annotation_version)

# Load in GWAS summary statistics
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(args.gwas_sumstat_file, chrom_dicti)
gwas_E_beta_sq = np.square(gwas_beta) - np.square(gwas_beta_se)


# Load in GWAS variant ld scores
gwas_variant_ld_scores, gwas_rsids_tmp, non_mediated_annotation_names = get_gwas_variant_ld_scores(args.variant_ldscore_filestem, args.variant_ldscore_filesuffix, chrom_arr, args.non_mediated_annotation_version)

# QUick error checking
if np.array_equal(gwas_rsids, gwas_rsids_tmp) == False:
	print('assumption eroror')
	pdb.set_trace()

# Create mapping from rsid to position
rsid_to_position = create_mapping_from_rsid_to_position(gwas_rsids)

# Load in eqtl data
eqtl_dataset_names, eqtl_dataset_Ns, eqtl_dataset_files = load_in_eqtl_dataset_summary_file(args.eqtl_summary_file)
genes, gene_info = load_in_eqtl_data(rsid_to_position, args.gene_ldscore_filestem, args.gene_ldscore_filesuffix, eqtl_dataset_names, chrom_arr)
gene_info = fill_in_eqtl_sumstats(gene_info, eqtl_dataset_files, eqtl_dataset_names, genes, chrom_dicti)

'''
#################
# Temp saving
f = open('gene_info.pickle', 'wb')
pickle.dump(gene_info, f)
f.close()

np.save('genes.npy', genes)
np.save('eqtl_dataset_names.npy', eqtl_dataset_names)
np.save('gwas_variant_ld_scores.npy', gwas_variant_ld_scores)
np.save('gwas_E_beta_sq.npy', gwas_E_beta_sq)
np.save('n_reference_snps.npy', n_reference_snps)
'''
'''
#################
# Temp loading
f = open('gene_info.pickle','rb')
gene_info = pickle.load(f)
f.close()

genes = np.load('genes.npy')
eqtl_dataset_names = np.load('eqtl_dataset_names.npy')
gwas_variant_ld_scores = np.load('gwas_variant_ld_scores.npy')
gwas_E_beta_sq = np.load('gwas_E_beta_sq.npy')
n_reference_snps = np.load('n_reference_snps.npy')
'''

#####################
# Open output file handle
#####################
output_file = args.output_stem + '_joint_ldsc_results.txt'
t = open(output_file,'w')
t.write('method\ttest_med_h2_\test_med_h2_per_tissue\test_nm_h2\test_mean_eqtl_h2\n')


#####################
# Run optimization
#####################
# Run two step joint ldsc
#obj = joint_ldsc.med_h2(version='two_step')
#obj.fit(genes, gene_info, gwas_variant_ld_scores, gwas_E_beta_sq, eqtl_dataset_names, n_reference_snps)
#t.write('two_step_joint_ldsc' + '\t' + str(np.sum(obj.med_h2)) + '\t' + ','.join(obj.med_h2.astype(str)) + '\t' + str(obj.nm_h2) + '\t' + str(obj.avg_eqtl_h2) + '\n')

# Run joint ldsc
obj = joint_ldsc.med_h2(max_iter=500,convergence_thresh=1e-10, fixed_variance=args.fixed_variance)
obj.fit(genes, gene_info, gwas_variant_ld_scores,gwas_E_beta_sq, eqtl_dataset_names, n_reference_snps)
t.write('joint_ldsc' + '\t' + str(np.sum(obj.med_h2)) + '\t' + ','.join(obj.med_h2.astype(str)) + '\t' + str(obj.nm_h2) + '\t' + str(obj.avg_eqtl_h2) + '\n')

t.close()
print(output_file)









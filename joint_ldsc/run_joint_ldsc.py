import numpy as np
import os
import sys
import pdb
import time
import pickle
import joint_ldsc
import joint_ldsc_gibbs
import argparse


def linear_regression(XX, YY, intercept=False):
	if intercept:
		print('not yet implemented. add column of 1s to X')
		pdb.set_trace()

	return np.dot(np.dot(np.linalg.pinv(np.dot(np.transpose(XX), XX)), np.transpose(XX)), YY)

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


def update_gene_info_to_include_eqtl_dataset_names(gene_info):
	# Get list of genes
	genes = np.asarray([*gene_info])
	# loop through genes
	for gene in genes:
		gene_info[gene]['eqtl_dataset_names'] = np.copy(gene_info[gene]['eqtl_category_names'])

	return gene_info

def bin_assignment(gene_cis_h2, dataset_bins, bin_names):
	n_bins = len(bin_names)
	bin_index_assignment = None
	for bin_iter in range(len(dataset_bins)-1):
		if gene_cis_h2 >= dataset_bins[bin_iter] and gene_cis_h2 < dataset_bins[(bin_iter+1)]:
			bin_index_assignment = bin_iter
	if bin_index_assignment is None:
		print('assumption error: gene not assigned to bin')
		pdb.set_trace()

	return bin_names[bin_index_assignment]

def update_gene_info_to_bin_genes_by_random(gene_info, eqtl_dataset_names, n_expr_cis_h2_bins):
	##############################
	# First estimate gene cis h2
	##############################
	# Get list of genes
	genes = np.asarray([*gene_info])

	bin_names = []
	for bin_iter in range(n_expr_cis_h2_bins):
		bin_names.append('bin' + str(bin_iter))
	bin_names = np.asarray(bin_names)

	# create list of new eqtl categories
	new_eqtl_categories = []
	for eqtl_dataset_name in eqtl_dataset_names:
		for bin_name in bin_names:
			new_eqtl_categories.append(eqtl_dataset_name + ':' + bin_name)
	new_eqtl_categories = np.asarray(new_eqtl_categories)


	##############################
	# Update gene_info file to now contain bin_info
	##############################
	# loop through genes
	for gene in genes:
		info = gene_info[gene]

		# Load in gene eqtl categories
		gene_eqtl_categories = info['eqtl_category_names']


		# Initialize vector to keep track of update gene_eqtl_categories
		updated_gene_eqtl_categories = []
		# Loop through eQTL categories
		for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):

			# Bin assignment
			bin_name = np.random.choice(bin_names)

			# Update gnee eqtl category
			updated_gene_eqtl_categories.append(eqtl_category_name + ':' + bin_name)
		# Put into nice array
		updated_gene_eqtl_categories = np.asarray(updated_gene_eqtl_categories)

		# Update gene info file
		gene_info[gene]['eqtl_dataset_names'] = np.copy(gene_info[gene]['eqtl_category_names'])
		gene_info[gene]['eqtl_category_names'] = updated_gene_eqtl_categories

	return gene_info, new_eqtl_categories



def update_gene_info_to_bin_genes_by_est_cis_h2(gene_info, eqtl_dataset_names, n_expr_cis_h2_bins):
	##############################
	# First estimate gene cis h2
	##############################
	# Get list of genes
	genes = np.asarray([*gene_info])
	# Mapping from dataset name to index
	dataset_name_to_index = {}
	dataset_name_to_cis_h2s = {}
	for ii, val in enumerate(eqtl_dataset_names):
		dataset_name_to_index[val] = ii
		dataset_name_to_cis_h2s[val] = []

	# Keep track of estimated heritability
	gene_cis_h2 = {}

	# loop through genes
	for gene in genes:
		info = gene_info[gene]
		# Load in squared ld 
		squared_ld = np.load(info['squared_ld_file'])

		# Initialize gene cis_h2 vec
		gene_cis_h2_vec = np.zeros(len(info['eqtl_category_names']))

		# Load in eqtl sumstats
		eqtl_sq_sumstats = info['squared_sumstats']

		# Load in gene eqtl categories
		gene_eqtl_categories = info['eqtl_category_names']

		# Loop through eQTL categories
		# Initialize seperately for each category
		for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):
			weights = linear_regression(squared_ld, eqtl_sq_sumstats[:, eqtl_category_iter])
			gene_est_cis_h2 = np.sum(weights*info['n_low_dimensional_snps'])
			#weights = linear_regression(np.sum(squared_ld,axis=1).reshape(-1,1), eqtl_sq_sumstats[:, eqtl_category_iter])
			#gene_est_cis_h2 = weights[0]*np.sum(info['n_low_dimensional_snps'])
			
			gene_cis_h2_vec[eqtl_category_iter] = gene_est_cis_h2
			dataset_name_to_cis_h2s[eqtl_category_name].append(gene_est_cis_h2)
		# Update dictionary
		gene_cis_h2[gene] = gene_cis_h2_vec

	##############################
	# Bin cis_h2s in each eqtl dataset
	##############################
	dataset_to_bins = {}
	quantiles = np.linspace(0, 1, n_expr_cis_h2_bins + 1)
	for eqtl_dataset_name in eqtl_dataset_names:
		bins = np.quantile(dataset_name_to_cis_h2s[eqtl_dataset_name], quantiles)
		bins[-1] = bins[-1] + 1  # Done to make upper bound <= instead of <
		bins[0] = bins[0] - 1  # Done to make lower bound >= instead of >
		dataset_to_bins[eqtl_dataset_name] = bins
	# Get names of bins
	bin_names = []
	for bin_iter in range(n_expr_cis_h2_bins):
		bin_names.append('bin' + str(bin_iter))
	bin_names = np.asarray(bin_names)

	# create list of new eqtl categories
	new_eqtl_categories = []
	for eqtl_dataset_name in eqtl_dataset_names:
		for bin_name in bin_names:
			new_eqtl_categories.append(eqtl_dataset_name + ':' + bin_name)
	new_eqtl_categories = np.asarray(new_eqtl_categories)

	##############################
	# Update gene_info file to now contain bin_info
	##############################
	# loop through genes
	for gene in genes:
		info = gene_info[gene]

		# Get gene cis_h2 vec
		gene_cis_h2_vec = gene_cis_h2[gene]

		# Load in gene eqtl categories
		gene_eqtl_categories = info['eqtl_category_names']

		# quick error check
		if len(gene_cis_h2_vec) != len(gene_eqtl_categories):
			print('assumptionerororo')
			pdb.set_trace()

		# Initialize vector to keep track of update gene_eqtl_categories
		updated_gene_eqtl_categories = []
		# Loop through eQTL categories
		for eqtl_category_iter, eqtl_category_name in enumerate(gene_eqtl_categories):

			# Bin assignment
			bin_name = bin_assignment(gene_cis_h2_vec[eqtl_category_iter], dataset_to_bins[eqtl_category_name], bin_names)
			
			# Update gnee eqtl category
			updated_gene_eqtl_categories.append(eqtl_category_name + ':' + bin_name)
		# Put into nice array
		updated_gene_eqtl_categories = np.asarray(updated_gene_eqtl_categories)

		# Update gene info file
		gene_info[gene]['eqtl_dataset_names'] = np.copy(gene_info[gene]['eqtl_category_names'])
		gene_info[gene]['eqtl_category_names'] = updated_gene_eqtl_categories

	return gene_info, new_eqtl_categories






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
parser.add_argument('--n-expr-cis-h2-bins', default=4, type=int,
                    help='Only used if args.gene_trait_architecture=="stdExpr"')
parser.add_argument('--inference-approach', default='mle', type=str,
                    help='How to perform optimization')
parser.add_argument('--jacknife', default=False, action='store_true',
                    help='Boolean on whether or not to Jacknife the estimtes')
parser.add_argument('--fixed-variance', default=False, action='store_true',
                    help='Boolean on whether the variance parameters should be fixed')
parser.add_argument('--per-dataset-variance', default=False, action='store_true',
                    help='Boolean on whether the variance parameters should be fixed')
parser.add_argument('--output-stem', default=None, type=str,
                    help='Output file stem to save data to')
args = parser.parse_args()

np.random.seed(1)

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
####################

# Create cis h2 bins
if args.gene_trait_architecture == 'linear':
	eqtl_category_names = np.copy(eqtl_dataset_names)
	gene_info = update_gene_info_to_include_eqtl_dataset_names(gene_info)
elif args.gene_trait_architecture == 'stdExpr':
	gene_info, eqtl_category_names = update_gene_info_to_bin_genes_by_est_cis_h2(gene_info, eqtl_dataset_names, args.n_expr_cis_h2_bins)
elif args.gene_trait_architecture == 'random_bins':
	gene_info, eqtl_category_names = update_gene_info_to_bin_genes_by_random(gene_info, eqtl_dataset_names, args.n_expr_cis_h2_bins)
else:
	print('assumption error: ' + str(args.gene_trait_architecture) + ' not an implemented option')
	pdb.set_trace()


#####################
# Open output file handle
#####################
output_file = args.output_stem + '_joint_ldsc_results.txt'
t = open(output_file,'w')
t.write('method\test_med_h2\test_nm_h2\test_mean_eqtl_h2\tper_dataset_med_h2\tper_category_med_h2\n')


#####################
# Run optimization
#####################
# Run two step joint ldsc
obj = joint_ldsc.med_h2(version='two_step')
obj.fit(genes, gene_info, gwas_variant_ld_scores, gwas_E_beta_sq, eqtl_category_names, eqtl_dataset_names, n_reference_snps)
t.write('two_step_joint_ldsc' + '\t' + str(np.sum(obj.category_med_h2)) + '\t' + str(obj.nm_h2) + '\t' + str(obj.avg_eqtl_h2) + '\t' + ','.join(obj.dataset_med_h2.astype(str)) + '\t' + ','.join(obj.category_med_h2.astype(str)) + '\n')
obj = joint_ldsc.med_h2(max_iter=1000,convergence_thresh=1e-10, fixed_variance=args.fixed_variance, per_dataset_variance=args.per_dataset_variance)
if args.inference_approach == 'mle':
	# Run joint ldsc
	obj = joint_ldsc.med_h2(max_iter=1000,convergence_thresh=1e-10, fixed_variance=args.fixed_variance, per_dataset_variance=args.per_dataset_variance)
	obj.fit(genes, gene_info, gwas_variant_ld_scores,gwas_E_beta_sq, eqtl_category_names, eqtl_dataset_names, n_reference_snps)
	t.write('joint_ldsc' + '\t' + str(np.sum(obj.category_med_h2)) + '\t' + str(obj.nm_h2) + '\t' + str(obj.avg_eqtl_h2) + '\t' + ','.join(obj.dataset_med_h2.astype(str)) + '\t' + ','.join(obj.category_med_h2.astype(str)) + '\n')
elif args.inference_approach == 'gibbs':
	# Run joint ldsc
	obj = joint_ldsc_gibbs.med_h2(burn_in_iter=500, max_iter=800,fixed_variance=args.fixed_variance, per_dataset_variance=args.per_dataset_variance)
	obj.fit(genes, gene_info, gwas_variant_ld_scores,gwas_E_beta_sq, eqtl_category_names, eqtl_dataset_names, n_reference_snps)
	t.write('joint_ldsc_gibbs' + '\t' + str(np.sum(obj.category_med_h2)) + '\t' + str(obj.nm_h2) + '\t' + str(obj.avg_eqtl_h2) + '\t' + ','.join(obj.dataset_med_h2.astype(str)) + '\t' + ','.join(obj.category_med_h2.astype(str)) + '\n')



t.close()
print(output_file)









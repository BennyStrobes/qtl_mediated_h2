import numpy as np
import os
import sys
import pdb
import time
import argparse
from bgen import BgenReader
from sklearn.linear_model import Lasso


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


def extract_expression_file_sample_names(expr_file):
	f = open(expr_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		sample_names = np.asarray(data[3:])
		# Only need first line
		break
	f.close()
	return sample_names

def extract_indices_corresponding_to_genotyped_samples(genotype_sample_names, expr_sample_names):
	dicti1 = {}
	for ii,val in enumerate(genotype_sample_names):
		dicti1[val.split('_')[0]] = ii

	mapping = []
	for expr_sample_name in expr_sample_names:
		mapping.append(dicti1[expr_sample_name])

	return mapping


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


def create_mapping_from_rsid_to_snp_index(G_obj_rsids):
	dicti = {}
	for ii, val in enumerate(G_obj_rsids):
		dicti[val] = ii
	return dicti

def get_non_zero_beta(G_cis_std, expr_vec, alpha_init, ens_id):
	factor_grid = [.5, .1, .01, .001, .0001]

	booler = False
	for factor in factor_grid:
		cur_alpha = alpha_init*factor
		lasso_obj = Lasso(alpha=cur_alpha, max_iter=100000).fit(G_cis_std, expr_vec)
		std_beta_hat = lasso_obj.coef_
		non_sparsity = np.sum(std_beta_hat!=0)/len(std_beta_hat)	
		if non_sparsity > 0.0:
			booler = True
			break
	if booler == False:
		print('Zero-gene ' + ens_id)

	return std_beta_hat, cur_alpha


#####################
# Parse command line arguments
#####################
parser = argparse.ArgumentParser()

parser.add_argument('--expr', default='None', type=str,
					help='Expression file')
parser.add_argument('--bgen', default='None', type=str,
					help='Bgen file')
parser.add_argument('--chromosome-file', default='None', type=str,
					help='Filename containing list of chromosomes to run (no header). If None, then defaults to 22 chromosomes')
parser.add_argument('--alpha', default=0.01, type=float,
					help='Lasso regularization parameter')
parser.add_argument('--cis-window', default=500000.0, type=float,
					help='BP region around gene to call cis-eqtls')
parser.add_argument('--output', default=None, type=str,
					help='Output file stem to save data to')
args = parser.parse_args()
np.random.seed(1)


# Open output file handle
t = open(args.output,'w')
# Print header
t.write('GENE\tCHR\tSNP\tSNP_POS\tA1\tA2\tDOSAGE_EFFECT\tSPARSITY_PARAM\tNON_SPARSITY_LEVEL\n')


# Extract array of chromosomes
chromosome_arr = extract_array_of_chromosomes(args.chromosome_file)

# Extract sample names
expr_sample_names = extract_expression_file_sample_names(args.expr)


counter = 0
# Loop through chromosomes
for chrom_num in chromosome_arr:


	# Load in genotype data for this chromosome
	genotype_stem = args.bgen + str(chrom_num) 
	genotype_obj = BgenReader(genotype_stem + '.bgen')
	G_obj_pos = np.asarray(genotype_obj.positions())
	G_obj_rsids = np.asarray(genotype_obj.rsids())
	snp_integers = np.arange(len(G_obj_rsids))
	
	# Create mapping from rsid to snp index
	rsid_to_snp_index = create_mapping_from_rsid_to_snp_index(G_obj_rsids)


	# Load in ref-alt alleles
	ref_alt_alleles = load_in_ref_alt_allele_arr(genotype_stem + '.pvar')
	ref_alt_alleles = np.asarray(ref_alt_alleles)

	# Get indices corresponding to genotyped samples
	genotype_sample_indices = extract_indices_corresponding_to_genotyped_samples(genotype_obj.samples, expr_sample_names)


	# Now loop through expression file (one gene at a time)
	f = open(args.expr)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()

		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Extract relevent fields
		ens_id = data[0]
		line_chrom_num = int(data[1])
		gene_position = int(data[2])

		# Skip genes not on current chromosome
		if line_chrom_num != chrom_num:
			continue
		#print(counter)
		counter = counter + 1

		# Extract vector of gene expression 
		expr_vec = np.asarray(data[3:]).astype(float)

		# Extract cis snp matrix
		cis_start_pos = gene_position - args.cis_window
		cis_end_pos = gene_position + args.cis_window
		cis_snp_indices = (G_obj_pos >= cis_start_pos) & (G_obj_pos <= cis_end_pos)
		cis_rsids = G_obj_rsids[cis_snp_indices]
		cis_pos = G_obj_pos[cis_snp_indices]
		cis_ref_alt_alleles = ref_alt_alleles[cis_snp_indices, :]
		G_cis_dosage = load_in_alt_allele_genotype_dosage_mat(genotype_obj, snp_integers[cis_snp_indices])
		# Filter to relevent individuals
		G_cis_dosage = G_cis_dosage[genotype_sample_indices, :]

		# Center Genotype dosage
		col_means = G_cis_dosage.mean(axis=0)
		G_cis_centered = G_cis_dosage - col_means

		# Standardize genotype
		snp_sdevs = G_cis_centered.std(axis=0)
		# Deal with nans if they exist
		if np.sum(snp_sdevs == 0.0) > 0.0:
			valid_indices = snp_sdevs != 0.0
			G_cis_centered = G_cis_centered[:, valid_indices]
			G_cis_dosage = G_cis_dosage[:, valid_indices]
			cis_rsids = cis_rsids[valid_indices]
			cis_pos = cis_pos[valid_indices]
			cis_ref_alt_alleles = cis_ref_alt_alleles[valid_indices]
			snp_sdevs = snp_sdevs[valid_indices]
		# Do the standardization
		G_cis_std = G_cis_centered/snp_sdevs


		# Fit lasso model
		lasso_obj = Lasso(alpha=args.alpha, max_iter=100000).fit(G_cis_std, expr_vec)
		std_beta_hat = lasso_obj.coef_

		sparsity_param = np.copy(args.alpha)*1.0
		non_sparsity_level = np.sum(std_beta_hat!=0)/len(std_beta_hat)
		if non_sparsity_level == 0.0:
			std_beta_hat, sparsity_param = get_non_zero_beta(G_cis_std, expr_vec, args.alpha, ens_id)
			non_sparsity_level = np.sum(std_beta_hat!=0)/len(std_beta_hat)

		# Convert to dosage based beta-hats
		dosage_beta_hat = np.copy(std_beta_hat)/snp_sdevs

		# Print to output
		for snp_iter, snp_rsid in enumerate(cis_rsids):
			# Skip snps that have no effect
			tmp_beta = dosage_beta_hat[snp_iter]
			if tmp_beta == 0.0:
				continue
			# Extract relevent snp information
			snp_pos = cis_pos[snp_iter]
			snp_a1 = cis_ref_alt_alleles[snp_iter, 1]
			snp_a2 = cis_ref_alt_alleles[snp_iter, 0]
			# print
			t.write(ens_id + '\t' + str(chrom_num) + '\t' + snp_rsid + '\t' + str(snp_pos) + '\t' + snp_a1 + '\t' + snp_a2 + '\t' + str(tmp_beta) + '\t' + str(sparsity_param) + '\t' + str(non_sparsity_level) + '\n')
	f.close()
t.close()






















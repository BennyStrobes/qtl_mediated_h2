import numpy as np
import os
import sys
import pdb
import argparse
import statsmodels.api as sm


def extract_data_from_gecs_input_data_dir(gecs_input_data_dir, regress_in_gene_windows_bool, genotype_intercept_only_bool, bad_chroms={}):
	variant_pair_names = []
	variant_pair_alleles = []
	variant_pair_ld = []
	pairwise_LD_scores = []
	gene_ld_scores = []


	# Loop through autosomal chromosomes
	for chrom_num in range(1,23):
		if chrom_num in bad_chroms:
			continue
		print(chrom_num)

		input_file = gecs_input_data_dir + 'UKBB_pqtl_susie_inf_pmces_' + str(chrom_num) + '_gecs_scores.txt'
		f = open(input_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')

			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue

			# Skip variant pairs not in gene window
			if regress_in_gene_windows_bool and data[7] == 'False':
				continue

			# Update variant pair names
			variant_pair_names.append([data[1], data[2]])
			variant_pair_alleles.append(data[3:7])

			# Add LD
			variant_pair_ld.append(float(data[8]))

			# Add gene Ld scores
			gene_ld_scores.append(float(data[9]))

			# Add pairwise LD scores
			if genotype_intercept_only_bool:
				pairwise_LD_scores.append(float(data[10]))
			else:
				pairwise_LD_scores.append(np.asarray(data[10:]).astype(float))
		f.close()

	# Convert to organized arrays
	variant_pair_names = np.asarray(variant_pair_names)
	variant_pair_alleles = np.asarray(variant_pair_alleles)
	variant_pair_ld = np.asarray(variant_pair_ld)
	pairwise_LD_scores = np.asarray(pairwise_LD_scores)
	gene_ld_scores = np.asarray(gene_ld_scores)

	if genotype_intercept_only_bool:
		pairwise_LD_scores = pairwise_LD_scores.reshape((len(pairwise_LD_scores),1))

	return variant_pair_names,variant_pair_alleles, variant_pair_ld, pairwise_LD_scores, gene_ld_scores

def load_in_gwas_data(gwas_data_file, variant_pair_names, variant_pair_alleles):
	# Create mapping from rsid to gwas z scores
	rsid_to_gwas_z = {}
	f = open(gwas_data_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			if data[0] != 'SNP' or data[1] != 'A1' or data[2] != 'A2' or data[3] != 'N' or data[4] != 'Z':
				print('gwas data assumption error')
				pdb.set_trace()
			continue
		if data[0] + '_' + data[1] + '_' + data[2] in rsid_to_gwas_z or data[0] + '_' + data[2] + '_' + data[1] in rsid_to_gwas_z:
			print('assumption eroorr')
			pdb.set_trace()
		rsid_to_gwas_z[data[0] + '_' + data[1] + '_' + data[2]] = float(data[4])
		rsid_to_gwas_z[data[0] + '_' + data[2] + '_' + data[1]] = -float(data[4])
		gwas_ss = float(data[3])
	f.close()

	# Now create variant pair mapping
	pairwise_zs = []
	for variant_pair_index in range(variant_pair_names.shape[0]):
		variant1 = variant_pair_names[variant_pair_index,0] + '_' + variant_pair_alleles[variant_pair_index,0] + '_' + variant_pair_alleles[variant_pair_index,1]
		variant2 = variant_pair_names[variant_pair_index,1] + '_' + variant_pair_alleles[variant_pair_index,2] + '_' + variant_pair_alleles[variant_pair_index,3]
		if variant1 not in rsid_to_gwas_z or variant2 not in rsid_to_gwas_z:
			pairwise_zs.append(np.nan)
		else:
			z1 = rsid_to_gwas_z[variant1]
			z2 = rsid_to_gwas_z[variant2]
			pairwise_zs.append(z1*z2)

	return np.asarray(pairwise_zs), gwas_ss

def extract_number_of_genes(n_genes_dir):
	n_genes = 0
	for chrom_num in range(1,23):
		filer = n_genes_dir + 'tglr_expression_scores_susie_inf_pmces.' + str(chrom_num) + '.l2.M'
		f = open(filer)
		gene_count = None
		for line in f:
			line = line.rstrip()
			if gene_count is not None:
				print('assumption eroror')
				pdb.set_trace()
			gene_count = float(line)
		f.close()
		if gene_count is None:
			print('assumption eroror')
			pdb.set_trace()
		n_genes = n_genes + gene_count
	return n_genes

def extract_total_h2_est(total_h2_input_stem):
	med_est_file = total_h2_input_stem + '_med_5_50.txt'
	nonmed_est_file = total_h2_input_stem + '_nonmed_5_50.txt'

	med_h2 = float(np.loadtxt(med_est_file, delimiter='\t',dtype=str)[1,0])
	nonmed_h2 = float(np.loadtxt(nonmed_est_file, delimiter='\t',dtype=str)[1,0])
	return med_h2 + nonmed_h2



###########################
# Parse command line args
###########################
parser = argparse.ArgumentParser()
parser.add_argument('--gwas-data-file', default='None', type=str,
                    help='GWAS summary stat file')
parser.add_argument('--gecs-input-data-dir', default='None', type=str,
                    help='Directory containing gecs input data (pairwise ld scores, expression scores, etc)')
parser.add_argument('--n-genes-dir', default='None', type=str,
                    help='Directory containing number of genes used')
parser.add_argument('--total-h2-input-stem', default='None', type=str,
                    help='Directory containing estimate of total h2 for trait')
parser.add_argument('--genotype-intercept-only', default=False, action='store_true',
                    help='Partition non-mediated variant effects by only genotype intercept')
parser.add_argument('--regress-in-gene-windows', default=False, action='store_true',
                    help='Only consider regression snp pairs such that both regression snps lie in cis region of a gene')
parser.add_argument('--out', default='None', type=str,
                    help='Output stem to print eqtl fine-mapping results to')
args = parser.parse_args()


###############################################
# First extract data from gecs input data dir

# Extract:
## 1. names of variant pairs
## 2. ordered alleles of variant pairs
## 3. LD for variants pairs
## 4. Pairwise LD scores (potentially annotations tratified)
## 5. Gene LD scores
'''
variant_pair_names,variant_pair_alleles, variant_pair_LD, pairwise_LD_scores, gene_ld_scores = extract_data_from_gecs_input_data_dir(args.gecs_input_data_dir, args.regress_in_gene_windows, args.genotype_intercept_only)

## Temporary saving (for fast loading)
np.save(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_variant_pair_names.npy', variant_pair_names)
np.save(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_variant_pair_alleles.npy', variant_pair_alleles)
np.save(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_variant_pair_LD.npy', variant_pair_LD)
np.save(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_pairwise_LD_scores.npy', pairwise_LD_scores)
np.save(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_gene_ld_scores.npy', gene_ld_scores)
'''

# Temporary loading
variant_pair_names = np.load(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_variant_pair_names.npy')
variant_pair_alleles = np.load(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_variant_pair_alleles.npy')
variant_pair_LD = np.load(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_variant_pair_LD.npy')
pairwise_LD_scores = np.load(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_pairwise_LD_scores.npy')
gene_ld_scores = np.load(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_gene_ld_scores.npy')



###############################################
# Extract number of genes
n_genes = extract_number_of_genes(args.n_genes_dir)

###############################################
# Extract estimated total h2
total_h2_est = extract_total_h2_est(args.total_h2_input_stem)

###############################################
# Load in GWAS data
pairwise_gwas_z_scores, gwas_ss = load_in_gwas_data(args.gwas_data_file, variant_pair_names, variant_pair_alleles)
# Filter to those snp-pairs we have gwas data for
valid_snp_pairs = np.isnan(pairwise_gwas_z_scores) == False
pairwise_gwas_z_scores = pairwise_gwas_z_scores[valid_snp_pairs]
variant_pair_names = variant_pair_names[valid_snp_pairs]
variant_pair_alleles = variant_pair_alleles[valid_snp_pairs]
variant_pair_LD = variant_pair_LD[valid_snp_pairs]
pairwise_LD_scores = pairwise_LD_scores[valid_snp_pairs]
gene_ld_scores = gene_ld_scores[valid_snp_pairs]


###############################################
# Run GECS regression
Y = pairwise_gwas_z_scores - variant_pair_LD
X = np.hstack((gene_ld_scores.reshape((len(gene_ld_scores),1)), pairwise_LD_scores))

# Tmp saving
xtx_inv = np.linalg.inv(np.dot(np.transpose(X), X))
#np.save(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_xtx_inv.npy', xtx_inv)
# Tmp loading
#xtx_inv = np.load(args.gecs_input_data_dir + 'tmp_save_' + str(args.regress_in_gene_windows) + '_' + str(args.genotype_intercept_only) + '_xtx_inv.npy')
xt_y = np.dot(np.transpose(X),Y)
params = np.dot(xtx_inv, xt_y)
med_h2 = params[0]*n_genes/gwas_ss
frac_mediated = med_h2/total_h2_est

'''
model = sm.OLS(Y,X).fit()
med_h2 = model.params[0]*n_genes/gwas_ss
frac_mediated = med_h2/total_h2_est
'''
###############################################
# Print to output
output_file = args.out + '_mediated_h2_point_estimates.txt'
t = open(output_file,'w')
t.write('med_h2\t' + str(med_h2) + '\n')
t.write('frac_h2_med\t' + str(frac_mediated) + '\n')
t.close()
print(output_file)


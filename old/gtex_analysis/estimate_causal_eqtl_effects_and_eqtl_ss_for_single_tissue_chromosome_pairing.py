import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
from pandas_plink import read_plink1_bin
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')
import argparse
import statsmodels.api as sm



def extract_gwas_variants(gwas_variant_bim_file):
	dicti = {}
	f = open(gwas_variant_bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		rsid = data[1]
		a1 = data[4]
		a2 = data[5]
		if rsid in dicti:
			print('assumption eroror')
			pdb.set_trace()
		dicti[rsid] = (a1, a2)
	f.close()

	return dicti


def load_in_eqtl_genotype_data(genotype_stem, chrom_num, gwas_variants, filter_strand_ambiguous):
	# Load in genotype data across chromosome for eQTL data set
	#genotype_stem = genotype_stem + chrom_num
	G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)
	G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	# a0 is the effect allele
	# For case of plink package, a0 is the first column in the plink bim file
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	# RSids
	G_obj_rsids = np.asarray(G_obj.snp)
	# Sample names
	G_obj_sample_names = np.asarray(G_obj.iid)
	# Snp ids
	G_obj_snp_ids = 'chr' + G_obj_chrom + '_' + (G_obj_pos.astype(str)) + '_' + G_obj_a0 + '_' + G_obj_a1

	# Mean impute and standardize genotype
	G_obj_geno_stand = mean_impute_and_standardize_genotype(G_obj_geno)

	########################
	# Filter variants to those in gwas and (optionally to strand ambigious variants)
	valid_variants = [] # Boolean vector to keep track of which variants pass filters

	# Loop through variants
	for var_iter, var_rsid in enumerate(G_obj_rsids):
		# Initialize boolean for variant to True
		booler = True

		# Identify those variants not in the gwas
		if var_rsid not in gwas_variants:
			booler = False

		# Filter out those variants not on this chromosome
		if G_obj_chrom[var_iter] != chrom_num:
			booler = False

		# Identify strand ambigious variants
		if filter_strand_ambiguous:
			var_a0 = G_obj_a0[var_iter]
			var_a1 = G_obj_a1[var_iter]
			if var_a0 == 'A' and var_a1 == 'T':
				booler = False
			if var_a0 == 'T' and var_a1 == 'A':
				booler = False
			if var_a0 == 'C' and var_a1 == 'G':
				booler = False
			if var_a0 == 'G' and var_a1 == 'C':
				booler = False

		if np.sum(np.isnan(G_obj_geno_stand[:, var_iter])) > 0:
			booler = False

		var_a0 = G_obj_a0[var_iter]
		var_a1 = G_obj_a1[var_iter]

		if gwas_variants[var_rsid][0] == var_a0 and gwas_variants[var_rsid][1] == var_a1:
			tmper = 1
		elif gwas_variants[var_rsid][1] == var_a0 and gwas_variants[var_rsid][0] == var_a1:
			tmper = 1
		else:
			booler = False

		valid_variants.append(booler)
	# Numpy array of valid variants
	valid_variants = np.asarray(valid_variants)

	print('Number of variants before variant filtering: ' + str(len(valid_variants)))
	print('Number of variants after variant filtering: ' + str(int(np.sum(valid_variants))))

	# Now filter output data structures to these variants

	return G_obj_geno_stand[:, valid_variants], G_obj_rsids[valid_variants], G_obj_a0[valid_variants], G_obj_a1[valid_variants], G_obj_pos[valid_variants], G_obj_chrom[valid_variants], G_obj_sample_names


def mean_impute_and_standardize_genotype(G_obj_geno):
	# Fill in missing values
	G_obj_geno_stand = np.copy(G_obj_geno)
	ncol = G_obj_geno_stand.shape[1]
	for col_iter in range(ncol):
		nan_indices = np.isnan(G_obj_geno[:,col_iter])
		non_nan_mean = np.mean(G_obj_geno[nan_indices==False, col_iter])
		G_obj_geno_stand[nan_indices, col_iter] = non_nan_mean


	# Standardize genotype (Mean 0, variance 1)
	G_obj_geno_stand = (G_obj_geno_stand -np.mean(G_obj_geno_stand,axis=0))/np.std(G_obj_geno_stand,axis=0)

	return G_obj_geno_stand

def save_gene_variant_info(gene_variant_info_output_file, rsids, positions, a0s, a1s, chroms):
	# Save in plink format
	t2 = open(gene_variant_info_output_file,'w')

	# Loop through variants
	for var_iter, rsid in enumerate(rsids):
		t2.write(chroms[var_iter] + '\t' + rsid + '\t' + '0' + '\t' + str(positions[var_iter]) + '\t' + a0s[var_iter] + '\t' + a1s[var_iter] + '\n')
	t2.close()
	return

def create_eqtl_sumstat_output_file(eqtl_sumstat_output_file, expr_vec, gene_geno, window_rsids, window_pos, window_a0, window_a1, window_chrom_num):
	# Get number of cis snps in window
	n_cis_snps = gene_geno.shape[1]

	# Quick error check
	if n_cis_snps != len(window_rsids):
		print('assumption eroror')
		pdb.set_trace()

	# Open output eqtl sumstat file
	t = open(eqtl_sumstat_output_file,'w')
	t.write('rsid\tchrom_num\tposition\ta1\ta2\tbeta\tbeta_se\n')

	# Iterate through snps
	for snp_iter in range(n_cis_snps):
		model = sm.OLS(expr_vec, gene_geno[:, snp_iter])
		results = model.fit()
		eqtl_effect_size = results.params[0]
		eqtl_effect_size_se = results.bse[0]

		t.write(window_rsids[snp_iter] + '\t' + window_chrom_num[snp_iter] + '\t' + str(window_pos[snp_iter]) + '\t' + window_a0[snp_iter] + '\t' + window_a1[snp_iter] + '\t' + str(eqtl_effect_size) + '\t' + str(eqtl_effect_size_se) + '\n')

	t.close()
	return

def load_in_covariate_matrix(cov_file):
	raw = np.loadtxt(cov_file, dtype=str, delimiter='\t')
	cov_sample_names = raw[1:,1]
	cov_mat = raw[1:,2:].astype(float)
	return cov_mat, cov_sample_names

def regress_out_covariates_on_gene_expression(expr_vec, cov_mat):
	model = sm.OLS(expr_vec, sm.add_constant(cov_mat))
	results = model.fit()
	pred_expr = np.dot(sm.add_constant(cov_mat), results.params)
	resid_expr_vec = expr_vec - pred_expr
	resid_expr_vec = (resid_expr_vec - np.mean(resid_expr_vec))/np.std(resid_expr_vec)
	return resid_expr_vec


def run_susie_eqtl_fine_mapping_with_individual_data(gwas_variants, expression_file, chrom_num, genotype_stem, cis_window_size, filter_strand_ambiguous, min_cis_snps_per_gene,cov_file, output_stem):
	#############################
	# Load in covariates
	#############################	
	cov_mat, cov_sample_names = load_in_covariate_matrix(cov_file)


	#############################
	# Output summary file (to keep track of all genes)
	#############################
	output_summary_file = output_stem + '_chr' + str(chrom_num) + '_gene_summary.txt'
	t = open(output_summary_file,'w')
	t.write('Gene\tCHR\tGENE_COORD\tINFO\tvarint_info_file\tsusie_alpha_file\tsusie_mu_file\tsusie_mu_var_file\tsusie_pmces_file\tsumstat_file\n')
	

	#############################
	# Load in eQTL Genotype data
	#############################
	# Outputs: 
	# G_mat: matrix of standarized genotype of dimension number of samples by number of snps
	# G_rsids: Vector of rsids corresponding to G_mat
	# G_a0: Vector of allele 0 values corresponding ot G_mat
	# G_a1: Vector of allele 1 values corresponding to G_mat
	# G_pos: Vector variant positions corresponding to G_mat
	# G_chrom: Vector of variant chromosomes corresponding to G_mat
	# G_sample_names: Vector of sample names corresponding to G_mat
	print('Load in Genotype data')
	G_mat, G_rsids, G_a0, G_a1, G_pos, G_chrom, G_sample_names = load_in_eqtl_genotype_data(genotype_stem, chrom_num, gwas_variants, filter_strand_ambiguous)


	#############################
	# Fit gene model in each gene, independently
	#############################
	print('Fit SuSiE eQTL gene models')
	# Loop through genes
	f = open(expression_file)
	head_count = 0  # To identify header
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Header
		if head_count == 0:
			head_count = head_count + 1
			# Make sure expression sample names and genotype sample names match
			expr_sample_names = np.asarray(data[3:])
			if np.array_equal(G_sample_names, expr_sample_names) == False:
				print('ASSUMPTION ERROR: Sample names in genotype file do not match sample names in expression matrix')
				pdb.set_trace()
			if np.array_equal(G_sample_names, cov_sample_names) == False:
				print('ASSUMPTION ERROR: Sample names in genotype file do not match sample names in expression matrix')
				pdb.set_trace()
			continue

		# Standard line (corresponding to a gene)
		ensamble_id = data[0]
		gene_chrom_num = data[1]
		# Skip genes not on this chromosome
		if gene_chrom_num != chrom_num:
			continue
		gene_position = int(data[2])
		expr_vec = np.asarray(data[3:]).astype(float)

		# Regress out covariates on gene expression
		resid_expr_vec = regress_out_covariates_on_gene_expression(expr_vec, cov_mat)
		expr_vec = np.copy(resid_expr_vec)

		# Get indices of variants corresponding to cis window fo this gene
		cis_window_start = gene_position - cis_window_size
		cis_window_end = gene_position + cis_window_size
		cis_snp_indices = (G_pos >= cis_window_start) & (G_pos < cis_window_end)

		#  Ignore genes with no or very few cis snps
		if np.sum(cis_snp_indices) < min_cis_snps_per_gene:
			print('gene skipped because it contained 0 or very small number of cis snps')
			t.write(ensamble_id + '\t' + gene_chrom_num + '\t' + str(gene_position) + '\t' + 'Fail_too_few_snps\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')
			continue

		# Extract standardized matrix of cis snps around the gene
		gene_geno = G_mat[:, cis_snp_indices]

		# Save variant info to output
		gene_variant_info_output_file = output_stem + '_' + ensamble_id + '_gene_variant_info.txt'
		save_gene_variant_info(gene_variant_info_output_file, G_rsids[cis_snp_indices], G_pos[cis_snp_indices], G_a0[cis_snp_indices], G_a1[cis_snp_indices], G_chrom[cis_snp_indices])

		# Generate eQTL summary stats
		eqtl_sumstat_output_file = output_stem + '_' + ensamble_id + '_eqtl_sumstats.txt'
		create_eqtl_sumstat_output_file(eqtl_sumstat_output_file, expr_vec, gene_geno, G_rsids[cis_snp_indices], G_pos[cis_snp_indices], G_a0[cis_snp_indices], G_a1[cis_snp_indices], G_chrom[cis_snp_indices])

		# Run eQTL variant fine-mapping with SuSiE
		susie_fitted = susieR_pkg.susie(gene_geno, expr_vec, L=10)

		# Test whether are 0 identified susie components for this gene
		if type(susie_fitted.rx2('sets').rx2('cs_index')) == rpy2.rinterface_lib.sexp.NULLType:
			t.write(ensamble_id+ '\t' + gene_chrom_num + '\t' + str(gene_position) + '\t' + 'Fail_purity_filter\t' + gene_variant_info_output_file + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + eqtl_sumstat_output_file + '\n')
			continue

		# This gene has passed purity filter
		susie_components = np.asarray(susie_fitted.rx2('sets').rx2('cs_index')) - 1
		pmces = np.sum(susie_fitted.rx2('alpha')*susie_fitted.rx2('mu'),axis=0)
		
		# Need to save individual SuSiE posterior objects
		# alpha
		alpha_model_output_file = output_stem + '_' + ensamble_id + '_gene_model_susie_alpha.npy'
		np.save(alpha_model_output_file, susie_fitted.rx2('alpha'))
		# mu
		mu_model_output_file = output_stem + '_' + ensamble_id + '_gene_model_susie_mu.npy'
		np.save(mu_model_output_file, susie_fitted.rx2('mu'))
		# mu_var
		mu_var = susie_fitted.rx2('mu2') - np.square(susie_fitted.rx2('mu'))
		mu_var_model_output_file = output_stem + '_' + ensamble_id + '_gene_model_susie_mu_var.npy'
		np.save(mu_var_model_output_file, mu_var)
		# PMCES
		pmces_model_output_file = output_stem + '_' + ensamble_id + '_gene_model_susie_pmces.npy'
		np.save(pmces_model_output_file, pmces)

		# Print filenames to summary file
		t.write(ensamble_id + '\t' + gene_chrom_num + '\t' + str(gene_position) + '\t' + 'Pass' + '\t' + gene_variant_info_output_file + '\t' + alpha_model_output_file + '\t' + mu_model_output_file + '\t' + mu_var_model_output_file + '\t' + pmces_model_output_file + '\t' + eqtl_sumstat_output_file + '\n')
		t.flush()
	f.close()
	t.close()
	return


tmp_plink_stem = sys.argv[1]
expr_file = sys.argv[2]
cov_file = sys.argv[3]
gwas_genotype_stem = sys.argv[4]
chrom_num = sys.argv[5]
eqtl_output_root = sys.argv[6]

cis_window_size=500000
filter_strand_ambiguous=True
min_cis_snps_per_gene= 10

# Extract gwas variants
gwas_variant_bim_file = gwas_genotype_stem + chrom_num + '.bim'
gwas_variants = extract_gwas_variants(gwas_variant_bim_file)


run_susie_eqtl_fine_mapping_with_individual_data(gwas_variants, expr_file, chrom_num, tmp_plink_stem, cis_window_size, filter_strand_ambiguous, min_cis_snps_per_gene, cov_file, eqtl_output_root)


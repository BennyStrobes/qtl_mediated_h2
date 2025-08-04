import numpy as np
import os
import sys
import pdb
from bgen import BgenReader
import gzip
import time

def update_rsid_info_with_pvar_file(pvar_file, rsid_to_alleles, rsid_to_position, rsid_to_chrom):
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
		rsid = data[2]
		snp_pos = data[1]
		snp_chrom = data[0]

		if rsid in rsid_to_alleles or rsid in rsid_to_position or rsid in rsid_to_chrom:
			print('assumptino eororor')
			pdb.set_trace()

		rsid_to_alleles[rsid] = (alt_allele, ref_allele)
		rsid_to_position[rsid] = snp_pos
		rsid_to_chrom[rsid] = snp_chrom

	f.close()
	return rsid_to_alleles, rsid_to_position, rsid_to_chrom

def extract_rsid_info(bgen_file_stem):
	rsid_to_chrom = {}
	rsid_to_position = {}
	rsid_to_alleles = {}

	for chrom_num in ['1', '2']:
		pvar_file = bgen_file_stem + chrom_num + '.pvar'
		rsid_to_alleles, rsid_to_position, rsid_to_chrom = update_rsid_info_with_pvar_file(pvar_file, rsid_to_alleles, rsid_to_position, rsid_to_chrom)


	return rsid_to_chrom, rsid_to_position, rsid_to_alleles


causal_eqtl_summary_file = sys.argv[1]
bgen_file_stem = sys.argv[2]
simulated_gene_expression_dir = sys.argv[3]
simulation_name_string = sys.argv[4]
study_names_file = sys.argv[5]

n_tissues=5

rsid_to_chrom, rsid_to_position, rsid_to_alleles = extract_rsid_info(bgen_file_stem)


# Open output files
tt = {}
for tissue_number in range(n_tissues):
	output_file = simulated_gene_expression_dir + simulation_name_string + '_tissue' + str(tissue_number) + '_true_causal_standardized_effects.txt'
	tt[tissue_number] = open(output_file,'w')
	tt[tissue_number].write('GENE\tCHR\tSNP\tSNP_POS\tA1\tA2\tSTANDARDIZED_EFFECT\tSPARSITY_PARAM\tNON_SPARSITY_LEVEL\n')



f = open(causal_eqtl_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract relevent fields
	ens_id = data[0]
	causal_eqtl_effect_file = data[3]
	rsid_file = data[4]

	causal_eqtl_effects = np.load(causal_eqtl_effect_file)
	rsids = np.load(rsid_file)

	if causal_eqtl_effects.shape[1] != n_tissues:
		print('assumption eroror')
		pdb.set_trace()

	n_snps = len(rsids)
	for tissue_iter in range(n_tissues):
		for snp_iter in range(n_snps):
			if causal_eqtl_effects[snp_iter, tissue_iter] == 0.0:
				continue

			rsid = rsids[snp_iter]
			snp_chrom = rsid_to_chrom[rsid]
			snp_position = rsid_to_position[rsid]
			snp_allele = rsid_to_alleles[rsid]

			tt[tissue_iter].write(ens_id + '\t' + snp_chrom + '\t' + rsid + '\t' + snp_position + '\t' + snp_allele[0] + '\t' + snp_allele[1] + '\t' + str(causal_eqtl_effects[snp_iter,tissue_iter]) + '\t' + 'NA\tNA\n')

f.close()


t_global = open(study_names_file,'w')
t_global.write('study_name\tstudy_causal_eqtl_effect_file\n')




# Close file handles
for tissue_number in range(n_tissues):
	tt[tissue_number].close()

	output_file = simulated_gene_expression_dir + simulation_name_string + '_tissue' + str(tissue_number) + '_true_causal_standardized_effects.txt'
	study_name = simulation_name_string + '_tissue' + str(tissue_number) + '_true'
	t_global.write(study_name + '\t' + output_file + '\n')
t_global.close()




import sys
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from bgen import BgenReader


def load_in_alt_allele_genotype_dosage_mat(bfile, window_indices, ref_alt_alleles):
	dosages = []

	for window_index in window_indices:
		var = bfile[window_index]
		dosage = var.minor_allele_dosage
		ma = var.minor_allele

		index_ref_alt_allele = ref_alt_alleles[window_index]

		# Flip dosage if alt-allele is not equal to minor allele
		if index_ref_alt_allele[1] != ma:
			# Quick error check
			if ma != index_ref_alt_allele[0]:
				print('assumptino eroror')
				pdb.set_trace()
			# Flip dosage
			dosage = 2.0 - dosage

		# Append snp dosage to global array
		dosages.append(dosage)

	# Convert to 2d matrix
	dosages = np.asarray(dosages)

	return np.transpose(dosages)


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

def standardize_genotype_dosage_matrix(genotype_dosage):
	# Quick error checking to make sure there do not exist missing entries
	n_missing = np.sum(np.isnan(genotype_dosage))
	if n_missing != 0:
		print('assumption eroror')
		pdb.set_trace()

	# Now standardize genotype of each snp
	n_snps = genotype_dosage.shape[1]
	# Initialize standardize genotype dosage matrix
	std_genotype_dosage = np.copy(genotype_dosage)
	for snp_iter in range(n_snps):
		# Standardize
		std_genotype_dosage[:, snp_iter] = (genotype_dosage[:,snp_iter] - np.mean(genotype_dosage[:,snp_iter]))/np.std(genotype_dosage[:,snp_iter])

	return std_genotype_dosage



def construct_variant_ld_scores_based_on_distance(genotype_stem, output_stem, window_size=1):
	# Load in ref-alt alleles
	ref_alt_alleles = load_in_ref_alt_allele_arr(genotype_stem + '.pvar')

	# Load in genotype object
	bfile = BgenReader(genotype_stem + '.bgen')

	snp_pos = np.asarray(genotype_obj['pos'])
	snp_rsids = np.asarray(genotype_obj['snp'])

	inner_window_start_pos = np.min(snp_pos) - 1

	mb = 1000000.0

	used_snps = {}

	output_file = output_stem + '.txt'
	t = open(output_file,'w')
	output_file2 = output_stem + '_bias_corrected.txt'
	t2 = open(output_file2,'w')	
	t.write('variant_id\tld_score\n')
	t2.write('variant_id\tld_score\n')

	window_counter = 0

	while inner_window_start_pos < np.max(snp_pos):

		print(window_counter)
		window_counter = window_counter + 1

		# Define windows
		inner_window_end_pos = inner_window_start_pos + window_size*mb
		outer_window_start_pos = inner_window_start_pos - (window_size*mb)
		outer_window_end_pos = inner_window_end_pos + (window_size*mb)

		window_position_indices = (snp_pos >= outer_window_start_pos) & (snp_pos < outer_window_end_pos)

		if np.sum(window_position_indices) == 0:
			inner_window_start_pos = inner_window_start_pos + (window_size*mb)
			continue

		# Extract genotype data for this window
		window_indices = np.arange(len(window_position_indices))[window_position_indices].astype(int)
		genotype_dosage = load_in_alt_allele_genotype_dosage_mat(bfile, window_indices, ref_alt_alleles)

		# Standardize genotype dosage matrix
		window_std_variant_genotype = standardize_genotype_dosage_matrix(genotype_dosage)

		pdb.set_trace()


		# Compute LD matrix for window
		ref_ss =window_std_variant_genotype.shape[0]
		squared_ld_mat = np.square(np.corrcoef(np.transpose(window_std_variant_genotype)))
		bias_corrected_squared_ld_mat = squared_ld_mat - ((1.0-squared_ld_mat)/(ref_ss-2.0))

		# Get into window orientation
		window_snp_pos = snp_pos[window_position_indices]
		window_rsids = snp_rsids[window_position_indices]

		# Loop through snps until we find a middle snp
		for tmp_snp_iter, tmp_snp_pos in enumerate(window_snp_pos):
			if tmp_snp_pos >= inner_window_start_pos and tmp_snp_pos < inner_window_end_pos:
				# This is a middle snp
				tmp_snp_rsid = window_rsids[tmp_snp_iter]

				if tmp_snp_rsid in used_snps:
					print('assumption eroror')
					pdb.set_trace()
				used_snps[tmp_snp_rsid] = 1

				tmp_snp_indices = (window_snp_pos >= tmp_snp_pos - (window_size*mb)) & (window_snp_pos < tmp_snp_pos + (window_size*mb))
				
				snp_ld_score = np.sum(squared_ld_mat[tmp_snp_iter, tmp_snp_indices])
				bias_corrected_snp_ld_score = np.sum(bias_corrected_squared_ld_mat[tmp_snp_iter, tmp_snp_indices])

				t.write(tmp_snp_rsid + '\t' + str(snp_ld_score) + '\n')
				t2.write(tmp_snp_rsid + '\t' + str(bias_corrected_snp_ld_score) + '\n')

		inner_window_start_pos = inner_window_start_pos + (window_size*mb)

	t.close()
	t2.close()


	return















##############################
# Command line argumemnts
##############################
simulation_number = int(sys.argv[1])
chrom_num = sys.argv[2]
simulation_name_string = sys.argv[3]
processed_genotype_data_dir = sys.argv[4]
simulated_gwas_dir = sys.argv[5]





ref_genotype_stem = processed_genotype_data_dir +'simulated_eqtl_1000_data_1'
output_stem = simulated_gwas_dir + simulation_name_string + '_variant_ref_geno_window_3_mb_ld_scores'
construct_variant_ld_scores_based_on_distance(ref_genotype_stem, output_stem, window_size=3)





in_sample_genotype_stem = processed_genotype_data_dir + 'simulated_gwas_data_1'


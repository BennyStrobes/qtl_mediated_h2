import sys
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from bgen import BgenReader

def load_in_alt_allele_genotype_dosage_mat(bfile, window_indices, ref_alt_alleles):
	dosages = []

	for window_index in window_indices:
		print(window_index)
		var = bfile[window_index]
		dosage = var.alt_dosage
		#ma = var.minor_allele

		index_ref_alt_allele = ref_alt_alleles[window_index]

		if var.alleles[0] != index_ref_alt_allele[1]:
			print('assumption errror')
			pdb.set_trace()
		if var.alleles[1] != index_ref_alt_allele[0]:
			print('assumption eroror')
			pdb.set_trace()
		'''
		# Flip dosage if alt-allele is not equal to minor allele
		if index_ref_alt_allele[1] != ma:
			# Quick error check
			if ma != index_ref_alt_allele[0]:
				print('assumptino eroror')
				pdb.set_trace()
			# Flip dosage
			dosage = 2.0 - dosage
		'''

		# Append snp dosage to global array
		dosages.append(dosage)

	# Convert to 2d matrix
	dosages = np.asarray(dosages)

	return np.transpose(dosages)


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



simulation_genotype_dir = sys.argv[1]
ld_dir = sys.argv[2]

chrom_num=1
gwas_plink_stem = simulation_genotype_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype directory

# Load in genotype object
genotype_obj = BgenReader(gwas_plink_stem + '.bgen')

print('loaded')
# Load in ref-alt alleles
ref_alt_alleles = load_in_ref_alt_allele_arr(gwas_plink_stem + '.pvar')
chrom_rsids = np.asarray(genotype_obj.rsids())
gene_indices = np.arange(200000, 215000)


genotype_dosage = load_in_alt_allele_genotype_dosage_mat(genotype_obj, gene_indices, ref_alt_alleles)
stand_eqtl_genotype = standardize_genotype_dosage_matrix(genotype_dosage)

LD = np.corrcoef(np.transpose(stand_eqtl_genotype))

np.save(ld_dir + 'ld.npy', LD)
np.save(ld_dir + 'std_geno.npy', stand_eqtl_genotype)



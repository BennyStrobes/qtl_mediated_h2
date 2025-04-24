import sys
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from bgen import BgenReader



# Extract dictionary list of hapmap3 rsids (ie regression snps in sldsc)
def extract_dictionary_list_of_hapmap3_rsids(hm3_rsid_file):
	# Initialize dictionary
	dicti = {}

	# Stream rsid file
	f = open(hm3_rsid_file)
	for line in f:
		line = line.rstrip()
		if line in dicti:
			print('assumption eroror')
			pdb.set_trace()
		dicti[line] = 1
	f.close()
	return dicti

def load_in_ordered_array_of_rsids_from_bim_file(genotype_bim):
	f = open(genotype_bim)
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		arr.append(data[1])
	f.close()
	arr = np.asarray(arr)

	# Quick error check 
	if len(np.unique(arr)) != len(arr):
		print('assumption eroror')
		pdb.set_trace()

	return arr


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


# Run GWAS on all snps
def run_gwas_on_all_snps(trait_values_file, gwas_plink_stem, gwas_output_file, batch_size=1000):
	# Load in trait vector
	trait_vector = np.loadtxt(trait_values_file)

	# Load in genotype object
	genotype_obj = BgenReader(gwas_plink_stem + '.bgen')

	# Load in ref-alt alleles
	ref_alt_alleles = load_in_ref_alt_allele_arr(gwas_plink_stem + '.pvar')

	snp_pos = np.asarray(genotype_obj.positions())
	ordered_rsids = np.asarray(genotype_obj.rsids())


	# Open output file handle
	t = open(gwas_output_file,'w')
	t.write('rsid\tbeta\tbeta_se\tz\n')


	snp_counter = 0
	# Loop through snps and run gwas for all snps
	for snp_iter, snp_rsid in enumerate(ordered_rsids):
		snp_counter = snp_counter + 1

		var = genotype_obj[snp_iter]
		dosage = var.minor_allele_dosage
		ma = var.minor_allele
		index_ref_alt_allele = ref_alt_alleles[snp_iter]
		if index_ref_alt_allele[1] != ma:
			# Quick error check
			if ma != index_ref_alt_allele[0]:
				print('assumptino errror')
			# flip dosage
			dosage = 2.0 - dosage
		std_variant_genotype = (dosage - np.mean(dosage))/np.std(dosage)

		# Now get effect size, standard erorr and z-score for association between standardized genotype and trait
		# Fit model using statsmodels
		mod = sm.OLS(trait_vector, sm.add_constant(std_variant_genotype))
		res = mod.fit()
		# Extract results
		effect_size = res.params[1]
		effect_size_se = res.bse[1]
		effect_size_z = effect_size/effect_size_se

		# Print to output file
		t.write(snp_rsid + '\t' + str(effect_size) + '\t' + str(effect_size_se) + '\t' + str(effect_size_z) + '\n')

	t.close()

	return




##############################
# Command line argumemnts
##############################
simulation_number = int(sys.argv[1])
chrom_num = sys.argv[2]
simulation_name_string = sys.argv[3]
processed_genotype_data_dir = sys.argv[4]
simulated_trait_dir = sys.argv[5]
simulated_gwas_dir = sys.argv[6]




####################################################
# Run GWAS on hapmap3 rsids (ie regression snps in sldsc)
####################################################
trait_values_file = simulated_trait_dir + simulation_name_string + '_trait_values.txt'  # Trait vector
gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype files
gwas_output_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results.txt'
run_gwas_on_all_snps(trait_values_file, gwas_plink_stem, gwas_output_file)










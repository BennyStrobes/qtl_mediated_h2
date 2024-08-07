import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os
import pdb
from pandas_plink import read_plink1_bin
import statsmodels.api as sm



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

# Run GWAS on all snps
def run_gwas_on_all_snps(trait_values_file, gwas_plink_stem, gwas_output_file, batch_size=1000):
	# Load in trait vector
	trait_vector = np.loadtxt(trait_values_file)

	# Load in ordered array of rsids
	genotype_bim = gwas_plink_stem + '.bim'
	ordered_rsids = load_in_ordered_array_of_rsids_from_bim_file(genotype_bim)

	# Load in genotype object
	genotype_obj = read_plink1_bin(gwas_plink_stem + '.bed', gwas_plink_stem + '.bim', gwas_plink_stem + '.fam', verbose=False)

	# Open output file handle
	t = open(gwas_output_file,'w')
	t.write('rsid\tbeta\tbeta_se\tz\n')

	# Initialize arrays to perform batch based analysis
	batch_rsids = []
	batch_variant_names = []

	snp_counter = 0
	# Loop through snps and run gwas for all snps
	for snp_iter, snp_rsid in enumerate(ordered_rsids):
		snp_counter = snp_counter + 1
		batch_rsids.append(snp_rsid)
		batch_variant_names.append('variant' + str(snp_iter))

		# Only extract genotypes for variants in bach
		if np.mod(len(batch_rsids), batch_size) == 0:
			print(snp_counter)
			# Extract genotype data for this batch of snps
			batch_variant_genotype = np.asarray(genotype_obj.sel(variant=batch_variant_names))

			# Loop through variants in batch
			for batch_iter, batch_rsid in enumerate(np.asarray(batch_rsids)):
				# Extract genotype data for this snp
				variant_genotype = batch_variant_genotype[:,batch_iter]
				std_variant_genotype = np.copy(variant_genotype)
				# Mean impute genotype
				nan_indices = np.isnan(variant_genotype)
				non_nan_mean = np.mean(variant_genotype[nan_indices==False])
				std_variant_genotype[nan_indices] = non_nan_mean
				# Standardize genotype
				std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)


				# Now get effect size, standard erorr and z-score for association between standardized genotype and trait
				# Fit model using statsmodels
				mod = sm.OLS(trait_vector, sm.add_constant(std_variant_genotype))
				res = mod.fit()
				# Extract results
				effect_size = res.params[1]
				effect_size_se = res.bse[1]
				effect_size_z = effect_size/effect_size_se

				# Print to output file
				t.write(batch_rsid + '\t' + str(effect_size) + '\t' + str(effect_size_se) + '\t' + str(effect_size_z) + '\n')

			# Reset batch_rsids and batch_variant_names
			batch_rsids = []
			batch_variant_names = []

	# Print stragglers
	if len(batch_rsids) > 0:
		# Extract genotype data for this batch of snps
		batch_variant_genotype = np.asarray(genotype_obj.sel(variant=batch_variant_names))

		# Loop through variants in batch
		for batch_iter, batch_rsid in enumerate(np.asarray(batch_rsids)):
			# Extract genotype data for this snp
			variant_genotype = batch_variant_genotype[:,batch_iter]
			std_variant_genotype = np.copy(variant_genotype)
			# Mean impute genotype
			nan_indices = np.isnan(variant_genotype)
			non_nan_mean = np.mean(variant_genotype[nan_indices==False])
			std_variant_genotype[nan_indices] = non_nan_mean
			# Standardize genotype
			std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)


			# Now get effect size, standard erorr and z-score for association between standardized genotype and trait
			# Fit model using statsmodels
			mod = sm.OLS(trait_vector, sm.add_constant(std_variant_genotype))
			res = mod.fit()
			# Extract results
			effect_size = res.params[1]
			effect_size_se = res.bse[1]
			effect_size_z = effect_size/effect_size_se

			# Print to output file
			t.write(batch_rsid + '\t' + str(effect_size) + '\t' + str(effect_size_se) + '\t' + str(effect_size_z) + '\n')


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










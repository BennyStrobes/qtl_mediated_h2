import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np
import os
import pdb
import statsmodels.api as sm



def load_in_trait_vec(trait_values_file):
	f = open(trait_values_file)
	arr = []
	for line in f:
		line = line.rstrip()
		arr.append(float(line))
	return np.asarray(arr)






#####################
# Command line args
#####################
simulation_name_string = sys.argv[1]
processed_genotype_data_dir = sys.argv[2]
simulated_gwas_data_dir = sys.argv[3]
selection_type = sys.argv[4]


# Load in gwas genotype data
#gwas_genotype = np.load(processed_genotype_data_dir + 'gwas_genotype_1.npy')
#n_snps = gwas_genotype.shape[1]


# Load in trait data
trait_values_file = simulated_gwas_data_dir + simulation_name_string + selection_type + '_selection_' + 'simulated_trait.txt'
trait_vector = load_in_trait_vec(trait_values_file)


# Compute GWAS summary statistics
gwas_output_file = simulated_gwas_data_dir + simulation_name_string + selection_type + '_selection_' + 'simulated_gwas_summary_stats.txt'
t = open(gwas_output_file,'w')
t.write('snp\tbeta\tbeta_se\tz\n')

n_bins=100


total_snp_iter = 0
bin_start = 0
for bin_iter in range(n_bins):
	bin_genotype = np.load(processed_genotype_data_dir + 'gwas_genotype_1_' + str(bin_iter) + '.npy')
	n_bin_snps = bin_genotype.shape[1]

	for snp_iter in range(n_bin_snps):
		std_variant_genotype = bin_genotype[:, snp_iter]

		mod = sm.OLS(trait_vector, sm.add_constant(std_variant_genotype))
		res = mod.fit()
		# Extract results
		effect_size = res.params[1]
		effect_size_se = res.bse[1]
		effect_size_z = effect_size/effect_size_se

		t.write('variant_' + str(total_snp_iter) + '\t' + str(effect_size) + '\t' + str(effect_size_se) + '\t' + str(effect_size_z) + '\n')
		total_snp_iter = total_snp_iter + 1



'''
for snp_iter in range(n_snps):
	std_variant_genotype = gwas_genotype[:, snp_iter]

	mod = sm.OLS(trait_vector, std_variant_genotype)
	res = mod.fit()
	# Extract results
	effect_size = res.params[0]
	effect_size_se = res.bse[0]
	effect_size_z = effect_size/effect_size_se

	t.write('variant_' + str(snp_iter) + '\t' + str(effect_size) + '\t' + str(effect_size_se) + '\t' + str(effect_size_z) + '\n')
'''
t.close()

print(gwas_output_file)







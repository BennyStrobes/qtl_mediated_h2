import numpy as np 
import os
import sys
import pdb
import rss_vi_variant_only


def load_in_gwas_z_scores(gwas_sum_stats_file):
	f = open(gwas_sum_stats_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(float(data[3]))
	f.close()
	arr = np.asarray(arr)
	return arr







#########################
# Command line args
#########################
simulation_name_string = sys.argv[1]
simulated_gwas_data_dir = sys.argv[2]
simulated_gene_models_dir = sys.argv[3]
eqtl_ss = sys.argv[4]
mediated_h2_results_dir = sys.argv[5]
processed_genotype_data_dir = sys.argv[6]


# File summarizing inferred gene models
gene_summary_file = simulated_gene_models_dir + simulation_name_string + 'model_summaries_' + eqtl_ss + '.txt' 

# File containing gwas summary statistics
gwas_sum_stats_file = simulated_gwas_data_dir + simulation_name_string + 'simulated_gwas_summary_stats.txt'
z_scores = load_in_gwas_z_scores(gwas_sum_stats_file)

# LD file
ld_file = processed_genotype_data_dir + 'gwas_genotype_LD_1.npy'
LD = np.load(ld_file)


# RSS VI variantly only
rss_vi_var_only = rss_vi_variant_only.RSS_VI_VARIANT_ONLY(LD=LD, z=z_scores)
rss_vi_var_only.fit()




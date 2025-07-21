import numpy as np
import os
import sys
import pdb



def generate_xt_eqtl_summary_file(mesc_expression_score_dir, simulation_name_string, method_name, output_file):
	t = open(output_file,'w')
	t.write('tissue_name\tN_rep1\tN_rep2\tN_full\teqtl_sumstat_rep1_file\teqtl_sumstat_rep2_file\teqtl_sumstat_full_file\n')
	tissue_names = ['tissue0', 'tissue1', 'tissue2', 'tissue3', 'tissue4']

	for tissue_name in tissue_names:

		rep1_file = mesc_expression_score_dir + simulation_name_string + '_' + tissue_name + '_' +  method_name + '_pred_eqtl_sumstats.txt'
		rep2_file = mesc_expression_score_dir + simulation_name_string + '_' + tissue_name + '_' +  method_name + '_pred_eqtl_sumstats.txt'
		t.write(tissue_name + '\t' + 'INF' + '\t' + 'INF' + '\t' + 'INF' + '\t' + rep1_file + '\t' + rep2_file + '\t' + 'NA' + '\n')

	t.close()

	return




#####################
# Command line args
#####################
simulation_number = sys.argv[1]
simulation_name_string = sys.argv[2]
mesc_expression_score_dir = sys.argv[3]





output_file = mesc_expression_score_dir + simulation_name_string + '_' + 'true_eqtls_corr_replicate_eqtl_sumstats_xt_summary.txt'
generate_xt_eqtl_summary_file(mesc_expression_score_dir,simulation_name_string, 'true', output_file)
output_file = mesc_expression_score_dir + simulation_name_string + '_' + 'true_eqtls_corr_standardized_replicate_eqtl_sumstats_xt_summary.txt'
generate_xt_eqtl_summary_file(mesc_expression_score_dir,simulation_name_string, 'true_standardized', output_file)

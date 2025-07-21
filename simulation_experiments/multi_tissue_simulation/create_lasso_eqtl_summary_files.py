import numpy as np
import os
import sys
import pdb



def generate_xt_eqtl_summary_file(mesc_expression_score_dir, simulation_name_string, eqtl_ss, method_name, output_file):
	t = open(output_file,'w')
	t.write('tissue_name\tN_rep1\tN_rep2\tN_full\teqtl_sumstat_rep1_file\teqtl_sumstat_rep2_file\teqtl_sumstat_full_file\n')
	tissue_names = ['tissue0', 'tissue1', 'tissue2', 'tissue3', 'tissue4']

	for tissue_name in tissue_names:

		if eqtl_ss == '100-1000':
			if tissue_name == 'tissue0':
				tissue_qtl_ss = '100'
			else:
				tissue_qtl_ss = '1000'
		else:
			tissue_qtl_ss = eqtl_ss
		sub_eqtl_ss = str(int(float(tissue_qtl_ss)/2.0))

		rep1_file = mesc_expression_score_dir + simulation_name_string + '_' + tissue_name + '_' + str(tissue_qtl_ss) + '_replicate1_' + method_name + '_pred_eqtl_sumstats.txt'
		rep2_file = mesc_expression_score_dir + simulation_name_string + '_' + tissue_name + '_' + str(tissue_qtl_ss) + '_replicate2_' + method_name + '_pred_eqtl_sumstats.txt'

		t.write(tissue_name + '\t' + sub_eqtl_ss + '\t' + sub_eqtl_ss + '\t' + tissue_qtl_ss + '\t' + rep1_file + '\t' + rep2_file + '\t' + 'NA' + '\n')



	t.close()

	return




#####################
# Command line args
#####################
simulation_number = sys.argv[1]
simulation_name_string = sys.argv[2]
mesc_expression_score_dir = sys.argv[3]
alpha_0 = sys.argv[4]


eqtl_sss = ['100', '200', '300', '1000', '100-1000']


for eqtl_ss in eqtl_sss:

	output_file = mesc_expression_score_dir + simulation_name_string + '_' + eqtl_ss + '_' + 'lasso_' + alpha_0 + '_corr_replicate_eqtl_sumstats_xt_summary.txt'
	generate_xt_eqtl_summary_file(mesc_expression_score_dir,simulation_name_string, eqtl_ss, 'lasso_' + str(alpha_0), output_file)

	output_file = mesc_expression_score_dir + simulation_name_string + '_' + eqtl_ss + '_' + 'lasso_' + alpha_0 + '_corr_standardized_replicate_eqtl_sumstats_xt_summary.txt'
	generate_xt_eqtl_summary_file(mesc_expression_score_dir,simulation_name_string, eqtl_ss, 'lasso_' + str(alpha_0) + '_standardized', output_file)


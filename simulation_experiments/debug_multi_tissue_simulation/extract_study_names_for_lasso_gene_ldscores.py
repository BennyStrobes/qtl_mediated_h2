import numpy as np
import os
import sys
import pdb










simulation_name_string = sys.argv[1]
alpha_0 = sys.argv[2]
lasso_gene_models_dir = sys.argv[3]
study_names_file = sys.argv[4]


# Open output file handle
t = open(study_names_file,'w')
t.write('study_name\tstudy_lasso_file\n')

eqtl_ss_arr = ['100', '200', '300', '1000']


for tissue_iter in range(5):
	for replicate_name in ['replicate1', 'replicate2']:
		for eqtl_ss in eqtl_ss_arr:
			study_name = simulation_name_string + '_tissue' + str(tissue_iter) + '_' + eqtl_ss + '_' + replicate_name + '_lasso_' + alpha_0
			lasso_file = lasso_gene_models_dir + study_name + '_est_causal_effects.txt'
			t.write(study_name + '\t' + lasso_file + '\n')
t.close()
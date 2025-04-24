import numpy as np
import os
import sys
import pdb





def get_tissue_names(tissue_info_file):
	arr = []
	f = open(tissue_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue = data[0]
		arr.append(tissue)

	f.close()


	return np.asarray(arr)





def create_chrom_output_list(chrom_output_list, tissue_names, mesc_run_name, expression_score_dir, chrom_num):
	t = open(chrom_output_list,'w')

	for tissue_name in tissue_names:
		t.write(expression_score_dir + 'overall_' + mesc_run_name + '_' + tissue_name + '_' + str(chrom_num) + '\n')

	t.close()

	return



########################
# Command line args
########################
tissue_info_file = sys.argv[1]
mesc_run_name = sys.argv[2]
expression_score_dir = sys.argv[3]



# Extract ordered list of tissue names
tissue_names = get_tissue_names(tissue_info_file)

# Make seperate list for each chromosome
for chrom_num in range(1,23):
	chrom_output_list = expression_score_dir + mesc_run_name + '_meta_analysis_list_' + str(chrom_num) + '.txt'
	create_chrom_output_list(chrom_output_list, tissue_names, mesc_run_name, expression_score_dir, chrom_num)



import numpy as np 
import os
import sys
import pdb




def get_expr_file_sample_names(expr_file):
	f = open(expr_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		sample_names = data[3:]
		break
	f.close()

	return sample_names







expr_file = sys.argv[1]
plink_fam_file = sys.argv[2]


# Ordered sample names in expression file
expr_file_sample_names = get_expr_file_sample_names(expr_file)

# Ordered sample names in plink file
tmp_plink = np.loadtxt(plink_fam_file, dtype=str, delimiter='\t')
plink_sample_names = tmp_plink[:,1]

# Check if they don't match
if np.array_equal(plink_sample_names, expr_file_sample_names) == False:
	print('ERROR: expr and genotype file sample names do not match')
	pdb.set_trace()
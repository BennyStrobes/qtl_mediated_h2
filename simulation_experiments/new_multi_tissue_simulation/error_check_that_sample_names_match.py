import numpy as np
import os
import sys
import pdb






# Command line args
sample_names_expr_file = sys.argv[1]
plink_fam_file = sys.argv[2]


aa = []
f = open(sample_names_expr_file)
for line in f:
	line = line.rstrip()
	data = line.split()
	if data[0] != data[1]:
		print('assumptioneororor')
		pdb.set_trace()
	aa.append(data[0])
f.close()
aa = np.asarray(aa)

counter = 0
f = open(plink_fam_file)
for line in f:
	line = line.rstrip()
	data = line.split()
	if data[0] != data[1]:
		print('assumption eroorro')
		pdb.set_trace()
	if data[0] != aa[counter]:
		print('assumpriont oeroro')
		pdb.set_trace()

	counter = counter + 1

f.close()
import numpy as np
import os
import sys
import pdb








########
# Command line args
########
expr_file = sys.argv[1]
output_file = sys.argv[2]



f = open(expr_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	sample_names = np.asarray(data[3:])

	# only care about first line
	break

f.close()


t = open(output_file,'w')
for sample_name in sample_names:
	t.write(sample_name + '\t' + sample_name + '\n')
t.close()
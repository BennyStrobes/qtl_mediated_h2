import numpy as np
import os
import sys
import pdb










orig_sumstats_file = sys.argv[1]
new_sumstats_file = sys.argv[2]


t = open(new_sumstats_file,'w')
t.write('SNP\tN\tZ\tA1\tA2\n')

f = open(orig_sumstats_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	rsid = data[0]
	NN = 100000
	Z = data[7]
	A1 = data[3]
	A2 = data[4]
	t.write(rsid + '\t' + str(NN) + '\t' + Z + '\t' + A1 + '\t' + A2 + '\n')
t.close()
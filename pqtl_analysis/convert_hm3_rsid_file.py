import numpy as np
import os
import sys
import pdb








orig_hm3_file = sys.argv[1]
new_hm3_file = sys.argv[2]

f = open(orig_hm3_file)
t = open(new_hm3_file,'w')

head_count = 0

for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	t.write(data[0] + '\n')
f.close()
t.close()
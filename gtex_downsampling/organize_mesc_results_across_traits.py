import numpy as np
import os
import sys
import pdb











########################
# Command line args
########################
mesc_run_identifier = sys.argv[1]
non_redundent_summary_statistics_file = sys.argv[2]
mesc_results_dir = sys.argv[3]
output_file = sys.argv[4]


t = open(output_file,'w')
t.write('trait_name\tmed_h2\tmed_h2_se\n')

f = open(non_redundent_summary_statistics_file)
head_count = 0
means = []
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	trait_name = data[0]

	mesc_file = mesc_results_dir + trait_name + '_' + mesc_run_identifier + '.all.h2med'
	tmp = np.loadtxt(mesc_file,dtype=str,delimiter='\t')

	t.write(trait_name + '\t' + tmp[1,3] + '\t' + tmp[1,4] + '\n')


	means.append(float(tmp[1,3]))
f.close()
t.close()

means = np.asarray(means)
print(np.mean(means))




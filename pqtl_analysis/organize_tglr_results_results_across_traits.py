import numpy as np
import os
import sys
import pdb


def extract_trait_names(non_redundent_summary_statistics_file):
	trait_names = []
	f = open(non_redundent_summary_statistics_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait_names.append(data[0])
	f.close()

	return np.unique(np.asarray(trait_names))






#####################
# Command line args
#####################
run_identifier = sys.argv[1]
non_redundent_summary_statistics_file = sys.argv[2]
tglr_results_dir = sys.argv[3]
xt_output_file = sys.argv[4]

trait_names = extract_trait_names(non_redundent_summary_statistics_file)


# Open output file handle
t = open(xt_output_file,'w')
t.write('trait_name\tfrac_h2_med\tfrac_h2_med_lb\tfrac_h2_med_ub\tmed_h2\tmed_h2_lb\tmed_h2_ub\tnm_h2\tnm_h2_lb\tnm_h2_ub\n')
#t.write('trait_name\tfrac_h2_med\tfrac_h2_med_lb\tfrac_h2_med_ub\n')

for trait_name in trait_names:

	filer = tglr_results_dir + trait_name + '_' + run_identifier + '_tglr_frac_h2_med_5_50.txt'
	aa = np.loadtxt(filer,dtype=str,delimiter='\t')
	frac = float(aa[1,0])
	frac_se = float(aa[1,1])
	frac_lb = frac - (1.96*frac_se)
	frac_ub = frac + (1.96*frac_se)

	filer = tglr_results_dir + trait_name + '_' + run_identifier + '_tglr_h2_med_5_50.txt'
	aa = np.loadtxt(filer,dtype=str,delimiter='\t')
	med = float(aa[1,0])
	med_se = float(aa[1,1])
	med_lb = med - (1.96*med_se)
	med_ub = med + (1.96*med_se)


	filer = tglr_results_dir + trait_name + '_' + run_identifier + '_tglr_h2_nonmed_5_50.txt'
	aa = np.loadtxt(filer,dtype=str,delimiter='\t')
	nonmed = float(aa[1,0])
	nonmed_se = float(aa[1,1])
	nonmed_lb = nonmed - (1.96*nonmed_se)
	nonmed_ub = nonmed + (1.96*nonmed_se)

	t.write(trait_name + '\t' + str(frac) + '\t' + str(frac_lb) + '\t' + str(frac_ub) + '\t' + str(med) + '\t' + str(med_lb) + '\t' + str(med_ub) + '\t' + str(nonmed) + '\t' + str(nonmed_lb) + '\t' + str(nonmed_ub) + '\n')
t.close()
print(xt_output_file)

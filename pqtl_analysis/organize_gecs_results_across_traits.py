import numpy as np
import os
import sys
import pdb





def summarize_arr(arr, method, anno):
	meaner = np.mean(arr)
	print(method + '\t' + anno + '\t' + str(meaner))

	return










##########################
# Command line args
##########################
tglr_results_dir = sys.argv[1]
gecs_results_dir = sys.argv[2]
non_redundent_summary_statistics_file = sys.argv[3]

tglr_baselineLD_summary_file = tglr_results_dir + 'cross_trait_med_h2_summary_susie_inf_pmces_baselineLD_no_qtl.txt'
tglr_genotype_intercept_summary_file = tglr_results_dir + 'cross_trait_med_h2_summary_susie_inf_pmces_genotype_intercept.txt'


output_file = gecs_results_dir + 'gecs_results_summary.txt'
t = open(output_file,'w')
t.write('trait_name\tmethod\tnon_mediated_annotations\tfrac_h2_med\tfrac_h2_med_lb\tfrac_h2_med_ub\n')

traits = []
head_count = 0
method='tglr'
anno='baselineLD'
aa = []
f = open(tglr_baselineLD_summary_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	trait_name = data[0]
	if trait_name == 'UKB_460K.biochemistry_Creatinine':
		continue
	frac_h2_med = data[1]
	frac_h2_med_lb = data[2]
	frac_h2_med_ub = data[3]
	t.write(trait_name + '\t' + method + '\t' + anno + '\t' + frac_h2_med + '\t' + frac_h2_med_lb + '\t' + frac_h2_med_ub + '\n')
	traits.append(trait_name)
	aa.append(float(frac_h2_med))
f.close()
aa = np.asarray(aa)
summarize_arr(aa, method, anno)
traits = np.asarray(traits)

head_count = 0
method='tglr'
anno='genotype_intercept'
bb = []
ii=0
f = open(tglr_genotype_intercept_summary_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	trait_name = data[0]
	if trait_name == 'UKB_460K.biochemistry_Creatinine':
		continue
	if trait_name != traits[ii]:
		print('assumption erororo')
		pdb.set_trace()
	frac_h2_med = data[1]
	frac_h2_med_lb = data[2]
	frac_h2_med_ub = data[3]
	t.write(trait_name + '\t' + method + '\t' + anno + '\t' + frac_h2_med + '\t' + frac_h2_med_lb + '\t' + frac_h2_med_ub + '\n')
	#traits[trait_name] = 1
	bb.append(float(frac_h2_med))
	ii = ii + 1
f.close()
bb = np.asarray(bb)
summarize_arr(bb, method, anno)

trait_arr = np.asarray(traits)
method= 'gecs'
full_anno='baselineLD_no_qtl'
anno = 'baselineLD'
cc = []
for trait_name in trait_arr:
	if trait_name == 'UKB_460K.biochemistry_Creatinine':
		continue
	filer = gecs_results_dir + trait_name + '_regress_in_gene_windows_' + full_anno + '_gecs_mediated_h2_point_estimates.txt'
	tmper = np.loadtxt(filer,dtype=str,delimiter='\t')
	frac_med = tmper[1,1]
	t.write(trait_name + '\t' + method + '\t' + anno + '\t' + frac_med + '\t' + frac_med + '\t' + frac_med + '\n')
	cc.append(float(frac_med))

cc = np.asarray(cc)
summarize_arr(cc, method, anno)


method= 'gecs'
full_anno='genotype_intercept'
anno = 'genotype_intercept'
dd = []
for trait_name in trait_arr:
	if trait_name == 'UKB_460K.biochemistry_Creatinine':
		continue
	filer = gecs_results_dir + trait_name + '_regress_in_gene_windows_' + full_anno + '_gecs_mediated_h2_point_estimates.txt'
	tmper = np.loadtxt(filer,dtype=str,delimiter='\t')
	frac_med = tmper[1,1]
	t.write(trait_name + '\t' + method + '\t' + anno + '\t' + frac_med + '\t' + frac_med + '\t' + frac_med + '\n')
	dd.append(float(frac_med))
dd = np.asarray(dd)
summarize_arr(dd, method, anno)

t.close()
print(output_file)

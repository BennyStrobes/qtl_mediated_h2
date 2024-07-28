import numpy as np 
import os
import sys
import pdb








mediated_h2_results_dir = sys.argv[1]
visualize_med_h2_results_dir = sys.argv[2]
simulation_name_string = sys.argv[3]



###################### 
# Aggregate results
#####################
eqtl_sss = ['100','300','1000']

method_names = ['true_gene', 'pred_gene', 'ridge_regr_pred_gene', 'blup_pred_gene', 'marginal_gene']

t = open(visualize_med_h2_results_dir+ "med_h2_summary.txt",'w')
t.write('method\teqtl_ss\tmed_h2\tmed_h2_lb\tmed_h2_ub\n')

for eqtl_ss in eqtl_sss:
	arr = []
	for simulation_iter in range(1,200):
		results_file = mediated_h2_results_dir + 'simulation_' + str(simulation_iter) + '_' + simulation_name_string + 'eqtl_ss_' + str(eqtl_ss) + '_med_h2_ldsc_style_summary.txt'
		aa = np.loadtxt(results_file, dtype=str,delimiter= '\t')
		res = aa[1,3:13].astype(float)
		med_h2s = []
		med_h2s.append(float(aa[1,4]))
		med_h2s.append(float(aa[1,6]))
		med_h2s.append(float(aa[1,8]))
		med_h2s.append(float(aa[1,10]))
		med_h2s.append(float(aa[1,14]))
		arr.append(np.asarray(med_h2s))
	arr = np.asarray(arr)

	for ii,method_name in enumerate(method_names):
		meany = np.mean(arr[:,ii])
		sdev = np.std(arr[:,ii])/np.sqrt(len(arr[:,ii]))
		ub = meany + 1.96*sdev
		lb = meany - 1.96*sdev
		t.write(method_name + '\t' + str(eqtl_ss) + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')
t.close()




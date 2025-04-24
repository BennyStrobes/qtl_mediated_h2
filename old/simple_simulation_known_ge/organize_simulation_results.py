import numpy as np 
import os
import sys
import pdb








mediated_h2_results_dir = sys.argv[1]
visualize_med_h2_results_dir = sys.argv[2]
simulation_name_string = sys.argv[3]



###################### 
# Aggregate results for selection analysis
#####################
output_file = visualize_med_h2_results_dir+ 'med_h2_summary_selection_differences.txt'
t = open(output_file,'w')
t.write('method\tsimulated_selection\tmed_h2\tmed_h2_lb\tmed_h2_ub\n')

method_names = ['ldsc_standardize_ge', 'ldsc_raw_ge']
selection_types = ['neutral', 'negative']
for selection_type in selection_types:
	arr = []
	for simulation_iter in range(1,201):
		for sub_iter in range(1,11):
			results_file = mediated_h2_results_dir + 'simulation_' + str(simulation_iter) + '_' + str(sub_iter) + '_' + simulation_name_string + selection_type + '_selection_med_h2_ldsc_style_summary.txt'
			aa = np.loadtxt(results_file, dtype=str,delimiter= '\t')
			med_h2s = []
			med_h2s.append(float(aa[1,4]))
			med_h2s.append(float(aa[1,6]))
			arr.append(np.asarray(med_h2s))
	arr = np.asarray(arr)

	for ii,method_name in enumerate(method_names):
		meany = np.mean(arr[:,ii])
		sdev = np.std(arr[:,ii])/np.sqrt(len(arr[:,ii]))
		ub = meany + 1.96*sdev
		lb = meany - 1.96*sdev
		t.write(method_name + '\t' + selection_type + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')
t.close()
print(output_file)

###################### 
# Aggregate results for attenuation analysis
#####################
noises = ['0.0', '0.0001', '0.001', '0.01', '0.1']
arr_standard = []
arr_correction = []
arr_correction2 = []
selection_type='neutral'
for simulation_iter in range(2,201):
	for sub_iter in range(1,11):
		results_file = mediated_h2_results_dir + 'simulation_' + str(simulation_iter) + '_' + str(sub_iter) + '_' + simulation_name_string + selection_type + '_selection_with_simulated_noise_med_h2_ldsc_style_summary.txt'
		aa = np.loadtxt(results_file, dtype=str,delimiter= '\t')


		arr_standard.append(aa[1,3:][::3].astype(float))
		arr_correction.append(aa[1,4:][::3].astype(float))
		arr_correction2.append(aa[1,5:][::3].astype(float))
arr_standard = np.asarray(arr_standard)
arr_correction = np.asarray(arr_correction)
arr_correction2 = np.asarray(arr_correction2)

output_file = visualize_med_h2_results_dir+ 'med_h2_summary_simulated_noise_correction.txt'
t = open(output_file,'w')
t.write('method\tsimulated_noise\tmed_h2\tmed_h2_lb\tmed_h2_ub\n')
agg_arr = [arr_standard, arr_correction, arr_correction2]
method_names = ['no_correction', 'correction_1', 'correction2']
for ii,method_name in enumerate(method_names):
	arr = agg_arr[ii]
	for jj, noise in enumerate(noises):
		meany = np.mean(arr[:,jj])
		sdev = np.std(arr[:,jj])/np.sqrt(len(arr[:,jj]))
		ub = meany + 1.96*sdev
		lb = meany - 1.96*sdev
		t.write(method_name + '\t' + noise + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')	

t.close()
print(output_file)









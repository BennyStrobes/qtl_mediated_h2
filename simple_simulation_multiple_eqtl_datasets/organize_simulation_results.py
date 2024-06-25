import numpy as np 
import os
import sys
import pdb








mediated_h2_results_dir = sys.argv[1]
visualize_med_h2_results_dir = sys.argv[2]
simulation_name_string = sys.argv[3]
n_eqtl_datasets = int(sys.argv[4])





###################### 
# Aggregate results for med h2 analysis with various methods
#####################
eqtl_sss = ['100', '300', '1000']

output_file = visualize_med_h2_results_dir+ 'med_h2_summary.txt'
t = open(output_file,'w')
t.write('method\teqtl_ss\tmed_h2\tmed_h2_lb\tmed_h2_ub\n')
sim_med_h2s = []

for eqtl_ss in eqtl_sss:
	arr = []
	selection_type='neutral'
	for simulation_iter in range(1,101):
		for sub_iter in range(1,6):
			results_file = mediated_h2_results_dir + 'simulation_' + str(simulation_iter) + '_' + str(sub_iter) + '_' + simulation_name_string + selection_type + '_selection_eqtl_ss_' + str(eqtl_ss) + 'med_h2_ldsc_style_summary.txt'
			aa = np.loadtxt(results_file, dtype=str,delimiter= '\t')
			tmp = []
			for ele in aa[1,3:]:
				tmp.append(float(ele.split(',')[0]) + float(ele.split(',')[1]))
			arr.append(np.asarray(tmp))
			method_names = aa[0,3:]
	arr = np.asarray(arr)

	for ii, method_name in enumerate(method_names):
		meany = np.mean(arr[:,ii])
		sdev = np.std(arr[:,ii])/np.sqrt(len(arr[:,ii]))
		ub = meany + 1.96*sdev
		lb = meany - 1.96*sdev
		t.write(method_name + '\t' + eqtl_ss + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')

t.close()
print(output_file)




###################### 
# Aggregate results for med h2 analysis with various methods
#####################
eqtl_sss = ['100', '300', '1000']

output_file = visualize_med_h2_results_dir+ 'med_h2_summary_individual_eqtl_datasets.txt'
t = open(output_file,'w')
t.write('eqtl_dataset_index\tmethod\teqtl_ss\tmed_h2\tmed_h2_lb\tmed_h2_ub\n')
sim_med_h2s = []

for eqtl_dataset_iter in range(n_eqtl_datasets): 
	for eqtl_ss in eqtl_sss:
		arr = []
		selection_type='neutral'
		for simulation_iter in range(1,101):
			for sub_iter in range(1,6):
				results_file = mediated_h2_results_dir + 'simulation_' + str(simulation_iter) + '_' + str(sub_iter) + '_' + simulation_name_string + selection_type + '_selection_eqtl_ss_' + str(eqtl_ss) + 'med_h2_ldsc_style_summary.txt'
				aa = np.loadtxt(results_file, dtype=str,delimiter= '\t')
				tmp = []
				for ele in aa[1,3:]:
					tmp.append(float(ele.split(',')[eqtl_dataset_iter]))
				arr.append(np.asarray(tmp))
				method_names = aa[0,3:]
		arr = np.asarray(arr)

		for ii, method_name in enumerate(method_names):
			meany = np.mean(arr[:,ii])
			sdev = np.std(arr[:,ii])/np.sqrt(len(arr[:,ii]))
			ub = meany + 1.96*sdev
			lb = meany - 1.96*sdev
			t.write(str(eqtl_dataset_iter + 1) + '\t' + method_name + '\t' + eqtl_ss + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')

t.close()
print(output_file)
sim_med_h2s = np.asarray(sim_med_h2s)


'''
###################### 
# Aggregate results for noise estimation
#####################
eqtl_sss = ['100', '300', '1000']

output_file = visualize_med_h2_results_dir+ 'noise_estimates.txt'
t = open(output_file,'w')
t.write('noise_type\teqtl_ss\tnoise_mean\tnoise_lb\tnoise_ub\n')


for eqtl_ss in eqtl_sss:
	arr_true_noise = []
	arr_est_noise = []
	selection_type='neutral'
	for simulation_iter in range(401,501):
		for sub_iter in range(1,5):
			results_file = mediated_h2_results_dir + 'simulation_' + str(simulation_iter) + '_' + str(sub_iter) + '_' + simulation_name_string + selection_type + '_selection_eqtl_ss_' + str(eqtl_ss) + 'med_h2_ldsc_style_summary.txt'
			aa = np.loadtxt(results_file, dtype=str,delimiter= '\t')
			arr_true_noise.append(aa[1,3].astype(float))
			arr_est_noise.append(aa[1,4].astype(float))
	arr_true_noise = np.asarray(arr_true_noise)
	arr_est_noise = np.asarray(arr_est_noise)

	method_names = ['true_noise', 'est_noise']
	meany = np.mean(arr_true_noise)
	sdev = np.std(arr_true_noise)/np.sqrt(len(arr_true_noise))
	ub = meany + 1.96*sdev
	lb = meany - 1.96*sdev
	t.write(method_names[0] + '\t' + eqtl_ss + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')

	meany = np.mean(arr_est_noise)
	sdev = np.std(arr_est_noise)/np.sqrt(len(arr_est_noise))
	ub = meany + 1.96*sdev
	lb = meany - 1.96*sdev
	t.write(method_names[1] + '\t' + eqtl_ss + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')	

t.close()
print(output_file)
'''










import numpy as np
import os
import sys
import pdb




def organize_results_across_sims(no_ld_simulation_res_dir, n_sims, simulation_name, simulation_suffix, output_file):
	t = open(output_file,'w')

	for sim_iter in range(2,n_sims+1):
		sim_file = no_ld_simulation_res_dir + 'sim_' + str(sim_iter) + '_' + simulation_name + '_' + simulation_suffix
		f = open(sim_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				if sim_iter == 1:
					t.write(line + '\n')
				continue
			t.write(line + '\n')
	t.close()
	return

def get_mean_param_estimates_1_var(output_file, output_file2):
	method_names = ['linear_reg_unobsered_true_effect_sizes', 'scipy_odr_effect_est', 'scipy_odr_effect_est_stringent', 'scipy_odr_effect_est_per_snp_var', 'scipy_odr_effect_est_stringent_per_snp_var']
	raw_data = np.loadtxt(output_file,dtype=str,delimiter='\t')
	t = open(output_file2,'w')
	t.write('method\tparam\tsim\test\test_95_lb\test_95_ub\n')

	for method_name in method_names:

		indices=raw_data[:,0]==method_name

		param_name ='alpha_1'
		sim_alpha = np.unique(raw_data[indices,2])[0]
		est_alpha = raw_data[indices,3].astype(float)
		meany = np.mean(est_alpha)
		se = np.std(est_alpha)/np.sqrt(len(est_alpha))
		ub = meany + (1.96*se)
		lb = meany - (1.96*se)

		t.write(method_name + '\t' + param_name + '\t' + str(sim_alpha) + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')

	t.close()
	print(output_file2)
	return


def get_mean_param_estimates(output_file, output_file2):
	method_names = ['linear_reg_unobsered_true_effect_sizes', 'scipy_odr_effect_est', 'scipy_odr_effect_est_stringent', 'scipy_odr_effect_est_per_snp_var', 'scipy_odr_effect_est_stringent_per_snp_var']
	raw_data = np.loadtxt(output_file,dtype=str,delimiter='\t')
	t = open(output_file2,'w')
	t.write('method\tparam\tsim\test\test_95_lb\test_95_ub\n')

	for method_name in method_names:

		indices=raw_data[:,0]==method_name

		param_name ='alpha_1'
		sim_alpha = np.unique(raw_data[indices,2])[0]
		est_alpha = raw_data[indices,3].astype(float)
		meany = np.mean(est_alpha)
		se = np.std(est_alpha)/np.sqrt(len(est_alpha))
		ub = meany + (1.96*se)
		lb = meany - (1.96*se)

		t.write(method_name + '\t' + param_name + '\t' + str(sim_alpha) + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')


		param_name ='alpha_2'
		sim_alpha = np.unique(raw_data[indices,4])[0]
		est_alpha = raw_data[indices,5].astype(float)
		meany = np.mean(est_alpha)
		se = np.std(est_alpha)/np.sqrt(len(est_alpha))
		ub = meany + (1.96*se)
		lb = meany - (1.96*se)

		t.write(method_name + '\t' + param_name + '\t' + str(sim_alpha) + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')


		param_name ='alpha_3'
		sim_alpha = np.unique(raw_data[indices,6])[0]
		est_alpha = raw_data[indices,7].astype(float)
		meany = np.mean(est_alpha)
		se = np.std(est_alpha)/np.sqrt(len(est_alpha))
		ub = meany + (1.96*se)
		lb = meany - (1.96*se)

		t.write(method_name + '\t' + param_name + '\t' + str(sim_alpha) + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')

	t.close()
	print(output_file2)
	return


def organize_polygenic_simulation_results(n_sims, eqtl_sss, no_ld_simulation_res_dir, output_file):
	t = open(output_file,'w')
	t.write('eqtl_ss\teffect_type\testimation_method\tbias\tbias_lb\tbias_ub\n')
	for eqtl_ss in eqtl_sss:
		med = []
		nm = []
		mean_med = []
		mean_nm = []
		for sim_iter in range(1, n_sims+1):
			filer = no_ld_simulation_res_dir + 'sim_' + str(sim_iter) + '_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_' + str(eqtl_ss) + '_eqtl_architecture_polygenic_posterior_eqtl_effect_est_res_summary.txt'
			data = np.loadtxt(filer,dtype=str,delimiter='\t')
			sim_nm = float(data[1,2])
			sim_med = float(data[1,3])
			est_nm = float(data[1,4])
			est_med = float(data[1,5])
			mean_est_nm = float(data[2,4])
			mean_est_med = float(data[2,5])
			nm.append(sim_nm -est_nm)
			med.append(sim_med-est_med)

			mean_nm.append(sim_nm-mean_est_nm)
			mean_med.append(sim_med-mean_est_med)
		nm = -np.asarray(nm)
		med = -np.asarray(med)

		nm2 = -np.asarray(mean_nm)
		med2 = -np.asarray(mean_med)


		meany = np.mean(nm)
		ub = meany + 1.96*np.std(nm)/np.sqrt(len(nm))
		lb = meany - 1.96*np.std(nm)/np.sqrt(len(nm))
		t.write(str(eqtl_ss) + '\t' + 'non-mediated\tposterior_distr\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')


		meany = np.mean(med)
		ub = meany + 1.96*np.std(med)/np.sqrt(len(med))
		lb = meany - 1.96*np.std(med)/np.sqrt(len(med))
		t.write(str(eqtl_ss) + '\t' + 'mediated\tposterior_distr\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')	


		meany = np.mean(nm2)
		ub = meany + 1.96*np.std(nm2)/np.sqrt(len(nm2))
		lb = meany - 1.96*np.std(nm2)/np.sqrt(len(nm2))
		t.write(str(eqtl_ss) + '\t' + 'non-mediated\tpmces\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')


		meany = np.mean(med2)
		ub = meany + 1.96*np.std(med2)/np.sqrt(len(med2))
		lb = meany - 1.96*np.std(med2)/np.sqrt(len(med2))
		t.write(str(eqtl_ss) + '\t' + 'mediated\tpmces\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')	




	t.close()

	return 



def organize_sparse_simulation_results(n_sims, eqtl_sss, no_ld_simulation_res_dir, output_file):
	t = open(output_file,'w')
	t.write('eqtl_ss\teffect_type\testimation_method\tbias\tbias_lb\tbias_ub\n')
	for eqtl_ss in eqtl_sss:
		med = []
		nm = []
		mean_med = []
		mean_nm = []
		for sim_iter in range(200, 200+n_sims+1):
			filer = no_ld_simulation_res_dir + 'sim_' + str(sim_iter) + '_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_' + str(eqtl_ss) + '_eqtl_architecture_sparse_posterior_eqtl_effect_est_res_summary.txt'
			data = np.loadtxt(filer,dtype=str,delimiter='\t')
			sim_nm = float(data[2,2])
			sim_med = float(data[2,3])
			est_nm = float(data[2,4])
			est_med = float(data[2,5])
			mean_est_nm = float(data[3,4])
			mean_est_med = float(data[3,5])
			nm.append(sim_nm -est_nm)
			med.append(sim_med-est_med)

			mean_nm.append(sim_nm-mean_est_nm)
			mean_med.append(sim_med-mean_est_med)
		nm = -np.asarray(nm)
		med = -np.asarray(med)

		nm2 = -np.asarray(mean_nm)
		med2 = -np.asarray(mean_med)


		meany = np.mean(nm)
		ub = meany + 1.96*np.std(nm)/np.sqrt(len(nm))
		lb = meany - 1.96*np.std(nm)/np.sqrt(len(nm))
		t.write(str(eqtl_ss) + '\t' + 'non-mediated\tposterior_distr\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')


		meany = np.mean(med)
		ub = meany + 1.96*np.std(med)/np.sqrt(len(med))
		lb = meany - 1.96*np.std(med)/np.sqrt(len(med))
		t.write(str(eqtl_ss) + '\t' + 'mediated\tposterior_distr\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')	


		meany = np.mean(nm2)
		ub = meany + 1.96*np.std(nm2)/np.sqrt(len(nm2))
		lb = meany - 1.96*np.std(nm2)/np.sqrt(len(nm2))
		t.write(str(eqtl_ss) + '\t' + 'non-mediated\tpmces\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')


		meany = np.mean(med2)
		ub = meany + 1.96*np.std(med2)/np.sqrt(len(med2))
		lb = meany - 1.96*np.std(med2)/np.sqrt(len(med2))
		t.write(str(eqtl_ss) + '\t' + 'mediated\tpmces\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')	




	t.close()


####################
# Command line args
####################
no_ld_simulation_res_dir = sys.argv[1]


n_sims=100

eqtl_sss = ['100', '300']
output_file = no_ld_simulation_res_dir + 'organized_results_polygenic_simulation.txt'
#organize_polygenic_simulation_results(n_sims, eqtl_sss, no_ld_simulation_res_dir, output_file)

eqtl_sss = ['300']
output_file = no_ld_simulation_res_dir + 'organized_results_sparse_simulation.txt'
organize_sparse_simulation_results(n_sims, eqtl_sss, no_ld_simulation_res_dir, output_file)

'''
eqtl_sss = ['300', '1000']
for eqtl_ss in eqtl_sss:
	simulation_name='no_ld_no_nm_variants_1_caus_tiss_eqtlss_' + eqtl_ss
	simulation_suffix = 'effect_est_res_summary.txt'
	output_file = no_ld_simulation_res_dir + 'merged_' + simulation_name + '_' + simulation_suffix
	organize_results_across_sims(no_ld_simulation_res_dir, n_sims, simulation_name, simulation_suffix, output_file)
	output_file2 = no_ld_simulation_res_dir + 'merged_' + simulation_name + '_' + 'param_summary.txt'
	get_mean_param_estimates_1_var(output_file, output_file2)



eqtl_sss = ['300', '1000']
sim_corrs = ['0.0', '0.7']
for eqtl_ss in eqtl_sss:
	for sim_corr in sim_corrs:
		simulation_name='no_ld_no_nm_variants_3_caus_tiss_eqtlss_' + eqtl_ss + '_sim_corr_' + sim_corr
		simulation_suffix = 'effect_est_res_summary.txt'
		output_file = no_ld_simulation_res_dir + 'merged_' + simulation_name + '_' + simulation_suffix
		organize_results_across_sims(no_ld_simulation_res_dir, n_sims, simulation_name, simulation_suffix, output_file)
		output_file2 = no_ld_simulation_res_dir + 'merged_' + simulation_name + '_' + 'param_summary.txt'
		get_mean_param_estimates(output_file, output_file2)
'''


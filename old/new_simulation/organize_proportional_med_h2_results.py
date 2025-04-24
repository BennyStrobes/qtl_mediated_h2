import numpy as np
import os
import sys
import pdb






def organize_simulation_results(summary_file, gene_frac, med_h2_simulation_dir, cis_h2s, eqtl_sss):
	t = open(summary_file,'w')

	cis_h2 = cis_h2s[0]
	t.write('method_name\teqtl_ss\theritability_type\ttrue_h2\test_h2\test_h2_lb\test_h2_ub\n')

	for eqtl_ss in eqtl_sss:


		ss_reml = []
		ss_reml_nm = []
		ss_reml_2_step = []
		ss_reml_2_step_nm = []
		ss_reml_eqtl_h2 = []
		true_med = []
		true_nm = []
		true_eqtl_h2 = []

		for sim_iter in range(1,100):
			filer = med_h2_simulation_dir +'proportional_med_h2_simulation_' + str(sim_iter) + '_gwas_ss_20000_eqtl_ss_' + str(eqtl_ss) + '_med_h2_0.1_nm_h2_0.3_mean_cis_h2_' + cis_h2 + '_gene_frac_' + gene_frac + '_arch_sparse_est_summary.txt'
			f = open(filer)
			for line in f:
				line = line.rstrip()
				data = line.split('\t')

				if data[1] == 'sumstat_reml':
					med_bias = float(data[6])
					ss_reml.append(med_bias)
					nm_bias = float(data[5])
					ss_reml_nm.append(nm_bias)
					true_med.append(float(data[3]))
					true_nm.append(float(data[2]))
					ss_reml_eqtl_h2.append(float(data[7]))
					true_eqtl_h2.append(float(data[4]))

				if data[1] == 'sumstat_reml_two_step':
					med_bias = float(data[6])
					ss_reml_2_step.append(med_bias)
					nm_bias = float(data[5])
					ss_reml_2_step_nm.append(nm_bias)


			f.close()
		ss_reml = np.asarray(ss_reml)
		ss_reml_2_step = np.asarray(ss_reml_2_step)
		ss_reml_nm = np.asarray(ss_reml_nm)
		ss_reml_2_step_nm = np.asarray(ss_reml_2_step_nm)
		true_med = np.asarray(true_med)
		true_nm = np.asarray(true_nm)
		true_eqtl_h2 = np.asarray(true_eqtl_h2)
		ss_reml_eqtl_h2 = np.asarray(ss_reml_eqtl_h2)


		t.write('ss_reml\t' + eqtl_ss + '\t' + 'med' + '\t' + '0.1' + '\t' + str(np.mean(ss_reml)) + '\t' + str(np.mean(ss_reml) - (1.96*np.std(ss_reml)/np.sqrt(len(ss_reml))))  + '\t' + str(np.mean(ss_reml) + (1.96*np.std(ss_reml)/np.sqrt(len(ss_reml)))) + '\n')
		#t.write('ss_reml_2_step\t' + cis_h2 + '\t' + 'med' + '\t' '0.1' + '\t' + str(np.mean(ss_reml_2_step)) + '\t' + str(np.mean(ss_reml_2_step) - (1.96*np.std(ss_reml_2_step)/np.sqrt(len(ss_reml_2_step))))  + '\t' + str(np.mean(ss_reml_2_step) + (1.96*np.std(ss_reml_2_step)/np.sqrt(len(ss_reml_2_step)))) + '\n')

		t.write('ss_reml\t' + eqtl_ss + '\t' + 'non-med' + '\t' + '0.3' + '\t' + str(np.mean(ss_reml_nm)) + '\t' + str(np.mean(ss_reml_nm) - (1.96*np.std(ss_reml_nm)/np.sqrt(len(ss_reml_nm))))  + '\t' + str(np.mean(ss_reml_nm) + (1.96*np.std(ss_reml_nm)/np.sqrt(len(ss_reml_nm)))) + '\n')
		#t.write('ss_reml_2_step\t' + cis_h2 + '\t' + 'non-med' + '\t' '0.3' + '\t' + str(np.mean(ss_reml_2_step_nm)) + '\t' + str(np.mean(ss_reml_2_step_nm) - (1.96*np.std(ss_reml_2_step_nm)/np.sqrt(len(ss_reml_2_step_nm))))  + '\t' + str(np.mean(ss_reml_2_step_nm) + (1.96*np.std(ss_reml_2_step_nm)/np.sqrt(len(ss_reml_2_step_nm)))) + '\n')
		

		t.write('ss_reml\t' + eqtl_ss + '\t' + 'eqtl-h2' + '\t' + str(np.mean(true_eqtl_h2)) + '\t' + str(np.mean(ss_reml_eqtl_h2)) + '\t' + str(np.mean(ss_reml_eqtl_h2) - (1.96*np.std(ss_reml_eqtl_h2)/np.sqrt(len(ss_reml_eqtl_h2))))  + '\t' + str(np.mean(ss_reml_eqtl_h2) + (1.96*np.std(ss_reml_eqtl_h2)/np.sqrt(len(ss_reml_eqtl_h2)))) + '\n')

		'''
		ss_bayes = []
		ss_bayes_nm = []
		true_med = []
		true_nm = []
		for sim_iter in range(11,100):
			filer = med_h2_simulation_dir +'proportional_med_h2_simulation_' + str(sim_iter) + '_gwas_ss_10000_eqtl_ss_500_med_h2_0.1_nm_h2_0.3_mean_cis_h2_' + cis_h2 + '_gene_frac_' + gene_frac + '_arch_sparse_est_summary.txt'
			if os.path.isfile(filer) == False:
				continue
			f = open(filer)
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if data[1] == 'sumstat_bayesian':
					med_bias = float(data[6])
					ss_bayes.append(med_bias)
					nm_bias = float(data[5])
					ss_bayes_nm.append(nm_bias)
					true_med.append(float(data[3]))
					true_nm.append(float(data[2]))					
			f.close()
		ss_bayes = np.asarray(ss_bayes)
		ss_bayes_nm = np.asarray(ss_bayes_nm)
		true_med = np.asarray(true_med)
		true_nm = np.asarray(true_nm)

		t.write('ss_bayes_lmm\t' + cis_h2 + '\t' + 'med' + '\t' '0.1' + '\t' + str(np.mean(ss_bayes)) + '\t' + str(np.mean(ss_bayes) - (1.96*np.std(ss_bayes)/np.sqrt(len(ss_bayes))))  + '\t' + str(np.mean(ss_bayes) + (1.96*np.std(ss_bayes)/np.sqrt(len(ss_bayes)))) + '\n')
		t.write('ss_bayes_lmm\t' + cis_h2 + '\t' + 'non-med' + '\t' '0.3' + '\t' + str(np.mean(ss_bayes_nm)) + '\t' + str(np.mean(ss_bayes_nm) - (1.96*np.std(ss_bayes_nm)/np.sqrt(len(ss_bayes_nm))))  + '\t' + str(np.mean(ss_bayes_nm) + (1.96*np.std(ss_bayes_nm)/np.sqrt(len(ss_bayes_nm)))) + '\n')
		'''
	t.close()
	return





########################
# Command line args
########################
med_h2_simulation_dir = sys.argv[1]



gene_frac = '0.5'
cis_h2s = ['0.05']
eqtl_sss = ['500', '1000', '5000']
summary_file = med_h2_simulation_dir + 'mediation_h2_summary.txt'
organize_simulation_results(summary_file, gene_frac, med_h2_simulation_dir, cis_h2s, eqtl_sss)
print(summary_file)

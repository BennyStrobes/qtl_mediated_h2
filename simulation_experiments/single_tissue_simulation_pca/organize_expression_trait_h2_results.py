import numpy as np
import os
import sys
import pdb







#####################
# Command line args
#####################
trait_med_h2_inference_dir = sys.argv[1]
visualize_expression_trait_h2 = sys.argv[2]

# Results summary 
output_file = visualize_expression_trait_h2 + 'expression_trait_h2_results_summary.txt'
t = open(output_file,'w')
t.write('method_name\th2_version\teQTL_SS\ttrue_mean\test_mean\test_mean_lb\test_mean_ub\n')

eqtl_sample_sizes = ['100','300', '1000', '10000']
n_sims = 30


broad_methods = ['VI_rss_lmm_multivariate_ig_prior', 'VI_rss_lmm_univariate_ig_prior', 'gibbs_rss_lmm_multivariate_ig_prior', 'gibbs_rss_lmm_univariate_ig_prior']
clean_broad_methods = ['VI_lmm_multivariate_ig_prior', 'VI_lmm_univariate_ig_prior', 'gibbs_lmm_multivariate_ig_prior', 'gibbs_lmm_univariate_ig_prior']

methods = []
clean_method_names = []
priors = ['0.0', '1e-10', '1e-06', '0.001']
for prior in priors:
	for ii, broad_method in enumerate(broad_methods):
		clean_broad_method = clean_broad_methods[ii]
		methods.append(broad_method + '_' + prior)
		clean_method_names.append(clean_broad_method + '_' + prior)
methods = np.asarray(methods)
clean_method_names = np.asarray(clean_method_names)


for ii, method in enumerate(methods):
	clean_method_name = clean_method_names[ii]
	for eqtl_sample_size in eqtl_sample_sizes:
		est_h2 = []
		est_full_h2 = []
		obs_h2 = []
		for sim_iter in range(1, n_sims+1):
			filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_iter) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_eqtl_' + str(eqtl_sample_size) + '_small_univariate_expression_trait_h2_inference.txt'
			f = open(filer)
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if data[0] != method:
					continue
				if float(data[2]) > .1:
					continue
				if data[3] == 'NA':
					continue
				obs_h2.append(float(data[2]))
				est_h2.append(float(data[3]))
				est_full_h2.append(float(data[4]))
			f.close()

		mean = np.mean(est_h2)
		mean_lb = mean - (1.96*np.std(est_h2)/np.sqrt(len(est_h2)))
		mean_ub = mean + (1.96*np.std(est_h2)/np.sqrt(len(est_h2)))
		t.write(clean_method_name + '\t' + 'iid_h2' + '\t' + eqtl_sample_size + '\t' + str(np.mean(obs_h2)) + '\t' + str(mean) + '\t' + str(mean_lb) + '\t' + str(mean_ub) + '\n')

		mean = np.mean(est_full_h2)
		mean_lb = mean - (1.96*np.std(est_full_h2)/np.sqrt(len(est_full_h2)))
		mean_ub = mean + (1.96*np.std(est_full_h2)/np.sqrt(len(est_full_h2)))
		t.write(clean_method_name + '\t' + 'joint_h2'+ '\t' + eqtl_sample_size + '\t' + str(np.mean(obs_h2)) + '\t' + str(mean) + '\t' + str(mean_lb) + '\t' + str(mean_ub) + '\n')
t.close()
print(output_file)



'''
methods = ['custom_gibbs_rss_lmm_ig_prior_0', 'vi_individual_const_resid_var_prior_0', 'vi_individual_const_resid_var_prior_1e-3', 'vi_individual_const_resid_var_prior_1e-6', 'vi_individual_const_resid_var_prior_1e-9', 'vi_individual_const_resid_var_prior_1e-14']
clean_method_names = ['gibbs_lmm_ig_prior_0', 'vi_lmm_ig_prior_0', 'vi_lmm_ig_prior_1e-3', 'vi_lmm_ig_prior_1e-6', 'vi_lmm_ig_prior_1e-9', 'vi_lmm_ig_prior_1e-14']




for ii, method in enumerate(np.asarray(methods)):
	clean_method_name = clean_method_names[ii]
	for eqtl_sample_size in eqtl_sample_sizes:
		est_h2 = []
		est_full_h2 = []
		obs_h2 = []
		for sim_iter in range(1, n_sims+1):
			filer = trait_med_h2_inference_dir +'simulation_' + str(sim_iter) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_eqtl_' + eqtl_sample_size + '_small_univariate_expression_trait_h2_inference_v2.txt'
			f = open(filer)
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if data[0] != method:
					continue
				if float(data[2]) > .1:
					continue
				obs_h2.append(float(data[2]))
				est_h2.append(float(data[3]))
				est_full_h2.append(float(data[4]))
			f.close()

		mean = np.mean(est_h2)
		mean_lb = mean - (1.96*np.std(est_h2)/np.sqrt(len(est_h2)))
		mean_ub = mean + (1.96*np.std(est_h2)/np.sqrt(len(est_h2)))
		t.write(clean_method_name + '\t' + 'iid_h2' + '\t' + eqtl_sample_size + '\t' + str(np.mean(obs_h2)) + '\t' + str(mean) + '\t' + str(mean_lb) + '\t' + str(mean_ub) + '\n')

		mean = np.mean(est_full_h2)
		mean_lb = mean - (1.96*np.std(est_full_h2)/np.sqrt(len(est_full_h2)))
		mean_ub = mean + (1.96*np.std(est_full_h2)/np.sqrt(len(est_full_h2)))
		t.write(clean_method_name + '\t' + 'joint_h2'+ '\t' + eqtl_sample_size + '\t' + str(np.mean(obs_h2)) + '\t' + str(mean) + '\t' + str(mean_lb) + '\t' + str(mean_ub) + '\n')

t.close()

print(output_file)
'''
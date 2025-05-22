import numpy as np
import os
import sys
import pdb







def create_mapping_from_simulation_number_to_simulated_h2_parameters(trait_med_h2_inference_dir, sim_nums):
	mapping = {}

	for sim_num in sim_nums:
		filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_300_joint_ldsc_multimethod4.txt'
	
		if os.path.isfile(filer) == False:
			mapping[sim_num] = (total_h2, nm_h2, med_h2)
			continue

		f = open(filer)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			total_h2 = float(data[2])
			nm_h2 = float(data[4])
			med_h2 = float(data[3])
			mapping[sim_num] = (total_h2, nm_h2, med_h2)
			break
		f.close()
	return mapping


def extract_sampled_params_from_marginal_gibbs_sampler_file(marginal_gibbs_sampler_file, burn_in_lb=15000, burn_in_ub=5000000000000):
	est_nm_h2 = []
	est_med_h2 = []
	est_total_h2 = []
	est_med_h2_per_tissue = []

	f = open(marginal_gibbs_sampler_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		itera = float(data[3])

		if itera <= burn_in_lb or itera > burn_in_ub:
			continue

		# Extract relevent params
		nm_h2 = float(data[4])
		med_h2 = float(data[5])
		total_h2 = float(data[6])
		med_h2_per_tissue = np.asarray(data[8].split(';')).astype(float)

		est_nm_h2.append(nm_h2)
		est_med_h2.append(med_h2)
		est_total_h2.append(total_h2)
		est_med_h2_per_tissue.append(med_h2_per_tissue)

	f.close()

	return np.asarray(est_nm_h2), np.asarray(est_med_h2), np.asarray(est_total_h2), np.asarray(est_med_h2_per_tissue)

def extract_sampled_params_from_rss_gibbs_sampler_file(rss_gibbs_sampler_file, burn_in_thresh=15000):
	est_nm_h2 = []
	alt_est_nm_h2 = []
	est_med_h2 = []
	alt_est_med_h2 = []
	est_total_h2 = []

	f = open(rss_gibbs_sampler_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		if float(data[0]) < burn_in_thresh:
			continue

		# Extract relevent params
		nm_h2 = float(data[1])
		alt_nm_h2 = float(data[2])
		med_h2 = float(data[3])
		alt_med_h2 = float(data[5])
		total_h2 = alt_nm_h2 + med_h2

		est_nm_h2.append(nm_h2)
		alt_est_nm_h2.append(alt_nm_h2)
		est_med_h2.append(med_h2)
		alt_est_med_h2.append(alt_med_h2)
		est_total_h2.append(total_h2)

	f.close()

	return np.asarray(est_nm_h2), np.asarray(alt_est_nm_h2), np.asarray(est_med_h2), np.asarray(alt_est_med_h2), np.asarray(est_total_h2)


def get_95_percent_ci_based_on_sampled_distribution(arr):
	lower_index = int(np.floor(len(arr)*.025))
	upper_index = int(np.floor(len(arr)*.975))

	sorted_arr = np.sort(arr)
	lower_bound = sorted_arr[lower_index]
	upper_bound = sorted_arr[upper_index]

	return lower_bound, upper_bound


def generate_per_iteration_results_summary_file(sim_nums, trait_med_h2_inference_dir, results_summary_file,eqtl_sample_sizes, sim_num_to_sim_h2_params, iteration_lb, iteration_ub):
	t = open(results_summary_file,'w')
	t.write('iteration\teqtl_sample_size\tmethod\tsim_total_h2\tsim_nm_h2\tsim_med_h2\test_nm_h2\test_med_h2\test_total_h2')
	for tissue_num in range(5):
		t.write('\t' + 'est_med_h2_tissue' + str(tissue_num))
	t.write('\n')

	# Loop through simulations
	for sim_iter in sim_nums:
		sim_total_h2, sim_nm_h2, sim_med_h2 = sim_num_to_sim_h2_params[sim_iter]

		# Loop through eqtl sample sizes
		for eqtl_sample_size in eqtl_sample_sizes:

			###########################################
			# First do method of marginal PCA gibbs sampler
			marginal_gibbs_sampler_file = trait_med_h2_inference_dir + 'simulation_' + str(sim_iter) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_xt_cov_resid_var_True_cc_0.0.txt'
			
			if os.path.isfile(marginal_gibbs_sampler_file) == False:
				continue

			# Extract distributions
			est_nm_h2, est_med_h2, est_total_h2, est_med_h2_per_tissue = extract_sampled_params_from_marginal_gibbs_sampler_file(marginal_gibbs_sampler_file, burn_in_lb=iteration_lb, burn_in_ub=iteration_ub)

			# Get 95 percent cis
			#alt_est_nm_lb, alt_est_nm_ub = get_95_percent_ci_based_on_sampled_distribution(alt_est_nm_h2)
			#alt_est_med_lb, alt_est_med_ub = get_95_percent_ci_based_on_sampled_distribution(alt_est_med_h2)
			#total_est_lb, total_est_ub = get_95_percent_ci_based_on_sampled_distribution(est_total_h2)

			#print(str(sim_total_h2) + '\t' + '[' + str(total_est_lb) + ', ' + str(total_est_ub) + ']')
			#print(str(sim_nm_h2) + '\t' + '[' + str(alt_est_nm_lb) + ', ' + str(alt_est_nm_ub) + ']')
			#print(str(sim_med_h2) + '\t' + '[' + str(alt_est_med_lb) + ', ' + str(alt_est_med_ub) + ']')


			# print to output
			t.write(str(sim_iter) + '\t' + str(eqtl_sample_size) + '\t' + 'marginal_pca_sumstat_gibbs_0.0' + '\t' + str(sim_total_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t')
			t.write(str(np.mean(est_nm_h2)) + '\t' + str(np.mean(est_med_h2)) + '\t' + str(np.mean(est_total_h2)) + '\t' + '\t'.join(np.mean(est_med_h2_per_tissue,axis=0).astype(str)) + '\n')

			'''
			###########################################
			# First do method of marginal PCA gibbs sampler
			marginal_gibbs_sampler_file = trait_med_h2_inference_dir + 'simulation_' + str(sim_iter) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_resid_var_True_cc_1e-06.txt'
			
			# Extract distributions
			est_nm_h2, alt_est_nm_h2, est_med_h2, alt_est_med_h2, est_total_h2 = extract_sampled_params_from_marginal_gibbs_sampler_file(marginal_gibbs_sampler_file, burn_in_thresh=burn_in_thresh)

			# Get 95 percent cis
			alt_est_nm_lb, alt_est_nm_ub = get_95_percent_ci_based_on_sampled_distribution(alt_est_nm_h2)
			alt_est_med_lb, alt_est_med_ub = get_95_percent_ci_based_on_sampled_distribution(alt_est_med_h2)
			total_est_lb, total_est_ub = get_95_percent_ci_based_on_sampled_distribution(est_total_h2)

			#print('')
			#print(str(sim_total_h2) + '\t' + '[' + str(total_est_lb) + ', ' + str(total_est_ub) + ']')
			#print(str(sim_nm_h2) + '\t' + '[' + str(alt_est_nm_lb) + ', ' + str(alt_est_nm_ub) + ']')
			#print(str(sim_med_h2) + '\t' + '[' + str(alt_est_med_lb) + ', ' + str(alt_est_med_ub) + ']')


			# print to output
			t.write(str(sim_iter) + '\t' + str(eqtl_sample_size) + '\t' + 'marginal_pca_sumstat_gibbs_1e-06' + '\t' + str(sim_total_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t')
			t.write(str(np.mean(est_nm_h2)) + '\t' + str(np.mean(alt_est_nm_h2)) + '\t' + str(np.mean(est_med_h2)) + '\t' + str(np.mean(alt_est_med_h2)) + '\t' + str(np.mean(est_total_h2)) + '\n')
			'''
			'''
			###########################################
			# First do method of marginal PCA gibbs sampler
			marginal_gibbs_sampler_file = trait_med_h2_inference_dir + 'simulation_' + str(sim_iter) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_resid_var_True_cc_0.001.txt'
			
			# Extract distributions
			est_nm_h2, alt_est_nm_h2, est_med_h2, alt_est_med_h2, est_total_h2 = extract_sampled_params_from_marginal_gibbs_sampler_file(marginal_gibbs_sampler_file, burn_in_thresh=15000)

			# Get 95 percent cis
			alt_est_nm_lb, alt_est_nm_ub = get_95_percent_ci_based_on_sampled_distribution(alt_est_nm_h2)
			alt_est_med_lb, alt_est_med_ub = get_95_percent_ci_based_on_sampled_distribution(alt_est_med_h2)
			total_est_lb, total_est_ub = get_95_percent_ci_based_on_sampled_distribution(est_total_h2)

			#print('')
			#print(str(sim_total_h2) + '\t' + '[' + str(total_est_lb) + ', ' + str(total_est_ub) + ']')
			#print(str(sim_nm_h2) + '\t' + '[' + str(alt_est_nm_lb) + ', ' + str(alt_est_nm_ub) + ']')
			#print(str(sim_med_h2) + '\t' + '[' + str(alt_est_med_lb) + ', ' + str(alt_est_med_ub) + ']')


			# print to output
			t.write(str(sim_iter) + '\t' + str(eqtl_sample_size) + '\t' + 'marginal_pca_sumstat_gibbs_0.001' + '\t' + str(sim_total_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t')
			t.write(str(np.mean(est_nm_h2)) + '\t' + str(np.mean(alt_est_nm_h2)) + '\t' + str(np.mean(est_med_h2)) + '\t' + str(np.mean(alt_est_med_h2)) + '\t' + str(np.mean(est_total_h2)) + '\n')
			'''
			'''
			###########################################
			# Second do method of rss gibbs sampler
			rss_gibbs_sampler_file = trait_med_h2_inference_dir + 'simulation_' + str(sim_iter) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_eqtl_' + str(eqtl_sample_size) + '_small_resid_var_variable_rss_tmp_res_big_windows_no_pca.txt'
			
			# Extract distributions
			est_nm_h2, alt_est_nm_h2, est_med_h2, alt_est_med_h2, est_total_h2 = extract_sampled_params_from_rss_gibbs_sampler_file(rss_gibbs_sampler_file, burn_in_thresh=15000)

			t.write(str(sim_iter) + '\t' + str(eqtl_sample_size) + '\t' + 'rss_gibbs' + '\t' + str(sim_total_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t')
			t.write(str(np.mean(est_nm_h2)) + '\t' + str(np.mean(alt_est_nm_h2)) + '\t' + str(np.mean(est_med_h2)) + '\t' + str(np.mean(alt_est_med_h2)) + '\t' + str(np.mean(est_total_h2)) + '\n')
			'''

	t.close()

	return

def print_temp_line(method, eqtl_sample_size, t, sim_total_h2, est_total_h2, heritability_type):
	sim_mean = np.mean(sim_total_h2)
	est_mean = np.mean(est_total_h2)
	est_mean_lb = est_mean - 1.96*(np.std(est_total_h2)/np.sqrt(len(est_total_h2)))
	est_mean_ub = est_mean + 1.96*(np.std(est_total_h2)/np.sqrt(len(est_total_h2)))


	t.write(method + '\t' + str(eqtl_sample_size) + '\t' + heritability_type + '\t' + str(sim_mean) + '\t' + str(est_mean) + '\t' + str(est_mean_lb) + '\t' + str(est_mean_ub) + '\n')

	return t


def average_results_across_simulations_single_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names):
	t = open(avg_results_summary_file,'w')
	t.write('method\teqtl_sample_size\theritability_type\tsim_h2\test_h2\test_h2_lb\test_h2_ub\n')
	for ii,method in enumerate(methods):
		for eqtl_sample_size in eqtl_sample_sizes:
			sim_total_h2 = []
			sim_med_h2 = []
			sim_nm_h2 = []
			est_nm_h2 = []
			est_med_h2 = []
			est_total_h2 = []
			for sim_num in sim_nums:
				filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod4.txt'
	
				if os.path.isfile(filer) == False:
					continue

				aa = np.loadtxt(filer,dtype=str,delimiter='\t')

				index = np.where(aa[:,0] == method)[0][0]
				vec = aa[index,:]

				sim_total_h2.append(float(vec[2]))
				sim_med_h2.append(float(vec[3]))
				sim_nm_h2.append(float(vec[4]))
				est_med_h2.append(float(vec[5]))
				est_nm_h2.append(float(vec[7]))
				est_total_h2.append(float(vec[5]) + float(vec[7]))
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_total_h2, est_total_h2, 'total_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_nm_h2, est_nm_h2, 'nm_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_med_h2, est_med_h2, 'med_h2')
	t.close()
	return


def average_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, permuted_eqtls=False):
	t = open(avg_results_summary_file,'w')
	t.write('method\teqtl_sample_size\theritability_type\tsim_h2\test_h2\test_h2_lb\test_h2_ub\n')
	for ii,method in enumerate(methods):
		for eqtl_sample_size in eqtl_sample_sizes:
			sim_total_h2 = []
			sim_med_h2 = []
			sim_nm_h2 = []
			est_nm_h2 = []
			est_med_h2 = []
			est_total_h2 = []
			est_per_tissue_h2 =[]
			est_per_category_h2 = []
			est_eqtl_h2 = []
			for sim_num in sim_nums:
				filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_2_cis_window_500000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + run_string + '_' + str(eqtl_sample_size) + '_calibrated_mesc_results.txt'
				#filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod13.txt'
				if permuted_eqtls:
					filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod13_permuted_eqtls.txt'

				if sim_num in invalid_sims:
					continue
				if os.path.isfile(filer) == False:
					print(sim_num)
					print('miss')
					continue

				#aa = np.loadtxt(filer,dtype=str,delimiter='\t')
				f = open(filer)
				arr = []
				for line in f:
					line = line.rstrip()
					data = line.split()
					if data[0].endswith(str(eqtl_sample_size)):
						new_data = []
						new_data.append(data[0].split(str(eqtl_sample_size))[0])
						new_data.append(eqtl_sample_size)
						for ele in data[1:]:
							new_data.append(ele)
						arr.append(np.asarray(new_data))
					else:
						arr.append(np.asarray(data))
				f.close()
				aa = np.asarray(arr)

				if aa.shape[0] != 7:
					print(sim_num)
					print(aa.shape)
					continue



				index = np.where(aa[:,0] == method)[0][0]
				vec = aa[index,:]


				per_tissue_est = np.asarray(vec[4].split(',')).astype(float)
				per_category_est = np.asarray(vec[5].split(',')).astype(float)


				sim_total_h2.append(float(.3))
				sim_med_h2.append(float(.03))
				sim_nm_h2.append(float(0.0))
				est_med_h2.append(float(vec[1]))
				est_nm_h2.append(float(vec[2]))
				est_total_h2.append(float(vec[1]) + float(vec[2]))
				est_per_tissue_h2.append(per_tissue_est)
				est_per_category_h2.append(per_category_est)
				est_eqtl_h2.append(float(vec[3]))



			est_per_tissue_h2 = np.asarray(est_per_tissue_h2)
			est_per_category_h2 = np.asarray(est_per_category_h2)

			if len(sim_total_h2) == 0:
				continue

			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_total_h2, est_total_h2, 'total_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_nm_h2, est_nm_h2, 'nm_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_med_h2, est_med_h2, 'med_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_med_h2, est_per_tissue_h2[:,0], 'causal_tissue_med_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, 0.0, np.sum(est_per_tissue_h2[:,1:],axis=1), 'non_causal_tissue_med_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, (.025 + .05 + .1)/3.0, est_eqtl_h2, 'eqtl_h2')
		

			for tissue_num in range(est_per_tissue_h2.shape[1]):
				est_single_tissue = est_per_tissue_h2[:, tissue_num]
				if tissue_num == 0:
					simmer = np.copy(sim_med_h2)
				else:
					simmer = np.copy(sim_med_h2)*0.0
				t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, simmer, est_single_tissue, 'med_h2_tissue' + str(tissue_num))

	t.close()
	print(avg_results_summary_file)
	return	


def average_results_across_simulations_1_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, permuted_eqtls=False):
	t = open(avg_results_summary_file,'w')
	t.write('method\teqtl_sample_size\theritability_type\tsim_h2\test_h2\test_h2_lb\test_h2_ub\n')
	for ii,method in enumerate(methods):
		for eqtl_sample_size in eqtl_sample_sizes:
			sim_total_h2 = []
			sim_med_h2 = []
			sim_nm_h2 = []
			est_nm_h2 = []
			est_med_h2 = []
			est_total_h2 = []
			est_per_tissue_h2 =[]
			for sim_num in sim_nums:
				filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + str(eqtl_sample_size) + '_joint_ldsc_single_tissue_multimethod12.txt'
				if permuted_eqtls:
					filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + str(eqtl_sample_size) + '_joint_ldsc_single_tissue_multimethod12_permuted_eqtls.txt'

				if os.path.isfile(filer) == False:
					continue

				#aa = np.loadtxt(filer,dtype=str,delimiter='\t')
				f = open(filer)
				arr = []
				for line in f:
					line = line.rstrip()
					data = line.split('\t')
					if data[0].endswith(str(eqtl_sample_size)):
						new_data = []
						new_data.append(data[0].split(str(eqtl_sample_size))[0])
						new_data.append(eqtl_sample_size)
						for ele in data[1:]:
							new_data.append(ele)
						arr.append(np.asarray(new_data))
					else:
						arr.append(np.asarray(data))
				f.close()
				aa = np.asarray(arr)


				if aa.shape[0] != 3:
					continue


				index = np.where(aa[:,0] == method)[0][0]
				vec = aa[index,:]

				per_tissue_est = np.asarray(vec[6].split(',')).astype(float)



				sim_total_h2.append(float(vec[2]))
				sim_med_h2.append(float(vec[3]))
				sim_nm_h2.append(float(vec[4]))
				est_med_h2.append(float(vec[5]))
				est_nm_h2.append(float(vec[7]))
				est_total_h2.append(float(vec[5]) + float(vec[7]))
				est_per_tissue_h2.append(per_tissue_est)

			est_per_tissue_h2 = np.asarray(est_per_tissue_h2)
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_total_h2, est_total_h2, 'total_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_nm_h2, est_nm_h2, 'nm_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_med_h2, est_med_h2, 'med_h2')
			

	t.close()
	return	

def average_mesc_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, permuted_eqtls=False):
	t = open(avg_results_summary_file,'w')
	t.write('method\teqtl_sample_size\theritability_type\tsim_h2\test_h2\test_h2_lb\test_h2_ub\n')
	for ii,method in enumerate(methods):
		for eqtl_sample_size in eqtl_sample_sizes:
			sim_total_h2 = []
			sim_med_h2 = []
			sim_nm_h2 = []
			est_nm_h2 = []
			est_med_h2 = []
			est_total_h2 = []
			est_per_tissue_h2 =[]
			for sim_num in sim_nums:
				filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_eqtl_SS_' + str(eqtl_sample_size) + '_tglr_estimates.txt'
				if permuted_eqtls:
					filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_eqtl_SS_' + str(eqtl_sample_size) + '_tglr_estimates_permuted_eqtls.txt'

				if os.path.isfile(filer) == False:
					continue

				#aa = np.loadtxt(filer,dtype=str,delimiter='\t')
				f = open(filer)
				arr = []
				for line in f:
					line = line.rstrip()
					data = line.split('\t')
					if data[0].endswith(str(eqtl_sample_size)):
						new_data = []
						new_data.append(data[0].split(str(eqtl_sample_size))[0])
						new_data.append(eqtl_sample_size)
						for ele in data[1:]:
							new_data.append(ele)
						arr.append(np.asarray(new_data))
					else:
						arr.append(np.asarray(data))
				f.close()
				aa = np.asarray(arr)

				if aa.shape[0] != 4:
					continue


				index = np.where(aa[:,0] == method)[0][0]
				vec = aa[index,:]


				per_tissue_est = np.asarray(vec[7].split(','))
				per_tissue_est[per_tissue_est == 'nan'] = '0'
				per_tissue_est = per_tissue_est.astype(float)
				vec[6] = str(np.sum(per_tissue_est))


				sim_total_h2.append(float(vec[1]))
				sim_med_h2.append(float(vec[3]))
				sim_nm_h2.append(float(vec[2]))
				est_med_h2.append(float(vec[6]))
				est_nm_h2.append(float(vec[5]))
				est_total_h2.append(float(vec[5]) + float(vec[6]))
				est_per_tissue_h2.append(per_tissue_est)

			est_per_tissue_h2 = np.asarray(est_per_tissue_h2)
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_total_h2, est_total_h2, 'total_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_nm_h2, est_nm_h2, 'nm_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_med_h2, est_med_h2, 'med_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_med_h2, est_per_tissue_h2[:,0], 'causal_tissue_med_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, 0.0, np.sum(est_per_tissue_h2[:,1:],axis=1), 'non_causal_tissue_med_h2')
			
			for tissue_num in range(est_per_tissue_h2.shape[1]):
				est_single_tissue = est_per_tissue_h2[:, tissue_num]
				if tissue_num == 0:
					simmer = np.copy(sim_med_h2)
				else:
					simmer = np.copy(sim_med_h2)*0.0
				t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, simmer, est_single_tissue, 'med_h2_tissue' + str(tissue_num))

	t.close()
	return	




def concatenate_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names):
	t = open(avg_results_summary_file,'w')
	t.write('eqtl_sample_size\tsim_iter\test_tissue0_joint\test_tissue0_two_step\tcorrelation_joint\tcorrelation_two_step\tcondition_number_joint\tcondition_number_two_step\n')
	for eqtl_sample_size in eqtl_sample_sizes:

		for sim_num in sim_nums:
			filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod12.txt'
			X_int_file = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod8_X_init.npy'
			if os.path.isfile(filer) == False or os.path.isfile(X_int_file) == False:
				continue

			#aa = np.loadtxt(filer,dtype=str,delimiter='\t')
			f = open(filer)
			arr = []
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if data[0].endswith(str(eqtl_sample_size)):
					new_data = []
					new_data.append(data[0].split(str(eqtl_sample_size))[0])
					new_data.append(eqtl_sample_size)
					for ele in data[1:]:
						new_data.append(ele)
					arr.append(np.asarray(new_data))
				else:
					arr.append(np.asarray(data))
			f.close()
			aa = np.asarray(arr)

			if aa.shape[0] != 4 :
				continue

			vec = aa[1,:]
			two_step_per_tissue_est = np.asarray(vec[6].split(',')).astype(float)
			vec = aa[2,:]
			joint_per_tissue_est = np.asarray(vec[6].split(',')).astype(float)


			X_file = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod12_X.npy'
			X_int_file = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + str(eqtl_sample_size) + '_joint_ldsc_multimethod8_X_init.npy'
			X = np.load(X_file)
			X_init = np.load(X_int_file)

			X_corr = np.corrcoef(np.transpose(X[:,1:]))
			X_init_corr = np.corrcoef(np.transpose(X_init[:,1:]))
			max_X_corr = np.max(X_corr[0,1:])
			max_X_init_corr = np.max(X_init_corr[0,1:])

			X_cond = np.linalg.cond(X[:,1:])
			X_init_cond = np.linalg.cond(X_init[:,1:])


			t.write(str(eqtl_sample_size) + '\t' + str(sim_num) + '\t' + str(joint_per_tissue_est[0]) + '\t' + str(two_step_per_tissue_est[0]) + '\t' + str(max_X_corr) + '\t' + str(max_X_init_corr) + '\t' + str(X_cond) + '\t' + str(X_init_cond) + '\n')
	

	t.close()
	return	

def average_results_across_simulations(results_summary_file, avg_results_summary_file, eqtl_sample_sizes, methods):
	t = open(avg_results_summary_file,'w')
	t.write('method\teqtl_sample_size\theritability_type\tsim_h2\test_h2\test_h2_lb\test_h2_ub\n')
	for method in methods:
		for eqtl_sample_size in eqtl_sample_sizes:
			f = open(results_summary_file)
			head_count=0
			sim_total_h2 = []
			sim_med_h2 = []
			sim_nm_h2 = []
			est_nm_h2 = []
			est_alt_nm_h2 = []
			est_med_h2 = []
			est_alt_med_h2 = []
			est_total_h2 = []
			est_per_tissue_h2 =[]
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count ==0:
					head_count = head_count + 1
					continue
				if data[1] != str(eqtl_sample_size):
					continue
				if data[2] != method:
					continue

				sim_total_h2.append(float(data[3]))
				sim_nm_h2.append(float(data[4]))
				sim_med_h2.append(float(data[5]))
				est_nm_h2.append(float(data[6]))
				est_med_h2.append(float(data[7]))
				est_total_h2.append(float(data[8]))
				per_tissue_h2 = np.asarray(data[9:]).astype(float)
				est_per_tissue_h2.append(per_tissue_h2)
			f.close()

			est_per_tissue_h2 = np.asarray(est_per_tissue_h2)
			t = print_temp_line(method, eqtl_sample_size, t, sim_total_h2, est_total_h2, 'total_h2')
			t = print_temp_line(method, eqtl_sample_size, t, sim_nm_h2, est_nm_h2, 'nm_h2')
			t = print_temp_line(method, eqtl_sample_size, t, sim_med_h2, est_med_h2, 'med_h2')

			for tissue_num in range(est_per_tissue_h2.shape[1]):
				est_single_tissue = est_per_tissue_h2[:, tissue_num]
				if tissue_num == 0:
					simmer = np.copy(sim_med_h2)
				else:
					simmer = np.copy(sim_med_h2)*0.0
				t = print_temp_line(method, eqtl_sample_size, t, simmer, est_single_tissue, 'med_h2_tissue' + str(tissue_num))

	t.close()
	print(avg_results_summary_file)
	return

def average_results_across_tglr_simulations(avg_results_summary_file, eqtl_sample_sizes,trait_med_h2_inference_dir, sim_nums):
	t = open(avg_results_summary_file,'w')
	t.write('method\teqtl_sample_size\theritability_type\tsim_h2\test_h2\test_h2_lb\test_h2_ub\n')

	methods = ['tglr_unscaled', 'tglr_scaled']
	for method in methods:
		for eqtl_sample_size in eqtl_sample_sizes:
			sim_total_h2 = []
			sim_med_h2 = []
			sim_nm_h2 = []
			est_nm_h2 = []
			est_med_h2 = []
			est_total_h2 = []
			for sim_iter in sim_nums:
				filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_iter) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_eqtl_SS_' + str(eqtl_sample_size) + '_tglr_estimates.txt'
				data = np.loadtxt(filer,dtype=str,delimiter='\t')

				if method == 'tglr_unscaled':
					sim_total_h2.append(data[1,1])
					sim_med_h2.append(data[1,3])
					sim_nm_h2.append(data[1,2])
					est_total_h2.append(data[1,4])
					est_nm_h2.append(data[1,5])
					est_med_h2.append(data[1,6])
				else:
					sim_total_h2.append(data[2,1])
					sim_med_h2.append(data[2,3])
					sim_nm_h2.append(data[2,2])
					est_total_h2.append(data[2,4])
					est_nm_h2.append(data[2,5])
					est_med_h2.append(data[2,6])					

			t = print_temp_line(method, eqtl_sample_size, t, np.asarray(sim_total_h2).astype(float), np.asarray(est_total_h2).astype(float), 'total_h2')
			t = print_temp_line(method, eqtl_sample_size, t, np.asarray(sim_nm_h2).astype(float), np.asarray(est_nm_h2).astype(float), 'nm_h2')
			t = print_temp_line(method, eqtl_sample_size, t, np.asarray(sim_med_h2).astype(float), np.asarray(est_med_h2).astype(float), 'med_h2')

	t.close()
	print(avg_results_summary_file)
	return



#####################
# Command line args
#####################
trait_med_h2_inference_dir = sys.argv[1]
visualize_trait_med_h2_dir = sys.argv[2]


# Simulation params
sim_nums = np.arange(1,201)
eqtl_sample_sizes = np.asarray(['100','200','300','1000'])
eqtl_sample_sizes = np.asarray(['100', '200'])
eqtl_sample_sizes = np.asarray(['100','200','300','1000', '100-1000'])



##################################
# Multi-tissue simulation
##################################


#####################
# Joint LDSC analysis
#sim_nums = np.arange(1,501)
invalid_sims = {}
invalid_sims[76] = 1 # This should be only bad one
invalid_sims[162] = 1 # These will get re-run



eqtl_snp_representation = 'pca_95'
eqtl_snp_representation = 'pca_90'
eqtl_snp_representation = 'bins_20'
beta_squared_thresh = '100.0'

non_med_anno = 'full_anno'
non_med_anno = 'genotype_intercept'

simulated_gt_architecture = 'linear'
inference_gt_architecture = 'linear'
gene_ld_score_type = 'ldsc_style_pred'
gene_ld_score_type = 'squared_marginal_sumstats'


run_string = eqtl_snp_representation + '_' + non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + gene_ld_score_type + '_' + beta_squared_thresh

avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_' + eqtl_snp_representation + '_' + non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + gene_ld_score_type + '_' + beta_squared_thresh + '_sim_results_calibrated_ldsc_summary_averaged.txt'
clean_method_names = ['two_step_ldsc', 'calibrated_two_step_ldsc']
methods = ['two_step_joint_ldsc', 'calibrated_two_step_joint_ldsc']
average_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string)

'''
# concatenate results
concat_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_sim_results_joint_ldsc_binned_summary_concatenated.txt'
methods = ['eqtl_5_binned_no_intercept_two_step', 'eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance', 'eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance_cis_window', 'eqtl_pced_no_intercept_two_step', 'eqtl_pced_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance', 'eqtl_pced_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance_cis_window']
clean_method_names = ['two_step_ldsc', 'joint_ldsc', 'joint_ldsc_cis_var', 'two_step_pca_ldsc', 'joint_pca_ldsc', 'joint_pca_ldsc_cis_var']
concatenate_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, concat_results_summary_file, methods, clean_method_names)
'''

'''
#####################
# Permuted Joint LDSC analysis
avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_permuted_eqtls_sim_results_joint_ldsc_binned_summary_averaged.txt'
methods = ['eqtl_5_binned_no_intercept_two_step', 'eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance']
clean_method_names = ['two_step_ldsc', 'joint_ldsc']
average_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names,permuted_eqtls=True)

#####################
# MESC analysis
sim_nums = np.arange(1,501)

avg_results_summary_file = visualize_trait_med_h2_dir + 'med_h2_5_causal_tissue_sim_results_mesc_summary_averaged.txt'
methods = ['tglr']
clean_method_names = ['mesc']
average_mesc_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names)


avg_results_summary_file = visualize_trait_med_h2_dir + 'med_h2_5_causal_tissue_permuted_sim_results_mesc_summary_averaged.txt'
methods = ['tglr']
clean_method_names = ['mesc']
average_mesc_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names,permuted_eqtls=True)
'''


#####################
# Only include single causal tissue
'''
avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_1_causal_tissue_sim_results_joint_ldsc_binned_summary_averaged.txt'
methods = ['eqtl_5_binned_no_intercept_two_step', 'eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance']
clean_method_names = ['two_step_ldsc', 'joint_ldsc']
average_results_across_simulations_1_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names)

#####################
# Permuted Joint LDSC analysis
avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_1_causal_tissue_permuted_eqtls_sim_results_joint_ldsc_binned_summary_averaged.txt'
methods = ['eqtl_5_binned_no_intercept_two_step', 'eqtl_5_binned_no_intercept_bayesian_gibbs_resid_var_multivariate_per_data_set_variance']
clean_method_names = ['two_step_ldsc', 'joint_ldsc']
average_results_across_simulations_1_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names,permuted_eqtls=True)
print(avg_results_summary_file)
'''


import numpy as np
import os
import sys
import pdb
from scipy.stats import norm




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

def print_temp_line_accounting_for_variance(method, eqtl_sample_size, t, sim_total_h2, est_total_h2, est_total_h2_se, heritability_type):

	se   = np.array(est_total_h2_se)       # their standard errors
	ww   = 1.0 / se**2     

	sim_mean = np.mean(sim_total_h2)
	

	est_mean = np.sum(est_total_h2*ww)/np.sum(ww)
	se_mu = np.sqrt(1.0/np.sum(ww))
	est_mean_lb = est_mean - 1.96*se_mu
	est_mean_ub = est_mean + 1.96*se_mu


	t.write(method + '\t' + str(eqtl_sample_size) + '\t' + heritability_type + '\t' + str(sim_mean) + '\t' + str(est_mean) + '\t' + str(est_mean_lb) + '\t' + str(est_mean_ub) + '\n')

	return t


def ci(mean, se, conf=0.95):
	alpha = 1 - conf
	z = norm.ppf(1 - alpha/2)
	return mean - z*se, mean + z*se

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


def compute_coverage(sim_values, est_values, est_value_ses, coverage):
	est_lb, est_ub = ci(est_values, est_value_ses, conf=coverage)


	n_covered = np.sum((sim_values >= est_lb) & (sim_values <= est_ub))
	n_tot = len(sim_values)
	coverage = n_covered/n_tot

	coverage_se = np.sqrt(((coverage)*(1.0-coverage))/n_tot)

	coverage_lb = coverage - (1.96*coverage_se)
	coverage_ub = coverage + (1.96*coverage_se)

	return coverage, coverage_lb, coverage_ub


def compute_power(sim_values, est_values, est_value_ses, coverage):
	est_lb, est_ub = ci(est_values, est_value_ses, conf=coverage)

	n_covered = np.sum(est_lb >= 0.0)
	n_tot = len(sim_values)
	coverage = n_covered/n_tot

	coverage_se = np.sqrt(((coverage)*(1.0-coverage))/n_tot)

	coverage_lb = coverage - (1.96*coverage_se)
	coverage_ub = coverage + (1.96*coverage_se)

	return coverage, coverage_lb, coverage_ub


def compute_t1e(est_values, est_value_ses, pvalue_threshold):
	est_lb, est_ub = ci(est_values, est_value_ses, conf=1.0-pvalue_threshold)


	n_covered = np.sum((est_lb <= 0) & (est_ub >= 0))
	n_tot = len(est_values)
	t1e = (n_tot-n_covered)/n_tot

	t1e_se = np.sqrt(((t1e)*(1.0-t1e))/n_tot)

	t1e_lb = t1e - (1.96*t1e_se)
	t1e_ub = t1e + (1.96*t1e_se)


	return t1e, t1e_lb, t1e_ub





def fstat_summary(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, permuted_eqtls=False, variance_weighting=False):

	t = open(avg_results_summary_file,'w')
	t.write('method_name\teQTL_ss\tfstat\n')


	for ii,method in enumerate(methods):
		for eqtl_sample_size in eqtl_sample_sizes:
			for sim_num in sim_nums:
				filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_2_cis_window_500000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + run_string + '_' + str(eqtl_sample_size) + '_' + weighting + '_step_1_f_stats.txt'

				if sim_num in invalid_sims:
					continue
				if os.path.isfile(filer) == False:
					print(sim_num)
					print('miss')
					continue
				f = open(filer)
				head_count = 0
				for line in f:
					line = line.rstrip()
					data = line.split()
					if head_count == 0:
						head_count = head_count + 1
						continue
					if data[1] != '1' or data[2] != 'bootstrapped':
						continue
					t.write(data[0] + '\t' + str(eqtl_sample_size) + '\t' + data[3] + '\n')
				f.close()
	t.close()
	print(avg_results_summary_file)
	return


def power_summary(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, permuted_eqtls=False, variance_weighting=False):
	pvalue_thresholds  = [.01, .05, .1, .2]

	t = open(avg_results_summary_file,'w')
	t.write('method_name\teQTL_ss\tgenetic_element\tpvalue_threshold\tpower\tpower_lb\tpower_ub\n')


	for ii,method in enumerate(methods):
		for eqtl_sample_size in eqtl_sample_sizes:
			sim_total_h2 = []
			sim_med_h2 = []
			sim_nm_h2 = []
			est_nm_h2 = []
			est_nm_h2_se = []
			est_med_h2 = []
			est_med_h2_se = []
			est_total_h2 = []
			est_total_h2_se = []
			est_per_tissue_h2 =[]
			est_per_tissue_h2_se = []
			est_per_category_h2 = []
			est_per_category_h2_se = []
			est_eqtl_h2 = []
			for sim_num in sim_nums:
				filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_2_cis_window_500000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + run_string + '_' + str(eqtl_sample_size) + '_' + weighting + '_' + method + '_jk_results.txt'

				if sim_num in invalid_sims:
					continue
				if os.path.isfile(filer) == False:
					print(sim_num)
					print('miss')
					continue
				f = open(filer)
				arr_est = []
				arr_names = []
				arr_est_se = []
				per_tissue_est = []
				per_tissue_est_se = []
				per_category_est = []
				per_category_est_se = []
				head_count = 0
				for line in f:
					line = line.rstrip()
					data = line.split()
					if head_count == 0:
						head_count = head_count + 1
						continue
					arr_est.append(float(data[1]))
					arr_names.append(data[0])
					arr_est_se.append(float(data[3]))
					if data[0].startswith('category_med'):
						per_category_est.append(float(data[1]))
						per_category_est_se.append(float(data[3]))
					if data[0].startswith('dataset_med'):
						per_tissue_est.append(float(data[1]))
						per_tissue_est_se.append(float(data[3]))
				f.close()
				arr_est = np.asarray(arr_est)
				arr_names = np.asarray(arr_names)
				arr_est_se = np.asarray(arr_est_se)
				per_tissue_est = np.asarray(per_tissue_est)
				per_tissue_est_se = np.asarray(per_tissue_est_se)
				per_category_est = np.asarray(per_category_est)
				per_category_est_se = np.asarray(per_category_est_se)


				tmp_total_h2_sim, tmp_total_nm_h2_sim, tmp_total_med_h2_sim = sim_heritabilities[sim_num]

				sim_total_h2.append(tmp_total_h2_sim)
				sim_med_h2.append(tmp_total_med_h2_sim)
				sim_nm_h2.append(tmp_total_nm_h2_sim)
				est_med_h2.append(float(arr_est[2]))
				est_med_h2_se.append(float(arr_est_se[2]))
				est_nm_h2.append(float(arr_est[1]))
				est_nm_h2_se.append(float(arr_est_se[1]))
				est_total_h2.append(float(arr_est[0]))
				est_total_h2_se.append(float(arr_est_se[0]))
				est_per_tissue_h2.append(per_tissue_est)
				est_per_tissue_h2_se.append(per_tissue_est_se)
				est_per_category_h2.append(per_category_est)
				est_per_category_h2_se.append(per_category_est_se)
				est_eqtl_h2.append(float(arr_est[-1]))	

			est_per_tissue_h2 = np.asarray(est_per_tissue_h2)
			est_per_category_h2 = np.asarray(est_per_category_h2)
			est_per_tissue_h2_se = np.asarray(est_per_tissue_h2_se)
			est_per_category_h2_se = np.asarray(est_per_category_h2_se)
			sim_med_h2 = np.asarray(sim_med_h2)
			est_med_h2 = np.asarray(est_med_h2)
			est_med_h2_se = np.asarray(est_med_h2_se)
			sim_nm_h2 = np.asarray(sim_nm_h2)
			est_nm_h2 = np.asarray(est_nm_h2)
			est_nm_h2_se = np.asarray(est_nm_h2_se)
			est_total_h2 = np.asarray(est_total_h2)
			est_total_h2_se = np.asarray(est_total_h2_se)


			# Calibration for med h2
			for pvalue_threshold in pvalue_thresholds:

				coverage = 1.0 - pvalue_threshold

				# NM H2
				observed_power, observed_power_lb, observed_power_ub = compute_power(sim_nm_h2, est_nm_h2, est_nm_h2_se, coverage)
				t.write(clean_method_names[ii] + '\t' + str(eqtl_sample_size) + '\t' + 'total_nm_h2' + '\t' + str(pvalue_threshold) + '\t' + str(observed_power) + '\t' + str(observed_power_lb) + '\t' + str(observed_power_ub) + '\n')

				# Total mediated h2
				observed_power, observed_power_lb, observed_power_ub = compute_power(sim_med_h2, est_med_h2, est_med_h2_se, coverage)
				t.write(clean_method_names[ii] + '\t' + str(eqtl_sample_size) + '\t' + 'total_med_h2' + '\t' + str(pvalue_threshold) + '\t' + str(observed_power) + '\t' + str(observed_power_lb) + '\t' + str(observed_power_ub) + '\n')

				# causal tissue mediated h2
				observed_power, observed_power_lb, observed_power_ub = compute_power(sim_med_h2, est_per_tissue_h2[:,0], est_per_tissue_h2_se[:,0], coverage)
				t.write(clean_method_names[ii] + '\t' + str(eqtl_sample_size) + '\t' + 'causal_tissue_med_h2' + '\t' + str(pvalue_threshold) + '\t' + str(observed_power) + '\t' + str(observed_power_lb) + '\t' + str(observed_power_ub) + '\n')




	t.close()

	print(avg_results_summary_file)

	return



def t1e_summary(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, t1e_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, permuted_eqtls=False, variance_weighting=False):
	pvalue_thresholds  = [.01, .05, .1, .2]

	t = open(t1e_results_summary_file,'w')
	t.write('method_name\teQTL_ss\tgenetic_element\tpvalue_threshold\tt1e\tt1e_lb\tt1e_ub\n')


	for ii,method in enumerate(methods):
		for eqtl_sample_size in eqtl_sample_sizes:
			sim_total_h2 = []
			sim_med_h2 = []
			sim_nm_h2 = []
			est_nm_h2 = []
			est_nm_h2_se = []
			est_med_h2 = []
			est_med_h2_se = []
			est_total_h2 = []
			est_total_h2_se = []
			est_per_tissue_h2 =[]
			est_per_tissue_h2_se = []
			est_per_category_h2 = []
			est_per_category_h2_se = []
			est_eqtl_h2 = []
			for sim_num in sim_nums:
				filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_2_cis_window_500000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + run_string + '_' + str(eqtl_sample_size) + '_' + weighting + '_' + method + '_jk_results.txt'

				if sim_num in invalid_sims:
					continue
				if os.path.isfile(filer) == False:
					print(sim_num)
					print('miss')
					continue
				f = open(filer)
				arr_est = []
				arr_names = []
				arr_est_se = []
				per_tissue_est = []
				per_tissue_est_se = []
				per_category_est = []
				per_category_est_se = []
				head_count = 0
				for line in f:
					line = line.rstrip()
					data = line.split()
					if head_count == 0:
						head_count = head_count + 1
						continue
					arr_est.append(float(data[1]))
					arr_names.append(data[0])
					arr_est_se.append(float(data[3]))
					if data[0].startswith('category_med'):
						per_category_est.append(float(data[1]))
						per_category_est_se.append(float(data[3]))
					if data[0].startswith('dataset_med'):
						per_tissue_est.append(float(data[1]))
						per_tissue_est_se.append(float(data[3]))
				f.close()
				arr_est = np.asarray(arr_est)
				arr_names = np.asarray(arr_names)
				arr_est_se = np.asarray(arr_est_se)
				per_tissue_est = np.asarray(per_tissue_est)
				per_tissue_est_se = np.asarray(per_tissue_est_se)
				per_category_est = np.asarray(per_category_est)
				per_category_est_se = np.asarray(per_category_est_se)


				tmp_total_h2_sim, tmp_total_nm_h2_sim, tmp_total_med_h2_sim = sim_heritabilities[sim_num]

				sim_total_h2.append(tmp_total_h2_sim)
				sim_med_h2.append(tmp_total_med_h2_sim)
				sim_nm_h2.append(tmp_total_nm_h2_sim)
				est_med_h2.append(float(arr_est[2]))
				est_med_h2_se.append(float(arr_est_se[2]))
				est_nm_h2.append(float(arr_est[1]))
				est_nm_h2_se.append(float(arr_est_se[1]))
				est_total_h2.append(float(arr_est[0]))
				est_total_h2_se.append(float(arr_est_se[0]))
				est_per_tissue_h2.append(per_tissue_est)
				est_per_tissue_h2_se.append(per_tissue_est_se)
				est_per_category_h2.append(per_category_est)
				est_per_category_h2_se.append(per_category_est_se)
				est_eqtl_h2.append(float(arr_est[-1]))	

			est_per_tissue_h2 = np.asarray(est_per_tissue_h2)
			est_per_category_h2 = np.asarray(est_per_category_h2)
			est_per_tissue_h2_se = np.asarray(est_per_tissue_h2_se)
			est_per_category_h2_se = np.asarray(est_per_category_h2_se)
			sim_med_h2 = np.asarray(sim_med_h2)
			est_med_h2 = np.asarray(est_med_h2)
			est_med_h2_se = np.asarray(est_med_h2_se)
			sim_nm_h2 = np.asarray(sim_nm_h2)
			est_nm_h2 = np.asarray(est_nm_h2)
			est_nm_h2_se = np.asarray(est_nm_h2_se)
			est_total_h2 = np.asarray(est_total_h2)
			est_total_h2_se = np.asarray(est_total_h2_se)


			# Calibration for med h2
			for pvalue_threshold in pvalue_thresholds:

				# NM H2
				observed_t1e, observed_t1e_lb, observed_t1e_ub = compute_t1e(est_per_tissue_h2[:,1], est_per_tissue_h2_se[:,1], pvalue_threshold)
				t.write(clean_method_names[ii] + '\t' + str(eqtl_sample_size) + '\t' + 'non_causal_tissue' + '\t' + str(pvalue_threshold) + '\t' + str(observed_t1e) + '\t' + str(observed_t1e_lb) + '\t' + str(observed_t1e_ub) + '\n')





	t.close()

	print(t1e_results_summary_file)

	return



def confidence_interval_calibration(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, permuted_eqtls=False, variance_weighting=False):
	coverages = [.5, .7, .9, .95]

	t = open(avg_results_summary_file,'w')
	t.write('method_name\teQTL_ss\tgenetic_element\texpected_coverage\tobserved_coverage\tobserved_coverage_lb\tobserved_coverage_ub\n')

	for ii,method in enumerate(methods):
		for eqtl_sample_size in eqtl_sample_sizes:
			sim_total_h2 = []
			sim_med_h2 = []
			sim_nm_h2 = []
			est_nm_h2 = []
			est_nm_h2_se = []
			est_med_h2 = []
			est_med_h2_se = []
			est_total_h2 = []
			est_total_h2_se = []
			est_per_tissue_h2 =[]
			est_per_tissue_h2_se = []
			est_per_category_h2 = []
			est_per_category_h2_se = []
			est_eqtl_h2 = []
			for sim_num in sim_nums:
				filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_2_cis_window_500000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + run_string + '_' + str(eqtl_sample_size) + '_' + weighting + '_' + method + '_jk_results.txt'

				if sim_num in invalid_sims:
					continue
				if os.path.isfile(filer) == False:
					print(sim_num)
					print('miss')
					continue

				f = open(filer)
				arr_est = []
				arr_names = []
				arr_est_se = []
				per_tissue_est = []
				per_tissue_est_se = []
				per_category_est = []
				per_category_est_se = []
				head_count = 0
				for line in f:
					line = line.rstrip()
					data = line.split()
					if head_count == 0:
						head_count = head_count + 1
						continue
					arr_est.append(float(data[1]))
					arr_names.append(data[0])
					arr_est_se.append(float(data[3]))
					if data[0].startswith('category_med'):
						per_category_est.append(float(data[1]))
						per_category_est_se.append(float(data[3]))
					if data[0].startswith('dataset_med'):
						per_tissue_est.append(float(data[1]))
						per_tissue_est_se.append(float(data[3]))
				f.close()
				arr_est = np.asarray(arr_est)
				arr_names = np.asarray(arr_names)
				arr_est_se = np.asarray(arr_est_se)
				per_tissue_est = np.asarray(per_tissue_est)
				per_tissue_est_se = np.asarray(per_tissue_est_se)
				per_category_est = np.asarray(per_category_est)
				per_category_est_se = np.asarray(per_category_est_se)


				tmp_total_h2_sim, tmp_total_nm_h2_sim, tmp_total_med_h2_sim = sim_heritabilities[sim_num]

				sim_total_h2.append(tmp_total_h2_sim)
				sim_med_h2.append(tmp_total_med_h2_sim)
				sim_nm_h2.append(tmp_total_nm_h2_sim)
				est_med_h2.append(float(arr_est[2]))
				est_med_h2_se.append(float(arr_est_se[2]))
				est_nm_h2.append(float(arr_est[1]))
				est_nm_h2_se.append(float(arr_est_se[1]))
				est_total_h2.append(float(arr_est[0]))
				est_total_h2_se.append(float(arr_est_se[0]))
				est_per_tissue_h2.append(per_tissue_est)
				est_per_tissue_h2_se.append(per_tissue_est_se)
				est_per_category_h2.append(per_category_est)
				est_per_category_h2_se.append(per_category_est_se)
				est_eqtl_h2.append(float(arr_est[-1]))	

			est_per_tissue_h2 = np.asarray(est_per_tissue_h2)
			est_per_category_h2 = np.asarray(est_per_category_h2)
			est_per_tissue_h2_se = np.asarray(est_per_tissue_h2_se)
			est_per_category_h2_se = np.asarray(est_per_category_h2_se)
			sim_med_h2 = np.asarray(sim_med_h2)
			est_med_h2 = np.asarray(est_med_h2)
			est_med_h2_se = np.asarray(est_med_h2_se)
			sim_nm_h2 = np.asarray(sim_nm_h2)
			est_nm_h2 = np.asarray(est_nm_h2)
			est_nm_h2_se = np.asarray(est_nm_h2_se)
			est_total_h2 = np.asarray(est_total_h2)
			est_total_h2_se = np.asarray(est_total_h2_se)
			sim_total_h2 = np.asarray(sim_total_h2)


			# Calibration for med h2
			for coverage in coverages:

				observed_coverage, observed_coverage_lb, observed_coverage_ub = compute_coverage(sim_total_h2, est_total_h2, est_total_h2_se, coverage)
				t.write(clean_method_names[ii] + '\t' + str(eqtl_sample_size) + '\t' + 'total_h2' + '\t' + str(coverage) + '\t' + str(observed_coverage) + '\t' + str(observed_coverage_lb) + '\t' + str(observed_coverage_ub) + '\n')

				observed_coverage, observed_coverage_lb, observed_coverage_ub = compute_coverage(sim_nm_h2, est_nm_h2, est_nm_h2_se, coverage)
				t.write(clean_method_names[ii] + '\t' + str(eqtl_sample_size) + '\t' + 'total_nm_h2' + '\t' + str(coverage) + '\t' + str(observed_coverage) + '\t' + str(observed_coverage_lb) + '\t' + str(observed_coverage_ub) + '\n')

				# Total mediated h2
				observed_coverage, observed_coverage_lb, observed_coverage_ub = compute_coverage(sim_med_h2, est_med_h2, est_med_h2_se, coverage)
				t.write(clean_method_names[ii] + '\t' + str(eqtl_sample_size) + '\t' + 'total_med_h2' + '\t' + str(coverage) + '\t' + str(observed_coverage) + '\t' + str(observed_coverage_lb) + '\t' + str(observed_coverage_ub) + '\n')

				# causal tissue mediated h2
				observed_coverage, observed_coverage_lb, observed_coverage_ub = compute_coverage(sim_med_h2, est_per_tissue_h2[:,0], est_per_tissue_h2_se[:,0], coverage)
				t.write(clean_method_names[ii] + '\t' + str(eqtl_sample_size) + '\t' + 'causal_tissue_med_h2' + '\t' + str(coverage) + '\t' + str(observed_coverage) + '\t' + str(observed_coverage_lb) + '\t' + str(observed_coverage_ub) + '\n')

				# Non-causal tissue mediated h2
				observed_coverage, observed_coverage_lb, observed_coverage_ub = compute_coverage(sim_med_h2*0.0, est_per_tissue_h2[:,1], est_per_tissue_h2_se[:,1], coverage)
				t.write(clean_method_names[ii] + '\t' + str(eqtl_sample_size) + '\t' + 'non_causal_tissue_med_h2' + '\t' + str(coverage) + '\t' + str(observed_coverage) + '\t' + str(observed_coverage_lb) + '\t' + str(observed_coverage_ub) + '\n')



	t.close()

	print(avg_results_summary_file)
	return




def average_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, permuted_eqtls=False, variance_weighting=False, weighted=True):
	t = open(avg_results_summary_file,'w')
	t.write('method\teqtl_sample_size\theritability_type\tsim_h2\test_h2\test_h2_lb\test_h2_ub\n')
	for ii,method in enumerate(methods):
		for eqtl_sample_size in eqtl_sample_sizes:
			sim_total_h2 = []
			sim_med_h2 = []
			sim_nm_h2 = []
			est_nm_h2 = []
			est_nm_h2_se = []
			est_med_h2 = []
			est_med_h2_se = []
			est_total_h2 = []
			est_total_h2_se = []
			est_per_tissue_h2 =[]
			est_per_tissue_h2_se = []
			est_per_category_h2 = []
			est_per_category_h2_se = []
			est_eqtl_h2 = []
			sim_counter = []
			for sim_num in sim_nums:
				if weighted:
					filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_2_cis_window_500000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + run_string + '_' + str(eqtl_sample_size) + '_weighted' + '_' + method + '_jk_results.txt'
				else:
					filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_2_cis_window_500000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + run_string + '_' + str(eqtl_sample_size) + '_unweighted' + '_' + method + '_jk_results.txt'
				if sim_num in invalid_sims:
					continue
				if os.path.isfile(filer) == False:
					print(sim_num)
					print('miss')
					continue

				f = open(filer)
				arr_est = []
				arr_names = []
				arr_est_se = []
				per_tissue_est = []
				per_tissue_est_se = []
				per_category_est = []
				per_category_est_se = []
				head_count = 0
				for line in f:
					line = line.rstrip()
					data = line.split()
					if head_count == 0:
						head_count = head_count + 1
						continue
					arr_est.append(float(data[1]))
					arr_names.append(data[0])
					#arr_est_se.append(float(data[3]))
					if data[0].startswith('category_med'):
						per_category_est.append(float(data[1]))
						#per_category_est_se.append(float(data[3]))
					if data[0].startswith('dataset_med'):
						per_tissue_est.append(float(data[1]))
						#per_tissue_est_se.append(float(data[3]))
				f.close()
				arr_est = np.asarray(arr_est)
				arr_names = np.asarray(arr_names)
				arr_est_se = np.asarray(arr_est_se)
				per_tissue_est = np.asarray(per_tissue_est)
				per_tissue_est_se = np.asarray(per_tissue_est_se)
				per_category_est = np.asarray(per_category_est)
				per_category_est_se = np.asarray(per_category_est_se)


				tmp_total_h2_sim, tmp_total_nm_h2_sim, tmp_total_med_h2_sim = sim_heritabilities[sim_num]

				sim_total_h2.append(tmp_total_h2_sim)
				sim_med_h2.append(tmp_total_med_h2_sim)
				sim_nm_h2.append(tmp_total_nm_h2_sim)
				est_med_h2.append(float(arr_est[2]))
				#est_med_h2_se.append(float(arr_est_se[2]))
				est_nm_h2.append(float(arr_est[1]))
				#est_nm_h2_se.append(float(arr_est_se[1]))
				est_total_h2.append(float(arr_est[0]))
				#est_total_h2_se.append(float(arr_est_se[0]))
				est_per_tissue_h2.append(per_tissue_est)
				#est_per_tissue_h2_se.append(per_tissue_est_se)
				est_per_category_h2.append(per_category_est)
				#est_per_category_h2_se.append(per_category_est_se)
				est_eqtl_h2.append(float(arr_est[-1]))
				sim_counter.append(sim_num)

			est_per_tissue_h2 = np.asarray(est_per_tissue_h2)
			est_per_category_h2 = np.asarray(est_per_category_h2)
			est_per_tissue_h2_se = np.asarray(est_per_tissue_h2_se)
			est_per_category_h2_se = np.asarray(est_per_category_h2_se)
			print(len(est_total_h2))

			if len(sim_total_h2) == 0:
				continue

			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_total_h2, est_total_h2, 'total_h2')
			t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_nm_h2, est_nm_h2, 'nm_h2')

			if clean_method_names[ii] == 'uncalibrated_mesc':
				t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_med_h2, np.asarray(est_med_h2), 'med_h2')
				t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, sim_med_h2, est_per_tissue_h2[:,0], 'causal_tissue_med_h2')
				t = print_temp_line(clean_method_names[ii], eqtl_sample_size, t, 0.0, np.sum(est_per_tissue_h2[:,1:],axis=1), 'non_causal_tissue_med_h2')
			else:
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




def average_results_across_simulations_5_causal_tissue_with_bootstrap_bias_correction(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, permuted_eqtls=False, variance_weighting=False):
	t = open(avg_results_summary_file,'w')
	t.write('method\teqtl_sample_size\theritability_type\tsim_h2\test_h2\test_h2_lb\test_h2_ub\n')
	for ii,method in enumerate(methods):
		for eqtl_sample_size in eqtl_sample_sizes:
			sim_total_h2 = []
			sim_med_h2 = []
			sim_nm_h2 = []
			est_nm_h2 = []
			est_nm_h2_se = []
			est_med_h2 = []
			est_med_h2_se = []
			est_total_h2 = []
			est_total_h2_se = []
			est_per_tissue_h2 =[]
			est_per_tissue_h2_se = []
			est_per_category_h2 = []
			est_per_category_h2_se = []
			est_eqtl_h2 = []
			for sim_num in sim_nums:
				filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_2_cis_window_500000_ss_100000_ge_h2_05_qtl_arch_default_n_tiss_5_' + run_string + '_' + str(eqtl_sample_size) + '_' + weighting + '_' + method + '_bs_results.txt'

				if sim_num in invalid_sims:
					continue
				if os.path.isfile(filer) == False:
					print(sim_num)
					print('miss')
					continue

				f = open(filer)
				arr_est = []
				arr_names = []
				arr_est_se = []
				per_tissue_est = []
				per_tissue_est_se = []
				per_category_est = []
				per_category_est_se = []
				head_count = 0
				for line in f:
					line = line.rstrip()
					data = line.split()
					if head_count == 0:
						head_count = head_count + 1
						continue
					arr_est.append(2.0*float(data[1]) - float(data[2]))
					arr_names.append(data[0])
					arr_est_se.append(float(data[3]))
					if data[0].startswith('category_med'):
						per_category_est.append(2.0*float(data[1]) - float(data[2]))
						per_category_est_se.append(float(data[3]))
					if data[0].startswith('dataset_med'):
						per_tissue_est.append(2.0*float(data[1]) - float(data[2]))
						per_tissue_est_se.append(float(data[3]))
				f.close()
				arr_est = np.asarray(arr_est)
				arr_names = np.asarray(arr_names)
				arr_est_se = np.asarray(arr_est_se)
				per_tissue_est = np.asarray(per_tissue_est)
				per_tissue_est_se = np.asarray(per_tissue_est_se)
				per_category_est = np.asarray(per_category_est)
				per_category_est_se = np.asarray(per_category_est_se)


				sim_total_h2.append(float(.3))
				sim_med_h2.append(float(.03))
				sim_nm_h2.append(float(0.0))
				est_med_h2.append(float(arr_est[1]))
				est_med_h2_se.append(float(arr_est_se[1]))
				est_nm_h2.append(float(arr_est[0]))
				est_nm_h2_se.append(float(arr_est_se[0]))
				est_total_h2.append(float(arr_est[1]) + float(arr_est[0]))
				est_per_tissue_h2.append(per_tissue_est)
				est_per_tissue_h2_se.append(per_tissue_est_se)
				est_per_category_h2.append(per_category_est)
				est_per_category_h2_se.append(per_category_est_se)
				est_eqtl_h2.append(float(arr_est[-1]))



			est_per_tissue_h2 = np.asarray(est_per_tissue_h2)
			est_per_category_h2 = np.asarray(est_per_category_h2)
			est_per_tissue_h2_se = np.asarray(est_per_tissue_h2_se)
			est_per_category_h2_se = np.asarray(est_per_category_h2_se)

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
				if variance_weighting:
					t = print_temp_line_accounting_for_variance(clean_method_names[ii], eqtl_sample_size, t, simmer, est_single_tissue, est_per_tissue_h2_se[:,tissue_num], 'med_h2_tissue' + str(tissue_num))
				else:
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



def create_mapping_from_simulation_number_to_simulated_heritabilities(simulated_trait_dir, tmp_simulation_name_string, sim_nums, inference_gt_architecture):
	mapping = {}
	for sim_number in sim_nums:
		genetic_trait_expr_med_file = simulated_trait_dir + 'simulation_' + str(sim_number) + tmp_simulation_name_string + '_gt_arch_' + inference_gt_architecture + '_expression_mediated_trait_values.txt'
		genetic_trait_nm_file = simulated_trait_dir + 'simulation_' + str(sim_number) + tmp_simulation_name_string + '_gt_arch_' + inference_gt_architecture + '_non_mediated_variant_mediated_trait_values.txt'
		sim_med_h2 = np.var(np.loadtxt(genetic_trait_expr_med_file))
		sim_nm_h2 = np.var(np.loadtxt(genetic_trait_nm_file))
		sim_h2 = np.var(np.loadtxt(genetic_trait_nm_file) + np.loadtxt(genetic_trait_expr_med_file))
		if sim_number in mapping:
			print('assumption erororo')
		mapping[sim_number] = (sim_h2, sim_nm_h2, sim_med_h2)


	return mapping


#####################
# Command line args
#####################
trait_med_h2_inference_dir = sys.argv[1]
visualize_trait_med_h2_dir = sys.argv[2]
simulated_trait_dir = sys.argv[3]
tmp_simulation_name_string = sys.argv[4]


# Simulation params
sim_nums = np.arange(1,201)
eqtl_sample_sizes = np.asarray(['100','200','300','1000'])
eqtl_sample_sizes = np.asarray(['100', '200'])
eqtl_sample_sizes = np.asarray(['100', '200', '300', '1000', '100-1000'])
eqtl_sample_sizes = np.asarray(['100', '200', '300', '1000', '100-1000'])



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
beta_squared_thresh = '1000.0'

non_med_anno = 'full_anno'
non_med_anno = 'genotype_intercept'


inference_gt_architecture = 'linear'
simulated_gt_architecture = 'linear'

inference_gt_architecture = 'linear'
simulated_gt_architecture = 'stdExpr'


gene_ld_score_type = 'ldsc_style_pred'
gene_ld_score_type = 'ashr_style_pred'
gene_ld_score_type = 'squared_marginal_sumstats'

weighting="weighted"

step1_gene_ldscores='MarginalSS'
step1_gene_ldscores='mescLassoPlusMarginalSS'
step1_gene_ldscores='mescLasso'



inference_approach="JIVE"
inference_approach="binned_2SLS"


cis_snp_h2_method="ldsc"
cis_snp_h2_method="ldsc"
cis_snp_h2_method="greml"
cis_snp_h2_method="avgChisq"
cis_snp_h2_method="true"

n_bins="6"


training_data_gene_ldscores_type="lasso_CV_corr_standardized"
validation_data_gene_ldscores_type="MarginalSS"

# Create mapping from simulation number to simulated heritabilities
sim_heritabilities = create_mapping_from_simulation_number_to_simulated_heritabilities(simulated_trait_dir, tmp_simulation_name_string, sim_nums, simulated_gt_architecture)

'''

#run_string = non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + beta_squared_thresh  + '_' + cis_snp_h2_method + '_' + training_data_gene_ldscores_type + '_' + inference_approach + '_' + cis_h2_thresh
#run_string = non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + beta_squared_thresh  + '_' + cis_snp_h2_method + '_' + training_data_gene_ldscores_type + '_'+ validation_data_gene_ldscores_type + '_' + inference_approach 
run_string = non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + beta_squared_thresh  + '_' + cis_snp_h2_method + '_' + training_data_gene_ldscores_type + '_'+ validation_data_gene_ldscores_type + '_' + inference_approach + '_' + n_bins

avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_' + run_string + '_sim_results_calibrated_ldsc_summary_averaged.txt'
clean_method_names = ['uncalibrated_mesc', 'calibrated_mesc']
methods = ['uncalibrated_mesc', 'calibrated_mesc']
average_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, variance_weighting=False)
'''




# Inf sample size
inference_approach="binned_2SLS"
inference_gt_architecture = 'linear'
simulated_gt_architecture = 'stdExpr'
beta_squared_thresh = '1000.0'
cis_snp_h2_method="true"
n_bins_arr = ['1', '2', '3', '4', '5', '6', '7']

for n_bins in n_bins_arr:
	eqtl_sample_sizes = np.asarray(['INF_SS'])


	run_string = non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + beta_squared_thresh  + '_' + cis_snp_h2_method  + '_' + inference_approach + '_' + n_bins

	avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_' + run_string + '_INF_SS_weighted_sim_results_calibrated_ldsc_summary_averaged.txt'
	clean_method_names = ['uncalibrated_mesc', 'calibrated_mesc']
	methods = ['uncalibrated_mesc', 'calibrated_mesc']
	average_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, variance_weighting=False, weighted=True)
	print(avg_results_summary_file)

inference_approach="2SLS"
inference_gt_architecture = 'stdExpr'
simulated_gt_architecture = 'stdExpr'
beta_squared_thresh = '1000.0'
cis_snp_h2_method="true"

eqtl_sample_sizes = np.asarray(['INF_SS'])

run_string = non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + beta_squared_thresh   + '_' + inference_approach 

avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_' + run_string + '_INF_SS_weighted_sim_results_calibrated_ldsc_summary_averaged.txt'
clean_method_names = ['uncalibrated_mesc', 'calibrated_mesc']
methods = ['uncalibrated_mesc', 'calibrated_mesc']
average_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, variance_weighting=False, weighted=True)
print(avg_results_summary_file)


'''
inference_approach="2SLS"
inference_gt_architecture = 'linear'
simulated_gt_architecture = 'stdExpr'
beta_squared_thresh = '1000.0'
cis_snp_h2_method="true"

eqtl_sample_sizes = np.asarray(['INF_SS'])

run_string = non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + beta_squared_thresh   + '_' + inference_approach 

avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_' + run_string + '_INF_SS_sim_results_calibrated_ldsc_summary_averaged.txt'
clean_method_names = ['uncalibrated_mesc', 'calibrated_mesc']
methods = ['uncalibrated_mesc', 'calibrated_mesc']
average_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, variance_weighting=False)
print(avg_results_summary_file)


# Create mapping from simulation number to simulated heritabilities
simulated_gt_architecture = 'linear'
sim_heritabilities = create_mapping_from_simulation_number_to_simulated_heritabilities(simulated_trait_dir, tmp_simulation_name_string, sim_nums, simulated_gt_architecture)

inference_approach="2SLS"
inference_gt_architecture = 'linear'
simulated_gt_architecture = 'linear'
beta_squared_thresh = '1000.0'
cis_snp_h2_method="true"

eqtl_sample_sizes = np.asarray(['INF_SS'])

run_string = non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + beta_squared_thresh   + '_' + inference_approach 

avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_' + run_string + '_INF_SS_sim_results_calibrated_ldsc_summary_averaged.txt'
clean_method_names = ['uncalibrated_mesc', 'calibrated_mesc']
methods = ['uncalibrated_mesc', 'calibrated_mesc']
average_results_across_simulations_5_causal_tissue(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, variance_weighting=False)
print(avg_results_summary_file)

'''









'''
calibration_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_' + eqtl_snp_representation + '_' + non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + gene_ld_score_type + '_' + beta_squared_thresh+ '_' + step1_gene_ldscores+ '_' + inference_approach + '_' + weighting + '_sim_results_calibration_summary.txt'
clean_method_names = ['calibrated_mesc']
methods = ['calibrated_mesc']
confidence_interval_calibration(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, calibration_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, variance_weighting=False)

t1e_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_' + eqtl_snp_representation + '_' + non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + gene_ld_score_type + '_' + beta_squared_thresh+ '_' + step1_gene_ldscores+ '_' + inference_approach + '_' + weighting + '_sim_results_t1e_summary.txt'
clean_method_names = ['calibrated_mesc']
methods = ['calibrated_mesc']
t1e_summary(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, t1e_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, variance_weighting=False)




power_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_' + eqtl_snp_representation + '_' + non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + gene_ld_score_type + '_' + beta_squared_thresh+ '_' + step1_gene_ldscores+ '_' + inference_approach + '_' + weighting + '_sim_results_power_summary.txt'
clean_method_names = ['uncalibrated_mesc','calibrated_mesc']
methods = ['uncalibrated_mesc','calibrated_mesc']
power_summary(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, power_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, variance_weighting=False)



fstat_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_' + eqtl_snp_representation + '_' + non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + gene_ld_score_type + '_' + beta_squared_thresh+ '_' + step1_gene_ldscores+ '_' + inference_approach + '_' + weighting + '_sim_results_fstat_summary.txt'
clean_method_names = ['calibrated_mesc']
methods = ['calibrated_mesc']
fstat_summary(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, fstat_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, sim_heritabilities, variance_weighting=False)

'''







avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_5_causal_tissue_' + eqtl_snp_representation + '_' + non_med_anno + '_' + simulated_gt_architecture + '_' +inference_gt_architecture + '_' + gene_ld_score_type + '_' + beta_squared_thresh + '_' + weighting + '_sim_results_calibrated_bs_bias_correction_ldsc_summary_averaged.txt'
clean_method_names = ['uncalibrated_mesc', 'calibrated_mesc']
methods = ['uncalibrated_mesc', 'calibrated_mesc']
#average_results_across_simulations_5_causal_tissue_with_bootstrap_bias_correction(sim_nums, eqtl_sample_sizes, trait_med_h2_inference_dir, avg_results_summary_file, methods, clean_method_names, invalid_sims, run_string, weighting, variance_weighting=False)





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


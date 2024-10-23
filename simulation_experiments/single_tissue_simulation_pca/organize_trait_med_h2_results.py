import numpy as np
import os
import sys
import pdb







def create_mapping_from_simulation_number_to_simulated_h2_parameters(trait_med_h2_inference_dir, sim_nums):
	mapping = {}

	for sim_num in sim_nums:
		filer = trait_med_h2_inference_dir + 'simulation_' + str(sim_num) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_eqtl_SS_100_med_h2_marginal_gibbs_sampler_resid_var_True_cc_0.0.txt'
		f = open(filer)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			total_h2 = float(data[0])
			nm_h2 = float(data[1])
			med_h2 = float(data[2])
			mapping[sim_num] = (total_h2, nm_h2, med_h2)
			break
		f.close()
	return mapping


def extract_sampled_params_from_marginal_gibbs_sampler_file(marginal_gibbs_sampler_file, burn_in_lb=15000, burn_in_ub=5000000000000):
	est_nm_h2 = []
	alt_est_nm_h2 = []
	est_med_h2 = []
	alt_est_med_h2 = []
	est_total_h2 = []

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
		alt_nm_h2 = float(data[5])
		med_h2 = float(data[6])
		alt_med_h2 = float(data[7])
		total_h2 = float(data[8])

		est_nm_h2.append(nm_h2)
		alt_est_nm_h2.append(alt_nm_h2)
		est_med_h2.append(med_h2)
		alt_est_med_h2.append(alt_med_h2)
		est_total_h2.append(total_h2)

	f.close()

	return np.asarray(est_nm_h2), np.asarray(alt_est_nm_h2), np.asarray(est_med_h2), np.asarray(alt_est_med_h2), np.asarray(est_total_h2)

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
	t.write('iteration\teqtl_sample_size\tmethod\tsim_total_h2\tsim_nm_h2\tsim_med_h2\test_nm_h2\talt_est_nm_h2\test_med_h2\talt_est_med_h2\test_total_h2\n')

	# Loop through simulations
	for sim_iter in sim_nums:
		sim_total_h2, sim_nm_h2, sim_med_h2 = sim_num_to_sim_h2_params[sim_iter]

		# Loop through eqtl sample sizes
		for eqtl_sample_size in eqtl_sample_sizes:

			###########################################
			# First do method of marginal PCA gibbs sampler
			marginal_gibbs_sampler_file = trait_med_h2_inference_dir + 'simulation_' + str(sim_iter) + '_chrom1_cis_window_100000_ss_100000_ge_h2_05_qtl_arch_default_eqtl_SS_' + str(eqtl_sample_size) + '_med_h2_marginal_gibbs_sampler_resid_var_True_cc_0.0_alt_.txt'
			
			# Extract distributions
			est_nm_h2, alt_est_nm_h2, est_med_h2, alt_est_med_h2, est_total_h2 = extract_sampled_params_from_marginal_gibbs_sampler_file(marginal_gibbs_sampler_file, burn_in_lb=iteration_lb, burn_in_ub=iteration_ub)

			# Get 95 percent cis
			alt_est_nm_lb, alt_est_nm_ub = get_95_percent_ci_based_on_sampled_distribution(alt_est_nm_h2)
			alt_est_med_lb, alt_est_med_ub = get_95_percent_ci_based_on_sampled_distribution(alt_est_med_h2)
			total_est_lb, total_est_ub = get_95_percent_ci_based_on_sampled_distribution(est_total_h2)

			#print(str(sim_total_h2) + '\t' + '[' + str(total_est_lb) + ', ' + str(total_est_ub) + ']')
			#print(str(sim_nm_h2) + '\t' + '[' + str(alt_est_nm_lb) + ', ' + str(alt_est_nm_ub) + ']')
			#print(str(sim_med_h2) + '\t' + '[' + str(alt_est_med_lb) + ', ' + str(alt_est_med_ub) + ']')


			# print to output
			t.write(str(sim_iter) + '\t' + str(eqtl_sample_size) + '\t' + 'marginal_pca_sumstat_gibbs_0.0' + '\t' + str(sim_total_h2) + '\t' + str(sim_nm_h2) + '\t' + str(sim_med_h2) + '\t')
			t.write(str(np.mean(est_nm_h2)) + '\t' + str(np.mean(alt_est_nm_h2)) + '\t' + str(np.mean(est_med_h2)) + '\t' + str(np.mean(alt_est_med_h2)) + '\t' + str(np.mean(est_total_h2)) + '\n')

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
				est_alt_nm_h2.append(float(data[7]))
				est_med_h2.append(float(data[8]))
				est_alt_med_h2.append(float(data[9]))
				est_total_h2.append(float(data[10]))
			f.close()


			t = print_temp_line(method, eqtl_sample_size, t, sim_total_h2, est_total_h2, 'total_h2')
			t = print_temp_line(method, eqtl_sample_size, t, sim_nm_h2, est_nm_h2, 'nm_h2')
			t = print_temp_line(method, eqtl_sample_size, t, sim_nm_h2, est_alt_nm_h2, 'alt_nm_h2')
			t = print_temp_line(method, eqtl_sample_size, t, sim_med_h2, est_med_h2, 'med_h2')
			t = print_temp_line(method, eqtl_sample_size, t, sim_med_h2, est_alt_med_h2, 'alt_med_h2')
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
raw_sim_nums = np.arange(1,101)
eqtl_sample_sizes = np.asarray([100,300,1000,10000])

sim_nums = []
for sim_num in raw_sim_nums:
	if sim_num == 45 or sim_num == 60:
		continue
	sim_nums.append(sim_num)
sim_nums = np.asarray(sim_nums)


# First create mapping from simulation number to simulated h2 parameters
sim_num_to_sim_h2_params = create_mapping_from_simulation_number_to_simulated_h2_parameters(trait_med_h2_inference_dir, sim_nums)


# Standard averaging
iteration_lb = 950000
iteration_ub = 1020000
# Organize results to per iteration summary
results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_results_summary.txt'
generate_per_iteration_results_summary_file(sim_nums, trait_med_h2_inference_dir, results_summary_file,eqtl_sample_sizes, sim_num_to_sim_h2_params, iteration_lb, iteration_ub)

# average results across per iterations
avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_results_summary_averaged.txt'
methods = np.asarray(['marginal_pca_sumstat_gibbs_0.0'])
average_results_across_simulations(results_summary_file, avg_results_summary_file, eqtl_sample_sizes,methods)

# Averaging across various windows
windows = np.arange(0, 1002000, 25000)

for ii,window_ub in enumerate(windows):
	if ii == 0:
		continue

	window_lb = windows[(ii-1)]

	print(str(window_lb) + ' : ' + str(window_ub))

	results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_results_summary_' + str(window_lb) + ':' + str(window_ub) + '.txt'
	generate_per_iteration_results_summary_file(sim_nums, trait_med_h2_inference_dir, results_summary_file,eqtl_sample_sizes, sim_num_to_sim_h2_params, window_lb, window_ub)

	# average results across per iterations
	avg_results_summary_file = visualize_trait_med_h2_dir+ 'med_h2_results_summary_averaged' + str(window_lb) + ':' + str(window_ub) + '.txt'
	methods = np.asarray(['marginal_pca_sumstat_gibbs_0.0'])
	average_results_across_simulations(results_summary_file, avg_results_summary_file, eqtl_sample_sizes,methods)




'''
sim_nums = np.arange(101,199)
eqtl_sample_sizes = np.asarray([100,300,1000,10000])
avg_results_summary_file = visualize_trait_med_h2_dir+ 'tglr_med_h2_results_summary_averaged.txt'

average_results_across_tglr_simulations(avg_results_summary_file, eqtl_sample_sizes,trait_med_h2_inference_dir, sim_nums)

'''



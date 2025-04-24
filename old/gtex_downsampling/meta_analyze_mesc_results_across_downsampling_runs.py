import numpy as np
import os
import sys
import pdb



def extract_med_h2_and_med_h2_se_from_mesc_res_file(mesc_res_file):
	g = open(mesc_res_file)
	med_h2 = []
	med_h2_se = []

	head_count = 0
	for line in g:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		med_h2.append(float(data[1]))
		med_h2_se.append(float(data[2]))


	g.close()

	return np.asarray(med_h2), np.asarray(med_h2_se)

def meta_analysis(means, std_errors):
	"""
	Perform meta-analysis of means given means and standard errors.
	
	Parameters:
	- means: list or numpy array of mean values from each study
	- std_errors: list or numpy array of standard errors corresponding to each mean
	
	Returns:
	- meta_mean: meta-analyzed mean
	- meta_std_error: standard error of the meta-analyzed mean
	"""
	# Convert inputs to numpy arrays to ensure calculations work smoothly
	means = np.array(means)
	std_errors = np.array(std_errors)
	
	# Calculate weights for each study based on inverse variance (1/se^2)
	weights = 1 / std_errors**2
	
	# Calculate weighted mean
	meta_mean = np.sum(weights * means) / np.sum(weights)
	
	# Calculate standard error of the meta-analyzed mean
	meta_std_error = np.sqrt(1 / np.sum(weights))

	return meta_mean, meta_std_error



###########################
# Command line args
###########################
downsample_run_summary_file = sys.argv[1]
mesc_results_dir = sys.argv[2]
variant_anno = sys.argv[3]


# Open file handle for meta-analyzed results
meta_analyzed_results_output_file = mesc_results_dir + 'downsampled_med_h2_meta_analyzed_summary_' + variant_anno + '.txt'
t = open(meta_analyzed_results_output_file,'w')
t.write('downsample_run_name\tsample_size\tmed_h2\tmed_h2_se\tmed_h2_lb\tmed_h2_ub\tunweighted_mean\n')


# Loop through downsampling runs
f = open(downsample_run_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue

	downsample_run_name = data[0]
	downsample_ss = data[1]

	mesc_res_file = mesc_results_dir + 'cross_trait_med_h2_summary_' + downsample_run_name + '_' + variant_anno + '.txt'
	
	# Extract vector of mediated h2 and mediated h2 se from mesc_res_file
	med_h2, med_h2_se = extract_med_h2_and_med_h2_se_from_mesc_res_file(mesc_res_file)

	# Run fixed effect meta-analysis
	mean_med_h2, mean_med_h2_se = meta_analysis(med_h2, med_h2_se)

	# Get 95% CI
	mean_med_h2_lb = mean_med_h2 - (1.96*mean_med_h2_se)
	mean_med_h2_ub = mean_med_h2 + (1.96*mean_med_h2_se)

	# Write to output
	t.write(downsample_run_name + '\t' + downsample_ss + '\t' + str(mean_med_h2) + '\t' + str(mean_med_h2_se) + '\t' + str(mean_med_h2_lb) + '\t' + str(mean_med_h2_ub)+ '\t' + str(np.mean(med_h2)) + '\n')


# Add full data
downsample_run_name = 'full_data'

downsample_ss = '13072'

mesc_res_file = mesc_results_dir + 'cross_trait_med_h2_summary_' + downsample_run_name + '_' + variant_anno + '.txt'
	
# Extract vector of mediated h2 and mediated h2 se from mesc_res_file
med_h2, med_h2_se = extract_med_h2_and_med_h2_se_from_mesc_res_file(mesc_res_file)

# Run fixed effect meta-analysis
mean_med_h2, mean_med_h2_se = meta_analysis(med_h2, med_h2_se)

# Get 95% CI
mean_med_h2_lb = mean_med_h2 - (1.96*mean_med_h2_se)
mean_med_h2_ub = mean_med_h2 + (1.96*mean_med_h2_se)

# Write to output
t.write(downsample_run_name + '\t' + downsample_ss + '\t' + str(mean_med_h2) + '\t' + str(mean_med_h2_se) + '\t' + str(mean_med_h2_lb) + '\t' + str(mean_med_h2_ub) + '\t' + str(np.mean(med_h2)) + '\n')



f.close()
t.close()


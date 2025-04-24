import numpy as np
import os
import sys
import pdb



def extract_trait_names(non_redundent_summary_statistics_file):
	f = open(non_redundent_summary_statistics_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)


def aggregate_frac_mediated_H2_estimates_cross_traits(mediated_h2_estimates_summary_file, trait_names, tglr_results_dir, qtl_version):
	t = open(mediated_h2_estimates_summary_file,'w')
	t.write('trait_name\tfrac_med_h2\tfrac_med_h2_se\tfrac_med_h2_lb\tfrac_med_h2_ub\n')
	for trait_name in trait_names:
		med_h2_summary = tglr_results_dir + trait_name + '_' + qtl_version + '_baselineLD_no_qtl_tglr_frac_h2_med_5_50.txt'
		aa = np.loadtxt(med_h2_summary,dtype=str,delimiter='\t')
		h2_med = float(aa[1,0])
		h2_med_se = float(aa[1,1])
		h2_med_lb = h2_med - (1.96*h2_med_se)
		h2_med_ub = h2_med + (1.96*h2_med_se)
		t.write(trait_name + '\t' + str(h2_med) + '\t' + str(h2_med_se) + '\t' + str(h2_med_lb) + '\t' + str(h2_med_ub) + '\n')
	t.close()
	return

def meta_analyze_normalized_per_gene_h2_estimates(meta_analyzed_normalized_per_gene_h2_file, trait_names, tglr_results_dir, qtl_version):
	per_gene_h2s = []
	per_gene_h2_ses = []
	for trait_name in trait_names:
		# Extract total h2
		total_h2_file = tglr_results_dir + trait_name + '_' + qtl_version + '_baselineLD_no_qtl_tglr_mediated_h2_5_50.txt'
		aa = np.loadtxt(total_h2_file,dtype=str,delimiter='\t')
		total_h2 = np.sum(aa[1:,1].astype(float))

		# Extract per gene h2
		per_gene_h2_file = tglr_results_dir + trait_name + '_' + qtl_version + '_baselineLD_no_qtl_tglr_organized_coef.txt'
		raw_coefs = np.loadtxt(per_gene_h2_file,dtype=str,delimiter='\t')
		coefs = raw_coefs[94:,1].astype(float)
		coefs_se = raw_coefs[94:,2].astype(float)
		coef_names = raw_coefs[94:,0]

		per_gene_h2s.append(coefs/total_h2)
		per_gene_h2_ses.append(coefs_se/total_h2)
	per_gene_h2s = np.asarray(per_gene_h2s)
	per_gene_h2_ses = np.asarray(per_gene_h2_ses)
	t = open(meta_analyzed_normalized_per_gene_h2_file,'w')
	t.write('eqtl_category\tmean\tmean_se\tlb\tub\n')
	for ii, coef_name in enumerate(coef_names):
		tmp_per_gene_h2 = per_gene_h2s[:,ii]
		tmp_per_gene_h2_se = per_gene_h2_ses[:,ii]
		meta_mean, meta_mean_se = meta_analysis(tmp_per_gene_h2, tmp_per_gene_h2_se)
		lb = meta_mean - (1.96*meta_mean_se)
		ub = meta_mean + (1.96*meta_mean_se)
		t.write(coef_name + '\t' + str(meta_mean) + '\t' + str(meta_mean_se) + '\t' + str(lb) + '\t' + str(ub) + '\n')
	t.close()
	return


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

def meta_analyze_mediated_h2_estimates(frac_mediated_h2_estimates_summary_file, meta_analyze_frac_med_h2_estimates):
	tmp = np.loadtxt(frac_mediated_h2_estimates_summary_file, dtype=str,delimiter='\t')
	med_h2_arr = np.asarray(tmp[1:,1]).astype(float)
	med_h2_se_arr = np.asarray(tmp[1:,2]).astype(float)

	meta_mean, meta_mean_se = meta_analysis(med_h2_arr, med_h2_se_arr)
	lb = meta_mean - (1.96*meta_mean_se)
	ub = meta_mean + (1.96*meta_mean_se)
	t = open(meta_analyze_frac_med_h2_estimates,'w')
	t.write('mean\tmean_se\tlb\tub\n')
	t.write(str(meta_mean) + '\t' + str(meta_mean_se) + '\t' + str(lb) + '\t' + str(ub) + '\n')
	t.close()
	return


####################
# Command line args
####################
non_redundent_summary_statistics_file = sys.argv[1]
tglr_results_dir = sys.argv[2]
qtl_version = sys.argv[3]

trait_names = extract_trait_names(non_redundent_summary_statistics_file)


# Mediated heritability estimates summary file
frac_mediated_h2_estimates_summary_file = tglr_results_dir + 'summary_frac_mediated_h2_across_traits_' + qtl_version + '.txt'
aggregate_frac_mediated_H2_estimates_cross_traits(frac_mediated_h2_estimates_summary_file, trait_names, tglr_results_dir, qtl_version)

# Meta analyze mediated h2 estimates
meta_analyze_frac_med_h2_estimates = tglr_results_dir + 'summary_frac_mediated_h2_meta_analyzed_cross_traits_' + qtl_version + '.txt'
meta_analyze_mediated_h2_estimates(frac_mediated_h2_estimates_summary_file, meta_analyze_frac_med_h2_estimates)



# Meta analyze normalized per gene heritabilities
meta_analyzed_normalized_per_gene_h2_file = tglr_results_dir + 'normalize_per_gene_heritability_meta_analyzed_cross_traits_' + qtl_version + '.txt'
meta_analyze_normalized_per_gene_h2_estimates(meta_analyzed_normalized_per_gene_h2_file, trait_names, tglr_results_dir, qtl_version)

print(meta_analyzed_normalized_per_gene_h2_file)







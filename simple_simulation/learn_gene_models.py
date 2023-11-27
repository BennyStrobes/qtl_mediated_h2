import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os
import pdb
import statsmodels.api as sm



def extract_non_colinear_predictors(corry, abs_correlation_threshold):
	valid_predictors = []
	valid_predictors.append(0)
	for predictor in range(1,corry.shape[0]):
		if np.max(np.abs(corry[predictor, valid_predictors])) < abs_correlation_threshold:
			valid_predictors.append(predictor)
	valid_predictors = np.asarray(valid_predictors)
	return valid_predictors


def estimate_causal_eqtl_effects_for_a_single_gene(expression_vec, genotype_mat):
	n_predictors = genotype_mat.shape[1]
	corry =  np.corrcoef(np.transpose(genotype_mat))
	non_colinear_predictors = extract_non_colinear_predictors(corry, .75)

	ols = sm.OLS(expression_vec, genotype_mat[:, non_colinear_predictors])
	ols_result = ols.fit()

	# Initialize matrix
	beta = np.zeros(n_predictors)
	beta_varcov = np.zeros((n_predictors, n_predictors))
	# Fill in matrix
	beta[non_colinear_predictors] = ols_result.params
	beta_varcov[non_colinear_predictors[:,None], non_colinear_predictors[None,:]] = ols_result.cov_params()

	return beta, beta_varcov, ols_result.rsquared_adj







######################
# Command line args
######################
simulation_name_string = sys.argv[1]
processed_genotype_data_dir = sys.argv[2]
simulated_eqtl_data_dir = sys.argv[3]
simulated_expression_data_dir = sys.argv[4]
simulated_gene_models_dir = sys.argv[5]
eqtl_ss = sys.argv[6]


# Load in genotype matrix
genotype_file = processed_genotype_data_dir + 'eqtl_' + str(eqtl_ss) + '_genotype_1.npy'
genotype = np.load(genotype_file)
# Load in gene expression
expression_file = simulated_expression_data_dir + simulation_name_string + 'eqtl_ss_' + str(eqtl_ss) + '.npy'
gene_expression = np.load(expression_file)


# Open output file handle
output_file = simulated_gene_models_dir + simulation_name_string + 'model_summaries_' + str(eqtl_ss) + '.txt'
t = open(output_file,'w')
t.write('gene_name\tgene_snp_indices\tgene_snp_causal_eqtl_effects\testimated_effect_file\testimated_effect_varcov_file\tgene_h2_est\n')
# Loop through genes
simulated_causal_eqtl_effect_file = simulated_eqtl_data_dir + simulation_name_string + 'causal_eqtl_effect_summary.txt'
f = open(simulated_causal_eqtl_effect_file)
head_count = 0
gene_index = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract relevent fields
	gene_name = data[0]
	gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
	gene_snp_causal_effects = np.asarray(data[2].split(',')).astype(float)
	if int(gene_name.split('_')[1]) != gene_index:
		print('assumption eroror')

	expression_vec = gene_expression[gene_index,:]

	beta, beta_varcov, gene_h2_est = estimate_causal_eqtl_effects_for_a_single_gene(expression_vec, genotype[:, gene_snp_indices])

	# Save beta file
	beta_file = simulated_gene_models_dir + simulation_name_string + 'estimated_eqtl_beta_' + str(eqtl_ss) + '_' + gene_name + '.npy'
	np.save(beta_file, beta)
	# save beta-varcov file
	beta_varcov_file = simulated_gene_models_dir + simulation_name_string + 'estimated_eqtl_beta_varcov_' + str(eqtl_ss) + '_' + gene_name + '.npy'
	np.save(beta_varcov_file, beta_varcov)


	t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + beta_file + '\t' + beta_varcov_file + '\t' + str(gene_h2_est) + '\n')

	gene_index = gene_index + 1
f.close()
t.close()

print(output_file)


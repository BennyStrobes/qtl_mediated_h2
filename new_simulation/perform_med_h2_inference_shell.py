import numpy as np 
import os
import sys
import pdb
import rss_vi_variant_only
from sklearn.linear_model import LinearRegression

def load_in_gwas_z_scores(gwas_sum_stats_file):
	f = open(gwas_sum_stats_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(float(data[3]))
	f.close()
	arr = np.asarray(arr)
	return arr

def load_in_gwas_se(gwas_sum_stats_file):
	f = open(gwas_sum_stats_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(float(data[2]))
	f.close()
	arr = np.asarray(arr)
	return arr	



def run_ld_score_regression(gwas_z_scores, LD, N_gwas, intercept_bool=True):
	chi_sq = np.square(gwas_z_scores)

	ld_scores = np.sum(np.square(LD),axis=0)

	n_snps = len(gwas_z_scores)

	if intercept_bool:
		model = LinearRegression(fit_intercept=True, normalize=False)
	else:
		model = LinearRegression(fit_intercept=False, normalize=False)
		chi_sq = chi_sq - 1


	ldsc_fit = model.fit(ld_scores.reshape((len(ld_scores),1)), chi_sq)

	ldsc_h2 = ldsc_fit.coef_[0]*n_snps/N_gwas

	return ldsc_h2

def run_mpldsc_with_marginal_iid_eqtl_effects(gwas_z_scores, LD, N_gwas,gene_summary_file, intercept_bool=True):
	chi_sq = np.square(gwas_z_scores)

	squared_ld = np.square(LD)

	ld_scores = np.sum(squared_ld,axis=0)

	n_snps = len(gwas_z_scores)

	gene_ld_scores = np.copy(ld_scores)*0.0

	f = open(gene_summary_file)
	head_count = 0
	n_genes = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_snp_indices = np.asarray(data[1].split(',')).astype(int)

		gene_ld_scores = gene_ld_scores + np.sum(squared_ld[:, gene_snp_indices],axis=1)/len(gene_snp_indices)
		n_genes = n_genes + 1
	f.close()

	if intercept_bool:
		model = LinearRegression(fit_intercept=True, normalize=False)
	else:
		model = LinearRegression(fit_intercept=False, normalize=False)
		chi_sq = chi_sq - 1

	joint_ld_scores = np.transpose(np.vstack((ld_scores,gene_ld_scores)))

	mpldsc_fit = model.fit(joint_ld_scores, chi_sq)


	mpldsc_var_h2 = mpldsc_fit.coef_[0]*n_snps/N_gwas
	mpldsc_gene_h2 = mpldsc_fit.coef_[1]*n_genes/N_gwas

	return mpldsc_var_h2, mpldsc_gene_h2

def run_mpldsc_with_marginal_eqtl_effects_using_eqtl_ld(gwas_z_scores, LD, N_gwas,gene_summary_file, eqtl_ss, intercept_bool=True):
	chi_sq = np.square(gwas_z_scores)

	squared_ld = np.square(LD)

	adj_squared_ld = squared_ld - ((1-squared_ld)/(eqtl_ss-2))

	ld_scores = np.sum(adj_squared_ld,axis=0)

	n_snps = len(gwas_z_scores)

	gene_ld_scores = np.copy(ld_scores)*0.0

	f = open(gene_summary_file)
	head_count = 0
	n_genes = 0
	gene_h2s = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
		marginal_beta_file = data[9]
		marginal_betas = np.load(marginal_beta_file)
		sampling_variance = float(data[7])
		gene_h2s.append(float(data[5]))

		#gene_ld_scores = gene_ld_scores + np.square(marginal_betas) + sampling_variance
		squared_marginal_betas = np.square(marginal_betas)
		adj_squared_marginal_betas = squared_marginal_betas - ((1-squared_marginal_betas)/(eqtl_ss-2))
		gene_ld_scores = gene_ld_scores + adj_squared_marginal_betas

		n_genes = n_genes + 1
	f.close()

	if intercept_bool:
		model = LinearRegression(fit_intercept=True, normalize=False)
	else:
		model = LinearRegression(fit_intercept=False, normalize=False)
		chi_sq = chi_sq - 1

	joint_ld_scores = np.transpose(np.vstack((ld_scores,gene_ld_scores)))

	mpldsc_fit = model.fit(joint_ld_scores, chi_sq)


	mpldsc_var_h2 = mpldsc_fit.coef_[0]*n_snps/N_gwas
	mpldsc_gene_h2 = mpldsc_fit.coef_[1]*n_genes/N_gwas


	return mpldsc_var_h2, mpldsc_gene_h2

def run_mpldsc_with_marginal_eqtl_effects(gwas_z_scores, LD, N_gwas,gene_summary_file, eqtl_ss, intercept_bool=True):
	chi_sq = np.square(gwas_z_scores)

	squared_ld = np.square(LD)

	ld_scores = np.sum(squared_ld,axis=0)

	n_snps = len(gwas_z_scores)

	gene_ld_scores = np.copy(ld_scores)*0.0

	f = open(gene_summary_file)
	head_count = 0
	n_genes = 0
	gene_h2s = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
		marginal_beta_file = data[6]
		marginal_betas = np.load(marginal_beta_file)
		sampling_variance = float(data[7])
		gene_h2s.append(float(data[5]))

		#gene_ld_scores = gene_ld_scores + np.square(marginal_betas) + sampling_variance
		squared_marginal_betas = np.square(marginal_betas)
		adj_squared_marginal_betas = squared_marginal_betas - ((1-squared_marginal_betas)/(eqtl_ss-2))
		gene_ld_scores = gene_ld_scores + adj_squared_marginal_betas

		n_genes = n_genes + 1
	f.close()

	if intercept_bool:
		model = LinearRegression(fit_intercept=True, normalize=False)
	else:
		model = LinearRegression(fit_intercept=False, normalize=False)
		chi_sq = chi_sq - 1

	joint_ld_scores = np.transpose(np.vstack((ld_scores,gene_ld_scores)))

	mpldsc_fit = model.fit(joint_ld_scores, chi_sq)


	mpldsc_var_h2 = mpldsc_fit.coef_[0]*n_snps/N_gwas
	mpldsc_gene_h2 = mpldsc_fit.coef_[1]*n_genes*.025/N_gwas


	return mpldsc_var_h2, mpldsc_gene_h2

def run_mpldsc_with_marginal_eqtl_effects_w_uncertainty(gwas_z_scores, LD, N_gwas,gene_summary_file, eqtl_ss, intercept_bool=True):
	chi_sq = np.square(gwas_z_scores)

	squared_ld = np.square(LD)

	ld_scores = np.sum(squared_ld,axis=0)

	n_snps = len(gwas_z_scores)

	gene_ld_scores = np.copy(ld_scores)*0.0

	f = open(gene_summary_file)
	head_count = 0
	n_genes = 0
	gene_h2s = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
		marginal_beta_file = data[6]
		marginal_betas = np.load(marginal_beta_file)
		sampling_variance = float(data[7])
		gene_h2s.append(float(data[5]))

		#gene_ld_scores = gene_ld_scores + np.square(marginal_betas) + sampling_variance

		marginal_betas = marginal_betas
		squared_marginal_betas = np.square(marginal_betas)
		adj_squared_marginal_betas = squared_marginal_betas #- ((1-squared_marginal_betas)/(eqtl_ss-2))

		#adj_squared_marginal_betas[gene_snp_indices] = adj_squared_marginal_betas[gene_snp_indices] - sampling_variance
		gene_ld_scores = gene_ld_scores + adj_squared_marginal_betas

		n_genes = n_genes + 1
	f.close()

	if intercept_bool:
		model = LinearRegression(fit_intercept=True, normalize=False)
	else:
		model = LinearRegression(fit_intercept=False, normalize=False)
		chi_sq = chi_sq - 1

	joint_ld_scores = np.transpose(np.vstack((ld_scores,gene_ld_scores)))

	mpldsc_fit = model.fit(joint_ld_scores, chi_sq)


	mpldsc_var_h2 = mpldsc_fit.coef_[0]*n_snps/N_gwas
	mpldsc_gene_h2 = mpldsc_fit.coef_[1]*n_genes*.025/N_gwas


	return mpldsc_var_h2, mpldsc_gene_h2



def run_mpldsc_he_with_marginal_eqtl_effects(gwas_z_scores, LD, LD_LD_t, N_gwas,gene_summary_file, eqtl_ss):
	chi_sq = np.square(gwas_z_scores)

	squared_ld = np.square(LD)

	n_snps = len(gwas_z_scores)


	z_mat = np.dot(gwas_z_scores.reshape(len(gwas_z_scores),1), gwas_z_scores.reshape(1,len(gwas_z_scores)))

	ld_scores = np.copy(LD_LD_t)
	gene_ld_scores = np.copy(z_mat)*0.0
	#np.fill_diagonal(ld_scores, np.sum(squared_ld,axis=0))


	gene_h2s = []

	f = open(gene_summary_file)
	head_count = 0
	n_genes = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_snp_indices = np.asarray(data[1].split(',')).astype(int)
		marginal_beta_file = data[6]
		marginal_betas = np.load(marginal_beta_file)
		sampling_variance = float(data[7])
		gene_h2s.append(float(data[5]))

		marginal_beta_transpose_beta = np.dot(marginal_betas.reshape(len(marginal_betas),1), marginal_betas.reshape(1,len(marginal_betas)))
		diags = np.diag(marginal_beta_transpose_beta)
		adj_diags = diags - ((1.0 - diags)/(eqtl_ss-2))
		np.fill_diagonal(marginal_beta_transpose_beta, adj_diags)

		gene_ld_scores = gene_ld_scores + marginal_beta_transpose_beta 

		n_genes = n_genes + 1
	f.close()


	model = LinearRegression(fit_intercept=False, normalize=False)
	z_mat = z_mat - LD

	joint_ld_scores = np.transpose(np.vstack((np.matrix.flatten(ld_scores),np.matrix.flatten(gene_ld_scores))))

	mpldsc_fit = model.fit(joint_ld_scores, np.matrix.flatten(z_mat))


	mpldsc_var_h2 = mpldsc_fit.coef_[0]*n_snps/N_gwas
	mpldsc_gene_h2 = mpldsc_fit.coef_[1]*n_genes*.025/N_gwas

	return mpldsc_var_h2, mpldsc_gene_h2




#########################
# Command line args
#########################
simulation_name_string = sys.argv[1]
simulated_gwas_data_dir = sys.argv[2]
simulated_gene_models_dir = sys.argv[3]
eqtl_ss = sys.argv[4]
mediated_h2_results_dir = sys.argv[5]
processed_genotype_data_dir = sys.argv[6]
print('start')


# Output file 
output_file = mediated_h2_results_dir + simulation_name_string + 'eqtl_ss_' + str(eqtl_ss) + '_h2_estimates.txt'
t = open(output_file,'w')
t.write('simulated_h2\tldsc_h2\tldsc_no_intercept_h2\tmpldsc_marginal_iid_eqtl_var_h2\tmpldsc_marginal_iid_eqtl_gene_h2\tmpldsc_marginal_eqtl_var_h2\tmpldsc_marginal_eqtl_gene_h2\tmpldsc_he_marginal_eqtl_var_h2\tmpldsc_he_marginal_eqtl_gene_h2\n')


# File summarizing inferred gene models
gene_summary_file = simulated_gene_models_dir + simulation_name_string + 'model_summaries_' + eqtl_ss + '.txt' 

# File containing gwas summary statistics
gwas_sum_stats_file = simulated_gwas_data_dir + simulation_name_string + 'simulated_gwas_summary_stats.txt'
gwas_z_scores = load_in_gwas_z_scores(gwas_sum_stats_file)
gwas_se = load_in_gwas_se(gwas_sum_stats_file)

# Simulated total genetic var
simulated_genetic_var_file = simulated_gwas_data_dir + simulation_name_string + 'simulated_genetic_var.npy'
sim_genetic_var = np.load(simulated_genetic_var_file) + 0.0



# LD file
ld_file = processed_genotype_data_dir + 'gwas_genotype_LD_1.npy'
LD = np.load(ld_file)

# eQTL LD file
eqtl_ld_file = processed_genotype_data_dir + 'eqtl_' + str(eqtl_ss) + '_genotype_LD_1.npy'
eqtl_LD = np.load(eqtl_ld_file)


'''
# LD_LD_T file
ld_ld_t_file = processed_genotype_data_dir + 'gwas_genotype_LD_LD_t_1.npy'
LD_LD_t = np.load(ld_ld_t_file)
'''

# Gwas sample size
N_gwas = 50000.0

# Quick LD Score regression fit
ldsc_h2 = run_ld_score_regression(gwas_z_scores, LD, N_gwas, intercept_bool=True)
ldsc_no_intercept_h2 = run_ld_score_regression(gwas_z_scores, LD, N_gwas, intercept_bool=False)

# MPLDSC regression with marginal_iid_eqtl_effect
mpldsc_marginal_iid_eqtl_var_h2, mpldsc_marginal_iid_eqtl_gene_h2  = run_mpldsc_with_marginal_iid_eqtl_effects(gwas_z_scores, LD, N_gwas,gene_summary_file, intercept_bool=True)


# MPLDSC regression with marginal_eqtl_effect
mpldsc_marginal_eqtl_var_h2, mpldsc_marginal_eqtl_gene_h2  = run_mpldsc_with_marginal_eqtl_effects(gwas_z_scores, LD, N_gwas,gene_summary_file, float(eqtl_ss), intercept_bool=True)



# MPLDSC regression with marginal_eqtl_effect using eqtl ld
mpldsc_marginal_eqtl_var_h2_eqtl, mpldsc_marginal_eqtl_gene_h2_eqtl  = run_mpldsc_with_marginal_eqtl_effects_using_eqtl_ld(gwas_z_scores, eqtl_LD, N_gwas,gene_summary_file, float(eqtl_ss), intercept_bool=True)




'''
# MPLDSC-he regression with marginal_eqtl_effect
mpldsc_he_marginal_eqtl_var_h2, mpldsc_he_marginal_eqtl_gene_h2  = run_mpldsc_he_with_marginal_eqtl_effects(gwas_z_scores, LD, LD_LD_t, N_gwas,gene_summary_file,float(eqtl_ss))

'''

t.write(str(sim_genetic_var) + '\t' + str(ldsc_h2) + '\t' + str(ldsc_no_intercept_h2) + '\t' + str(mpldsc_marginal_iid_eqtl_var_h2) + '\t' + str(mpldsc_marginal_iid_eqtl_gene_h2) + '\t' + str(mpldsc_marginal_eqtl_var_h2) + '\t' + str(mpldsc_marginal_eqtl_gene_h2) + '\t' + str(mpldsc_marginal_eqtl_var_h2_eqtl) + '\t' + str(mpldsc_marginal_eqtl_gene_h2_eqtl) + '\n')

t.close()

print(output_file)


# RSS VI variantly only
# Not working for whatever reason
#rss_vi_var_only = rss_vi_variant_only.RSS_VI_VARIANT_ONLY(LD=LD, z=gwas_z_scores, se=gwas_se)
#rss_vi_var_only.fit()



'''
# Just goes to show that two ways to compute gene LD scores are equivalent when using summary statistics
temp_R = LD[:100,:][:,:100]
temp_beta = np.random.normal(size=100)
sigma_2_p = .3
temp_R_inv = np.linalg.inv(temp_R)

mean = np.dot(np.dot(temp_R_inv, temp_beta).reshape(100,1), np.dot(temp_R_inv, temp_beta).reshape(1,100))
cov = sigma_2_p*temp_R_inv

value = np.sum((np.dot(np.transpose(temp_R[:1]),temp_R[:1]))*(mean + cov))


pdb.set_trace()
'''

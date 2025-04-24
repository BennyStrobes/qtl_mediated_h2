import numpy as np
import os
import sys
import pdb
import statsmodels.api as sm


def print_95_ci(arr):
	meany = np.mean(arr)
	se = np.std(arr)/np.sqrt(len(arr))
	ub = meany + (1.96*se)
	lb = meany - (1.96*se)
	print(str(meany) + ':   ' + '[' + str(lb) + ', ' + str(ub) + ']')
	return


def formal_print_95_ci(t, arr, method_name, sim_identifier):

	meany = np.mean(arr)
	se = np.std(arr)/np.sqrt(len(arr))
	ub = meany + (1.96*se)
	lb = meany - (1.96*se)


	t.write(method_name + '\t' + sim_identifier + '\t' + str(meany) + '\t' + str(lb) + '\t' + str(ub) + '\n')

	return t


#######################
# Input data
#######################
raw_ld_file = sys.argv[1]
output_dir = sys.argv[2]
n_causal_eqtl_snps_string = sys.argv[3]


###########################
# Simulation parameters
###########################
nm_h2 = 0.25
med_h2 = 0.05
eqtl_h2 = 0.05
alpha = med_h2/eqtl_h2
n_sims = 10000

n_causal_nm_snps = 200


############################
# Step 1: Load in LD data
############################
LD = np.loadtxt(raw_ld_file)
# tmp saving
#np.save('tmp_ld.npy', LD)
#LD = np.load('tmp_ld.npy')

padded_inv = np.linalg.inv(LD + (np.eye(LD.shape[0])*.01))
#np.save('padded_inv.npy', padded_inv)
#padded_inv = np.load('padded_inv.npy')


a_terms = np.diag(np.dot(LD, padded_inv))

n_snps = LD.shape[0]
n_window_eqtl_snps = int(np.math.floor(n_snps/2))
n_causal_eqtl_snps = 50
squared_LD = np.square(LD)

LD_t_LD = np.dot(LD, LD)
upper_triangle_indices = np.triu_indices(LD.shape[0], k=1)
low_ld_t_ld_indices_10 = np.abs(LD_t_LD[upper_triangle_indices]) < 10.0
low_ld_t_ld_indices_2 = np.abs(LD_t_LD[upper_triangle_indices]) < 2.0
low_ld_t_ld_indices_half = np.abs(LD_t_LD[upper_triangle_indices]) < .5



# LD Scores
ld_second_moments = np.sum(squared_LD,axis=0)
gene_window_ld_second_moments = np.sum(squared_LD[:n_window_eqtl_snps,:],axis=0)
total_possible_nm_snps = int(np.math.floor(n_snps/2)) - 2



output_file = output_dir + 'single_gene_sparse_eqtl_unshared_causal_snps_extra_nm_var_in_gene_' + str(n_causal_eqtl_snps_string) + '_eqtl_snps_v2.txt'
t = open(output_file,'w')
t.write('method\tsim_identifier\test\test_lb\test_ub\n')

#n_causal_eqtl_snps_arr = [1, 5, 50, 100]
n_causal_eqtl_snps_arr = [int(n_causal_eqtl_snps_string)]


for n_causal_eqtl_snps in n_causal_eqtl_snps_arr:

	################
	# Initialize vector to keep track of estimated parmams
	ldsc_h2_est = []
	mesc_alpha_sq_est = []
	mesc_alpha_sq_prime_est = []
	gene_window_mesc_alpha_sq_est = []
	gene_window_mesc_alpha_sq_prime_est = []
	gecs_alpha_sq_est = []
	gecs_alpha_sq_prime_est = []

	corr_causal_sq_effects = []
	corr_marginal_sq_effects = []
	sq_corr_causal_effects = []
	sq_corr_marginal_effects = []
	corr_cross_term_marginal_effects = []
	no_ld_mesc_alpha_sq_est = []
	no_ld_mesc_alpha_sq_prime_est = []


	for sim_iter in range(n_sims):
		print(sim_iter)
		############################
		# Step 2: simulate data
		############################
		# Simulate causal nm effects
		sim_causal_nm_effects = np.zeros(n_snps)


		causal_nm_indices1 = np.random.choice(np.arange(n_window_eqtl_snps, n_snps),int(n_causal_nm_snps*(5/25)), replace=False)
		causal_nm_indices2 = np.random.choice(np.arange(n_window_eqtl_snps),int(n_causal_nm_snps*(20/25)), replace=False)

		causal_nm_indices = np.hstack((causal_nm_indices1, causal_nm_indices2))
		sim_causal_nm_effects[causal_nm_indices] = np.random.normal(loc=0, scale=np.sqrt(nm_h2/n_causal_nm_snps), size=n_causal_nm_snps)


		# Simulate causal eqtl effects (deltas)
		sim_causal_eqtl_effects = np.zeros(n_snps)
		sim_alt_tissue_causal_eqtl_effects = np.zeros(n_snps)
		causal_eqtl_indices1 = np.random.choice(np.arange(n_window_eqtl_snps),n_causal_eqtl_snps, replace=False)
		causal_eqtl_indices2 = np.random.choice(np.arange(n_window_eqtl_snps),n_causal_eqtl_snps, replace=False)

		sim_causal_eqtl_effects[causal_eqtl_indices1] = np.random.normal(loc=0, scale=np.sqrt(eqtl_h2/n_causal_eqtl_snps), size=n_causal_eqtl_snps)
		sim_alt_tissue_causal_eqtl_effects[causal_eqtl_indices2] = np.random.normal(loc=0, scale=np.sqrt(eqtl_h2/n_causal_eqtl_snps), size=n_causal_eqtl_snps)


		sim_causal_med_effects = sim_causal_eqtl_effects*alpha
		# Calculate marginal effects
		beta_marginal = np.dot(LD, sim_causal_nm_effects + sim_causal_med_effects)
		delta_marginal = np.dot(LD, sim_causal_eqtl_effects)
		delta_prime_marginal = np.dot(LD, sim_alt_tissue_causal_eqtl_effects)

		beta_removed_ld = np.dot(padded_inv, beta_marginal)
		delta_removed_ld = np.dot(padded_inv, delta_marginal)
		delta_prime_removed_ld = np.dot(padded_inv, delta_prime_marginal)


		############################
		# Step 3: fit models
		############################	
		# LD score regression
		ldsc_model = sm.OLS(np.square(beta_marginal),ld_second_moments).fit()
		ldsc_h2_est.append(ldsc_model.params[0]*n_snps)

		# MESC regression
		X = np.transpose(np.vstack((ld_second_moments, np.square(delta_marginal))))
		mesc_model = sm.OLS(np.square(beta_marginal),X).fit()
		mesc_alpha_sq_est.append(mesc_model.params[1])

		# MESC regression (alt tissue)
		X = np.transpose(np.vstack((ld_second_moments, np.square(delta_prime_marginal))))
		mesc_alt_tiss_model = sm.OLS(np.square(beta_marginal),X).fit()
		mesc_alpha_sq_prime_est.append(mesc_alt_tiss_model.params[1])

		# MESC gene window regression
		X = np.transpose(np.vstack((ld_second_moments, gene_window_ld_second_moments, np.square(delta_marginal))))
		mesc_model = sm.OLS(np.square(beta_marginal),X).fit()
		gene_window_mesc_alpha_sq_est.append(mesc_model.params[2])	

		# MESC-prime gene window regression
		X = np.transpose(np.vstack((ld_second_moments, gene_window_ld_second_moments, np.square(delta_prime_marginal))))
		mesc_model = sm.OLS(np.square(beta_marginal),X).fit()
		gene_window_mesc_alpha_sq_prime_est.append(mesc_model.params[2])	


		# No LD MESC regression
		X = np.transpose(np.vstack((a_terms, np.square(delta_removed_ld))))
		mesc_model = sm.OLS(np.square(beta_removed_ld),X).fit()
		no_ld_mesc_alpha_sq_est.append(mesc_model.params[1])

		# No LD MESC regression (alt tissue)
		X = np.transpose(np.vstack((a_terms, np.square(delta_prime_removed_ld))))
		mesc_model = sm.OLS(np.square(beta_removed_ld),X).fit()
		no_ld_mesc_alpha_sq_prime_est.append(mesc_model.params[1])


		'''
		low_ld_t_ld_indices_10 = np.abs(LD_t_LD[upper_triangle_indices]) < 10.0
		low_ld_t_ld_indices_2 = np.abs(LD_t_LD[upper_triangle_indices]) < 2.0
		low_ld_t_ld_indices_half = np.abs(LD_t_LD[upper_triangle_indices]) < .5

		'''



		# low_ld_t_ld_indices_2
		# Version of GE co-score regression
		beta_beta_transpose = np.dot(beta_marginal.reshape(-1,1), beta_marginal.reshape(1,-1))
		delta_delta_transpose = np.dot(delta_marginal.reshape(-1,1), delta_marginal.reshape(1,-1))
		Y = beta_beta_transpose[upper_triangle_indices][low_ld_t_ld_indices_2]
		cross_ld_scores = LD_t_LD[upper_triangle_indices][low_ld_t_ld_indices_2]
		cross_delta_terms = delta_delta_transpose[upper_triangle_indices][low_ld_t_ld_indices_2]
		X = np.transpose(np.vstack((cross_ld_scores, cross_delta_terms)))
		gecs_model = sm.OLS(Y,X).fit()
		gecs_alpha_sq_est.append(gecs_model.params[1])

		# Alt tissue
		delta_prime_delta_prime_transpose = np.dot(delta_prime_marginal.reshape(-1,1), delta_prime_marginal.reshape(1,-1))
		cross_delta_prime_terms = delta_prime_delta_prime_transpose[upper_triangle_indices][low_ld_t_ld_indices_2]
		X = np.transpose(np.vstack((cross_ld_scores, cross_delta_prime_terms)))
		gecs_model = sm.OLS(Y,X).fit()
		gecs_alpha_sq_prime_est.append(gecs_model.params[1])






		corr_causal_sq_effects.append(np.corrcoef(np.square(sim_causal_eqtl_effects), np.square(sim_alt_tissue_causal_eqtl_effects))[0,1])
		corr_marginal_sq_effects.append(np.corrcoef(np.square(delta_prime_marginal), np.square(delta_marginal))[0,1])


		sq_corr_causal_effects.append(np.square(np.corrcoef(sim_causal_eqtl_effects, sim_alt_tissue_causal_eqtl_effects)[0,1]))
		sq_corr_marginal_effects.append(np.square(np.corrcoef(delta_prime_marginal, delta_marginal)[0,1]))

		corr_cross_term_marginal_effects.append(np.corrcoef(cross_delta_terms, cross_delta_prime_terms)[0,1])

	# Organize results
	ldsc_h2_est = np.asarray(ldsc_h2_est)
	t = formal_print_95_ci(t, ldsc_h2_est, 'ldsc', str(n_causal_eqtl_snps))
	mesc_alpha_sq_est = np.asarray(mesc_alpha_sq_est)
	t = formal_print_95_ci(t, mesc_alpha_sq_est, 'mesc', str(n_causal_eqtl_snps))
	mesc_alpha_sq_prime_est = np.asarray(mesc_alpha_sq_prime_est)
	t = formal_print_95_ci(t, mesc_alpha_sq_prime_est, 'mesc_prime', str(n_causal_eqtl_snps))

	gene_window_mesc_alpha_sq_est = np.asarray(gene_window_mesc_alpha_sq_est)
	t = formal_print_95_ci(t, gene_window_mesc_alpha_sq_est, 'mesc_w_gene_window', str(n_causal_eqtl_snps))
	gene_window_mesc_alpha_sq_prime_est = np.asarray(gene_window_mesc_alpha_sq_prime_est)
	t = formal_print_95_ci(t, gene_window_mesc_alpha_sq_prime_est, 'mesc_prime_w_gene_window', str(n_causal_eqtl_snps))

	no_ld_mesc_alpha_sq_est = np.asarray(no_ld_mesc_alpha_sq_est)
	t = formal_print_95_ci(t, no_ld_mesc_alpha_sq_est, 'mesc_removed_ld', str(n_causal_eqtl_snps))
	no_ld_mesc_alpha_sq_prime_est = np.asarray(no_ld_mesc_alpha_sq_prime_est)
	t = formal_print_95_ci(t, no_ld_mesc_alpha_sq_prime_est, 'mesc_prime_removed_ld', str(n_causal_eqtl_snps))

	gecs_alpha_sq_est = np.asarray(gecs_alpha_sq_est)
	t = formal_print_95_ci(t, gecs_alpha_sq_est, 'gecs', str(n_causal_eqtl_snps))
	gecs_alpha_sq_prime_est = np.asarray(gecs_alpha_sq_prime_est)
	t = formal_print_95_ci(t, gecs_alpha_sq_prime_est, 'gecs_prime', str(n_causal_eqtl_snps))


	corr_cross_term_marginal_effects = np.asarray(corr_cross_term_marginal_effects)
	t = formal_print_95_ci(t, corr_cross_term_marginal_effects, 'corr_cross_terms_marginal_effects', str(n_causal_eqtl_snps))

	corr_causal_sq_effects = np.asarray(corr_causal_sq_effects)
	t = formal_print_95_ci(t, corr_causal_sq_effects, 'corr_causal_sq_effects', str(n_causal_eqtl_snps))
	corr_marginal_sq_effects = np.asarray(corr_marginal_sq_effects)
	t = formal_print_95_ci(t, corr_marginal_sq_effects, 'corr_marginal_sq_effects', str(n_causal_eqtl_snps))

	sq_corr_causal_effects = np.asarray(sq_corr_causal_effects)
	t = formal_print_95_ci(t, sq_corr_causal_effects, 'sq_corr_causal_effects', str(n_causal_eqtl_snps))
	sq_corr_marginal_effects = np.asarray(sq_corr_marginal_effects)
	t = formal_print_95_ci(t, sq_corr_marginal_effects, 'sq_corr_marginal_effects', str(n_causal_eqtl_snps))




t.close()

print(output_file)

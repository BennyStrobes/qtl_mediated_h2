import numpy as np
import os
import sys
import pdb
import statsmodels.api as sm
import time

def print_95_ci(arr):
	meany = np.mean(arr)
	se = np.std(arr)/np.sqrt(len(arr))
	ub = meany + (1.96*se)
	lb = meany - (1.96*se)

	#print(str(meany) + ':   ' + '[' + str(lb) + ', ' + str(ub) + ']')
	return meany, lb, ub

def simulate_data(ld_mat, gwas_ss, eqtl_ss, nm_h2, med_h2, eqtl_h2, sim_sparsity):
	n_snps = ld_mat.shape[0]
	delta = np.zeros(n_snps)
	sub_n_snps = int(n_snps/2)
	delta[:sub_n_snps] = np.random.normal(loc=0, scale=np.sqrt(eqtl_h2/sub_n_snps),size=sub_n_snps)


	gamma = np.random.normal(loc=0, scale=np.sqrt(nm_h2/n_snps),size=n_snps)


	alpha_sq = med_h2/eqtl_h2
	beta = gamma + (delta*np.sqrt(alpha_sq))




	beta_hat = np.random.normal(np.dot(ld_mat, beta), scale=np.sqrt(1.0/gwas_ss))
	#beta_hat = np.random.multivariate_normal(mean=np.dot(ld_mat, beta), cov=ld_mat/gwas_ss)

	delta_hat = np.random.normal(np.dot(ld_mat, delta), scale=np.sqrt(1.0/eqtl_ss))

	delta_hat1 = np.random.normal(np.dot(ld_mat, delta), scale=np.sqrt(1.0/(eqtl_ss/2)))
	delta_hat2 = np.random.normal(np.dot(ld_mat, delta), scale=np.sqrt(1.0/(eqtl_ss/2)))


	return beta_hat, delta_hat, np.dot(ld_mat, beta), np.dot(ld_mat, delta), delta_hat1, delta_hat2



def run_standard_ldsc(ldscores, beta_hat, sample_size):
	#ldscores = np.sum(np.square(ld_mat), axis=0)

	y_terms = np.square(beta_hat) - (1.0/sample_size)

	modeler = sm.OLS(y_terms, ldscores).fit()


	n_snps = ldscores.shape[0]
	est_h2 = modeler.params[0]*n_snps


	return est_h2

def run_two_step_ldsc_binned(ldscores, binned_ld_scores, snps_per_bin, beta_hat, gwas_ss, delta_hat, eqtl_ss, eqtl_h2):
	y_terms = np.square(delta_hat) - (1.0/eqtl_ss)
	modeler = sm.OLS(y_terms, binned_ld_scores).fit()

	est_delta_sq = modeler.params
	y_terms = np.square(beta_hat) - (1.0/gwas_ss)

	x_terms = np.transpose(np.vstack((ldscores, np.dot(binned_ld_scores, modeler.params))))

	modeler = sm.OLS(y_terms, x_terms).fit()

	n_snps = ldscores.shape[0]
	est_nm_h2 = modeler.params[0]*n_snps

	est_eqtl_h2 = np.sum(est_delta_sq*snps_per_bin)
	est_med_h2 = modeler.params[1]*est_eqtl_h2

	return est_nm_h2, est_med_h2

def run_two_step_corrected_ldsc_binned(ldscores, binned_ld_scores, snps_per_bin, beta_hat, gwas_ss, delta_hat, eqtl_ss, eqtl_h2, delta_sub1_hat, delta_sub2_hat):
	y_terms = np.square(delta_hat) - (1.0/eqtl_ss)
	modeler = sm.OLS(y_terms, binned_ld_scores).fit()
	est_delta_sq = modeler.params

	y_terms1 = np.square(delta_sub1_hat) - (1.0/(eqtl_ss/2))
	modeler = sm.OLS(y_terms1, binned_ld_scores).fit()
	est_delta_sq1 = modeler.params

	y_terms2 = np.square(delta_sub2_hat) - (1.0/(eqtl_ss/2))
	modeler = sm.OLS(y_terms2, binned_ld_scores).fit()
	est_delta_sq2 = modeler.params


	gene_ldscores = np.dot(binned_ld_scores, est_delta_sq)
	gene_ldscores1 = np.dot(binned_ld_scores, est_delta_sq1)
	gene_ldscores2 = np.dot(binned_ld_scores, est_delta_sq2)


	y_terms = np.square(beta_hat) - (1.0/gwas_ss)
	x_terms = np.transpose(np.vstack((ldscores, gene_ldscores)))
	x_terms1 = np.transpose(np.vstack((ldscores, gene_ldscores1)))
	x_terms2 = np.transpose(np.vstack((ldscores, gene_ldscores2)))
	params = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(x_terms1), x_terms2)), np.transpose(x_terms)), y_terms)


	#x_terms = np.transpose(np.vstack((ldscores, np.square(delta_hat) - (1.0/eqtl_ss))))
	#x_terms1 = np.transpose(np.vstack((ldscores, np.square(delta_sub1_hat) - (1.0/(eqtl_ss/2)))))
	#x_terms2 = np.transpose(np.vstack((ldscores, np.square(delta_sub2_hat) - (1.0/(eqtl_ss/2)))))
	#params = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(x_terms1), x_terms2)), np.transpose(x_terms)), y_terms)




	#tmp = sm.OLS(gene_ldscores1, sm.add_constant(x_terms2)).fit()
	#pred_gene_ldscores = tmp.params[0] + (tmp.params[1]*ldscores) + (tmp.params[2]*gene_ldscores2)
	#x_terms_tmp = np.transpose(np.vstack((ldscores, pred_gene_ldscores)))

	#x_terms3 = np.transpose(np.vstack((ldscores, y_terms2)))
	#tmp = sm.OLS(y_terms1, sm.add_constant(x_terms3)).fit()
	#pred_gene_ldscores = tmp.params[0] + (tmp.params[1]*ldscores) + (tmp.params[2]*y_terms2)
	#x_terms_tmp = np.transpose(np.vstack((ldscores, pred_gene_ldscores)))

	#params = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(x_terms_tmp), x_terms_tmp)), np.transpose(x_terms_tmp)), y_terms)


	n_snps = ldscores.shape[0]
	est_nm_h2 = params[0]*n_snps

	est_eqtl_h2 = np.sum(est_delta_sq*snps_per_bin)
	est_med_h2 = params[1]*.05

	return est_nm_h2, est_med_h2, est_eqtl_h2

def run_joint_ldsc_binned_true_resid(ldscores, binned_ld_scores, snps_per_bin, beta_hat, gwas_ss, delta_hat, eqtl_ss, eqtl_h2, max_iter=10000, conv_thresh=1e-12):

	eqtl_qq = np.square(delta_tilde)
	gwas_qq = np.square(beta_tilde)


	# Initialize parameters
	y_terms = np.square(delta_hat) - (1.0/eqtl_ss)
	modeler = sm.OLS(y_terms, binned_ld_scores).fit()
	est_delta_sq = modeler.params

	est_eqtl_residual_var = np.sum(np.square(y_terms-np.dot(binned_ld_scores, modeler.params)))/len(y_terms)
	est_eqtl_residual_var = np.mean(np.square(eqtl_qq - np.dot(binned_ld_scores, modeler.params)))

	y_terms = np.square(beta_hat) - (1.0/gwas_ss)
	x_terms = np.transpose(np.vstack((ldscores, np.dot(binned_ld_scores, est_delta_sq))))

	modeler = sm.OLS(y_terms, x_terms).fit()
	est_gamma_sq = modeler.params[0]
	est_alpha_sq = modeler.params[1]

	n_snps = ldscores.shape[0]
	est_nm_h2 = est_gamma_sq*n_snps
	est_med_h2 = est_alpha_sq*eqtl_h2

	est_gwas_residual_var = np.sum(np.square(y_terms-np.dot(x_terms, modeler.params)))/len(y_terms)
	est_gwas_residual_var = np.mean(np.square(gwas_qq - np.dot(x_terms, modeler.params)))

	old_nm_h2 = np.copy(est_nm_h2)

	for itera in range(max_iter):

		################################
		# Joint update of delta-sq
		cur_gwas_resid = np.square(beta_hat) - (1.0/gwas_ss) - ldscores*est_gamma_sq
		cur_eqtl_resid = np.square(delta_hat) - (1.0/eqtl_ss)

		big_resid = np.hstack((cur_gwas_resid, cur_eqtl_resid))
		big_ldscore = np.vstack((binned_ld_scores*est_alpha_sq, binned_ld_scores))
		gwas_weights = np.ones(len(cur_gwas_resid))*(1.0/est_gwas_residual_var)
		eqtl_weights = np.ones(len(cur_eqtl_resid))*(1.0/est_eqtl_residual_var)
		big_weights = np.hstack((gwas_weights, eqtl_weights))

		modeler = sm.WLS(big_resid, big_ldscore, weights=big_weights).fit()
		est_delta_sq = modeler.params

		est_eqtl_residual_var = np.sum(np.square(cur_eqtl_resid-np.dot(binned_ld_scores, est_delta_sq)))/len(cur_eqtl_resid)
		tmp = np.square(eqtl_qq - np.dot(binned_ld_scores, est_delta_sq))
		est_eqtl_residual_var = np.mean(tmp[:int(n_snps/2)])

		###############################
		# Update alpha_sq and gamma_sq
		cur_gwas_resid = np.square(beta_hat) - (1.0/gwas_ss)
		x_terms = np.transpose(np.vstack((ldscores, np.dot(binned_ld_scores, est_delta_sq))))
		modeler = sm.OLS(cur_gwas_resid, x_terms).fit()
		est_gamma_sq = modeler.params[0]
		est_alpha_sq = modeler.params[1]

		est_gwas_residual_var = np.sum(np.square(cur_gwas_resid-np.dot(x_terms, modeler.params)))/len(cur_gwas_resid)

		est_gwas_residual_var = np.mean(np.square(gwas_qq-np.dot(x_terms, modeler.params)))

		n_snps = ldscores.shape[0]
		est_nm_h2 = est_gamma_sq*n_snps
		est_eqtl_h2 = np.sum(est_delta_sq*snps_per_bin)
		est_med_h2 = est_alpha_sq*est_eqtl_h2


		diff = np.abs(est_nm_h2 - old_nm_h2)
		#print('############################')
		#print(diff)
		#print(est_nm_h2)
		old_nm_h2 = np.copy(est_nm_h2)
		if diff < conv_thresh:
			break

	print(itera)
	return est_nm_h2, est_med_h2, est_eqtl_h2

def run_joint_ldsc_binned(ldscores, binned_ld_scores, snps_per_bin, beta_hat, gwas_ss, delta_hat, eqtl_ss, eqtl_h2, max_iter=20000, conv_thresh=1e-8):

	# Initialize parameters
	y_terms = np.square(delta_hat) - (1.0/eqtl_ss)
	modeler = sm.OLS(y_terms, binned_ld_scores).fit()
	est_delta_sq = modeler.params

	est_eqtl_residual_var = np.sum(np.square(y_terms-np.dot(binned_ld_scores, modeler.params)))/len(y_terms)

	y_terms = np.square(beta_hat) - (1.0/gwas_ss)
	x_terms = np.transpose(np.vstack((ldscores, np.dot(binned_ld_scores, est_delta_sq))))

	modeler = sm.OLS(y_terms, x_terms).fit()
	est_gamma_sq = modeler.params[0]
	est_alpha_sq = modeler.params[1]

	n_snps = ldscores.shape[0]
	est_nm_h2 = est_gamma_sq*n_snps
	est_med_h2 = est_alpha_sq*eqtl_h2

	est_gwas_residual_var = np.sum(np.square(y_terms-np.dot(x_terms, modeler.params)))/len(y_terms)

	old_nm_h2 = np.copy(est_nm_h2)

	for itera in range(max_iter):

		################################
		# Joint update of delta-sq
		cur_gwas_resid = np.square(beta_hat) - (1.0/gwas_ss) - ldscores*est_gamma_sq
		cur_eqtl_resid = np.square(delta_hat) - (1.0/eqtl_ss)

		big_resid = np.hstack((cur_gwas_resid, cur_eqtl_resid))
		big_ldscore = np.vstack((binned_ld_scores*est_alpha_sq, binned_ld_scores))
		gwas_weights = np.ones(len(cur_gwas_resid))*(1.0/est_gwas_residual_var)
		eqtl_weights = np.ones(len(cur_eqtl_resid))*(1.0/est_eqtl_residual_var)
		big_weights = np.hstack((gwas_weights, eqtl_weights))

		modeler = sm.WLS(big_resid, big_ldscore, weights=big_weights).fit()
		est_delta_sq = modeler.params

		est_eqtl_residual_var = np.sum(np.square(cur_eqtl_resid-np.dot(binned_ld_scores, est_delta_sq)))/len(cur_eqtl_resid)

		###############################
		# Update alpha_sq and gamma_sq
		cur_gwas_resid = np.square(beta_hat) - (1.0/gwas_ss)
		x_terms = np.transpose(np.vstack((ldscores, np.dot(binned_ld_scores, est_delta_sq))))
		modeler = sm.OLS(cur_gwas_resid, x_terms).fit()
		est_gamma_sq = modeler.params[0]
		est_alpha_sq = modeler.params[1]

		est_gwas_residual_var = np.sum(np.square(cur_gwas_resid-np.dot(x_terms, modeler.params)))/len(cur_gwas_resid)


		n_snps = ldscores.shape[0]
		est_nm_h2 = est_gamma_sq*n_snps
		est_med_h2 = est_alpha_sq*eqtl_h2

		print('############################')
		diff = np.abs(est_nm_h2 - old_nm_h2)
		print(diff)
		#print(est_nm_h2)
		#print(est_gwas_residual_var)
		#print(est_eqtl_residual_var)
		old_nm_h2 = np.copy(est_nm_h2)
		if diff < conv_thresh:
			break
	print(itera)


	return est_nm_h2, est_med_h2



####################
# Command line args
####################
ld_dir = sys.argv[1]
output_dir = sys.argv[2]
gwas_ss = float(sys.argv[3])
eqtl_ss = float(sys.argv[4])
nm_h2 = float(sys.argv[5])
med_h2 = float(sys.argv[6])
eqtl_h2 = float(sys.argv[7])
sim_sparsity = sys.argv[8]


# Load in LD
ld_mat = np.load(ld_dir + 'ld.npy')
ldscores = np.sum(np.square(ld_mat), axis=0)

# Load in LD data based on bins
n_bins=30
n_snps = ld_mat.shape[0]
bins = np.array_split(np.arange(int(n_snps/2)), n_bins)
binned_ld_scores = []
snps_per_bin = []
for ii, biner in enumerate(bins):
	binned_ld_scores.append(np.sum((np.square(ld_mat))[:, biner],axis=1))
	snps_per_bin.append(len(biner))
binned_ld_scores = np.transpose(np.asarray(binned_ld_scores))
snps_per_bin = np.asarray(snps_per_bin)


gwas_h2_ests = []
eqtl_h2_ests = []

two_step_ldsc_nm_h2_ests = []
two_step_ldsc_med_h2_ests = []
joint_ldsc_nm_h2_ests = []
joint_ldsc_med_h2_ests = []
joint_ldsc_eqtl_h2_ests = []
corrected_ldsc_nm_h2_ests = []
corrected_ldsc_med_h2_ests = []
corrected_ldsc_eqtl_h2_ests = []


output_file = output_dir + 'simple_sim_' + str(gwas_ss) + '_' + str(eqtl_ss) + '_' + str(nm_h2) + '_' + str(med_h2) + '_' + str(eqtl_h2) + '_' + str(sim_sparsity) + '_results.txt'
t = open(output_file,'w')
t.write('param\ttruth\test\test_lb\test_ub\n')


for sim_iter in range(500):
	print(sim_iter)
	# Simulate data
	beta_hat, delta_hat, beta_tilde, delta_tilde, delta_sub1_hat, delta_sub2_hat = simulate_data(ld_mat, gwas_ss, eqtl_ss, nm_h2, med_h2, eqtl_h2, sim_sparsity)


	######################################
	# Run ldscore regression on GWAS data
	est_gwas_h2 = run_standard_ldsc(ldscores, beta_hat, gwas_ss)
	gwas_h2_ests.append(est_gwas_h2)

	######################################
	# Run ldscore regression on eQTL data
	est_eqtl_h2 = run_standard_ldsc(ldscores, delta_hat, eqtl_ss)
	eqtl_h2_ests.append(est_eqtl_h2)

	######################################
	two_step_ldsc_nm_h2_est, two_step_ldsc_med_h2_est = run_two_step_ldsc_binned(ldscores, binned_ld_scores, snps_per_bin, beta_hat, gwas_ss, delta_hat, eqtl_ss,eqtl_h2)
	two_step_ldsc_nm_h2_ests.append(two_step_ldsc_nm_h2_est)
	two_step_ldsc_med_h2_ests.append(two_step_ldsc_med_h2_est)


	corrected_ldsc_nm_h2_est, corrected_ldsc_med_h2_est, corrected_ldsc_eqtl_h2_est = run_two_step_corrected_ldsc_binned(ldscores, binned_ld_scores, snps_per_bin, beta_hat, gwas_ss, delta_hat, eqtl_ss,eqtl_h2, delta_sub1_hat, delta_sub2_hat)
	corrected_ldsc_nm_h2_ests.append(corrected_ldsc_nm_h2_est)
	corrected_ldsc_med_h2_ests.append(corrected_ldsc_med_h2_est)
	corrected_ldsc_eqtl_h2_ests.append(corrected_ldsc_eqtl_h2_est)

	######################################
	#two_step_ldsc_nm_h2_est, two_step_ldsc_med_h2_est = run_joint_ldsc_binned(ldscores, binned_ld_scores, snps_per_bin, beta_hat, gwas_ss, delta_hat, eqtl_ss,eqtl_h2)
	'''
	joint_ldsc_nm_h2_est, joint_ldsc_med_h2_est, joint_ldsc_eqtl_h2_est = run_joint_ldsc_binned_true_resid(ldscores, binned_ld_scores, snps_per_bin, beta_hat, gwas_ss, delta_hat, eqtl_ss,eqtl_h2)
	joint_ldsc_nm_h2_ests.append(joint_ldsc_nm_h2_est)
	joint_ldsc_med_h2_ests.append(joint_ldsc_med_h2_est)
	joint_ldsc_eqtl_h2_ests.append(joint_ldsc_eqtl_h2_est)
	'''

gwas_h2_ests = np.asarray(gwas_h2_ests)
eqtl_h2_ests = np.asarray(eqtl_h2_ests)
two_step_ldsc_nm_h2_ests = np.asarray(two_step_ldsc_nm_h2_ests)
two_step_ldsc_med_h2_ests = np.asarray(two_step_ldsc_med_h2_ests)
joint_ldsc_nm_h2_ests = np.asarray(joint_ldsc_nm_h2_ests)
joint_ldsc_med_h2_ests = np.asarray(joint_ldsc_med_h2_ests)
joint_ldsc_eqtl_h2_ests = np.asarray(joint_ldsc_eqtl_h2_ests)
corrected_ldsc_nm_h2_ests = np.asarray(corrected_ldsc_nm_h2_ests)
corrected_ldsc_med_h2_ests = np.asarray(corrected_ldsc_med_h2_ests)
corrected_ldsc_eqtl_h2_ests = np.asarray(corrected_ldsc_eqtl_h2_ests)


gwas_h2_est_mean, gwas_h2_est_lb, gwas_h2_est_ub = print_95_ci(gwas_h2_ests)
t.write('gwas_ldsc_h2\t' + str(med_h2+nm_h2) + '\t' + str(gwas_h2_est_mean) + '\t' + str(gwas_h2_est_lb) + '\t' + str(gwas_h2_est_ub) + '\n')
eqtl_h2_est_mean, eqtl_h2_est_lb, eqtl_h2_est_ub = print_95_ci(eqtl_h2_ests)
t.write('eqtl_ldsc_h2\t' + str(eqtl_h2) + '\t' + str(eqtl_h2_est_mean) + '\t' + str(eqtl_h2_est_lb) + '\t' + str(eqtl_h2_est_ub) + '\n')
two_step_ldsc_nm_h2_est_mean, two_step_ldsc_nm_h2_est_lb, two_step_ldsc_nm_h2_est_ub = print_95_ci(two_step_ldsc_nm_h2_ests)
t.write('two_step_ldsc_nm_h2\t' + str(nm_h2) + '\t' + str(two_step_ldsc_nm_h2_est_mean) + '\t' + str(two_step_ldsc_nm_h2_est_lb) + '\t' + str(two_step_ldsc_nm_h2_est_ub) + '\n')
two_step_ldsc_med_h2_est_mean, two_step_ldsc_med_h2_est_lb, two_step_ldsc_med_h2_est_ub = print_95_ci(two_step_ldsc_med_h2_ests)
t.write('two_step_ldsc_med_h2\t' + str(med_h2) + '\t' + str(two_step_ldsc_med_h2_est_mean) + '\t' + str(two_step_ldsc_med_h2_est_lb) + '\t' + str(two_step_ldsc_med_h2_est_ub) + '\n')


joint_ldsc_nm_h2_est_mean, joint_ldsc_nm_h2_est_lb, joint_ldsc_nm_h2_est_ub = print_95_ci(corrected_ldsc_nm_h2_ests)
t.write('two_step_corrected_ldsc_nm_h2\t' + str(nm_h2) + '\t' + str(joint_ldsc_nm_h2_est_mean) + '\t' + str(joint_ldsc_nm_h2_est_lb) + '\t' + str(joint_ldsc_nm_h2_est_ub) + '\n')
joint_ldsc_med_h2_est_mean, joint_ldsc_med_h2_est_lb, joint_ldsc_med_h2_est_ub = print_95_ci(corrected_ldsc_med_h2_ests)
t.write('two_step_corrected_ldsc_med_h2\t' + str(med_h2) + '\t' + str(joint_ldsc_med_h2_est_mean) + '\t' + str(joint_ldsc_med_h2_est_lb) + '\t' + str(joint_ldsc_med_h2_est_ub) + '\n')
joint_ldsc_eqtl_h2_est_mean, joint_ldsc_eqtl_h2_est_lb, joint_ldsc_eqtl_h2_est_ub = print_95_ci(corrected_ldsc_eqtl_h2_ests)
t.write('two_step_corrected_ldsc_eqtl_h2\t' + str(eqtl_h2) + '\t' + str(joint_ldsc_eqtl_h2_est_mean) + '\t' + str(joint_ldsc_eqtl_h2_est_lb) + '\t' + str(joint_ldsc_eqtl_h2_est_ub) + '\n')

'''
joint_ldsc_nm_h2_est_mean, joint_ldsc_nm_h2_est_lb, joint_ldsc_nm_h2_est_ub = print_95_ci(joint_ldsc_nm_h2_ests)
t.write('joint_ldsc_nm_h2\t' + str(nm_h2) + '\t' + str(joint_ldsc_nm_h2_est_mean) + '\t' + str(joint_ldsc_nm_h2_est_lb) + '\t' + str(joint_ldsc_nm_h2_est_ub) + '\n')
joint_ldsc_med_h2_est_mean, joint_ldsc_med_h2_est_lb, joint_ldsc_med_h2_est_ub = print_95_ci(joint_ldsc_med_h2_ests)
t.write('joint_ldsc_med_h2\t' + str(med_h2) + '\t' + str(joint_ldsc_med_h2_est_mean) + '\t' + str(joint_ldsc_med_h2_est_lb) + '\t' + str(joint_ldsc_med_h2_est_ub) + '\n')
joint_ldsc_eqtl_h2_est_mean, joint_ldsc_eqtl_h2_est_lb, joint_ldsc_eqtl_h2_est_ub = print_95_ci(joint_ldsc_eqtl_h2_ests)
t.write('joint_ldsc_eqtl_h2\t' + str(eqtl_h2) + '\t' + str(joint_ldsc_eqtl_h2_est_mean) + '\t' + str(joint_ldsc_eqtl_h2_est_lb) + '\t' + str(joint_ldsc_eqtl_h2_est_ub) + '\n')
'''

t.close()

print(output_file)


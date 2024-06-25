import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np
import os
import pdb
import statsmodels.api as sm
import joint_sldsc_no_ld


def posterior_distribution_on_causal_eqtl_effects_closed_form(expression_vec, genotype_mat, ge_h2, resid_var):
	n_snps = genotype_mat.shape[1]
	per_snp_h2 = ge_h2/n_snps

	tmp_prec = (np.eye(n_snps)/per_snp_h2) + (np.dot(np.transpose(genotype_mat), genotype_mat)/(resid_var))
	posterior_var = np.linalg.inv(tmp_prec)
	posterior_mean = (1.0/resid_var)*np.dot(np.dot(posterior_var, np.transpose(genotype_mat)), expression_vec)

	return posterior_mean, posterior_var

def posterior_distribution_on_causal_eqtl_effects_closed_form_iid(expression_vec, genotype_mat, ge_h2, resid_var):
	n_snps = genotype_mat.shape[1]
	per_snp_h2 = ge_h2/n_snps
	N_eqtl = genotype_mat.shape[0]

	tmp_prec = (1.0/per_snp_h2) + (N_eqtl/(resid_var))
	posterior_var = 1.0/tmp_prec
	posterior_mean = (1.0/resid_var)*posterior_var*np.dot(np.transpose(genotype_mat), expression_vec)

	return posterior_mean, posterior_var







#####################
# Command line args
#####################
gwas_ss = int(sys.argv[1])
n_snps = int(sys.argv[2])
nm_h2 = float(sys.argv[3])
n_genes = int(sys.argv[4])
med_h2 = float(sys.argv[5])
snps_per_gene = int(sys.argv[6])
ge_h2 = float(sys.argv[7])
eqtl_ss = int(sys.argv[8])
sim_iter = int(sys.argv[9])
output_dir = sys.argv[10]

avg_gene_h2s = []
med_h2s = []
nm_h2s = []

for sim_iter in range(100):
	print(sim_iter)
	#####################
	# SIMULATION
	#####################
	# Simulate causal variant- trait effects
	sim_gamma = np.random.normal(loc=0,scale=np.sqrt(nm_h2/n_snps),size=n_snps)
	# Simulate causal gene-trait effects
	sim_alpha = np.random.normal(loc=0,scale=np.sqrt(med_h2/(ge_h2*n_genes)),size=n_genes)
	# Simulate causal variant_gene effects
	gene_info = {}
	snp_index = 0
	for gene_iter in range(n_genes):
		gene_snps = np.arange(snp_index, snp_index + snps_per_gene)
		gene_causal_eqtl_effects = np.random.normal(loc=0, scale=np.sqrt(ge_h2/snps_per_gene), size=snps_per_gene)
		gene_info['gene' + str(gene_iter)] = {}
		gene_info['gene' + str(gene_iter)]['indices'] = gene_snps
		gene_info['gene' + str(gene_iter)]['causal_eqtl_effects'] = gene_causal_eqtl_effects
		snp_index = snp_index + snps_per_gene

	# Estimate causal eqtl effects
	passer = True
	aa = []
	for gene_iter in range(n_genes):
		gene_causal_eqtl_effects = gene_info['gene' + str(gene_iter)]['causal_eqtl_effects']
		specific_n_snps = len(gene_causal_eqtl_effects)
		# Simulate eqtl genotype data
		eqtl_geno = np.random.normal(loc=0,scale=1.0,size=(eqtl_ss,specific_n_snps))
		for snp_iter in range(specific_n_snps):
			eqtl_geno[:,snp_iter] = (eqtl_geno[:,snp_iter] - np.mean(eqtl_geno[:,snp_iter]))/np.std(eqtl_geno[:,snp_iter])
		genetic_gene = np.dot(eqtl_geno, gene_causal_eqtl_effects)
		if np.var(genetic_gene,ddof=1) > 1:
			passer = False
			continue
		aa.append(np.var(genetic_gene,ddof=1))
		EE = np.random.normal(genetic_gene, scale=np.sqrt(1.0-np.var(genetic_gene,ddof=1)))
		posterior_mean, posterior_var = posterior_distribution_on_causal_eqtl_effects_closed_form(EE, eqtl_geno, .05, .95)
		gene_info['gene' + str(gene_iter)]['init'] = np.square(posterior_mean) + np.diagonal(posterior_var)
		# Get gwas summary stats
		eqtl_chi_sq = []
		for snp_iter in range(specific_n_snps):
			olser = sm.OLS(EE, sm.add_constant(eqtl_geno[:,snp_iter])).fit()
			beta = olser.params[1]
			beta_se = olser.bse[1]
			zed = beta/beta_se
			eqtl_chi_sq.append(np.square(zed))
		eqtl_chi_sq = np.asarray(eqtl_chi_sq)
		gene_info['gene' + str(gene_iter)]['chi_sq'] = eqtl_chi_sq
		gene_info['gene' + str(gene_iter)]['ss'] = len(EE)
	#print(np.mean(aa))
	#print('####')

	if passer == False:
		continue

	# Simulate GWAS genotype data
	gwas_geno = np.random.normal(loc=0,scale=1.0,size=(gwas_ss,n_snps))
	for snp_iter in range(n_snps):
		gwas_geno[:,snp_iter] = (gwas_geno[:,snp_iter] - np.mean(gwas_geno[:,snp_iter]))/np.std(gwas_geno[:,snp_iter])

	# Simulate GWAS trait value
	genetic_trait = np.dot(gwas_geno, sim_gamma)
	for gene_iter in range(n_genes):
		gene_genetic_trait = sim_alpha[gene_iter]*np.dot(gwas_geno[:, gene_info['gene' + str(gene_iter)]['indices']], gene_info['gene' + str(gene_iter)]['causal_eqtl_effects'])
		genetic_trait = genetic_trait + gene_genetic_trait
	Y = np.random.normal(loc=genetic_trait, scale=np.sqrt(1.0-np.var(genetic_trait)))
	Y = (Y - np.mean(Y))/np.std(Y)

	# Get gwas summary stats
	gwas_chi_sq = []
	for snp_iter in range(n_snps):
		olser = sm.OLS(Y, gwas_geno[:,snp_iter]).fit()
		beta = olser.params[0]
		beta_se = olser.bse[0]
		zed = beta/beta_se
		gwas_chi_sq.append(np.square(zed))
	gwas_chi_sq = np.asarray(gwas_chi_sq)

	gene_names = np.sort([*gene_info])



	joint_sldsc = joint_sldsc_no_ld.JOINT_SLDSC(max_iter=30)
	joint_sldsc.fit(gwas_chi_sq, gwas_ss, gene_info)

	avg_gene_h2 = np.mean(joint_sldsc.gene_h2s)

	avg_gene_h2s.append(avg_gene_h2)
	med_h2s.append(joint_sldsc.gene_med_h2)
	nm_h2s.append(joint_sldsc.nm_h2)

avg_gene_h2s = np.asarray(avg_gene_h2s)
med_h2s = np.asarray(med_h2s)
nm_h2s = np.asarray(nm_h2s)

pdb.set_trace()


'''
# Get gene ld scores
gene_ld_scores = np.zeros(n_snps)
for gene_iter in range(n_genes):
	gene_snp_indices = gene_info['gene' + str(gene_iter)]['indices']
	gene_snp_est_causal_eqtl_effects =  gene_info['gene' + str(gene_iter)]['est_causal_eqtl_effects']
	gene_ld_scores[gene_snp_indices] = gene_ld_scores[gene_snp_indices] + np.square(gene_snp_est_causal_eqtl_effects)

var_ld_scores = np.ones((n_snps,1))
joint_ld_scores = np.hstack((var_ld_scores,gene_ld_scores.reshape(n_snps,1)))
ldsc_reg = sm.OLS(gwas_chi_sq - 1.0, joint_ld_scores).fit()
'''
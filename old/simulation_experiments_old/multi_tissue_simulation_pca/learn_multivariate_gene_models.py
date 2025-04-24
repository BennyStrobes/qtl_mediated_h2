import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')
import statsmodels.api as sm



def get_cis_eqtl_h2_est_w_ldsc(geno, expr):
	NN, GG = geno.shape
	marginal_z_scores = []
	for snp_iter in range(GG):
		model = sm.OLS(expr,sm.add_constant(geno[:, snp_iter])).fit()
		marginal_zed = model.params[1]/model.bse[1]
		marginal_z_scores.append(marginal_zed)
	marginal_z_scores = np.asarray(marginal_z_scores)

	ld_mat = np.corrcoef(np.transpose(geno))
	squared_ld_mat = np.square(ld_mat)

	squared_adj_ld_mat = squared_ld_mat - ((1.0 - squared_ld_mat)/(NN-2.0))

	ld_scores = np.sum(squared_adj_ld_mat,axis=0)


	# LD score regression
	model = sm.OLS(np.square(marginal_z_scores) - 1.0, ld_scores).fit()
	h2_est = model.params[0]*GG/NN

	return h2_est





def get_gene_names(sumstat_output_file):
	f = open(sumstat_output_file)
	head_count = 0
	genes = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		genes.append(data[0])
	f.close()

	return np.unique(np.asarray(genes))



#######################
# Command line args
#######################
simulation_number= sys.argv[1]
simulated_learned_gene_models_dir = sys.argv[2]
simulation_name_string = sys.argv[3]
eqtl_sample_size = sys.argv[4]
alt_simulated_learned_gene_models_dir = sys.argv[5]


version = 'small'

# Get names of genes
sumstat_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
gene_names = get_gene_names(sumstat_output_file)

multivariate_gene_model_output_file = alt_simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_multivariate_gene_model_output.txt'
t = open(multivariate_gene_model_output_file,'w')
t.write('gene\tgene_model_boolean\tldsc_cis_snp_h2\tcausal_effects_file\tsnp_names_file\tpredicted_genetic_variance\n')

# Loop through genes
for ensamble_id in gene_names:

	expr_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_' + ensamble_id + '_expression.npy' 
	genotype_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_' + ensamble_id.split(':')[0] + '_genotype.npy'
	cis_snp_names_output = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_' + ensamble_id.split(':')[0] + '_genotype_cis_snp_names.npy' 


	expr = np.load(expr_output_file)
	geno = np.load(genotype_output_file)
	geno_rsids = np.load(cis_snp_names_output)

	# Get cis-eqtl heritability estimate from S-LDSC
	cis_h2_est = get_cis_eqtl_h2_est_w_ldsc(geno, expr)

	# Run eQTL variant fine-mapping with SuSiE
	susie_fitted = susieR_pkg.susie(geno, expr,L=10)


	# Test whether there exist any identified susie components for this gene
	if type(susie_fitted.rx2('sets').rx2('cs_index')) != rpy2.rinterface_lib.sexp.NULLType:
		susie_components = np.asarray(susie_fitted.rx2('sets').rx2('cs_index')) - 1
		pmces = np.sum(susie_fitted.rx2('alpha')*susie_fitted.rx2('mu'),axis=0)
		booler = 'True'

		pred_var = np.var(np.dot(geno, pmces))

		susie_pmces_file = alt_simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_' + ensamble_id + '_susie_pmces.npy' 
		np.save(susie_pmces_file,pmces)
		susie_alpha_file = alt_simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_' + ensamble_id + '_susie_alpha.npy' 
		np.save(susie_alpha_file,susie_fitted.rx2('alpha'))
		susie_mu_file = alt_simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_' + ensamble_id + '_susie_mu.npy' 
		np.save(susie_mu_file,susie_fitted.rx2('mu'))
		susie_mu2_file = alt_simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_' + ensamble_id + '_susie_mu2.npy' 
		np.save(susie_mu2_file,susie_fitted.rx2('mu2'))

		t.write(ensamble_id + '\t' + booler + '\t' + str(cis_h2_est) + '\t' + susie_pmces_file + '\t' + cis_snp_names_output + '\t' + str(pred_var) + '\n')

	else:
		booler = 'False'
		t.write(ensamble_id + '\t' + booler + '\t' + str(cis_h2_est) + '\t' + 'NA' + '\t' + cis_snp_names_output + '\t' + '0.0' + '\n')

t.close()
print(multivariate_gene_model_output_file)

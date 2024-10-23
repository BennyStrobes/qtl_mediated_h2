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

version = 'small'

# Get names of genes
sumstat_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_small_window_eqtl_sumstats.txt'
gene_names = get_gene_names(sumstat_output_file)

multivariate_gene_model_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_multivariate_gene_model_output.txt'
t = open(multivariate_gene_model_output_file,'w')
t.write('gene\tgene_model_boolean\tcausal_effects_file\tsnp_names_file\n')


# Loop through genes
for ensamble_id in gene_names:

	expr_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_' + ensamble_id + '_expression.npy' 
	genotype_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_' + ensamble_id + '_genotype.npy'
	cis_snp_names_output = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_' + ensamble_id + '_genotype_cis_snp_names.npy' 

	expr = np.load(expr_output_file)
	geno = np.load(genotype_output_file)
	geno_rsids = np.load(cis_snp_names_output)


	# Run eQTL variant fine-mapping with SuSiE
	susie_fitted = susieR_pkg.susie(geno, expr,L=10)


	# Test whether there exist any identified susie components for this gene
	if type(susie_fitted.rx2('sets').rx2('cs_index')) != rpy2.rinterface_lib.sexp.NULLType:
		susie_components = np.asarray(susie_fitted.rx2('sets').rx2('cs_index')) - 1
		pmces = np.sum(susie_fitted.rx2('alpha')*susie_fitted.rx2('mu'),axis=0)
		booler = 'True'

		susie_pmces_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + version + '_' + ensamble_id + '_susie_pmces.npy' 
		np.save(susie_pmces_file,pmces)
		t.write(ensamble_id + '\t' + booler + '\t' + susie_pmces_file + '\t' + cis_snp_names_output + '\n')

	else:
		booler = 'False'
		t.write(ensamble_id + '\t' + booler + '\t' + 'NA' + '\t' + cis_snp_names_output + '\n')

t.close()
print(multivariate_gene_model_output_file)

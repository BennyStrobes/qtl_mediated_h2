import sys
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from bgen import BgenReader
import time




def load_in_ref_alt_allele_arr(pvar_file):
	ref_alt_alleles = []
	f = open(pvar_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 5:
			print('assumptino erooror')
			pdb.set_trace()
		ref_allele = data[3]
		alt_allele = data[4]
		if ref_allele == alt_allele:
			print('assumptino eororor')
			pdb.set_trace()
		ref_alt_alleles.append((ref_allele, alt_allele))
	f.close()
	return ref_alt_alleles

def load_in_alt_allele_genotype_dosage_mat(bfile, window_indices, ref_alt_alleles):
	dosages = []

	for window_index in window_indices:
		var = bfile[window_index]
		dosage = var.alt_dosage
		#ma = var.minor_allele

		index_ref_alt_allele = ref_alt_alleles[window_index]

		if var.alleles[0] != index_ref_alt_allele[1]:
			print('assumption errror')
			pdb.set_trace()
		if var.alleles[1] != index_ref_alt_allele[0]:
			print('assumption eroror')
			pdb.set_trace()

		# Append snp dosage to global array
		dosages.append(dosage)

	# Convert to 2d matrix
	dosages = np.asarray(dosages)

	return np.transpose(dosages)


def load_in_alt_allele_sdevs(bfile, window_indices, ref_alt_alleles):
	sdevs = []
	for window_index in window_indices:
		var = bfile[window_index]
		dosage = var.alt_dosage
		#ma = var.minor_allele

		index_ref_alt_allele = ref_alt_alleles[window_index]

		if var.alleles[0] != index_ref_alt_allele[1]:
			print('assumption errror')
			pdb.set_trace()
		if var.alleles[1] != index_ref_alt_allele[0]:
			print('assumption eroror')
			pdb.set_trace()

		# Compute standard deviation
		sdevs.append(np.std(dosage))

	return np.asarray(sdevs)


def get_genotype_standard_deviations(genotype_dosage, rep1_indices, rep2_indices):
	# Quick error checking to make sure there do not exist missing entries
	n_missing = np.sum(np.isnan(genotype_dosage))
	if n_missing != 0:
		print('assumption eroror')
		pdb.set_trace()
	# Now standardize genotype of each snp
	n_snps = genotype_dosage.shape[1]

	global_sdevs = []
	rep1_sdevs = []
	rep2_sdevs = []

	for snp_iter in range(n_snps):
		xx = genotype_dosage[:,snp_iter]

		global_sdevs.append(np.std(xx))
		rep1_sdevs.append(np.std(xx[rep1_indices]))
		rep2_sdevs.append(np.std(xx[rep2_indices]))


	return np.asarray(global_sdevs), np.asarray(rep1_sdevs), np.asarray(rep2_sdevs)



def standardize_genotype_dosage_matrix(genotype_dosage):
	# Quick error checking to make sure there do not exist missing entries
	n_missing = np.sum(np.isnan(genotype_dosage))
	if n_missing != 0:
		print('assumption eroror')
		pdb.set_trace()

	# Now standardize genotype of each snp
	n_snps = genotype_dosage.shape[1]
	# Initialize standardize genotype dosage matrix
	std_genotype_dosage = np.copy(genotype_dosage)
	for snp_iter in range(n_snps):
		# Standardize
		std_genotype_dosage[:, snp_iter] = (genotype_dosage[:,snp_iter] - np.mean(genotype_dosage[:,snp_iter]))/np.std(genotype_dosage[:,snp_iter])

	return std_genotype_dosage


def create_mapping_from_rsid_to_snp_index(G_obj_rsids):
	dicti = {}
	for ii, val in enumerate(G_obj_rsids):
		dicti[val] = ii
	return dicti

def extract_regression_snp_indices(regression_snp_summary_file, rsid_to_snp_index):
	f = open(regression_snp_summary_file)
	indices = []
	snp_names = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[1]
		indices.append(rsid_to_snp_index[rsid])
		snp_names.append(rsid)
	f.close()
	return np.asarray(indices), np.asarray(snp_names)


def simulate_gene_expression(G, sim_beta):
	#gene_heritability = np.sum(np.square(sim_beta))
	genetic_pred_ge_tmp = np.dot(G, sim_beta)
	genetic_pred_ge = genetic_pred_ge_tmp - np.mean(genetic_pred_ge_tmp)


	gene_heritability = np.var(genetic_pred_ge)
	if gene_heritability > 1:
		gene_heritability = .99
		print('assumption eroorr')
		pdb.set_trace()

	noise = np.random.normal(loc=0, scale=np.sqrt(1.0-gene_heritability), size=len(genetic_pred_ge))

	# Find right amount of noise such that variance is 1.0
	scale_factors = np.arange(0,3.0,.0001)
	diffs = []
	for scale_factor in scale_factors:
		diff = np.abs(1.0 - np.std(genetic_pred_ge + (scale_factor*noise)))
		diffs.append(diff)
	diffs = np.asarray(diffs)
	best_index = np.argmin(diffs)
	best_scale_factor = scale_factors[best_index]

	# Get gene expression
	# Should have variance really close to 1.0

	ge = genetic_pred_ge + (best_scale_factor*noise)

	if np.abs(np.std(ge) - 1.0) > .025:
		print('assumption eroror')
		pdb.set_trace()
	#ge = np.random.normal(loc=genetic_pred_ge, scale=np.sqrt(1.0-gene_heritability))
	return ge

def compute_marginal_regression_coefficients(sim_stand_expr, gene_geno):
	'''
	n_snps = gene_geno.shape[1]
	marginal_effects = np.zeros(n_snps)
	marginal_effects_se = np.zeros(n_snps)
	pvalue_arr = np.zeros(n_snps)
	for snp_iter in range(n_snps):
		# Now get effect size, standard erorr and z-score for association between standardized genotype and trait
		# Fit model using statsmodels
		mod = sm.OLS(sim_stand_expr, sm.add_constant(gene_geno[:, snp_iter]))
		res = mod.fit()
		# Extract results
		effect_size = res.params[1]
		effect_size_se = res.bse[1]
		pvalue = res.pvalues[1]
		effect_size_z = effect_size/effect_size_se

		#marginal_effects.append(effect_size)
		#marginal_effects_se.append(effect_size_se)
		#pvalue_arr.append(pvalue)
		marginal_effects[snp_iter] = effect_size
		marginal_effects_se[snp_iter] = effect_size_se
		pvalue_arr[snp_iter] = pvalue

	'''

	marginal_effects2 = np.dot(np.transpose(gene_geno), sim_stand_expr)/len(sim_stand_expr)
	marginal_effects_se2 = np.sqrt(np.var(np.transpose((gene_geno*marginal_effects2)) - sim_stand_expr, ddof=2,axis=1)/len(sim_stand_expr)) # 2 df needed for interecept model

	return marginal_effects2, marginal_effects_se2


def marginal_ols(Y, X):
	"""
	Run Y ~ 1 + X[:,k] for each k, returning arrays of slopes and SEs.
	Y: (n,) or (n,1)
	X: (n, K)
	"""
	Y = Y.ravel()
	n, K = X.shape

	# 1) compute means
	Y_mean = Y.mean()
	X_means = X.mean(axis=0)

	# 2) center
	Yc = Y - Y_mean               # shape (n,)
	Xc = X - X_means[None, :]     # shape (n, K)

	# 3) sums of squares/cross‐products
	Sxy = Xc.T.dot(Yc)            # shape (K,)
	Sxx = np.sum(Xc**2, axis=0)   # shape (K,)
	Syy = np.sum(Yc**2)           # scalar

	# 4) slope estimates
	beta = Sxy / Sxx              # shape (K,)

	# 5) residual sum of squares for each k
	RSS = Syy - beta * Sxy        # shape (K,)

	# 6) estimate of σ² (using df = n‐2)
	mse = RSS / (n - 2)           # shape (K,)

	# 7) standard errors
	se = np.sqrt(mse / Sxx)       # shape (K,)

	return beta, se

def fixed_effect_meta(beta, se):
    """
    beta : array‐like, shape (K,)
        Point estimates from each study.
    se : array‐like, shape (K,)
        Standard errors of the estimates.
    Returns: dict with keys 'beta_FE', 'se_FE', 'ci95', 'z', 'p'
    """
    beta = np.asarray(beta, dtype=float)
    se   = np.asarray(se,   dtype=float)
    # inverse‐variance weights
    w = 1.0 / (se**2)

    # fixed‐effect estimate
    beta_FE = np.sum(w * beta) / np.sum(w)

    # standard error of the pooled estimate
    se_FE = np.sqrt(1.0 / np.sum(w))

    return beta_FE, se_FE





def simulate_gene_expression_and_and_compute_eqtl_ss_all_genes_shell(simulated_causal_eqtl_effect_summary_file, eqtl_sample_size, simulation_name_string, processed_genotype_data_dir, simulated_learned_gene_models_dir, chrom_string, n_tissues):

	# Open output file handles
	t = {}
	for tissue_iter in range(n_tissues):
		tissue_sumstat1_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_tissue' + str(tissue_iter) + '_' + str(eqtl_sample_size) + '_replicate1' + '_eqtl_sumstats.txt'
		t['tissue' + str(tissue_iter) + 'replicate1'] = open(tissue_sumstat1_output_file,'w')
		t['tissue' + str(tissue_iter) + 'replicate1'].write('gene_id\trsid\tchr\tpos\ta1\ta2\tbeta\tbeta_se\tz\tin_sample_sdev\n')
		tissue_sumstat2_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_tissue' + str(tissue_iter) + '_' + str(eqtl_sample_size) + '_replicate2' + '_eqtl_sumstats.txt'
		t['tissue' + str(tissue_iter) + 'replicate2'] = open(tissue_sumstat2_output_file,'w')
		t['tissue' + str(tissue_iter) + 'replicate2'].write('gene_id\trsid\tchr\tpos\ta1\ta2\tbeta\tbeta_se\tz\tin_sample_sdev\n')
		tissue_sumstatfull_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_tissue' + str(tissue_iter) + '_' + str(eqtl_sample_size) + '_full' + '_eqtl_sumstats.txt'
		t['tissue' + str(tissue_iter) + 'full'] = open(tissue_sumstatfull_output_file,'w')
		t['tissue' + str(tissue_iter) + 'full'].write('gene_id\trsid\tchr\tpos\ta1\ta2\tbeta\tbeta_se\tz\tin_sample_sdev\n')



	# Get chromosomes corresponding to chrom_string
	if chrom_string == '1_2':
		chrom_arr = np.asarray([1,2])


	# Split samples into two replicates
	rep1_indices, rep2_indices = np.array_split(np.random.permutation(np.arange(eqtl_sample_size)),2)

	rep1_ss = len(rep1_indices)
	rep2_ss = len(rep2_indices)

	for chrom_num in chrom_arr:

		# Load in genotype object for gwas data to get population level standard deviations
		ref_genotype_stem = processed_genotype_data_dir + 'simulated_reference_genotype_data_' + str(chrom_num)
		ref_genotype_obj = BgenReader(ref_genotype_stem + '.bgen')
		ref_G_obj_rsids = np.asarray(ref_genotype_obj.rsids())
		ref_geno_ref_alt_alleles = load_in_ref_alt_allele_arr(ref_genotype_stem + '.pvar')
		ref_geno_sdevs = load_in_alt_allele_sdevs(ref_genotype_obj, np.arange(len(ref_G_obj_rsids)), ref_geno_ref_alt_alleles)
		del ref_genotype_obj

		# Load in genotype object
		genotype_stem = processed_genotype_data_dir + 'simulated_eqtl_' + str(eqtl_sample_size) + '_data_' + str(chrom_num)
		genotype_obj = BgenReader(genotype_stem + '.bgen')
		G_obj_pos = np.asarray(genotype_obj.positions())
		G_obj_rsids = np.asarray(genotype_obj.rsids())
		# Load in ref-alt alleles
		ref_alt_alleles = load_in_ref_alt_allele_arr(genotype_stem + '.pvar')
		genotype_dosage = load_in_alt_allele_genotype_dosage_mat(genotype_obj, np.arange(len(G_obj_rsids)), ref_alt_alleles)
		global_genotype_sdev, rep1_genotype_sdev, rep2_genotype_sdev = get_genotype_standard_deviations(genotype_dosage, rep1_indices, rep2_indices)

		# Quick error checking
		if np.array_equal(np.asarray(ref_alt_alleles), np.asarray(ref_geno_ref_alt_alleles)) == False:
			print('assumption eorooror')
			pdb.set_trace()
		if np.array_equal(ref_G_obj_rsids, G_obj_rsids) == False:
			print('assumption erororo')
			pdb.set_trace()

		# Create mapping from rsid to snp index
		rsid_to_snp_index = create_mapping_from_rsid_to_snp_index(G_obj_rsids)

		# Loop through genes
		head_count = 0
		f = open(simulated_causal_eqtl_effect_summary_file)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			ensamble_id = data[0]
			line_chrom_num = data[1]
			# Skip genes not on this chromosome
			if line_chrom_num != str(chrom_num):
				continue
			gene_tss = int(data[2])
			gene_causal_eqtl_effect_file = data[3]
			gene_cis_snp_ids = data[4]
			gene_cis_snp_indices_file = data[5]

			# Extract simulated causal effect sizes for this gene (cis_snpsXn_tissues)
			sim_causal_eqtl_effect_sizes = np.load(gene_causal_eqtl_effect_file)

			# Extract indices of cis_snps
			cis_snp_indices_raw = np.load(gene_cis_snp_indices_file)
			n_cis_snps = len(cis_snp_indices_raw)
			cis_rsids = G_obj_rsids[cis_snp_indices_raw]

			# Quick error checking
			if np.array_equal(cis_rsids, np.load(gene_cis_snp_ids)) == False:
				print('assumption eroror')
				pdb.set_trace()

			# Extract indices of regression snps
			regression_snp_summary_file = processed_genotype_data_dir + 'gene_level_ld_chr' + str(chrom_num) + '_bins_10_' + ensamble_id + '_regression_snp_summary.txt'
			regression_snp_indices_raw, regression_snps = extract_regression_snp_indices(regression_snp_summary_file, rsid_to_snp_index)
			# quick error check
			if np.array_equal(G_obj_rsids[regression_snp_indices_raw], regression_snps) == False:
				print('assumption eroror')
				pdb.set_trace()

			# Extract scaled matrix of cis snps around the gene (note: not mean centered)
			scaled_cis_geno = genotype_dosage[:,cis_snp_indices_raw]/ref_geno_sdevs[cis_snp_indices_raw]
			cis_regr_geno_rep1 = genotype_dosage[:,regression_snp_indices_raw][rep1_indices,:]
			cis_regr_geno_rep2 = genotype_dosage[:,regression_snp_indices_raw][rep2_indices,:]
			cis_regr_geno = genotype_dosage[:,regression_snp_indices_raw]

			regression_snp_alleles = np.asarray(ref_alt_alleles)[regression_snp_indices_raw,:]
			regression_snp_pos = G_obj_pos[regression_snp_indices_raw]
			regression_snp_rep1_sdev = rep1_genotype_sdev[regression_snp_indices_raw]
			regression_snp_rep2_sdev = rep2_genotype_sdev[regression_snp_indices_raw]
			regression_snp_global_sdev = global_genotype_sdev[regression_snp_indices_raw]



			# Now loop through tissues
			n_tiss = sim_causal_eqtl_effect_sizes.shape[1]

			for tiss_iter in range(n_tiss):
				tissue_sim_causal_eqtl_effect_sizes = sim_causal_eqtl_effect_sizes[:, tiss_iter]

				# Skip genes that are not cis heritable
				if np.max(np.abs(tissue_sim_causal_eqtl_effect_sizes)) == 0:
					continue
			
				# Simulate gene expression in this gene
				sim_expr = simulate_gene_expression(scaled_cis_geno, tissue_sim_causal_eqtl_effect_sizes)
				# Standardize simulated gene expression
				sim_stand_expr = (sim_expr - np.mean(sim_expr))/np.std(sim_expr)
				

				###############################
				# First replicate
				###############################
				# Note these are per-allele effect sizes. Can be converted to standardized by multiplying 
				# the marginal effects and the marginal_effects_se by the desired standard deviation
				marginal_effects_rep1, marginal_effects_se_rep1 = marginal_ols(sim_stand_expr[rep1_indices], cis_regr_geno_rep1)


				# Print to output
				for regression_snp_index, regression_snp in enumerate(regression_snps):
					t['tissue' + str(tiss_iter) + 'replicate1'].write(ensamble_id + '\t' + regression_snp + '\t' + str(chrom_num) + '\t' + str(regression_snp_pos[regression_snp_index]) + '\t' + str(regression_snp_alleles[regression_snp_index,1]) + '\t' + str(regression_snp_alleles[regression_snp_index,0]) + '\t' + str(marginal_effects_rep1[regression_snp_index]) + '\t' + str(marginal_effects_se_rep1[regression_snp_index]) + '\t' + str(marginal_effects_rep1[regression_snp_index]/marginal_effects_se_rep1[regression_snp_index]) + '\t' + str(regression_snp_rep1_sdev[regression_snp_index]) + '\n')
		
				###############################
				# Second replicate
				###############################
				# Note these are per-allele effect sizes. Can be converted to standardized by multiplying 
				# the marginal effects and the marginal_effects_se by the desired standard deviation
				marginal_effects_rep2, marginal_effects_se_rep2 = marginal_ols(sim_stand_expr[rep2_indices], cis_regr_geno_rep2)

				# Print to output
				for regression_snp_index, regression_snp in enumerate(regression_snps):
					t['tissue' + str(tiss_iter) + 'replicate2'].write(ensamble_id + '\t' + regression_snp + '\t' + str(chrom_num) + '\t' + str(regression_snp_pos[regression_snp_index]) + '\t' + str(regression_snp_alleles[regression_snp_index,1]) + '\t' + str(regression_snp_alleles[regression_snp_index,0]) + '\t' + str(marginal_effects_rep2[regression_snp_index]) + '\t' + str(marginal_effects_se_rep2[regression_snp_index]) + '\t' + str(marginal_effects_rep2[regression_snp_index]/marginal_effects_se_rep2[regression_snp_index]) + '\t' + str(regression_snp_rep2_sdev[regression_snp_index]) + '\n')

				###############################
				# Full data
				###############################
				marginal_effects, marginal_effects_se = marginal_ols(sim_stand_expr, cis_regr_geno)
				# Print to output
				for regression_snp_index, regression_snp in enumerate(regression_snps):
					t['tissue' + str(tiss_iter) + 'full'].write(ensamble_id + '\t' + regression_snp + '\t' + str(chrom_num) + '\t' + str(regression_snp_pos[regression_snp_index]) + '\t' + str(regression_snp_alleles[regression_snp_index,1]) + '\t' + str(regression_snp_alleles[regression_snp_index,0]) + '\t' + str(marginal_effects[regression_snp_index]) + '\t' + str(marginal_effects_se[regression_snp_index]) + '\t' + str(marginal_effects[regression_snp_index]/marginal_effects_se[regression_snp_index]) + '\t' + str(regression_snp_global_sdev[regression_snp_index]) + '\n')

		f.close()

	for tissue_iter in range(n_tissues):
		t['tissue' + str(tissue_iter)+ 'replicate1'].close()
		t['tissue' + str(tissue_iter)+ 'replicate2'].close()
		t['tissue' + str(tiss_iter) + 'full'].close()

	return rep1_ss, rep2_ss










#####################
# Command line args
#####################
simulation_number = sys.argv[1]
chrom_string = sys.argv[2]
simulated_gene_expression_dir = sys.argv[3]
simulated_learned_gene_models_dir = sys.argv[4]
simulation_name_string = sys.argv[5]
processed_genotype_data_dir = sys.argv[6]
eqtl_sample_size = sys.argv[7]
n_tissues = int(sys.argv[8])


# Set seed
np.random.seed(int(simulation_number))



############################
# Get previously simulated causal eQTL effect sizes across tissues
############################
# Created file to keep track of causal eqtl effect sizes across genes
simulated_causal_eqtl_effect_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'


############################
# Simulate Gene expression and fit gene models for each data-set (eqtl sample-size), tissue
############################
rep1_ss, rep2_ss = simulate_gene_expression_and_and_compute_eqtl_ss_all_genes_shell(simulated_causal_eqtl_effect_summary_file, int(eqtl_sample_size), simulation_name_string, processed_genotype_data_dir, simulated_learned_gene_models_dir, chrom_string, n_tissues)


global_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_replicate_eqtl_sumstats_xt_summary.txt'
t = open(global_output_file,'w')
t.write('tissue_name\tN_rep1\tN_rep2\tN_full\teqtl_sumstat_rep1_file\teqtl_sumstat_rep2_file\teqtl_sumstat_full_file\n')
for tissue_iter in range(n_tissues):
	tissue_sumstat1_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_tissue' + str(tissue_iter) + '_' + str(eqtl_sample_size) + '_replicate1_eqtl_sumstats.txt'
	tissue_sumstat2_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_tissue' + str(tissue_iter) + '_' + str(eqtl_sample_size) + '_replicate2_eqtl_sumstats.txt'
	tissue_sumstat_full_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_tissue' + str(tissue_iter) + '_' + str(eqtl_sample_size) + '_full_eqtl_sumstats.txt'

	t.write('tissue' + str(tissue_iter) + '\t' + str(rep1_ss) + '\t' + str(rep2_ss) + '\t' + str(rep1_ss + rep2_ss) + '\t' + tissue_sumstat1_output_file + '\t' + tissue_sumstat2_output_file + '\t' + tissue_sumstat_full_output_file + '\n')
t.close()
print(global_output_file)

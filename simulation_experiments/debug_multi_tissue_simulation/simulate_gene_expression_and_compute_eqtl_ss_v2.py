import sys
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from bgen import BgenReader





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
		'''
		# Flip dosage if alt-allele is not equal to minor allele
		if index_ref_alt_allele[1] != ma:
			# Quick error check
			if ma != index_ref_alt_allele[0]:
				print('assumptino eroror')
				pdb.set_trace()
			# Flip dosage
			dosage = 2.0 - dosage
		'''

		# Append snp dosage to global array
		dosages.append(dosage)

	# Convert to 2d matrix
	dosages = np.asarray(dosages)

	return np.transpose(dosages)

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
	genetic_pred_ge = np.dot(G, sim_beta)
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


def simulate_gene_expression_and_and_compute_eqtl_ss_all_genes_shell(simulated_causal_eqtl_effect_summary_file, eqtl_sample_size, simulation_name_string, processed_genotype_data_dir, simulated_learned_gene_models_dir, chrom_string, n_tissues):

	# Open output file handles
	t = {}
	for tissue_iter in range(n_tissues):
		tissue_sumstat_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_tissue' + str(tissue_iter) + '_' + str(eqtl_sample_size) + '_eqtl_sumstats_repeat.txt'
		t['tissue' + str(tissue_iter)] = open(tissue_sumstat_output_file,'w')
		t['tissue' + str(tissue_iter)].write('gene_id\trsid\tchr\tpos\ta1\ta2\tbeta\tbeta_se\tz\n')

	# Get chromosomes corresponding to chrom_string
	if chrom_string == '1_2':
		chrom_arr = np.asarray([1,2])

	for chrom_num in chrom_arr:

		# Load in genotype object
		genotype_stem = processed_genotype_data_dir + 'simulated_eqtl_' + str(eqtl_sample_size) + '_data_' + str(chrom_num)
		genotype_obj = BgenReader(genotype_stem + '.bgen')
		G_obj_pos = np.asarray(genotype_obj.positions())
		G_obj_rsids = np.asarray(genotype_obj.rsids())
		# Load in ref-alt alleles
		ref_alt_alleles = load_in_ref_alt_allele_arr(genotype_stem + '.pvar')
		genotype_dosage = load_in_alt_allele_genotype_dosage_mat(genotype_obj, np.arange(len(G_obj_rsids)), ref_alt_alleles)
		G_obj_geno_stand = standardize_genotype_dosage_matrix(genotype_dosage)
		del genotype_dosage
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

			# Extract standardized matrix of cis snps around the gene
			gene_cis_geno = G_obj_geno_stand[:, cis_snp_indices_raw]
			gene_regr_geno = G_obj_geno_stand[:, regression_snp_indices_raw]
			regression_snp_alleles = np.asarray(ref_alt_alleles)[regression_snp_indices_raw,:]
			regression_snp_pos = G_obj_pos[regression_snp_indices_raw]


			# Now loop through tissues
			n_tiss = sim_causal_eqtl_effect_sizes.shape[1]

			for tiss_iter in range(n_tiss):
				tissue_sim_causal_eqtl_effect_sizes = sim_causal_eqtl_effect_sizes[:, tiss_iter]

				# Skip genes that are not cis heritable
				if np.max(np.abs(tissue_sim_causal_eqtl_effect_sizes)) == 0:
					continue
			
				# Simulate gene expression in this gene
				sim_expr = simulate_gene_expression(gene_cis_geno, tissue_sim_causal_eqtl_effect_sizes)
				# Standardize simulated gene expression
				sim_stand_expr = (sim_expr - np.mean(sim_expr))/np.std(sim_expr)

				# Compute marginal association statistics for this gene
				marginal_effects, marginal_effects_se = compute_marginal_regression_coefficients(sim_stand_expr, gene_regr_geno)

				for regression_snp_index, regression_snp in enumerate(regression_snps):
					t['tissue' + str(tiss_iter)].write(ensamble_id + '\t' + regression_snp + '\t' + str(chrom_num) + '\t' + str(regression_snp_pos[regression_snp_index]) + '\t' + str(regression_snp_alleles[regression_snp_index,1]) + '\t' + str(regression_snp_alleles[regression_snp_index,0]) + '\t' + str(marginal_effects[regression_snp_index]) + '\t' + str(marginal_effects_se[regression_snp_index]) + '\t' + str(marginal_effects[regression_snp_index]/marginal_effects_se[regression_snp_index]) + '\n')
		f.close()


	for tissue_iter in range(n_tissues):
		t['tissue' + str(tissue_iter)].close()

	return










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
np.random.seed(int(simulation_number)+1)



############################
# Get previously simulated causal eQTL effect sizes across tissues
############################
# Created file to keep track of causal eqtl effect sizes across genes
simulated_causal_eqtl_effect_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'


############################
# Simulate Gene expression and fit gene models for each data-set (eqtl sample-size), tissue
############################
simulate_gene_expression_and_and_compute_eqtl_ss_all_genes_shell(simulated_causal_eqtl_effect_summary_file, int(eqtl_sample_size), simulation_name_string, processed_genotype_data_dir, simulated_learned_gene_models_dir, chrom_string, n_tissues)



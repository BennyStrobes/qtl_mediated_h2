import sys
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from bgen import BgenReader



def load_in_sldsc_tau_from_sldsc_results_file(real_sldsc_results_file):
	aa = np.loadtxt(real_sldsc_results_file,dtype=str)
	tau = aa[1:,-3].astype(float)
	return tau

def extract_expected_per_snp_heritabilities(sldsc_tau, variant_annotation_file):
	snp_h2 = []
	snp_names = []
	f = open(variant_annotation_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rs_id = data[2]
		snp_names.append(rs_id)
		snp_anno_vec = np.asarray(data[4:]).astype(float)
		snp_h2.append(np.dot(sldsc_tau, snp_anno_vec))
	f.close()
	return np.asarray(snp_h2), np.asarray(snp_names)


# Simulate non-mediated variant causal effect sizes
def simulate_non_mediated_variant_causal_effect_sizes(gwas_plink_stem, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, non_mediated_causal_effects_output):
	# Extract underlying heritability model
	#sldsc_tau = load_in_sldsc_tau_from_sldsc_results_file(real_sldsc_results_file)

	# Extract expected per snp heritability of snps
	#snp_expected_heritability, snp_names = extract_expected_per_snp_heritabilities(sldsc_tau, variant_annotation_file)
	# Set negative per snp heritabilities to zero
	#snp_expected_heritability[snp_expected_heritability < 0.0] = 0.0

	# Convert whole thing to probability
	#snp_probability = snp_expected_heritability/np.sum(snp_expected_heritability)


	snp_names = []
	bim_file = gwas_plink_stem + '.pvar'
	f = open(bim_file)
	head_count = 0
	for line in f:
		if head_count == 0:
			head_count = head_count + 1
			continue
		line = line.rstrip()
		data = line.split('\t')
		snp_names.append(data[2])
	f.close()

	snp_probability = np.ones(len(snp_names))/len(snp_names)

	# Calculate number of causal snps
	total_non_mediated_h2 = total_heritability - (fraction_expression_mediated_heritability*total_heritability)
	n_causal_snps = int(np.round(total_non_mediated_h2/per_element_heritability))

	# Total number of snps
	total_num_snps = len(snp_probability)

	# Extract indices of non-mediated causal snps
	non_med_causal_snp_indices = np.random.choice(np.arange(total_num_snps), size=n_causal_snps, replace=False, p=snp_probability)

	# Extract causal effect sizes
	non_mediated_causal_effect_sizes = np.zeros(total_num_snps)  # Initialization
	non_mediated_causal_effect_sizes[non_med_causal_snp_indices] = np.random.normal(loc=0.0, scale=np.sqrt(per_element_heritability),size=n_causal_snps)

	# Print to output file
	t = open(non_mediated_causal_effects_output,'w')
	for var_iter in range(total_num_snps):
		t.write(snp_names[var_iter] + '\t' + str(non_mediated_causal_effect_sizes[var_iter]) + '\n')
	t.close()

def extract_ordered_gene_names_from_gene_summary_file(simulated_expression_summary_file):
	f = open(simulated_expression_summary_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)


# Extract boolean matrix of total_genesXtissues that is a 1 if gene, tissue is cis_h2
# Causal eqtl effects can only be selected in gene, tissue pairs that are cis heritable
def extract_boolean_matrix_of_cis_heritable_genes(simulated_expression_summary_file):
	boolean_arr = []
	f = open(simulated_expression_summary_file)
	head_count = 0
	gene_counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# extract causal eqtl effect matrix (n_snpsXn_tiss)
		causal_eqtl_effect_size_mat = np.load(data[3])

		# Genes that are not cis heritable will always have 0 eqtl effect size and therefor zero variance
		cis_h2_genes = 1.0*(np.var(causal_eqtl_effect_size_mat,axis=0) != 0.0)
		boolean_arr.append(cis_h2_genes)
	f.close()

	# Convert from list to array
	boolean_mat = np.asarray(boolean_arr)

	return boolean_mat


def simulate_expression_mediated_gene_causal_effect_sizes(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output, ge_h2):
	simulate_expression_mediated_gene_causal_effect_sizes_one_causal_tissues(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output,ge_h2)

	return





# Simulate mediated gene-expression causal effect sizes
def simulate_expression_mediated_gene_causal_effect_sizes_one_causal_tissues(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output, ge_h2):
	# Extract list of ordered gene names
	ordered_gene_names = extract_ordered_gene_names_from_gene_summary_file(simulated_expression_summary_file)
	total_genes = len(ordered_gene_names)

	# Extract boolean matrix of total_genesXtissues that is a 1 if gene, tissue is cis_h2
	# Causal eqtl effects can only be selected in gene, tissue pairs that are cis heritable
	cis_h2_gene_boolean_matrix = extract_boolean_matrix_of_cis_heritable_genes(simulated_expression_summary_file)


	# Quick error checking
	if cis_h2_gene_boolean_matrix.shape[0] != total_genes:
		print('assumption eroror')
		pdb.set_trace()

	# Calculate number of causal genes
	total_mediated_h2 = (fraction_expression_mediated_heritability*total_heritability)
	total_n_causal_genes = int(np.round(total_mediated_h2/per_element_heritability))

	# We are going to assume 50% of gene-mediated h2 is coming from one tissue and 50% of gene-mediated h2 is coming from another tissue
	# And as a reminder there are 10 tissues
	n_causal_genes_per_causal_tissue = int(np.round(total_mediated_h2/per_element_heritability))


	# Get average ge_h2 
	avg_ge_h2 = ((ge_h2/2.0) + ge_h2 + (2.0*ge_h2))/3.0

	# Number of tissues
	n_tiss = cis_h2_gene_boolean_matrix.shape[1]

	# Initialize matrix of gene causal effect sizes
	gene_causal_effect_sizes = np.zeros((total_genes, n_tiss))


	# Extract indices of mediated causal genes for this causal tissue
	cis_h2_genes = np.where(cis_h2_gene_boolean_matrix[:,0] == 1.0)[0]

	tissue_med_causal_gene_indices = np.random.choice(cis_h2_genes, size=n_causal_genes_per_causal_tissue, replace=False)
	# Randomly sample gene causal effect sizes at causal indices
	gene_causal_effect_sizes[tissue_med_causal_gene_indices, 0] = np.random.normal(loc=0.0, scale=np.sqrt(per_element_heritability/avg_ge_h2),size=n_causal_genes_per_causal_tissue)

	# Quick error check
	if np.array_equal(cis_h2_gene_boolean_matrix[tissue_med_causal_gene_indices,0], np.ones(n_causal_genes_per_causal_tissue)) == False:
		print('assumption eroror')

	# Print to output file
	t = open(expression_mediated_causal_effects_output,'w')
	for gene_iter in range(total_genes):
		t.write(ordered_gene_names[gene_iter] + '\t' + '\t'.join(gene_causal_effect_sizes[gene_iter,:].astype(str)) + '\n')
	t.close()

	return



# Simulate mediated gene-expression causal effect sizes
def simulate_expression_mediated_gene_causal_effect_sizes_two_causal_tissues_with_selection(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output, simulated_gt_selection_summary_file):
	# Extract list of ordered gene names
	ordered_gene_names = extract_ordered_gene_names_from_gene_summary_file(simulated_expression_summary_file)
	total_genes = len(ordered_gene_names)

	# Extract boolean matrix of total_genesXtissues that is a 1 if gene, tissue is cis_h2
	# Causal eqtl effects can only be selected in gene, tissue pairs that are cis heritable
	cis_h2_gene_boolean_matrix = extract_boolean_matrix_of_cis_heritable_genes(simulated_expression_summary_file)

	selection_summary_mat = np.loadtxt(simulated_gt_selection_summary_file)

	# QUick error checking
	if np.array_equal(cis_h2_gene_boolean_matrix==1.0,selection_summary_mat>=0.0) == False:
		print('assumption erororo')
		pdb.set_trace()

	# Quick error checking
	if cis_h2_gene_boolean_matrix.shape[0] != total_genes:
		print('assumption eroror')
		pdb.set_trace()

	# Calculate number of causal genes
	total_mediated_h2 = (fraction_expression_mediated_heritability*total_heritability)
	total_n_causal_genes = int(np.round(total_mediated_h2/per_element_heritability))

	# We are going to assume 50% of gene-mediated h2 is coming from one tissue and 50% of gene-mediated h2 is coming from another tissue
	# And as a reminder there are 10 tissues
	n_causal_genes_per_causal_tissue = int(np.round(.5*total_mediated_h2/per_element_heritability))

	# Causal tissues
	causal_tissues = [0,3]

	# Initialize matrix of gene causal effect sizes
	gene_causal_effect_sizes = np.zeros((total_genes, 10))

	for causal_tissue in causal_tissues:
		# Extract indices of mediated causal genes for this causal tissue
		cis_h2_genes = np.where(cis_h2_gene_boolean_matrix[:, causal_tissue] == 1.0)[0]
		#tissue_med_causal_gene_indices = np.random.choice(cis_h2_genes, size=n_causal_genes_per_causal_tissue, replace=False)
		tissue_med_causal_gene_indices = np.where(selection_summary_mat[:, causal_tissue] == 1.0)[0]

		# Quick error check
		if len(tissue_med_causal_gene_indices) != n_causal_genes_per_causal_tissue:
			print('assumption eroro')
			pdb.set_trace()
		if np.sum(cis_h2_gene_boolean_matrix[tissue_med_causal_gene_indices,causal_tissue] != 1.0) > 0.0:
			print('assumption erororo')
			pdb.set_trace()

		# Randomly sample gene causal effect sizes at causal indices
		gene_causal_effect_sizes[tissue_med_causal_gene_indices, causal_tissue] = np.random.normal(loc=0.0, scale=np.sqrt(per_element_heritability),size=n_causal_genes_per_causal_tissue)

		# Quick error check
		if np.array_equal(cis_h2_gene_boolean_matrix[tissue_med_causal_gene_indices,causal_tissue], np.ones(n_causal_genes_per_causal_tissue)) == False:
			print('assumption eroror')

	# Print to output file
	t = open(expression_mediated_causal_effects_output,'w')
	for gene_iter in range(total_genes):
		t.write(ordered_gene_names[gene_iter] + '\t' + '\t'.join(gene_causal_effect_sizes[gene_iter,:].astype(str)) + '\n')
	t.close()

	return


# Simulate mediated gene-expression causal effect sizes
def simulate_expression_mediated_gene_causal_effect_sizes_two_causal_tissues(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output):
	# Extract list of ordered gene names
	ordered_gene_names = extract_ordered_gene_names_from_gene_summary_file(simulated_expression_summary_file)
	total_genes = len(ordered_gene_names)

	# Extract boolean matrix of total_genesXtissues that is a 1 if gene, tissue is cis_h2
	# Causal eqtl effects can only be selected in gene, tissue pairs that are cis heritable
	cis_h2_gene_boolean_matrix = extract_boolean_matrix_of_cis_heritable_genes(simulated_expression_summary_file)

	# Quick error checking
	if cis_h2_gene_boolean_matrix.shape[0] != total_genes:
		print('assumption eroror')
		pdb.set_trace()

	# Calculate number of causal genes
	total_mediated_h2 = (fraction_expression_mediated_heritability*total_heritability)
	total_n_causal_genes = int(np.round(total_mediated_h2/per_element_heritability))

	# We are going to assume 50% of gene-mediated h2 is coming from one tissue and 50% of gene-mediated h2 is coming from another tissue
	# And as a reminder there are 10 tissues
	n_causal_genes_per_causal_tissue = int(np.round(.5*total_mediated_h2/per_element_heritability))

	# Causal tissues
	causal_tissues = [0,3]

	# Initialize matrix of gene causal effect sizes
	gene_causal_effect_sizes = np.zeros((total_genes, 10))

	for causal_tissue in causal_tissues:
		# Extract indices of mediated causal genes for this causal tissue
		cis_h2_genes = np.where(cis_h2_gene_boolean_matrix[:, causal_tissue] == 1.0)[0]
		tissue_med_causal_gene_indices = np.random.choice(cis_h2_genes, size=n_causal_genes_per_causal_tissue, replace=False)
		# Randomly sample gene causal effect sizes at causal indices
		gene_causal_effect_sizes[tissue_med_causal_gene_indices, causal_tissue] = np.random.normal(loc=0.0, scale=np.sqrt(per_element_heritability),size=n_causal_genes_per_causal_tissue)

		# Quick error check
		if np.array_equal(cis_h2_gene_boolean_matrix[tissue_med_causal_gene_indices,causal_tissue], np.ones(n_causal_genes_per_causal_tissue)) == False:
			print('assumption eroror')

	# Print to output file
	t = open(expression_mediated_causal_effects_output,'w')
	for gene_iter in range(total_genes):
		t.write(ordered_gene_names[gene_iter] + '\t' + '\t'.join(gene_causal_effect_sizes[gene_iter,:].astype(str)) + '\n')
	t.close()

	return


def load_in_alt_allele_genotype_dosage_mat(bfile, window_indices, ref_alt_alleles):
	dosages = []

	for window_index in window_indices:
		var = bfile[window_index]
		dosage = var.minor_allele_dosage
		ma = var.minor_allele

		index_ref_alt_allele = ref_alt_alleles[window_index]

		# Flip dosage if alt-allele is not equal to minor allele
		if index_ref_alt_allele[1] != ma:
			# Quick error check
			if ma != index_ref_alt_allele[0]:
				print('assumptino eroror')
				pdb.set_trace()
			# Flip dosage
			dosage = 2.0 - dosage

		# Append snp dosage to global array
		dosages.append(dosage)

	# Convert to 2d matrix
	dosages = np.asarray(dosages)

	return np.transpose(dosages)


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




def compute_expression_mediated_trait_values_for_single_gene(genotype_obj, ref_alt_alleles, gene_trait_effect_size_vec, causal_eqtl_effect_sizes, eqtl_snp_indices, computation_version='v1'):
	# Total number of snps
	global_n_snps = len(eqtl_snp_indices)

	# Extract eqtl variant names
	eqtl_index_positions = np.arange(global_n_snps)[eqtl_snp_indices].astype(int)

	# Extract genotype matrix for these snps
	genotype_dosage = load_in_alt_allele_genotype_dosage_mat(genotype_obj, eqtl_index_positions, ref_alt_alleles)
	stand_eqtl_genotype = standardize_genotype_dosage_matrix(genotype_dosage)
	#eqtl_genotype = np.asarray(genotype_obj.sel(variant=eqtl_index_names))
	#stand_eqtl_genotype = mean_impute_and_standardize_genotype(eqtl_genotype)

	# NOTE: THE FOLLOWING TWO VERSIONS GIVE EQUIVALENT RESULTS
	# IT WAS MORE PLACED HERE FOR  A TEACHING EXCERCISE
	if computation_version == 'v1':
		# Get genetically predicted gene expression (in each tissue)
		genetically_predicted_gene_expression = np.dot(stand_eqtl_genotype, causal_eqtl_effect_sizes)  # Matrix of number of samplesXnumber of tissues

		'''
		# Compute variance (standard deviation) of genetically predicted expression
		genetically_predicted_gene_expression_sdev = np.std(genetically_predicted_gene_expression,axis=0)

		# Quick error check
		if np.sum((genetically_predicted_gene_expression_sdev==0) & (gene_trait_effect_size_vec > 0)) != 0:
			print('assumption error')
			pdb.set_trace()

		# hack to get rid of nans in tissues for this gene that are not genetically heritable
		genetically_predicted_gene_expression_sdev[genetically_predicted_gene_expression_sdev==0] = 1

		# Standardize genetically predicted gene expression
		standardized_genetically_predicted_gene_expression = genetically_predicted_gene_expression/genetically_predicted_gene_expression_sdev
		'''
		# Get expression mediated trait values
		expr_med_trait_for_current_gene = np.dot(genetically_predicted_gene_expression, gene_trait_effect_size_vec)
	else:
		print('assumption eroror')
		pdb.set_trace()

	return expr_med_trait_for_current_gene




def compute_non_mediated_variant_mediated_trait_values(non_mediated_variant_causal_effects_file, gwas_plink_stem, variant_mediated_trait_values_output, n_gwas_individuals):
	# Initialize variant mediated trait values
	variant_med_trait = np.zeros(n_gwas_individuals)

	# Load in genotype object
	# Load in genotype object
	genotype_obj = BgenReader(gwas_plink_stem + '.bgen')

	# Load in ref-alt alleles
	ref_alt_alleles = load_in_ref_alt_allele_arr(gwas_plink_stem + '.pvar')


	# Extract non-mediated variant-trait effect sizes
	var_trait_effect_sizes = []
	ordered_rsids = []
	f = open(non_mediated_variant_causal_effects_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		ordered_rsids.append(data[0])
		var_trait_effect_sizes.append(float(data[1]))
	f.close()
	ordered_rsids = np.asarray(ordered_rsids)
	var_trait_effect_sizes = np.asarray(var_trait_effect_sizes)



	# Loop through variants
	total_nvar = len(ordered_rsids)
	for var_iter in range(total_nvar):
		# Only consider variants with non-zero variant-trait effect
		var_trait_effect_size = var_trait_effect_sizes[var_iter]
		if var_trait_effect_size == 0.0:
			continue

		# Ok, we are at a variant that has non-zero effect on the trait

		# Extract genotype for this variant
		'''
		variant_genotype = np.asarray(genotype_obj.sel(variant='variant' + str(var_iter)))
		std_variant_genotype = np.copy(variant_genotype)
		# Mean impute genotype
		nan_indices = np.isnan(variant_genotype)
		non_nan_mean = np.mean(variant_genotype[nan_indices==False])
		std_variant_genotype[nan_indices] = non_nan_mean
		# Standardize genotype
		std_variant_genotype = (std_variant_genotype - np.mean(std_variant_genotype))/np.std(std_variant_genotype)
		'''
		var = genotype_obj[var_iter]
		dosage = var.minor_allele_dosage
		ma = var.minor_allele
		index_ref_alt_allele = ref_alt_alleles[var_iter]
		if index_ref_alt_allele[1] != ma:
			# Quick error check
			if ma != index_ref_alt_allele[0]:
				print('assumptino errror')
			# flip dosage
			dosage = 2.0 - dosage
		std_variant_genotype = (dosage - np.mean(dosage))/np.std(dosage)

		# Get predicted trait effects from this one variant
		trait_effects_from_single_variant = std_variant_genotype*var_trait_effect_size

		# Update global effects
		variant_med_trait = variant_med_trait + trait_effects_from_single_variant

	# Save to output
	np.savetxt(variant_mediated_trait_values_output, variant_med_trait, fmt="%s", delimiter='\n')

	print('Non-mediated h2: ' + str(np.var(variant_med_trait)))

	return




def compute_expression_mediated_trait_values(simulated_expression_summary_file, expression_mediated_causal_effects_file, gwas_plink_stem, expression_mediated_trait_values_output, n_gwas_individuals):
	# First extract gene-trait effect sizes
	#gene_trait_effect_sizes_raw = np.loadtxt(expression_mediated_causal_effects_file, dtype=str, delimiter='\t')
	#ordered_gene_names = gene_trait_effect_sizes_raw[:,0]
	#gene_trait_effect_sizes = gene_trait_effect_sizes_raw[:,1:].astype(float)
	ordered_gene_names = []
	gene_trait_effect_sizes = []
	f = open(expression_mediated_causal_effects_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		ordered_gene_names.append(data[0])
		gene_trait_effect_sizes.append(np.asarray(data[1:]).astype(float))
	f.close()
	ordered_gene_names = np.asarray(ordered_gene_names)
	gene_trait_effect_sizes = np.asarray(gene_trait_effect_sizes)

	# Initialize expression mediated trait values
	expr_med_trait = np.zeros(n_gwas_individuals)

	# Load in genotype object
	genotype_obj = BgenReader(gwas_plink_stem + '.bgen')

	# Load in ref-alt alleles
	ref_alt_alleles = load_in_ref_alt_allele_arr(gwas_plink_stem + '.pvar')


	#global_tissue_corr = np.zeros((10,10))
	#tissue_corr_counter = 0
	# Loop through genes in simulated_expression_summary_file (should be of same length as ordered_gene_names) and only stop at genes with non-zero gene-trait-effect sizes
	f = open(simulated_expression_summary_file)
	head_count = 0
	gene_counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		#print(gene_counter)
		# This corresponds to a single gene
		# Extract relevent fields for this gene
		gene_name = data[0]
		causal_eqtl_effect_file = data[3]
		eqtl_snp_id_file = data[4]
		eqtl_snp_indices_file = data[5]
		total_n_genome_snps = int(data[6])

		# Quick error check
		if gene_name != ordered_gene_names[gene_counter]:
			print('assumption error')
			pdb.set_trace()

		# We can skip genes that have no gene to trait effect
		gene_trait_effect_size_vec = gene_trait_effect_sizes[gene_counter, :]
		if np.array_equal(gene_trait_effect_size_vec, np.zeros(len(gene_trait_effect_size_vec))):
			gene_counter = gene_counter + 1
			continue

		# We are at a gene that has a non-zero effect on the trait
		
		# Load in causal eqtl effect sizes for this gene
		causal_eqtl_effect_sizes = np.load(causal_eqtl_effect_file)


		# Load in "global" indices of snps defining causal eqtl effect sizes
		eqtl_snp_indices_raw = np.load(eqtl_snp_indices_file)
		eqtl_snp_indices = np.asarray([False]*total_n_genome_snps)
		eqtl_snp_indices[eqtl_snp_indices_raw] = True


		# Compute expression-mediated trait value coming from this single gene (across tissues)
		expr_med_trait_for_current_gene = compute_expression_mediated_trait_values_for_single_gene(genotype_obj,ref_alt_alleles, gene_trait_effect_size_vec, causal_eqtl_effect_sizes, eqtl_snp_indices, computation_version='v1')
		#expr_med_trait_for_current_gene_v2 = compute_expression_mediated_trait_values_for_single_gene(genotype_obj, gene_trait_effect_size_vec, causal_eqtl_effect_sizes, eqtl_snp_indices, computation_version='v2')
		#global_tissue_corr = global_tissue_corr + tissue_corr
		#tissue_corr_counter = tissue_corr_counter + 1

		# Now add current gene to total
		expr_med_trait = expr_med_trait + expr_med_trait_for_current_gene

		# Update gene counter
		gene_counter = gene_counter + 1
	f.close()

	# Save to output
	np.savetxt(expression_mediated_trait_values_output, expr_med_trait, fmt="%s", delimiter='\n')

	print('mediated h2: ' + str(np.var(expr_med_trait)))

	return


def simulate_trait_values(expression_mediated_trait_values_file, variant_mediated_trait_values_file, trait_values_output):
	# Load in genetically predicted trait values
	expr_med_trait = np.loadtxt(expression_mediated_trait_values_file)
	non_med_trait = np.loadtxt(variant_mediated_trait_values_file)

	# Calculate total genetic variance
	expr_med_var = np.var(expr_med_trait)
	non_med_var = np.var(non_med_trait)
	total_genetic_var = np.var(expr_med_trait + non_med_trait)

	# Residual variance
	residual_var = 1.0 - total_genetic_var

	# Draw trait values
	trait_values = np.random.normal(loc=(expr_med_trait + non_med_trait), scale=np.sqrt(residual_var))

	# Standardize trait values
	standardized_trait_values = (trait_values - np.mean(trait_values))/np.std(trait_values)

	# Save to output
	np.savetxt(trait_values_output, standardized_trait_values, fmt="%s", delimiter='\n')

	return


# Command line args
simulation_number = int(sys.argv[1])
chrom_num = sys.argv[2]
cis_window = int(sys.argv[3])
simulated_gene_expression_dir = sys.argv[4]
simulation_name_string = sys.argv[5]
processed_genotype_data_dir = sys.argv[6]
per_element_heritability = float(sys.argv[7])
total_heritability = float(sys.argv[8])
fraction_expression_mediated_heritability = float(sys.argv[9])
simulated_trait_dir = sys.argv[10]  # Output dir
n_gwas_individuals = int(sys.argv[11])
eqtl_architecture = sys.argv[12]
ge_h2 = float('.' + sys.argv[13])


# Set seed
np.random.seed(simulation_number)


####################################################
# Simulate mediated gene-expression causal effect sizes
####################################################
simulated_expression_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'  # This file contains a line for each gene, and we will use it to select which genes in which tissue are used
expression_mediated_causal_effects_output = simulated_trait_dir + simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
simulate_expression_mediated_gene_causal_effect_sizes(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output, ge_h2)


####################################################
# Simulate non-mediated variant causal effect sizes
####################################################
# Need to remove dependency on annotations
gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype directory
#variant_annotation_file = processed_genotype_data_dir + 'baseline.' + chrom_num + '.annot'  # Matrix of annotations describing genetic variants
non_mediated_causal_effects_output = simulated_trait_dir + simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
simulate_non_mediated_variant_causal_effect_sizes(gwas_plink_stem, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, non_mediated_causal_effects_output)



####################################################
# Extract component of trait values coming from genetically predicted gene expression
####################################################
simulated_expression_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'  # This file contains a line for each gene, and also contains causal eqtl effect variant-to-gene effect sizes
expression_mediated_causal_effects_file = simulated_trait_dir + simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'  # File containing gene-to-trait effect sizes
gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype directory
expression_mediated_trait_values_output = simulated_trait_dir + simulation_name_string + '_expression_mediated_trait_values.txt'  # Output file
compute_expression_mediated_trait_values(simulated_expression_summary_file, expression_mediated_causal_effects_file, gwas_plink_stem, expression_mediated_trait_values_output, n_gwas_individuals)


####################################################
# Extract component of trait values coming from non-mediated variant effects
####################################################
non_mediated_variant_causal_effects_file = simulated_trait_dir + simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
variant_mediated_trait_values_output = simulated_trait_dir + simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'  # Output file
compute_non_mediated_variant_mediated_trait_values(non_mediated_variant_causal_effects_file, gwas_plink_stem, variant_mediated_trait_values_output, n_gwas_individuals)


####################################################
# Simulate trait values
####################################################
expression_mediated_trait_values_file = simulated_trait_dir + simulation_name_string + '_expression_mediated_trait_values.txt'  # Now an input file
variant_mediated_trait_values_file = simulated_trait_dir + simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'  # Now an input file
trait_values_output = simulated_trait_dir + simulation_name_string + '_trait_values.txt'  # Output file
simulate_trait_values(expression_mediated_trait_values_file, variant_mediated_trait_values_file, trait_values_output)

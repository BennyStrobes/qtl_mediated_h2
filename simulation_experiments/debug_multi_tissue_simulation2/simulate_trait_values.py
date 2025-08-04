import sys
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from bgen import BgenReader






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

def extract_cis_snp_h2_in_causal_tissues(cis_snp_h2_output):
	cis_snp_h2_arr = []
	gene_names_arr = []

	f = open(cis_snp_h2_output)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_names_arr.append(data[0])
		cis_snp_h2_arr.append(float(data[1]))
	f.close()


	return np.asarray(cis_snp_h2_arr), np.asarray(gene_names_arr)


def simulate_expression_mediated_gene_causal_effect_sizes(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output, ge_h2, gene_trait_architecture, cis_snp_h2_output):
	simulate_expression_mediated_gene_causal_effect_sizes_one_causal_tissues(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output,ge_h2, gene_trait_architecture, cis_snp_h2_output)

	return




# Simulate mediated gene-expression causal effect sizes
def simulate_expression_mediated_gene_causal_effect_sizes_one_causal_tissues(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output, ge_h2, gene_trait_architecture, cis_snp_h2_output, constant_h2=.00053):
	
	# Extract cis snp h2 in causal genes
	gene_cis_snp_h2s, ordered_gene_names2 = extract_cis_snp_h2_in_causal_tissues(cis_snp_h2_output)

	# Extract list of ordered gene names
	ordered_gene_names = extract_ordered_gene_names_from_gene_summary_file(simulated_expression_summary_file)
	total_genes = len(ordered_gene_names)

	# Extract boolean matrix of total_genesXtissues that is a 1 if gene, tissue is cis_h2
	# Causal eqtl effects can only be selected in gene, tissue pairs that are cis heritable
	cis_h2_gene_boolean_matrix = extract_boolean_matrix_of_cis_heritable_genes(simulated_expression_summary_file)


	# Quick error checking
	if np.array_equal(ordered_gene_names2, ordered_gene_names) == False:
		print('assumption errror')
		pdb.set_trace()
	if cis_h2_gene_boolean_matrix.shape[0] != total_genes:
		print('assumption eroror')
		pdb.set_trace()

	if np.array_equal(1.0*(gene_cis_snp_h2s!=0), cis_h2_gene_boolean_matrix[:,0]) == False:
		print('assumption erororor')
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
	if gene_trait_architecture == 'linear':
		gene_causal_effect_sizes[tissue_med_causal_gene_indices, 0] = np.random.normal(loc=0.0, scale=np.sqrt(per_element_heritability/avg_ge_h2),size=n_causal_genes_per_causal_tissue)
	elif gene_trait_architecture == 'stdExpr':
		gene_causal_effect_sizes[tissue_med_causal_gene_indices, 0] = np.random.normal(loc=0.0, scale=np.sqrt(per_element_heritability),size=n_causal_genes_per_causal_tissue)
	elif gene_trait_architecture == 'exprH2Linear':
		for gene_index in tissue_med_causal_gene_indices:
			gene_cis_snp_h2 = gene_cis_snp_h2s[gene_index]
			if gene_cis_snp_h2 == 0.0:
				print('assumption errror')
				pdb.set_trace()

			slope = (per_element_heritability - constant_h2)/(-avg_ge_h2)
			gene_heritability = constant_h2 - gene_cis_snp_h2*slope
			gene_causal_effect_sizes[gene_index, 0] = np.random.normal(loc=0.0, scale=np.sqrt(gene_heritability/avg_ge_h2))

	else:
		print('not yet implemented')
		pdb.set_trace()

	# Quick error check
	if np.array_equal(cis_h2_gene_boolean_matrix[tissue_med_causal_gene_indices,0], np.ones(n_causal_genes_per_causal_tissue)) == False:
		print('assumption eroror')


	# Print to output file
	t = open(expression_mediated_causal_effects_output,'w')
	for gene_iter in range(total_genes):
		t.write(ordered_gene_names[gene_iter] + '\t' + '\t'.join(gene_causal_effect_sizes[gene_iter,:].astype(str)) + '\n')
	t.close()

	return


# Simulate non-mediated variant causal effect sizes
def simulate_non_mediated_variant_causal_effect_sizes(processed_genotype_data_dir, chrom_string, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, non_mediated_causal_effects_output):
	# Get ordered array of all snps
	snp_names = []

	# Get chromosomes corresponding to chrom_string
	if chrom_string == '1_2':
		chrom_arr = np.asarray([1,2])

	chrom_nums = []
	for chrom_num in chrom_arr:
		gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  
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
			chrom_nums.append(str(chrom_num))
		f.close()
	snp_names = np.asarray(snp_names)
	chrom_nums = np.asarray(chrom_nums)

	if len(snp_names) != len(chrom_nums):
		print('assumptione rororo')
		pdb.set_trace()

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
		t.write(snp_names[var_iter] + '\t' + str(chrom_nums[var_iter]) + '\t' + str(non_mediated_causal_effect_sizes[var_iter]) + '\n')
	t.close()

	return

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


def compute_cis_h2_single_gene(genotype_obj, ref_alt_alleles, causal_eqtl_effect_sizes, eqtl_snp_indices):
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
	genetically_predicted_gene_expression = np.dot(stand_eqtl_genotype, causal_eqtl_effect_sizes)  # Matrix of number of samplesXnumber of tissues


	return np.var(genetically_predicted_gene_expression,axis=0)



def compute_expression_mediated_trait_values_for_single_gene(genotype_obj, ref_alt_alleles, gene_trait_effect_size_vec, causal_eqtl_effect_sizes, eqtl_snp_indices, gene_trait_architecture):
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
	if gene_trait_architecture == 'linear' or gene_trait_architecture == 'exprH2Linear':
		# Get genetically predicted gene expression (in each tissue)
		genetically_predicted_gene_expression = np.dot(stand_eqtl_genotype, causal_eqtl_effect_sizes)  # Matrix of number of samplesXnumber of tissues

		# Get expression mediated trait values
		expr_med_trait_for_current_gene = np.dot(genetically_predicted_gene_expression, gene_trait_effect_size_vec)
	elif gene_trait_architecture == 'stdExpr':
		genetically_predicted_gene_expression = np.dot(stand_eqtl_genotype, causal_eqtl_effect_sizes)  # Matrix of number of samplesXnumber of tissues
		
		genetically_predicted_gene_expression_sdev = np.std(genetically_predicted_gene_expression,axis=0)

		# Quick error check
		if np.sum((genetically_predicted_gene_expression_sdev==0) & (gene_trait_effect_size_vec > 0)) != 0:
			print('assumption error')
			pdb.set_trace()

		# hack to get rid of nans in tissues for this gene that are not genetically heritable
		genetically_predicted_gene_expression_sdev[genetically_predicted_gene_expression_sdev==0] = 1

		# Standardize genetically predicted gene expression
		standardized_genetically_predicted_gene_expression = genetically_predicted_gene_expression/genetically_predicted_gene_expression_sdev

		# Get expression mediated trait values
		expr_med_trait_for_current_gene = np.dot(standardized_genetically_predicted_gene_expression, gene_trait_effect_size_vec)


	else:
		print('assumption eroror')
		pdb.set_trace()

	return expr_med_trait_for_current_gene, np.var(genetically_predicted_gene_expression)



def compute_expression_mediated_trait_values_for_single_gene_v2(genotype_obj, ref_alt_alleles, gene_trait_effect_size_vec, causal_eqtl_effect_sizes, eqtl_snp_indices, gene_trait_architecture, snp_names, snp_gene_effects_dicti, gene_name):
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
	if gene_trait_architecture == 'linear' or gene_trait_architecture == 'exprH2Linear':
		# Get genetically predicted gene expression (in each tissue)
		genetically_predicted_gene_expression = np.dot(stand_eqtl_genotype, causal_eqtl_effect_sizes)  # Matrix of number of samplesXnumber of tissues

		# Get expression mediated trait values
		expr_med_trait_for_current_gene = np.dot(genetically_predicted_gene_expression, gene_trait_effect_size_vec)
	elif gene_trait_architecture == 'stdExpr':
		genetically_predicted_gene_expression = np.dot(stand_eqtl_genotype, causal_eqtl_effect_sizes)  # Matrix of number of samplesXnumber of tissues
		
		genetically_predicted_gene_expression_sdev = np.std(genetically_predicted_gene_expression,axis=0)

		# Quick error check
		if np.sum((genetically_predicted_gene_expression_sdev==0) & (gene_trait_effect_size_vec > 0)) != 0:
			print('assumption error')
			pdb.set_trace()

		# hack to get rid of nans in tissues for this gene that are not genetically heritable
		genetically_predicted_gene_expression_sdev[genetically_predicted_gene_expression_sdev==0] = 1

		# Standardize genetically predicted gene expression
		standardized_genetically_predicted_gene_expression = genetically_predicted_gene_expression/genetically_predicted_gene_expression_sdev

		# Get expression mediated trait values
		expr_med_trait_for_current_gene = np.dot(standardized_genetically_predicted_gene_expression, gene_trait_effect_size_vec)


	else:
		print('assumption eroror')
		pdb.set_trace()


	orig_effects = []
	new_effects = []
	for snp_iter, snp_id in enumerate(snp_names):
		gene_snp_name = gene_name + '_' + snp_id
		if gene_snp_name not in snp_gene_effects_dicti:
			continue
		orig_effect = snp_gene_effects_dicti[gene_snp_name]
		new_effect = sm.OLS(standardized_genetically_predicted_gene_expression[:,0], sm.add_constant(stand_eqtl_genotype[:,snp_iter])).fit().params[1]
		orig_effects.append(orig_effect)
		new_effects.append(new_effect)
	print(np.corrcoef(orig_effects, new_effects)[0,1])
	print(sm.OLS(orig_effects, new_effects).fit().params)


	return expr_med_trait_for_current_gene, np.var(genetically_predicted_gene_expression)



def create_gene_snp_effect_mapping(file_name):
	dicti = {}

	f = open(file_name)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')

		if head_count == 0:
			head_count = head_count + 1
			continue

		gene_id = data[0]
		rs_id = data[1]

		effect_size = float(data[6])
		gene_snp_pair = gene_id + '_' + rs_id
		if gene_snp_pair in dicti:
			print('assumptione rororor')
			pdb.set_trace()
		dicti[gene_snp_pair] = effect_size

	f.close()

	return dicti


def compute_expression_mediated_trait_values(simulated_expression_summary_file, expression_mediated_causal_effects_file, processed_genotype_data_dir, chrom_string, expression_mediated_trait_values_output, n_gwas_individuals, gene_trait_architecture):
	# First extract gene-trait effect sizes
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


	# Get chromosomes corresponding to chrom_string
	if chrom_string == '1_2':
		chrom_arr = np.asarray([1,2])


	gene_counter = 0
	for chrom_num in chrom_arr:
		gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype directory

		# Load in genotype object
		genotype_obj = BgenReader(gwas_plink_stem + '.bgen')

		# Load in ref-alt alleles
		ref_alt_alleles = load_in_ref_alt_allele_arr(gwas_plink_stem + '.pvar')
		chrom_rsids = np.asarray(genotype_obj.rsids())

		total_n_chrom_genome_snps = len(ref_alt_alleles)


		# Loop through genes in simulated_expression_summary_file (should be of same length as ordered_gene_names) and only stop at genes with non-zero gene-trait-effect sizes
		f = open(simulated_expression_summary_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			# This corresponds to a single gene
			# Extract relevent fields for this gene
			gene_name = data[0]
			line_chrom_num = data[1]
			causal_eqtl_effect_file = data[3]
			eqtl_snp_id_file = data[4]
			eqtl_snp_indices_file = data[5]
			

			if line_chrom_num != str(chrom_num):
				continue

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
			eqtl_snp_indices = np.asarray([False]*total_n_chrom_genome_snps)
			eqtl_snp_indices[eqtl_snp_indices_raw] = True

			# Quick error checking
			if np.array_equal(chrom_rsids[eqtl_snp_indices], np.load(eqtl_snp_id_file)) == False:
				print('assumption eroror')
				pdb.set_trace()

			# Compute expression-mediated trait value coming from this single gene (across tissues)
			expr_med_trait_for_current_gene, cis_snp_h2 = compute_expression_mediated_trait_values_for_single_gene(genotype_obj, ref_alt_alleles, gene_trait_effect_size_vec, causal_eqtl_effect_sizes, eqtl_snp_indices, gene_trait_architecture)


			# Now add current gene to total
			expr_med_trait = expr_med_trait + expr_med_trait_for_current_gene

			# Update gene counter
			gene_counter = gene_counter + 1
		f.close()

	# Save to output
	np.savetxt(expression_mediated_trait_values_output, expr_med_trait, fmt="%s", delimiter='\n')
	print('mediated h2: ' + str(np.var(expr_med_trait)))

	return


def compute_true_cis_snp_h2(simulated_expression_summary_file, processed_genotype_data_dir, chrom_string, n_gwas_individuals, gene_trait_architecture, cis_snp_h2_output):


	# Initialize expression mediated trait values
	expr_med_trait = np.zeros(n_gwas_individuals)


	# Get chromosomes corresponding to chrom_string
	if chrom_string == '1_2':
		chrom_arr = np.asarray([1,2])

	t = open(cis_snp_h2_output,'w')
	t.write('gene')
	for tissue_iter in range(5):
		t.write('\tcis_h2_tissue' + str(tissue_iter))
	t.write('\n')

	gene_counter = 0
	for chrom_num in chrom_arr:
		gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype directory

		# Load in genotype object
		genotype_obj = BgenReader(gwas_plink_stem + '.bgen')

		# Load in ref-alt alleles
		ref_alt_alleles = load_in_ref_alt_allele_arr(gwas_plink_stem + '.pvar')
		chrom_rsids = np.asarray(genotype_obj.rsids())

		total_n_chrom_genome_snps = len(ref_alt_alleles)

		# Loop through genes in simulated_expression_summary_file (should be of same length as ordered_gene_names) and only stop at genes with non-zero gene-trait-effect sizes
		f = open(simulated_expression_summary_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			# This corresponds to a single gene
			# Extract relevent fields for this gene
			gene_name = data[0]
			line_chrom_num = data[1]
			causal_eqtl_effect_file = data[3]
			eqtl_snp_id_file = data[4]
			eqtl_snp_indices_file = data[5]
			

			if line_chrom_num != str(chrom_num):
				continue

			print(gene_counter)


			# We are at a gene that has a non-zero effect on the trait
		
			# Load in causal eqtl effect sizes for this gene
			causal_eqtl_effect_sizes = np.load(causal_eqtl_effect_file)


			# Load in "global" indices of snps defining causal eqtl effect sizes
			eqtl_snp_indices_raw = np.load(eqtl_snp_indices_file)
			eqtl_snp_indices = np.asarray([False]*total_n_chrom_genome_snps)
			eqtl_snp_indices[eqtl_snp_indices_raw] = True

			# Quick error checking
			if np.array_equal(chrom_rsids[eqtl_snp_indices], np.load(eqtl_snp_id_file)) == False:
				print('assumption eroror')
				pdb.set_trace()


			# Compute expression-mediated trait value coming from this single gene (across tissues)
			per_tissue_cis_snp_h2 = compute_cis_h2_single_gene(genotype_obj, ref_alt_alleles, causal_eqtl_effect_sizes, eqtl_snp_indices)

			t.write(gene_name + '\t' + '\t'.join(per_tissue_cis_snp_h2.astype(str)) + '\n')

			# Update gene counter
			gene_counter = gene_counter + 1
		f.close()

	t.close()

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


def compute_non_mediated_variant_mediated_trait_values(non_mediated_variant_causal_effects_file, processed_genotype_data_dir, chrom_string, variant_mediated_trait_values_output, n_gwas_individuals):
	# Initialize variant mediated trait values
	variant_med_trait = np.zeros(n_gwas_individuals)

	# Get chromosomes corresponding to chrom_string
	if chrom_string == '1_2':
		chrom_arr = np.asarray([1,2])


	gene_counter = 0
	for chrom_num in chrom_arr:

		gwas_plink_stem = processed_genotype_data_dir + 'simulated_gwas_data_' + str(chrom_num)  # Genotype directory


		# Load in genotype object
		# Load in genotype object
		genotype_obj = BgenReader(gwas_plink_stem + '.bgen')

		# Load in ref-alt alleles
		ref_alt_alleles = load_in_ref_alt_allele_arr(gwas_plink_stem + '.pvar')

		chrom_rsids = np.asarray(genotype_obj.rsids())


		# Extract non-mediated variant-trait effect sizes
		var_trait_effect_sizes = []
		ordered_rsids = []
		f = open(non_mediated_variant_causal_effects_file)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if data[1] != str(chrom_num):
				continue
			ordered_rsids.append(data[0])
			var_trait_effect_sizes.append(float(data[2]))
		f.close()
		ordered_rsids = np.asarray(ordered_rsids)
		var_trait_effect_sizes = np.asarray(var_trait_effect_sizes)

		# Quick error check
		if np.array_equal(chrom_rsids, ordered_rsids) == False:
			print('assumptioneoror')
			pdb.set_trace()

		# Loop through variants
		total_nvar = len(ordered_rsids)
		for var_iter in range(total_nvar):
			# Only consider variants with non-zero variant-trait effect
			var_trait_effect_size = var_trait_effect_sizes[var_iter]
			if var_trait_effect_size == 0.0:
				continue

			# Ok, we are at a variant that has non-zero effect on the trait

			# Extract genotype for this variant
			var = genotype_obj[var_iter]

			if var.rsid != ordered_rsids[var_iter]:
				print('assumption eroor')
				pdb.set_trace()

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




# Command line args
simulation_number = int(sys.argv[1])
chrom_string = sys.argv[2]
cis_window = int(sys.argv[3])
simulated_gene_expression_dir = sys.argv[4]
simulation_name_string = sys.argv[5]
gwas_simulation_name_string = sys.argv[6]
processed_genotype_data_dir = sys.argv[7]
per_element_heritability = float(sys.argv[8])
total_heritability = float(sys.argv[9])
fraction_expression_mediated_heritability = float(sys.argv[10])
simulated_trait_dir = sys.argv[11]  # Output dir
n_gwas_individuals = int(sys.argv[12])
eqtl_architecture = sys.argv[13]
ge_h2 = float('.' + sys.argv[14])
gene_trait_architecture = sys.argv[15]

# Set seed
np.random.seed(simulation_number + 3000)


####################################################
# Extract true cis-snp heritabilities
####################################################
simulated_expression_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'  # This file contains a line for each gene, and we will use it to select which genes in which tissue are used
cis_snp_h2_output = simulated_trait_dir + gwas_simulation_name_string + '_true_expression_cis_snp_h2.txt'
compute_true_cis_snp_h2(simulated_expression_summary_file, processed_genotype_data_dir, chrom_string, n_gwas_individuals, gene_trait_architecture, cis_snp_h2_output)

####################################################
# Simulate mediated gene-expression causal effect sizes
####################################################
simulated_expression_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'  # This file contains a line for each gene, and we will use it to select which genes in which tissue are used
expression_mediated_causal_effects_output = simulated_trait_dir + gwas_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'
simulate_expression_mediated_gene_causal_effect_sizes(simulated_expression_summary_file, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, expression_mediated_causal_effects_output, ge_h2, gene_trait_architecture, cis_snp_h2_output)


####################################################
# Simulate non-mediated variant causal effect sizes
####################################################
# Need to remove dependency on annotations
non_mediated_causal_effects_output = simulated_trait_dir + gwas_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
simulate_non_mediated_variant_causal_effect_sizes(processed_genotype_data_dir, chrom_string, per_element_heritability, total_heritability, fraction_expression_mediated_heritability, non_mediated_causal_effects_output)



####################################################
# Extract component of trait values coming from genetically predicted gene expression
####################################################
simulated_expression_summary_file = simulated_gene_expression_dir + simulation_name_string + '_causal_eqtl_effect_summary.txt'  # This file contains a line for each gene, and also contains causal eqtl effect variant-to-gene effect sizes
expression_mediated_causal_effects_file = simulated_trait_dir + gwas_simulation_name_string + '_expression_mediated_gene_causal_effect_sizes.txt'  # File containing gene-to-trait effect sizes
expression_mediated_trait_values_output = simulated_trait_dir + gwas_simulation_name_string + '_expression_mediated_trait_values.txt'  # Output file
compute_expression_mediated_trait_values(simulated_expression_summary_file, expression_mediated_causal_effects_file, processed_genotype_data_dir, chrom_string, expression_mediated_trait_values_output, n_gwas_individuals, gene_trait_architecture)



####################################################
# Extract component of trait values coming from non-mediated variant effects
####################################################
non_mediated_variant_causal_effects_file = simulated_trait_dir + gwas_simulation_name_string + '_non_mediated_variant_causal_effect_sizes.txt'
variant_mediated_trait_values_output = simulated_trait_dir + gwas_simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'  # Output file
compute_non_mediated_variant_mediated_trait_values(non_mediated_variant_causal_effects_file, processed_genotype_data_dir, chrom_string, variant_mediated_trait_values_output, n_gwas_individuals)

####################################################
# Simulate trait values
####################################################
expression_mediated_trait_values_file = simulated_trait_dir + gwas_simulation_name_string + '_expression_mediated_trait_values.txt'  # Now an input file
variant_mediated_trait_values_file = simulated_trait_dir + gwas_simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'  # Now an input file
trait_values_output = simulated_trait_dir + gwas_simulation_name_string + '_trait_values.txt'  # Output file
simulate_trait_values(expression_mediated_trait_values_file, variant_mediated_trait_values_file, trait_values_output)

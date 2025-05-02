import numpy as np
import os
import sys
import pdb
import time
import pickle
import joint_ldsc




def load_in_gwas_data(gwas_summary_file):
	rsids = []
	betas = []
	beta_ses = []

	head_count = 0

	f = open(gwas_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count = head_count + 1
			continue

		rsids.append(data[0])
		betas.append(float(data[5]))
		beta_ses.append(float(data[6]))

	f.close()

	return np.asarray(rsids), np.asarray(betas), np.asarray(beta_ses)

def load_in_regression_snps_for_gene(regression_snp_file):
	f = open(regression_snp_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[1])
	f.close()
	return np.asarray(arr)

def load_in_eqtl_data(rsid_to_position, simulation_genotype_dir, chrom_string, eqtl_snp_representation, tissue_names):
	if chrom_string == '1_2':
		chrom_arr = np.asarray([1,2])
	genes = []
	gene_info = {}

	for chrom_num in chrom_arr:
		gene_summary_file = simulation_genotype_dir + 'gene_level_ld_chr' + str(chrom_num) + '_' + eqtl_snp_representation + '_summary_file.txt'
		f = open(gene_summary_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue

			# Each line is a gene. Extract relevent fields
			ensamble_id = data[1]
			genes.append(ensamble_id)
			squared_ld_mat_file = data[7]
			n_low_dim_snp_file = data[8]
			regression_snp_file = data[6]
			if ensamble_id in gene_info:
				print('assumption eroror')
				pdb.set_trace()
			# Load in regression snps
			gene_regression_snps = load_in_regression_snps_for_gene(regression_snp_file)
			# Get corresponding indices of regression snps
			regression_snp_indices = []
			for rsid in gene_regression_snps:
				regression_snp_indices.append(rsid_to_position[rsid])
			regression_snp_indices = np.asarray(regression_snp_indices)
			# Create mapping from rsid to gene position
			rsid_to_gene_position = {}
			for ii, val in enumerate(gene_regression_snps):
				rsid_to_gene_position[val] = ii


			gene_info[ensamble_id] = {}
			gene_info[ensamble_id]['squared_ld_file'] = squared_ld_mat_file
			gene_info[ensamble_id]['n_low_dimensional_snps'] = np.load(n_low_dim_snp_file)
			gene_info[ensamble_id]['regression_snp_indices'] = regression_snp_indices
			gene_info[ensamble_id]['rsid_to_gene_position'] = rsid_to_gene_position
			gene_info[ensamble_id]['n_gene_regression_snps'] = len(gene_regression_snps)
			gene_info[ensamble_id]['squared_sumstats'] = np.zeros((len(gene_regression_snps), len(tissue_names)))
			gene_info[ensamble_id]['eqtl_category_names'] = np.copy(tissue_names)
		f.close()


	return np.asarray(genes), gene_info


def get_gwas_variant_ld_scores(simulation_genotype_dir, chrom_string):
	if chrom_string == '1_2':
		chrom_arr = np.asarray([1,2])
	rsids = []
	ldscores = []

	for chrom_num in chrom_arr:
		filer = simulation_genotype_dir + 'variant_reference_genotype_data_ldscores_chrom' + str(chrom_num) + '.txt'
		f = open(filer)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			rsid = data[1]
			ldscore = float(data[6])
			rsids.append(rsid)
			ldscores.append(ldscore)
		f.close()


	return np.asarray(ldscores), np.asarray(rsids)


def create_mapping_from_rsid_to_position(rsids):
	rsid_to_position = {}

	for ii, rsid in enumerate(rsids):
		if rsid in rsid_to_position:
			print('assumption eroror')
			pdb.set_trace()
		rsid_to_position[rsid] = ii
	return rsid_to_position

def fill_in_eqtl_sumstats(gene_info, simulated_learned_gene_models_dir, simulation_name_string, eqtl_sample_size, tissue_names,genes):
	for tt,tissue_name in enumerate(tissue_names):
		eqtl_sumstat_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + tissue_name + '_' + str(eqtl_sample_size) + '_eqtl_sumstats.txt'
		f = open(eqtl_sumstat_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			ens_id = data[0]
			rsid = data[1]
			beta = float(data[6])
			beta_se = float(data[7])
			gene_row_index = gene_info[ens_id]['rsid_to_gene_position'][rsid]
			gene_info[ens_id]['squared_sumstats'][gene_row_index, tt] = np.square(beta) - np.square(beta_se)
		f.close()

	# Filter out missing eqtl categories for each gene
	for ens_id in genes:
		valid_categories = []
		for ii, category_name in enumerate(gene_info[ens_id]['eqtl_category_names']):
			gene_info[ens_id]['squared_sumstats'][:, ii]
			if np.array_equal(gene_info[ens_id]['squared_sumstats'][:, ii], np.zeros(len(gene_info[ens_id]['squared_sumstats'][:, ii]))):
				valid_categories.append(False)
			else:
				valid_categories.append(True)
		valid_categories = np.asarray(valid_categories)
		gene_info[ens_id]['squared_sumstats'] = gene_info[ens_id]['squared_sumstats'][:,valid_categories]
		gene_info[ens_id]['eqtl_category_names'] = gene_info[ens_id]['eqtl_category_names'][valid_categories]


	return gene_info

def extract_number_of_reference_snps(simulation_genotype_dir, chrom_string):
	if chrom_string == '1_2':
		chrom_arr = np.asarray([1,2])

	n_snps = 0
	for chrom_num in chrom_arr:
		filer = simulation_genotype_dir + 'simulated_reference_genotype_data_' + str(chrom_num) + '.pvar'
		f = open(filer)
		head_count = 0
		for line in f:
			line = line.rstrip()
			if head_count == 0:
				head_count = head_count + 1
				continue
			n_snps = n_snps + 1
		f.close()


	return n_snps



#######################
# Command line args
#######################
simulation_number= int(sys.argv[1])
simulation_name_string= sys.argv[2]
simulated_trait_dir=sys.argv[3]
simulated_gwas_dir=sys.argv[4]
simulation_genotype_dir=sys.argv[5]
simulated_learned_gene_models_dir=sys.argv[6]
n_gwas_individuals=int(sys.argv[7])
eqtl_sample_size=int(sys.argv[8])
trait_med_h2_inference_dir=sys.argv[9]
simulated_gene_expression_dir=sys.argv[10]
chrom_string = sys.argv[11]
eqtl_snp_representation = sys.argv[12]

N_gwas = n_gwas_individuals


####################
# Load in data
####################
# Load in true simulated data parameters
genetic_trait_expr_med_file = simulated_trait_dir + simulation_name_string +'_expression_mediated_trait_values.txt'
genetic_trait_nm_file = simulated_trait_dir +simulation_name_string + '_non_mediated_variant_mediated_trait_values.txt'
sim_med_h2 = np.var(np.loadtxt(genetic_trait_expr_med_file))
sim_nm_h2 = np.var(np.loadtxt(genetic_trait_nm_file))
sim_h2 = np.var(np.loadtxt(genetic_trait_nm_file) + np.loadtxt(genetic_trait_expr_med_file))

print(sim_h2)
print(sim_med_h2)
print(sim_nm_h2)

# Extract number of reference snps
n_reference_snps = extract_number_of_reference_snps(simulation_genotype_dir, chrom_string)
# Load in GWAS summary statistics
gwas_summary_file = simulated_gwas_dir + simulation_name_string + '_simualated_gwas_results.txt'
gwas_rsids, gwas_beta, gwas_beta_se = load_in_gwas_data(gwas_summary_file)
gwas_E_beta_sq = np.square(gwas_beta) - np.square(gwas_beta_se)

# Load in GWAS variant ld scores
gwas_variant_ld_scores, gwas_rsids_tmp = get_gwas_variant_ld_scores(simulation_genotype_dir, chrom_string)

# QUick error checking
if np.array_equal(gwas_rsids, gwas_rsids_tmp) == False:
	print('assumption eroror')
	pdb.set_trace()

# Create mapping from rsid to position
rsid_to_position = create_mapping_from_rsid_to_position(gwas_rsids)

# Load in eqtl data
tissue_names = []
for itera in range(5):
	tissue_names.append('tissue' + str(itera))
tissue_names = np.asarray(tissue_names)
genes, gene_info = load_in_eqtl_data(rsid_to_position, simulation_genotype_dir, chrom_string, eqtl_snp_representation, tissue_names)
gene_info = fill_in_eqtl_sumstats(gene_info, simulated_learned_gene_models_dir, simulation_name_string, eqtl_sample_size, tissue_names, genes)

'''
#################
# Temp saving
f = open('gene_info.pickle', 'wb')
pickle.dump(gene_info, f)
f.close()

np.save('genes.npy', genes)
np.save('tissue_names.npy', tissue_names)
np.save('gwas_variant_ld_scores.npy', gwas_variant_ld_scores)
np.save('gwas_beta.npy', gwas_beta)
np.save('gwas_beta_se.npy', gwas_beta_se)

#################
# Temp loading
f = open('gene_info.pickle','rb')
gene_info = pickle.load(f)
f.close()

genes = np.load('genes.npy')
tissue_names = np.load('tissue_names.npy')
gwas_variant_ld_scores = np.load('gwas_variant_ld_scores.npy')
gwas_beta = np.load('gwas_beta.npy')
gwas_beta_se = np.load('gwas_beta_se.npy')
gwas_E_beta_sq = np.square(gwas_beta) - np.square(gwas_beta_se)

tissue_names = []
for itera in range(5):
	tissue_names.append('tissue' + str(itera))
tissue_names = np.asarray(tissue_names)
'''

output_file = trait_med_h2_inference_dir + simulation_name_string+ '_' + eqtl_snp_representation + '_' + str(eqtl_sample_size) + '_joint_ldsc_results.txt'
t = open(output_file,'w')
t.write('method\teQTL_SS\tsim_h2\tsim_med_h2\tsim_nm_h2\test_med_h2_\test_med_h2_per_tissue\test_nm_h2_joint_reml\test_mean_eqtl_h2\n')


# Run two step joint ldsc
obj = joint_ldsc.med_h2(version='two_step')
obj.fit(genes, gene_info, gwas_variant_ld_scores,gwas_E_beta_sq, tissue_names, n_reference_snps)
t.write('two_step_joint_ldsc\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(np.sum(obj.med_h2)) + '\t' + ','.join(obj.med_h2.astype(str)) + '\t' + str(obj.nm_h2) + '\t' + str(obj.avg_eqtl_h2) + '\n')

# Run joint ldsc
obj = joint_ldsc.med_h2(max_iter=600)
obj.fit(genes, gene_info, gwas_variant_ld_scores,gwas_E_beta_sq, tissue_names, n_reference_snps)
t.write('joint_ldsc\t' + str(eqtl_sample_size) + '\t' + str(sim_h2) + '\t' + str(sim_med_h2) + '\t' + str(sim_nm_h2) + '\t' + str(np.sum(obj.med_h2)) + '\t' + ','.join(obj.med_h2.astype(str)) + '\t' + str(obj.nm_h2) + '\t' + str(obj.avg_eqtl_h2) + '\n')

t.close()
print(output_file)









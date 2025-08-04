import numpy as np
import os
import sys
import pdb

def linear_regression(XX, YY, intercept=False):
	if intercept:
		print('not yet implemented. add column of 1s to X')
		pdb.set_trace()

	return np.dot(np.dot(np.linalg.pinv(np.dot(np.transpose(XX), XX)), np.transpose(XX)), YY)

def load_in_eqtl_dataset_summary_file(eqtl_summary_file):
	names = []
	sample_sizes1 = []
	sample_sizes2 = []
	sample_sizesfull = []
	filers1 = []
	filers2 = []
	filersfull = []
	f = open(eqtl_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		names.append(data[0])
		sample_sizes1.append(float(data[1]))
		sample_sizes2.append(float(data[2]))
		sample_sizesfull.append(float(data[3]))
		filers1.append(data[4])
		filers2.append(data[5])
		filersfull.append(data[6])
	f.close()
	return np.asarray(names), np.asarray(sample_sizes1), np.asarray(sample_sizes2), np.asarray(filers1), np.asarray(filers2), np.asarray(sample_sizesfull), np.asarray(filersfull)


def extract_chromosome_names(chromosome_file):
	if chromosome_file is None:
		tmp_chrom_arr = np.arange(1,23)
		arr = []
		dicti = {}
		for ele in tmp_chrom_arr:
			arr.append(str(ele))
			dicti[str(ele)] = 1
	else:
		arr = []
		dicti = {}
		f = open(chromosome_file)
		for line in f:
			line = line.rstrip()
			arr.append(line)
			dicti[line] = 1
		f.close()
	return np.asarray(arr), dicti

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


def extract_gene_info(gene_ldscore_filestem, gene_ldscore_filesuffix, eqtl_dataset_names, chrom_arr):
	genes = []
	gene_info = {}

	for chrom_num in chrom_arr:
		gene_summary_file = gene_ldscore_filestem + str(chrom_num) + gene_ldscore_filesuffix
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
			squared_ld_mat_file = data[7]
			n_low_dim_snp_file = data[8]
			regression_snp_file = data[6]
			if ensamble_id in gene_info:
				print('assumption eroror')
				pdb.set_trace()
			# Load in regression snps
			gene_regression_snps = load_in_regression_snps_for_gene(regression_snp_file)
			rsid_to_gene_index = {}
			for ii,val in enumerate(gene_regression_snps):
				rsid_to_gene_index[val] = ii

			gene_info[ensamble_id] = {}
			gene_info[ensamble_id]['squared_ld_file'] = squared_ld_mat_file
			gene_info[ensamble_id]['n_low_dimensional_snps'] = np.load(n_low_dim_snp_file)
			gene_info[ensamble_id]['n_gene_regression_snps'] = len(gene_regression_snps)
			gene_info[ensamble_id]['rsid_to_regression_snp_index'] = rsid_to_gene_index
		f.close()
	return gene_info


def extract_cis_snp_h2_with_ldsc(eqtl_sumstats_file, gene_info, rsid_to_variant_stdev, dataset_cis_h2_file, squared_eqtl_effect_threshold=100.0):
	#######################################
	# First extract squared summary stats
	#######################################
	sq_sum_stats = {}
	valid_genes = {}
	for geneid in [*gene_info]:
		sq_sum_stats[geneid] = np.zeros(gene_info[geneid]['n_gene_regression_snps'])
	f = open(eqtl_sumstats_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		rsid = data[1]
		beta = float(data[6])
		beta_se = float(data[7])
		variant_sdev = rsid_to_variant_stdev[rsid]
		# Update betas and standard errors according to this
		beta = beta*variant_sdev
		beta_se = beta_se*variant_sdev

		reg_index = gene_info[gene_id]['rsid_to_regression_snp_index'][rsid]

		sq_sum_stats[gene_id][reg_index] = np.square(beta) - np.square(beta_se)
		valid_genes[gene_id] = 1
	f.close()

	#######################################
	# Loop through valid genes
	#######################################
	t = open(dataset_cis_h2_file,'w')
	t.write('gene_id\tcis_snp_h2\n')
	dicti = {}
	for gene_id in [*valid_genes]:
		eqtl_sq_sumstats = sq_sum_stats[gene_id]
		if np.sum(eqtl_sq_sumstats == 0) > 0:
			print('assumption eororor')
			pdb.set_trace()

		info = gene_info[gene_id]
		# Load in squared ld 
		squared_ld = np.load(info['squared_ld_file'])

		temp_Y = np.copy(eqtl_sq_sumstats)
		valid_indices = (np.abs(temp_Y) < squared_eqtl_effect_threshold) & (np.isnan(temp_Y) == False)
			
		weights = linear_regression(squared_ld[valid_indices,:], temp_Y[valid_indices])

		gene_est_cis_h2 = np.sum(weights*info['n_low_dimensional_snps'])
		t.write(gene_id + '\t' + str(gene_est_cis_h2) + '\n')

		if gene_id in dicti:
			print('assumption eorroro')
			pdb.set_trace()
		dicti[gene_id] = gene_est_cis_h2
	t.close()

	return dicti


def extract_cis_snp_h2_with_avgChisq(eqtl_sumstats_file, gene_info, rsid_to_variant_stdev, dataset_cis_h2_file, squared_eqtl_effect_threshold=100.0):
	#######################################
	# First extract squared summary stats
	#######################################
	sq_sum_stats = {}
	valid_genes = {}
	for geneid in [*gene_info]:
		sq_sum_stats[geneid] = np.zeros(gene_info[geneid]['n_gene_regression_snps'])
	f = open(eqtl_sumstats_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		rsid = data[1]
		beta = float(data[6])
		beta_se = float(data[7])
		variant_sdev = rsid_to_variant_stdev[rsid]
		# Update betas and standard errors according to this
		beta = beta*variant_sdev
		beta_se = beta_se*variant_sdev

		reg_index = gene_info[gene_id]['rsid_to_regression_snp_index'][rsid]

		sq_sum_stats[gene_id][reg_index] = np.square(beta) - np.square(beta_se)
		valid_genes[gene_id] = 1
	f.close()

	#######################################
	# Loop through valid genes
	#######################################
	t = open(dataset_cis_h2_file,'w')
	t.write('gene_id\tcis_snp_h2\tcis_snp_h2_se\n')
	dicti = {}
	for gene_id in [*valid_genes]:
		eqtl_sq_sumstats = sq_sum_stats[gene_id]
		if np.sum(eqtl_sq_sumstats == 0) > 0:
			print('assumption eororor')
			pdb.set_trace()

		info = gene_info[gene_id]
		# Load in squared ld 
		squared_ld = np.load(info['squared_ld_file'])

		temp_Y = np.copy(eqtl_sq_sumstats)
		valid_indices = (np.abs(temp_Y) < squared_eqtl_effect_threshold) & (np.isnan(temp_Y) == False)

		gene_est_cis_h2 = np.mean(temp_Y[valid_indices])*np.sum(info['n_low_dimensional_snps'])/np.mean(np.sum(squared_ld[valid_indices,:],axis=1))

		gene_est_cis_h2_se = (np.std(temp_Y[valid_indices])*np.sum(info['n_low_dimensional_snps'])/np.mean(np.sum(squared_ld[valid_indices,:],axis=1)))/np.sqrt(len(temp_Y[valid_indices])/np.mean(np.sum(squared_ld[valid_indices,:],axis=1)))

		t.write(gene_id + '\t' + str(gene_est_cis_h2) + '\t' + str(gene_est_cis_h2_se) + '\n')

		if gene_id in dicti:
			print('assumption eorroro')
			pdb.set_trace()
		dicti[gene_id] = gene_est_cis_h2
	t.close()

	return 




def create_mapping_from_rsid_to_variant_standard_deviation(variant_stdev_filestem, chrom_arr):
	dicti = {}

	for chrom_num in chrom_arr:
		filer = variant_stdev_filestem + chrom_num + '.txt'
		f = open(filer)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			rsid = data[0]
			stdev = float(data[5])
			if rsid in dicti:
				print('assumption eororro')
				pdb.set_trace()
			dicti[rsid] = stdev
		f.close()
	return dicti


def extract_cis_snp_h2_with_greml(mesc_file_stem, chrom_arr, ldsc_h2s, output_file):
	# Create mapping from gene to greml cis snp h2
	mapping = {}
	for chrom_num in chrom_arr:
		mesc_file = mesc_file_stem + '_' + str(chrom_num) + '.' + str(chrom_num) + '.hsq'
		f = open(mesc_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				continue
			ens_id = data[0]
			cis_h2_est = data[2]
			cis_h2_est_se = data[3]
			if ens_id in mapping:
				print('assumptikoneororor')
				pdb.set_trace()
			mapping[ens_id] = (cis_h2_est, cis_h2_est_se)
		f.close()

	if len(ldsc_h2s) != len(mapping):
		print('asssumption eroror')
		pdb.set_trace()

	t = open(output_file,'w')
	t.write('gene_id\tcis_snp_h2\tcis_snp_h2_se\n')

	for gene_id in [*ldsc_h2s]:
		if gene_id not in mapping:
			print('assumprion error')
			pdb.set_trace()
		greml_est, greml_est_se = mapping[gene_id]
		if greml_est == 'NA':
			final_est = ldsc_h2s[gene_id]
			final_est_se = '1.0'
		else:
			final_est = mapping[gene_id][0]
			final_est_se = mapping[gene_id][1]

		t.write(gene_id + '\t' + str(final_est) + '\t' + str(final_est_se) + '\n')

	t.close()
	return


def extract_posterior_cis_snp_h2(greml_dataset_cis_h2_file, posterior_greml_dataset_cis_h2_file, tmp_output_stem):
	os.system('Rscript ashr_for_cis_snp_h2_posterior.R ' + greml_dataset_cis_h2_file + ' ' + tmp_output_stem)
	aa = np.loadtxt(tmp_output_stem)

	f = open(greml_dataset_cis_h2_file)
	t = open(posterior_greml_dataset_cis_h2_file,'w')

	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		t.write(data[0] + '\t' + str(aa[counter]) + '\t' + 'NA\n')
		counter = counter + 1
	f.close()
	t.close()
	os.system('rm ' +tmp_output_stem)
	return

def print_true_cis_h2_file(true_cis_snp_h2_file, dataset_cis_h2_file, dataset_iter):
	f = open(true_cis_snp_h2_file)
	t = open(dataset_cis_h2_file,'w')
	t.write('gene_id\tcis_snp_h2\n')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[0]
		cis_snp_h2 = data[(1+dataset_iter)]
		if float(cis_snp_h2) == 0.0:
			continue
		t.write(ensamble_id + '\t' + cis_snp_h2 + '\n')


	f.close()
	t.close()
	return




gene_ldscore_filestem = sys.argv[1]
gene_ldscore_filesuffix = sys.argv[2]
eqtl_summary_file = sys.argv[3]
chrom_nums_file = sys.argv[4]
variant_stdev_filestem = sys.argv[5]
eqtl_sample_size = sys.argv[6]
est_cis_snp_h2_dir = sys.argv[7]
simulation_name_string = sys.argv[8]
mesc_expression_score_dir = sys.argv[9]
simulated_trait_dir = sys.argv[10]


'''
# Names of chromosomes to run analysis on
chrom_arr, chrom_dicti = extract_chromosome_names(chrom_nums_file)


# Create mapping from rsid to variant standard deviation
rsid_to_variant_stdev = create_mapping_from_rsid_to_variant_standard_deviation(variant_stdev_filestem, chrom_arr)

##############################
# Load in eqtl data
eqtl_dataset_names, eqtl_dataset_Ns_training, eqtl_dataset_Ns_validation, eqtl_dataset_files_training, eqtl_dataset_files_validation, eqtl_dataset_Ns_full, eqtl_dataset_files_full = load_in_eqtl_dataset_summary_file(eqtl_summary_file)

##############################
# Load in gene info (ldscores and what not)
gene_info = extract_gene_info(gene_ldscore_filestem, gene_ldscore_filesuffix, eqtl_dataset_names, chrom_arr)
'''


##############################
# Extract cis-snp h2s using unknown truth
##############################
# Create global output file
global_output_file = est_cis_snp_h2_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_true_est_cis_snp_h2_summary.txt'
tt = open(global_output_file,'w')
tt.write('tissue_name\tN_rep1\tN_rep2\tcis_h2_rep1_file\tcis_h2_rep2_file\n')

eqtl_dataset_names = ['tissue0', 'tissue1', 'tissue2', 'tissue3', 'tissue4']

# Loop through eqtl data sets
gene_to_ldsc_cis_h2s = {}
for dataset_iter, eqtl_dataset_name in enumerate(eqtl_dataset_names):

	true_cis_snp_h2_file = simulated_trait_dir + simulation_name_string + '_gt_arch_stdExpr_true_expression_cis_snp_h2.txt'
	dataset_cis_h2_file = est_cis_snp_h2_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_' + eqtl_dataset_name + '_true_cis_snp_h2.txt'

	print_true_cis_h2_file(true_cis_snp_h2_file, dataset_cis_h2_file, dataset_iter)

	###################
	# Replicate 1
	###################
	#N_training = eqtl_dataset_Ns_training[dataset_iter]
	#training_data_file = eqtl_dataset_files_training[dataset_iter]

	###################
	# Replicate 2
	###################
	#N_validation = eqtl_dataset_Ns_validation[dataset_iter]
	#validation_data_file = eqtl_dataset_files_validation[dataset_iter]


	tt.write(eqtl_dataset_name + '\t' + 'INF' + '\t' + 'INF' + '\t' + dataset_cis_h2_file + '\t' + dataset_cis_h2_file + '\n')
tt.close()
print(global_output_file)

'''
##############################
# Extract cis-snp h2s using LDSC
##############################
# Create global output file
global_output_file = est_cis_snp_h2_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_ldsc_est_cis_snp_h2_summary.txt'
tt = open(global_output_file,'w')
tt.write('tissue_name\tN_rep1\tN_rep2\tcis_h2_rep1_file\tcis_h2_rep2_file\n')

# Loop through eqtl data sets
gene_to_ldsc_cis_h2s = {}
for dataset_iter, eqtl_dataset_name in enumerate(eqtl_dataset_names):
	###################
	# Replicate 1
	###################
	N_training = eqtl_dataset_Ns_training[dataset_iter]
	training_data_file = eqtl_dataset_files_training[dataset_iter]

	file_stem = eqtl_dataset_files_training[dataset_iter].split('_eqtl_sumstats')[0].split('/')[-1]
	dataset_cis_h2_rep1_file = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_ldsc_est_cis_snp_h2.txt'
	gene_id_to_cis_snp_h2_mapping_rep1 = extract_cis_snp_h2_with_ldsc(training_data_file, gene_info, rsid_to_variant_stdev, dataset_cis_h2_rep1_file)

	gene_to_ldsc_cis_h2s[eqtl_dataset_name + ':' + 'replicate1'] = gene_id_to_cis_snp_h2_mapping_rep1

	###################
	# Replicate 2
	###################
	N_validation = eqtl_dataset_Ns_validation[dataset_iter]
	validation_data_file = eqtl_dataset_files_validation[dataset_iter]

	file_stem = eqtl_dataset_files_validation[dataset_iter].split('_eqtl_sumstats')[0].split('/')[-1]

	dataset_cis_h2_rep2_file = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_ldsc_est_cis_snp_h2.txt'
	gene_id_to_cis_snp_h2_mapping_rep2 = extract_cis_snp_h2_with_ldsc(validation_data_file, gene_info, rsid_to_variant_stdev, dataset_cis_h2_rep2_file)

	tt.write(eqtl_dataset_name + '\t' + str(N_training) + '\t' + str(N_validation) + '\t' + dataset_cis_h2_rep1_file + '\t' + dataset_cis_h2_rep2_file + '\n')
	gene_to_ldsc_cis_h2s[eqtl_dataset_name + ':' + 'replicate2'] = gene_id_to_cis_snp_h2_mapping_rep2
tt.close()
print(global_output_file)


##############################
# Extract cis-snp h2s using previously created greml estimates
##############################
# Create global output file
global_output_file = est_cis_snp_h2_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_greml_est_cis_snp_h2_summary.txt'
tt = open(global_output_file,'w')
tt.write('tissue_name\tN_rep1\tN_rep2\tcis_h2_rep1_file\tcis_h2_rep2_file\n')

# Loop through eqtl data sets
for dataset_iter, eqtl_dataset_name in enumerate(eqtl_dataset_names):

	###################
	# Replicate 1
	###################
	N_training = eqtl_dataset_Ns_training[dataset_iter]
	file_stem = eqtl_dataset_files_training[dataset_iter].split('_eqtl_sumstats')[0].split('/')[-1]
	mesc_file_stem = mesc_expression_score_dir + file_stem
	dataset_cis_h2_rep1_file = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_greml_est_cis_snp_h2.txt'
	extract_cis_snp_h2_with_greml(mesc_file_stem, chrom_arr, gene_to_ldsc_cis_h2s[eqtl_dataset_name + ':' + 'replicate1'], dataset_cis_h2_rep1_file)

	###################
	# Replicate 2
	###################
	N_validation = eqtl_dataset_Ns_validation[dataset_iter]
	file_stem = eqtl_dataset_files_validation[dataset_iter].split('_eqtl_sumstats')[0].split('/')[-1]
	mesc_file_stem = mesc_expression_score_dir + file_stem
	dataset_cis_h2_rep2_file = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_greml_est_cis_snp_h2.txt'
	extract_cis_snp_h2_with_greml(mesc_file_stem, chrom_arr, gene_to_ldsc_cis_h2s[eqtl_dataset_name + ':' + 'replicate2'], dataset_cis_h2_rep2_file)

	tt.write(eqtl_dataset_name + '\t' + str(N_training) + '\t' + str(N_validation) + '\t' + dataset_cis_h2_rep1_file + '\t' + dataset_cis_h2_rep2_file + '\n')
tt.close()
print(global_output_file)



##############################
# Extract cis-snp h2s using posterior greml estimates
##############################
# Create global output file
greml_output_file = est_cis_snp_h2_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_greml_est_cis_snp_h2_summary.txt'

global_output_file = est_cis_snp_h2_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_posterior_greml_est_cis_snp_h2_summary.txt'
tt = open(global_output_file,'w')
tt.write('tissue_name\tN_rep1\tN_rep2\tcis_h2_rep1_file\tcis_h2_rep2_file\n')

print(greml_output_file)
f = open(greml_output_file)
head_count = 0
dataset_iter = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	eqtl_dataset_name = eqtl_dataset_names[dataset_iter]

	###################
	# Replicate 1
	###################
	N_training = eqtl_dataset_Ns_training[dataset_iter]
	file_stem = eqtl_dataset_files_training[dataset_iter].split('_eqtl_sumstats')[0].split('/')[-1]
	greml_dataset_cis_h2_rep1_file = data[3]
	tmp_output_rep1_stem = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_posterior_greml_est_cis_snp_h2_tmp.txt'
	posterior_dataset_cis_h2_rep1_file = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_posterior_greml_est_cis_snp_h2.txt'
	extract_posterior_cis_snp_h2(greml_dataset_cis_h2_rep1_file, posterior_dataset_cis_h2_rep1_file, tmp_output_rep1_stem)

	###################
	# Replicate 2
	###################
	N_validation = eqtl_dataset_Ns_validation[dataset_iter]
	file_stem = eqtl_dataset_files_validation[dataset_iter].split('_eqtl_sumstats')[0].split('/')[-1]
	greml_dataset_cis_h2_rep2_file = data[4]
	tmp_output_rep2_stem = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_posterior_greml_est_cis_snp_h2_tmp.txt'
	posterior_dataset_cis_h2_rep2_file = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_posterior_greml_est_cis_snp_h2.txt'
	extract_posterior_cis_snp_h2(greml_dataset_cis_h2_rep2_file, posterior_dataset_cis_h2_rep2_file, tmp_output_rep2_stem)

	tt.write(eqtl_dataset_name + '\t' + str(N_training) + '\t' + str(N_validation) + '\t' + posterior_dataset_cis_h2_rep1_file + '\t' + posterior_dataset_cis_h2_rep2_file + '\n')


	dataset_iter = dataset_iter + 1
f.close()
tt.close()
print(global_output_file)
'''


##############################
# Extract cis-snp h2s using avg_chi_sq
##############################
'''
# Create global output file
global_output_file = est_cis_snp_h2_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_avgChisq_est_cis_snp_h2_summary.txt'
tt = open(global_output_file,'w')
tt.write('tissue_name\tN_rep1\tN_rep2\tcis_h2_rep1_file\tcis_h2_rep2_file\n')

# Loop through eqtl data sets
gene_to_ldsc_cis_h2s = {}
for dataset_iter, eqtl_dataset_name in enumerate(eqtl_dataset_names):
	###################
	# Replicate 1
	###################
	N_training = eqtl_dataset_Ns_training[dataset_iter]
	training_data_file = eqtl_dataset_files_training[dataset_iter]

	file_stem = eqtl_dataset_files_training[dataset_iter].split('_eqtl_sumstats')[0].split('/')[-1]
	dataset_cis_h2_rep1_file = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_avgChisq_est_cis_snp_h2.txt'
	extract_cis_snp_h2_with_avgChisq(training_data_file, gene_info, rsid_to_variant_stdev, dataset_cis_h2_rep1_file)

	###################
	# Replicate 2
	###################
	N_validation = eqtl_dataset_Ns_validation[dataset_iter]
	validation_data_file = eqtl_dataset_files_validation[dataset_iter]

	file_stem = eqtl_dataset_files_validation[dataset_iter].split('_eqtl_sumstats')[0].split('/')[-1]

	dataset_cis_h2_rep2_file = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_avgChisq_est_cis_snp_h2.txt'
	extract_cis_snp_h2_with_avgChisq(validation_data_file, gene_info, rsid_to_variant_stdev, dataset_cis_h2_rep2_file)

	tt.write(eqtl_dataset_name + '\t' + str(N_training) + '\t' + str(N_validation) + '\t' + dataset_cis_h2_rep1_file + '\t' + dataset_cis_h2_rep2_file + '\n')
tt.close()
print(global_output_file)
'''

'''
##############################
# Extract cis-snp h2s using posterior avg_chi_sq estimates
##############################
# Create global output file
greml_output_file = est_cis_snp_h2_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_avgChisq_est_cis_snp_h2_summary.txt'

global_output_file = est_cis_snp_h2_dir + simulation_name_string + '_' + str(eqtl_sample_size) + '_posterior_avgChisq_est_cis_snp_h2_summary.txt'
tt = open(global_output_file,'w')
tt.write('tissue_name\tN_rep1\tN_rep2\tcis_h2_rep1_file\tcis_h2_rep2_file\n')

print(greml_output_file)
f = open(greml_output_file)
head_count = 0
dataset_iter = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	eqtl_dataset_name = eqtl_dataset_names[dataset_iter]

	###################
	# Replicate 1
	###################
	N_training = eqtl_dataset_Ns_training[dataset_iter]
	file_stem = eqtl_dataset_files_training[dataset_iter].split('_eqtl_sumstats')[0].split('/')[-1]
	greml_dataset_cis_h2_rep1_file = data[3]
	tmp_output_rep1_stem = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_posterior_avgChisq_est_cis_snp_h2_tmp.txt'
	posterior_dataset_cis_h2_rep1_file = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_posterior_avgChisq_est_cis_snp_h2.txt'
	extract_posterior_cis_snp_h2(greml_dataset_cis_h2_rep1_file, posterior_dataset_cis_h2_rep1_file, tmp_output_rep1_stem)

	###################
	# Replicate 2
	###################
	N_validation = eqtl_dataset_Ns_validation[dataset_iter]
	file_stem = eqtl_dataset_files_validation[dataset_iter].split('_eqtl_sumstats')[0].split('/')[-1]
	greml_dataset_cis_h2_rep2_file = data[4]
	tmp_output_rep2_stem = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_posterior_avgChisq_est_cis_snp_h2_tmp.txt'
	posterior_dataset_cis_h2_rep2_file = est_cis_snp_h2_dir + str(eqtl_sample_size) + '_' + file_stem + '_posterior_avgChisq_est_cis_snp_h2.txt'
	extract_posterior_cis_snp_h2(greml_dataset_cis_h2_rep2_file, posterior_dataset_cis_h2_rep2_file, tmp_output_rep2_stem)

	tt.write(eqtl_dataset_name + '\t' + str(N_training) + '\t' + str(N_validation) + '\t' + posterior_dataset_cis_h2_rep1_file + '\t' + posterior_dataset_cis_h2_rep2_file + '\n')

	dataset_iter = dataset_iter + 1
f.close()
tt.close()
print(global_output_file)


'''








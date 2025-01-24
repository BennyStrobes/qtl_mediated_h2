import sys
import numpy as np 
import os
import pdb
import statsmodels.api as sm
from bgen import BgenReader

def mean_impute_and_standardize_genotype(G_obj_geno):
	# Fill in missing values
	G_obj_geno_stand = np.copy(G_obj_geno)
	ncol = G_obj_geno_stand.shape[1]
	n_missing = []
	for col_iter in range(ncol):
		nan_indices = np.isnan(G_obj_geno[:,col_iter])
		non_nan_mean = np.mean(G_obj_geno[nan_indices==False, col_iter])
		G_obj_geno_stand[nan_indices, col_iter] = non_nan_mean
		n_missing.append(np.sum(nan_indices))
	n_missing = np.asarray(n_missing)



	# Quick error check
	if np.sum(np.std(G_obj_geno_stand,axis=0) == 0) > 0:
		print('no variance genotype assumption error')
		pdb.set_trace()

	# Standardize genotype
	G_obj_geno_stand3 = (G_obj_geno_stand -np.mean(G_obj_geno_stand,axis=0))/np.std(G_obj_geno_stand,axis=0)


	return G_obj_geno_stand3


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


def construct_variant_ld_mat_based_on_distance(genotype_stem, output_stem, window_size=1):
	# Load in ref-alt alleles
	ref_alt_alleles = load_in_ref_alt_allele_arr(genotype_stem + '.pvar')

	# Load in genotype object
	bfile = BgenReader(genotype_stem + '.bgen')

	snp_pos = np.asarray(bfile.positions())
	snp_rsids = np.asarray(bfile.rsids())

	inner_window_start_pos = np.min(snp_pos) - 1

	mb = 1000000.0

	used_snps = {}

	output_file = output_stem + '_summary.txt'
	t = open(output_file,'w')
	t.write('window_name\twindow_start\twindow_end\tld_file\twindow_snp_position_file\trsid_file\tmiddle_boolean_file\n')

	window_counter = 0

	while inner_window_start_pos < np.max(snp_pos):

		print(window_counter)
		window_counter = window_counter + 1

		# Define windows
		inner_window_end_pos = inner_window_start_pos + window_size*mb
		outer_window_start_pos = inner_window_start_pos - (window_size*mb)
		outer_window_end_pos = inner_window_end_pos + (window_size*mb)

		window_position_indices = (snp_pos >= outer_window_start_pos) & (snp_pos < outer_window_end_pos)

		if np.sum(window_position_indices) == 0:
			inner_window_start_pos = inner_window_start_pos + (window_size*mb)
			continue

		# Extract genotype data for this window
		window_indices = np.arange(len(window_position_indices))[window_position_indices].astype(int)
		genotype_dosage = load_in_alt_allele_genotype_dosage_mat(bfile, window_indices, ref_alt_alleles)

		# Standardize genotype dosage matrix
		std_genotype_dosage = standardize_genotype_dosage_matrix(genotype_dosage)

		# Compute LD matrix for window
		ref_ss = std_genotype_dosage.shape[0]
		ld_mat = np.corrcoef(np.transpose(std_genotype_dosage))

		# Get into window orientation
		window_snp_pos = snp_pos[window_position_indices]
		window_rsids = snp_rsids[window_position_indices]

		middle_rsids_bool = []
		# Loop through snps until we find a middle snp
		for tmp_snp_iter, tmp_snp_pos in enumerate(window_snp_pos):
			if tmp_snp_pos >= inner_window_start_pos and tmp_snp_pos < inner_window_end_pos:
				# This is a middle snp
				middle_rsids_bool.append(1.0)
			else:
				middle_rsids_bool.append(0.0)
		middle_rsids_bool = np.asarray(middle_rsids_bool)
		print(len(middle_rsids_bool))

		window_global_positions = np.where(window_position_indices==True)[0]


		# Start saving things to output file
		window_name = 'window_' + str(int(outer_window_start_pos)) + '_' + str(int(outer_window_end_pos))

		# Save LD mat to output
		ld_mat_file = output_stem + '_' + window_name + '_ld.npy'
		np.save(ld_mat_file, ld_mat)

		# Save rsids to outtput file
		rsid_file = output_stem + '_' + window_name + '_rsids.npy'
		np.save(rsid_file, window_rsids)

		# Save middle snp boolean to output file
		middle_boolean_file = output_stem + '_' + window_name + '_middle_boolean.npy'
		np.save(middle_boolean_file, middle_rsids_bool)

		# Save window positions to output
		window_position_file = output_stem + '_' + window_name + '_window_positions.npy'
		np.save(window_position_file, window_global_positions)


		t.write(window_name + '\t' + str(int(outer_window_start_pos)) + '\t' + str(int(outer_window_end_pos)) + '\t' + ld_mat_file + '\t' + window_position_file + '\t' + rsid_file + '\t' + middle_boolean_file + '\n')	
		inner_window_start_pos = inner_window_start_pos + (window_size*mb)

	t.close()
	return

def compute_lambda_thresh(lambdas, rho_thresh):
	totaler = np.sum(lambdas)
	cur_total = 0
	lambda_thresh = -1
	for lambda_val in -np.sort(-lambdas):
		cur_total = cur_total + lambda_val
		if cur_total/totaler > rho_thresh:
			if lambda_thresh == -1:
				lambda_thresh = lambda_val


	if lambda_thresh == -1:
		print('assumption eroror')
		pdb.set_trace()

	return lambda_thresh



def construct_variant_ld_mat_based_on_quasi_windows(genotype_stem, output_stem, quasi_windows):
	# Load in ref-alt alleles
	ref_alt_alleles = load_in_ref_alt_allele_arr(genotype_stem + '.pvar')

	# Load in genotype object
	bfile = BgenReader(genotype_stem + '.bgen')

	snp_pos = np.asarray(bfile.positions())
	snp_rsids = np.asarray(bfile.rsids())


	used_snps = {}

	output_file = output_stem + '_summary.txt'
	t = open(output_file,'w')
	t.write('window_name\twindow_start\twindow_end\tld_file\twindow_snp_position_file\trsid_file\tmiddle_boolean_file\tQ_mat_file\tw_premult_file\n')

	window_counter = 0

	for quasi_window in quasi_windows:

		print(window_counter)
		window_counter = window_counter + 1

		# Define windows
		outer_window_start_pos = quasi_window[0]
		outer_window_end_pos = quasi_window[1]

		window_position_indices = (snp_pos >= outer_window_start_pos) & (snp_pos < outer_window_end_pos)

		if np.sum(window_position_indices) == 0:
			continue

		# Extract genotype data for this window
		window_indices = np.arange(len(window_position_indices))[window_position_indices].astype(int)
		genotype_dosage = load_in_alt_allele_genotype_dosage_mat(bfile, window_indices, ref_alt_alleles)

		# Standardize genotype dosage matrix
		std_genotype_dosage = standardize_genotype_dosage_matrix(genotype_dosage)

		# Compute LD matrix for window
		ref_ss = std_genotype_dosage.shape[0]
		print(std_genotype_dosage.shape)
		ld_mat = np.corrcoef(np.transpose(std_genotype_dosage))

		# Get into window orientation
		window_snp_pos = snp_pos[window_position_indices]
		window_rsids = snp_rsids[window_position_indices]

		middle_rsids_bool = []
		# Loop through snps until we find a middle snp
		for tmp_snp_iter, tmp_snp_pos in enumerate(window_snp_pos):
			middle_rsids_bool.append(1.0)
		middle_rsids_bool = np.asarray(middle_rsids_bool)

		window_global_positions = np.where(window_position_indices==True)[0]


		# EIG value decomp
		lambdas_full, U_full = np.linalg.eig(ld_mat)
		non_negative_components = lambdas_full > 0.0
		lambdas = lambdas_full[non_negative_components]
		U = U_full[:, non_negative_components]
		real_components = np.iscomplex(lambdas) == False
		lambdas = lambdas[real_components]
		U = U[:, real_components]
		if np.sum(np.iscomplex(lambdas)) > 0:
			print('assumption eroror')
			pdb.set_trace()
		lambdas = lambdas.astype(float)
		U = U.astype(float)

		rho_thresh = 0.99
		lambda_thresh = compute_lambda_thresh(lambdas, rho_thresh)
		thresh_components = lambdas >= lambda_thresh
		lambdas = lambdas[thresh_components]
		U = U[:, thresh_components]


		# Note that reconstruction of ld_mat is achieved with np.dot(np.dot(U, np.diag(lambdas)), np.transpose(U))

		# Compute some relevent quantities
		Q_mat = np.dot(np.diag(lambdas**(.5)), np.transpose(U))
		w_premult = np.dot(np.diag(lambdas**(-.5)), np.transpose(U))

		# Start saving things to output file
		window_name = 'window_' + str(int(outer_window_start_pos)) + '_' + str(int(outer_window_end_pos))

		# Save LD mat to output
		ld_mat_file = output_stem + '_' + window_name + '_ld.npy'
		np.save(ld_mat_file, ld_mat)

		# Save rsids to outtput file
		rsid_file = output_stem + '_' + window_name + '_rsids.npy'
		np.save(rsid_file, window_rsids)

		# Save middle snp boolean to output file
		middle_boolean_file = output_stem + '_' + window_name + '_middle_boolean.npy'
		np.save(middle_boolean_file, middle_rsids_bool)

		# Save window positions to output
		window_position_file = output_stem + '_' + window_name + '_window_positions.npy'
		np.save(window_position_file, window_global_positions)

		# Save Q_mat to output
		Q_mat_file = output_stem + '_' + window_name + '_Q_mat.npy'
		np.save(Q_mat_file, Q_mat)
		
		# Save w_premult to output
		w_premult_file = output_stem + '_' + window_name + '_w_premult_mat.npy'
		np.save(w_premult_file, w_premult)

		t.write(window_name + '\t' + str(int(outer_window_start_pos)) + '\t' + str(int(outer_window_end_pos)) + '\t' + ld_mat_file + '\t' + window_position_file + '\t' + rsid_file + '\t' + middle_boolean_file + '\t' + Q_mat_file + '\t' + w_premult_file + '\n')	

	t.close()
	return




def get_quasi_windows(chrom1_quasi_independent_ld_blocks_file):
	f = open(chrom1_quasi_independent_ld_blocks_file)
	windows = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		start = int(data[1])
		end = int(data[2])
		if end < start:
			print('assumption eororo')
			pdb.set_trace()
		windows.append((start, end))
	f.close()
	return windows


def get_big_quasi_windows(chrom1_quasi_independent_ld_blocks_file):
	f = open(chrom1_quasi_independent_ld_blocks_file)
	windows = []
	head_count = 0
	count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if count == 0:
			start = int(data[1])
		elif np.mod(count, 30) == 0:
			end = int(data[2])
			windows.append((start, end))
			start = int(data[2])
		tmp_end = int(data[2])
		count = count + 1
	if np.mod(count, 30) != 1:
		windows.append((start, tmp_end))
	f.close()
	return windows



##############################
# Command line argumemnts
##############################
processed_genotype_data_dir = sys.argv[1]
quasi_independent_dir = sys.argv[2]

window_size = 1

eqtl_sample_sizes = ['100','300', '500', '1000', '10000']

chrom1_quasi_independent_ld_blocks_file = quasi_independent_dir + 'EUR/' + 'fourier_ls-chr1.bed'
quasi_windows = get_quasi_windows(chrom1_quasi_independent_ld_blocks_file)


in_sample_genotype_stem = processed_genotype_data_dir +'simulated_gwas_data_1'
output_stem = processed_genotype_data_dir + 'variant_ref_geno_gwas_quasi_independent_windows_ld'
construct_variant_ld_mat_based_on_quasi_windows(in_sample_genotype_stem, output_stem, quasi_windows)


for eqtl_ss in eqtl_sample_sizes:
	in_sample_genotype_stem = processed_genotype_data_dir +'simulated_eqtl_' + eqtl_ss + '_data_1'
	output_stem = processed_genotype_data_dir + 'variant_ref_geno_eqtl_' + eqtl_ss + '_quasi_independent_windows_ld'
	construct_variant_ld_mat_based_on_quasi_windows(in_sample_genotype_stem, output_stem, quasi_windows)

# Make on full region
starter = quasi_windows[0][0]
ender = quasi_windows[-1][1]
biggest_quasi_windows = [(starter, ender)]

in_sample_genotype_stem = processed_genotype_data_dir +'simulated_gwas_data_1'
output_stem = processed_genotype_data_dir + 'variant_ref_geno_gwas_full_space_ld_ld'
construct_variant_ld_mat_based_on_quasi_windows(in_sample_genotype_stem, output_stem, biggest_quasi_windows)

print('done')

chrom1_quasi_independent_ld_blocks_file = quasi_independent_dir + 'EUR/' + 'fourier_ls-chr1.bed'
big_quasi_windows = get_big_quasi_windows(chrom1_quasi_independent_ld_blocks_file)

in_sample_genotype_stem = processed_genotype_data_dir +'simulated_gwas_data_1'
output_stem = processed_genotype_data_dir + 'variant_ref_geno_gwas_big_quasi_independent_windows_ld'
construct_variant_ld_mat_based_on_quasi_windows(in_sample_genotype_stem, output_stem, big_quasi_windows)

for eqtl_ss in eqtl_sample_sizes:
	in_sample_genotype_stem = processed_genotype_data_dir +'simulated_eqtl_' + eqtl_ss + '_data_1'
	output_stem = processed_genotype_data_dir + 'variant_ref_geno_eqtl_' + eqtl_ss + '_big_quasi_independent_windows_ld'
	construct_variant_ld_mat_based_on_quasi_windows(in_sample_genotype_stem, output_stem, big_quasi_windows)


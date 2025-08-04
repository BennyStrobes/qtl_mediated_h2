import numpy as np
import os
import sys
import pdb
import gzip


def make_linear_G_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_G_file_stem, per_tissue_sample_sizes):
	for chrom_num in chrom_nums:
		chrom_G_file = new_G_file_stem + '.' + str(chrom_num) + '.G'
		gene_counts = []
		for tissue_num in tissue_nums:
			orig_G_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_num) + '_' + str(per_tissue_sample_sizes[tissue_num]) + '_' + sample_split + '_' + str(chrom_num) + '.' + str(chrom_num) + '.G'
			f = open(orig_G_file)
			counter = 0
			for line in f:
				line = line.rstrip()
				data = line.split(' ')
				tissue_gene_counts = np.asarray(data).astype(int)
				if counter > 0:
					print('assumption eroror')
					pdb.set_trace()
				counter = counter + 1
			f.close()
			gene_counts.append(np.sum(tissue_gene_counts))
		gene_counts = np.asarray(gene_counts)

		# Print to new per-chromosome output
		t = open(chrom_G_file,'w')
		t.write(' '.join(gene_counts.astype(str)) + '\n')
		t.close()
	return


def make_stdExpr_G_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_G_file_stem, per_tissue_sample_sizes):
	for chrom_num in chrom_nums:
		chrom_G_file = new_G_file_stem + '.' + str(chrom_num) + '.G'
		gene_counts = []
		for tissue_num in tissue_nums:
			orig_G_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_num) + '_' + str(per_tissue_sample_sizes[tissue_num]) + '_' + sample_split + '_' + str(chrom_num) + '.' + str(chrom_num) + '.G'
			f = open(orig_G_file)
			counter = 0
			for line in f:
				line = line.rstrip()
				data = line.split(' ')
				tissue_gene_counts = np.asarray(data).astype(int)
				if counter > 0:
					print('assumption eroror')
					pdb.set_trace()
				counter = counter + 1
			f.close()
			gene_counts.append(tissue_gene_counts)
		gene_counts = np.hstack(gene_counts)

		# Print to new per-chromosome output
		t = open(chrom_G_file,'w')
		t.write(' '.join(gene_counts.astype(str)) + '\n')
		t.close()
	return



def make_linear_ave_h2cis_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_ave_h2_cis_file_stem, per_tissue_sample_sizes):
	for chrom_num in chrom_nums:
		h2_ciss = []
		for tissue_num in tissue_nums:
			# Get G in this tissue
			orig_G_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_num) + '_' + str(per_tissue_sample_sizes[tissue_num]) + '_' + sample_split + '_' + str(chrom_num) + '.' + str(chrom_num) + '.G'
			f = open(orig_G_file)
			counter = 0
			for line in f:
				line = line.rstrip()
				data = line.split(' ')
				tissue_gene_counts = np.asarray(data).astype(int)
				if counter > 0:
					print('assumption eroror')
					pdb.set_trace()
				counter = counter + 1
			f.close()
			# Get cis h2 in this tissue
			orig_h2_cis_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_num) + '_' + str(per_tissue_sample_sizes[tissue_num]) + '_' + sample_split + '_' + str(chrom_num) + '.' + str(chrom_num) + '.ave_h2cis'
			f = open(orig_h2_cis_file)
			counter = 0
			for line in f:
				line = line.rstrip()
				data = line.split(' ')
				tissue_cis_h2 = np.asarray(data).astype(float)
				if counter > 0:
					print('assumption eroror')
					pdb.set_trace()
				counter = counter + 1
			f.close()
			h2_ciss.append(np.dot(tissue_gene_counts, tissue_cis_h2)/np.sum(tissue_gene_counts))
		h2_ciss = np.asarray(h2_ciss)
		chrom_cish2_file = new_ave_h2_cis_file_stem + '.' + str(chrom_num) + '.ave_h2cis'
		t = open(chrom_cish2_file,'w')
		t.write(' '.join(h2_ciss.astype(str)) + '\n')
		t.close()		
	return

def make_stdExpr_ave_h2cis_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_ave_h2_cis_file_stem, per_tissue_sample_sizes):
	for chrom_num in chrom_nums:
		h2_ciss = []
		for tissue_num in tissue_nums:
			# Get G in this tissue
			orig_G_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_num) + '_' + str(per_tissue_sample_sizes[tissue_num]) + '_' + sample_split + '_' + str(chrom_num) + '.' + str(chrom_num) + '.G'
			f = open(orig_G_file)
			counter = 0
			for line in f:
				line = line.rstrip()
				data = line.split(' ')
				tissue_gene_counts = np.asarray(data).astype(int)
				if counter > 0:
					print('assumption eroror')
					pdb.set_trace()
				counter = counter + 1
			f.close()
			# Get cis h2 in this tissue
			orig_h2_cis_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_num) + '_' + str(per_tissue_sample_sizes[tissue_num]) + '_' + sample_split + '_' + str(chrom_num) + '.' + str(chrom_num) + '.ave_h2cis'
			f = open(orig_h2_cis_file)
			counter = 0
			for line in f:
				line = line.rstrip()
				data = line.split(' ')
				tissue_cis_h2 = np.asarray(data).astype(float)
				if counter > 0:
					print('assumption eroror')
					pdb.set_trace()
				counter = counter + 1
			f.close()
			#h2_ciss.append(np.dot(tissue_gene_counts, tissue_cis_h2)/np.sum(tissue_gene_counts))
			h2_ciss.append(tissue_cis_h2)
		h2_ciss = np.hstack(h2_ciss)
		chrom_cish2_file = new_ave_h2_cis_file_stem + '.' + str(chrom_num) + '.ave_h2cis'
		t = open(chrom_cish2_file,'w')
		t.write(' '.join(h2_ciss.astype(str)) + '\n')
		t.close()		
	return


def make_linear_gannot_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_gannot_file_stem, per_tissue_sample_sizes):
	for chrom_num in chrom_nums:

		chrom_gannot_file = new_gannot_file_stem + '.' + str(chrom_num) + '.gannot.gz'
		t = gzip.open(chrom_gannot_file, 'wt')
		t.write('Gene')
		for tissue_num in tissue_nums:
			t.write('\t' + 'Cis_herit_bin_' + str(tissue_num+1))
		t.write('\n')

		for tissue_num in tissue_nums:
			# Get G in this tissue
			orig_gannot_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_num) + '_' + str(per_tissue_sample_sizes[tissue_num]) + '_' + sample_split + '_' + str(chrom_num) + '.' + str(chrom_num) + '.gannot.gz'
			f = gzip.open(orig_gannot_file)
			head_count = 0
			for line in f:
				line = line.decode('utf-8').strip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				gene_id = data[0] + '_' + 'tiss' + str(tissue_num)
				bin_assignment_vec = np.zeros(len(tissue_nums)).astype(int)
				bin_assignment_vec[tissue_num] = 1
				t.write(gene_id + '\t' + '\t'.join(bin_assignment_vec.astype(str)) + '\n')
			f.close()
		t.close()
	return


def make_stdExpr_gannot_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_gannot_file_stem, per_tissue_sample_sizes):
	for chrom_num in chrom_nums:

		chrom_gannot_file = new_gannot_file_stem + '.' + str(chrom_num) + '.gannot.gz'
		t = gzip.open(chrom_gannot_file, 'wt')
		t.write('Gene')
		col_counter = 1
		for tissue_num in tissue_nums:
			for bin_iter in range(5):
				t.write('\t' + 'Cis_herit_bin_' + str(col_counter))
				col_counter = col_counter + 1
		t.write('\n')

		for tissue_num in tissue_nums:
			# Get G in this tissue
			orig_gannot_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_num) + '_' + str(per_tissue_sample_sizes[tissue_num]) + '_' + sample_split + '_' + str(chrom_num) + '.' + str(chrom_num) + '.gannot.gz'
			f = gzip.open(orig_gannot_file)
			head_count = 0
			for line in f:
				line = line.decode('utf-8').strip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				gene_id = data[0] + '_' + 'tiss' + str(tissue_num)
				bin_assignment_vec = np.zeros(5*len(tissue_nums)).astype(int)
				bin_assignment_vec[tissue_num*5:(tissue_num*5+5)] = np.asarray(data[1:]).astype(int)
				t.write(gene_id + '\t' + '\t'.join(bin_assignment_vec.astype(str)) + '\n')
			f.close()
		t.close()
	return


def make_linear_expscore_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_expscore_file_stem, per_tissue_sample_sizes):
	for chrom_num in chrom_nums:
		xt_exp_scores = []
		tissue_names = []
		n_snps = []
		for tissue_num in tissue_nums:
			tissue_exp_scores = []
			orig_expscore_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_num) + '_' + str(per_tissue_sample_sizes[tissue_num]) + '_' + sample_split + '_' + str(chrom_num) + '.' + str(chrom_num) + '.expscore.gz'
			f = gzip.open(orig_expscore_file)
			head_count = 0
			for line in f:
				line = line.decode('utf-8').rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				tissue_exp_scores.append(data)
			f.close()
			tissue_exp_scores = np.asarray(tissue_exp_scores)
			xt_exp_scores.append(tissue_exp_scores)
			tissue_names.append('Cis_herit_bin_' + str(tissue_num+1))
			n_snps.append(tissue_exp_scores.shape[0])
		tissue_names = np.asarray(tissue_names)
		n_snps = np.asarray(n_snps)
		if len(np.unique(n_snps)) != 1:
			print('assumptione rororor')
			pdb.set_trace()

		# Now print to global output
		chrom_expscores_file = new_expscore_file_stem + '.' + str(chrom_num) + '.expscore.gz'
		t = gzip.open(chrom_expscores_file,'wt')
		# Print header 
		t.write('CHR\tSNP\tBP\t' + '\t'.join(tissue_names) + '\n')
		for snp_iter in range(n_snps[0]):
			# Shared columns
			t.write('\t'.join(xt_exp_scores[0][snp_iter,:3]))
			for tissue_num in tissue_nums:
				tmp_exp_score_vec = xt_exp_scores[tissue_num][snp_iter][3:]
				tissue_score = np.sum(tmp_exp_score_vec.astype(float))
				t.write('\t' + str(tissue_score))
			t.write('\n')
		t.close()
	return

def make_stdExpr_expscore_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_expscore_file_stem, per_tissue_sample_sizes):
	for chrom_num in chrom_nums:
		xt_exp_scores = []
		tissue_names = []
		n_snps = []
		col_counter = 1
		for tissue_num in tissue_nums:
			tissue_exp_scores = []
			orig_expscore_file = mesc_expression_score_dir + simulation_name_string + '_tissue' + str(tissue_num) + '_' + str(per_tissue_sample_sizes[tissue_num]) + '_' + sample_split + '_' + str(chrom_num) + '.' + str(chrom_num) + '.expscore.gz'
			f = gzip.open(orig_expscore_file)
			head_count = 0
			for line in f:
				line = line.decode('utf-8').rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				tissue_exp_scores.append(data)
			f.close()
			tissue_exp_scores = np.asarray(tissue_exp_scores)
			xt_exp_scores.append(tissue_exp_scores)
			for ii in range(tissue_exp_scores[:,3:].shape[1]):
				tissue_names.append('Cis_herit_bin_' + str(col_counter))
				col_counter = col_counter + 1
			n_snps.append(tissue_exp_scores.shape[0])
		tissue_names = np.asarray(tissue_names)
		n_snps = np.asarray(n_snps)
		if len(np.unique(n_snps)) != 1:
			print('assumptione rororor')
			pdb.set_trace()

		# Now print to global output
		chrom_expscores_file = new_expscore_file_stem + '.' + str(chrom_num) + '.expscore.gz'
		t = gzip.open(chrom_expscores_file,'wt')
		# Print header 
		t.write('CHR\tSNP\tBP\t' + '\t'.join(tissue_names) + '\n')
		for snp_iter in range(n_snps[0]):
			# Shared columns
			t.write('\t'.join(xt_exp_scores[0][snp_iter,:3]))
			for tissue_num in tissue_nums:
				tmp_exp_score_vec = xt_exp_scores[tissue_num][snp_iter][3:]
				t.write('\t' + '\t'.join(tmp_exp_score_vec))
			t.write('\n')
		t.close()
	return



def preprocess_linear_expression_scores(simulation_number, simulation_name_string, simulation_genotype_dir, mesc_expression_score_dir, eqtl_sample_size, sample_split, mesc_processed_input_dir, chrom_nums, tissue_nums, per_tissue_sample_sizes):
	# Make G file
	new_G_file_stem = mesc_processed_input_dir + simulation_name_string + '_linear_' + eqtl_sample_size + '_' + sample_split 
	make_linear_G_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_G_file_stem, per_tissue_sample_sizes)

	# Make .ave_h2_cis file
	new_ave_h2_cis_file_stem = mesc_processed_input_dir + simulation_name_string + '_linear_' + eqtl_sample_size + '_' + sample_split 
	make_linear_ave_h2cis_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_ave_h2_cis_file_stem, per_tissue_sample_sizes)

	# Make .expscore.gz
	new_expscore_file_stem = mesc_processed_input_dir + simulation_name_string + '_linear_' + eqtl_sample_size + '_' + sample_split 
	make_linear_expscore_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_expscore_file_stem, per_tissue_sample_sizes)

	# Make .gannot file
	new_gannot_file_stem = mesc_processed_input_dir + simulation_name_string + '_linear_' + eqtl_sample_size + '_' + sample_split 
	make_linear_gannot_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_gannot_file_stem, per_tissue_sample_sizes)

	return

def preprocess_stdExpr_expression_scores(simulation_number, simulation_name_string, simulation_genotype_dir, mesc_expression_score_dir, eqtl_sample_size, sample_split, mesc_processed_input_dir, chrom_nums, tissue_nums, per_tissue_sample_sizes):
	# Make G file
	new_G_file_stem = mesc_processed_input_dir + simulation_name_string + '_stdExpr_' + eqtl_sample_size + '_' + sample_split 
	make_stdExpr_G_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_G_file_stem, per_tissue_sample_sizes)

	# Make .ave_h2_cis file
	new_ave_h2_cis_file_stem = mesc_processed_input_dir + simulation_name_string + '_stdExpr_' + eqtl_sample_size + '_' + sample_split 
	make_stdExpr_ave_h2cis_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_ave_h2_cis_file_stem, per_tissue_sample_sizes)

	# Make .expscore.gz
	new_expscore_file_stem = mesc_processed_input_dir + simulation_name_string + '_stdExpr_' + eqtl_sample_size + '_' + sample_split 
	make_stdExpr_expscore_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_expscore_file_stem, per_tissue_sample_sizes)


	# Make .gannot file
	new_gannot_file_stem = mesc_processed_input_dir + simulation_name_string + '_stdExpr_' + eqtl_sample_size + '_' + sample_split 
	make_stdExpr_gannot_file(mesc_expression_score_dir, simulation_name_string, eqtl_sample_size, sample_split, chrom_nums, tissue_nums, new_gannot_file_stem, per_tissue_sample_sizes)

	return



#####################
# Command line arguments
#####################
simulation_number = sys.argv[1]
chrom_string = sys.argv[2]
simulation_name_string = sys.argv[3]
simulation_genotype_dir = sys.argv[4]
mesc_expression_score_dir = sys.argv[5]
eqtl_sample_size = sys.argv[6]
sample_split = sys.argv[7]
mesc_processed_input_dir = sys.argv[8]
gt_arch = sys.argv[9]


if chrom_string == '1_2':
	chrom_nums = np.arange(1,3)

if eqtl_sample_size == '100-1000':
	per_tissue_sample_sizes = np.asarray(['1000']*5)
	per_tissue_sample_sizes[0] = '100'
else:
	per_tissue_sample_sizes = np.asarray([eqtl_sample_size]*5)


tissue_nums = np.arange(5)

if gt_arch == 'linear':
	preprocess_linear_expression_scores(simulation_number, simulation_name_string, simulation_genotype_dir, mesc_expression_score_dir, eqtl_sample_size, sample_split, mesc_processed_input_dir, chrom_nums, tissue_nums, per_tissue_sample_sizes)
elif gt_arch == 'stdExpr':
	preprocess_stdExpr_expression_scores(simulation_number, simulation_name_string, simulation_genotype_dir, mesc_expression_score_dir, eqtl_sample_size, sample_split, mesc_processed_input_dir, chrom_nums, tissue_nums, per_tissue_sample_sizes)


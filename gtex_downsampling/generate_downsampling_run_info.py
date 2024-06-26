import numpy as np
import os
import sys
import pdb




def extract_all_samples_and_tissue_names(tissue_info_file):
	f = open(tissue_info_file)
	head_count = 0
	tissue_names = []
	sample_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_names.append(data[0])
		tissue_sample_names = np.loadtxt(data[2], dtype=str)
		if len(tissue_sample_names) != int(data[1]):
			print('assumption eroror')
			pdb.set_trace()
		for sample_name in tissue_sample_names:
			sample_names.append(data[0] + ':' + sample_name)
	f.close()
	return np.asarray(sample_names), np.asarray(tissue_names)

def downsample_sample_names(downsampled_names_dicti, tissue_name, full_sample_names_file, downsampled_sample_names_file):
	f3 = open(full_sample_names_file)
	t3 = open(downsampled_sample_names_file,'w')
	n_tissue_samples = 0
	for line in f3:
		indi = line.rstrip()
		sample_name = tissue_name + ':' + indi
		if sample_name in downsampled_names_dicti:
			t3.write(indi + '\n')
			n_tissue_samples = n_tissue_samples + 1
	f3.close()
	t3.close()
	return n_tissue_samples

def make_tissue_info_file_for_downsampled_run(downsampled_names_arr, tissue_info_file, downsampled_tissue_info_file, run_name, downsampling_info_dir):
	# First convert downsampled_names_arr to downsampled_names_dicti
	downsampled_names_dicti = {}
	for sample_name in downsampled_names_arr:
		downsampled_names_dicti[sample_name] = 1

	# Quick error check
	if len(downsampled_names_dicti) != len(downsampled_names_arr):
		print('asssumption eroror')
		pdb.set_trace()


	f = open(tissue_info_file)
	t2 = open(downsampled_tissue_info_file,'w')

	# Stream full data tissue info file
	head_count = 0 # for header
	NN = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t2.write(line + '\n')
			continue

		# Extract relevent fields
		tissue_name = data[0]
		full_sample_names_file = data[2]

		# Create sample names file for this downsampled tissue
		downsampled_sample_names_file = downsampling_info_dir + run_name + '_' + tissue_name + '_EUR_sample_names.txt'
		n_tissue_samples = downsample_sample_names(downsampled_names_dicti, tissue_name, full_sample_names_file, downsampled_sample_names_file)
		
		# print to tissue info file
		t2.write(tissue_name + '\t' + str(n_tissue_samples) + '\t' + downsampled_sample_names_file + '\n')

		# Keep track of total samples for error checking
		NN = NN + n_tissue_samples

	f.close()
	t2.close()

	if NN != len(downsampled_names_dicti):
		print('assumption eroror')
		pdb.set_trace()

	return

def make_donor_id_file_for_downsampled_run(downsampled_names_arr, donor_id_file):
	donor_ids = []
	for name in downsampled_names_arr:
		donor_ids.append(name.split(':')[1])
	donor_ids = np.unique(donor_ids)
	t4 = open(donor_id_file,'w')
	for donor_id in donor_ids:
		t4.write(donor_id + '\n')
	t4.close()
	return




######################
# Command line args
######################
tissue_info_file = sys.argv[1]
downsampling_info_dir = sys.argv[2]
downsample_run_summary_file = sys.argv[3]

# Set random seed
np.random.seed(1)

# Extract array containing list of all samples (where a sample is a tissue:indi_id pair)
# also extract list of tissues
full_data_all_samples, tissue_names = extract_all_samples_and_tissue_names(tissue_info_file)

# Decide grid subsampling fractions
subsampling_fractions = np.arange(.5, .99,.02)

# Total number of samples
total_n_samples = len(full_data_all_samples)

# Open downsample_run_summary_file handle
# It has a row for each downsampling run (ie. subsampling fraction)
t = open(downsample_run_summary_file, 'w')
# Header
t.write('run_name\tdownsample_sample_size\ttissue_info_file\n')

# Loop through subsampling fractions
for subsampling_fraction in subsampling_fractions:

	# Get number of samples in this downsampling run
	n_downsampled_samples = int(np.floor(total_n_samples*subsampling_fraction))

	# Name of downsampling run
	run_name = 'downsample_' + str(n_downsampled_samples)

	# Randomly select samples
	downsampled_names_arr = np.sort(np.random.choice(full_data_all_samples, size=n_downsampled_samples, replace=False))

	# Quick error check
	if len(downsampled_names_arr) != n_downsampled_samples:
		print('assumption eroror')
		pdb.set_trace()


	# Make tissue info file for this downsampled run
	downsampled_tissue_info_file = downsampling_info_dir + run_name + '_tissue_info.txt'
	make_tissue_info_file_for_downsampled_run(downsampled_names_arr, tissue_info_file, downsampled_tissue_info_file, run_name,downsampling_info_dir)

	# Make Donor id file for this downsampled run
	donor_id_file = downsampling_info_dir + run_name + '_donor_ids.txt'
	#make_donor_id_file_for_downsampled_run(downsampled_names_arr, donor_id_file)

	# Write to run summary file
	t.write(run_name + '\t' + str(n_downsampled_samples) + '\t' + downsampled_tissue_info_file + '\n')
t.close()







import numpy as np
import os
import sys
import pdb
import gzip



def extract_tissue_names_by_looping_expression_matrices_dir(gtex_v8_normalized_expression_matrices_dir):
	tissue_names = []
	for file_name in os.listdir(gtex_v8_normalized_expression_matrices_dir):
		if file_name.endswith('.gz') == False:
			continue
		tissue_name = file_name.split('.v8.normalized')[0]
		tissue_names.append(tissue_name)
	return np.sort(np.asarray(tissue_names))

def extract_dictionary_list_of_european_ancestry_individuals(gtex_v8_european_list_file):
	dicti = {}
	f = open(gtex_v8_european_list_file)
	for line in f:
		line = line.rstrip()
		dicti[line] = 1
	f.close()

	return dicti




def extract_tissue_sample_names(tissue_expr_file):
	f = gzip.open(tissue_expr_file)

	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		sample_names = np.asarray(data[4:])
		break
	f.close()
	return np.asarray(sample_names)


def print_EUR_tissue_sample_names_to_output_file(tissue_sample_names_all, eur_indi_list, sample_names_output_file):
	t = open(sample_names_output_file,'w')

	n_samp = 0

	for sample_name in tissue_sample_names_all:
		if sample_name in eur_indi_list:
			n_samp = n_samp + 1
			t.write(sample_name + '\n')

	t.close()
	return n_samp



#######################
# Command line args
#######################
gtex_v8_normalized_expression_matrices_dir = sys.argv[1]  # used to get tissue names and sample names in each tissue
gtex_v8_european_list_file = sys.argv[2]  # Used to extract european individuals
tissue_sample_names_dir = sys.argv[3]  # Output directory
tissue_info_file = sys.argv[4]  # Output file with line for each tissue



# Get GTEx tissue names
tissue_names = extract_tissue_names_by_looping_expression_matrices_dir(gtex_v8_normalized_expression_matrices_dir)

# Extract dictionary list of European ancestry individuals
eur_indi_list = extract_dictionary_list_of_european_ancestry_individuals(gtex_v8_european_list_file)


# Create global output file handle for tissue_info_file
t_glob = open(tissue_info_file,'w')
# Header of tissue info file
t_glob.write('Tissue_name\tsample_size\tsample_names_file\n')


# Loop through tissues
for tissue_name in tissue_names:

	# Tissue expression file
	tissue_expr_file = gtex_v8_normalized_expression_matrices_dir + tissue_name + '.v8.normalized_expression.bed.gz'

	# Extract tissue sample names
	tissue_sample_names_all = extract_tissue_sample_names(tissue_expr_file)


	# Print European ancestry tissue sample names to output file
	# Return sample size (restricted to EUR ancestry)
	sample_names_output_file = tissue_sample_names_dir + tissue_name + '_EUR_sample_names.txt'
	n_samp = print_EUR_tissue_sample_names_to_output_file(tissue_sample_names_all, eur_indi_list, sample_names_output_file)


	# Print tissue to tissue_info file
	t_glob.write(tissue_name + '\t' + str(n_samp) + '\t' + sample_names_output_file + '\n')


# Close tissue info file handle
t_glob.close()










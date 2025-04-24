import numpy as np
import os
import sys
import pdb


def extract_names_of_tissues(tissue_info_file):
	f = open(tissue_info_file)
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

	return arr




#######################
# Command line args
#######################
tissue_info_file = sys.argv[1]
gtex_causal_eqtl_effects_dir = sys.argv[2]

# Extract names of tissues
tissue_names = extract_names_of_tissues(tissue_info_file)

# Open Output summary file
output_file = gtex_causal_eqtl_effects_dir + 'GTEx_v8_genotype_EUR_cross_tissue_gene_model_summary.txt'
t = open(output_file,'w')

# Loop through tissues
for tiss_iter, tissue_name in enumerate(tissue_names):
	# Loop through chromosomes
	for chrom_num in range(1,23):
		tiss_chrom_summary_file = gtex_causal_eqtl_effects_dir + 'GTEx_v8_genotype_EUR_' + tissue_name + '/' + 'GTEx_v8_eqtl_' + tissue_name + '_' + str(chrom_num) + '_chr' + str(chrom_num) + '_gene_summary.txt'
		# Error check
		if os.path.exists(tiss_chrom_summary_file) == False:
			print('assumption erororo')
			pdb.set_trace()
		f = open(tiss_chrom_summary_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				# Print header only once
				if tiss_iter == 0 and chrom_num == 1:
					t.write('QTL_class\t' + line + '\n')
				continue
			t.write(tissue_name + '\t' + line + '\n')
		f.close()
t.close()

print(output_file)
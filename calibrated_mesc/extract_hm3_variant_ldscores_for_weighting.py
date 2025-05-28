import numpy as np
import os
import sys
import pdb
from bgen import BgenReader
import argparse
import gzip





def create_dictionary_list_of_hm3_rsids(hm3_rs_id_file):
	dicti = {}
	f = open(hm3_rs_id_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		dicti[line] = 1
	f.close()
	return dicti


def create_dictionary_mapping_from_rsid_to_cm_position(sldsc_annotation_file):
	cm_dicti = {}
	anno_dicti = {}
	head_count = 0
	f = gzip.open(sldsc_annotation_file)
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[2]
		cm = float(data[3])
		if rsid in cm_dicti:
			print('assumption erororr')
			pdb.set_trace()
		cm_dicti[rsid] = cm
	f.close()
	return cm_dicti

def extract_snp_cm_pos(ordered_rsids, rsid_to_cm):
	cm_arr = []
	for rsid in ordered_rsids:
		cm_arr.append(rsid_to_cm[rsid])
	return np.asarray(cm_arr)

def extract_snp_anno(ordered_rsids, rsid_to_anno):
	anno_arr = []
	for rsid in ordered_rsids:
		anno_arr.append(rsid_to_anno[rsid])
	return np.asarray(anno_arr)	

def print_m_vec_to_output(anno_names, m_vec, variant_M_file):
	if len(anno_names) != len(m_vec):
		print('assumtpioneororr')
		pdb.set_trace()
	t = open(variant_M_file,'w')
	t.write('annotation\tn_snps\n')
	for ii, anno_name in enumerate(anno_names):
		t.write(anno_name + '\t' + str(m_vec[ii]) + '\n')
	t.close()
	return

#####################
# Parse command line arguments
#####################
parser = argparse.ArgumentParser()
parser.add_argument('--chrom', default='None', type=str,
                    help='Chromosome number (e.g. 1)')
parser.add_argument('--bgen-file', default='None', type=str,
                    help='Absolute path to reference bgen genotype file')
parser.add_argument('--hm3-rsid-file', default='None', type=str,
                    help='Absolute path to file containing list of hm3 rsids (note: file has no header. just one line per snp)')
parser.add_argument('--cm-position-file', default='None', type=str,
                    help='Absolute path to gziped sldsc annotation file')
parser.add_argument('--variant-ld-score-file', default='None', type=str,
                    help='Absolute path to LD output file')
parser.add_argument('--cm-window-size', default=1.0, type=float,
                    help='size of CM window around genetic elements to use')
args = parser.parse_args()


cm_window_size = args.cm_window_size

# Create dictionary list of hm3rsids
hm3_rsids = create_dictionary_list_of_hm3_rsids(args.hm3_rsid_file)

# Create dictionary mapping from rsid to CM position
rsid_to_cm = create_dictionary_mapping_from_rsid_to_cm_position(args.cm_position_file)

# Load in genotype object
genotype_obj = BgenReader(args.bgen_file)
snp_pos = np.asarray(genotype_obj.positions())
ordered_rsids = np.asarray(genotype_obj.rsids())
snp_cm = extract_snp_cm_pos(ordered_rsids, rsid_to_cm)
snp_integers = np.arange(len(ordered_rsids))



# Open output file handle
t = open(args.variant_ld_score_file,'w')
t.write('chr\trsid\tcm\tpos\tref_allele\talt_allele\t' + 'hm3_weights' + '\n')

mid_window_left = np.min(snp_cm)
mid_window_right = mid_window_left + cm_window_size
left_window_left = mid_window_left - cm_window_size
right_window_right = mid_window_right + cm_window_size

counter = 0
while mid_window_left <= np.max(snp_cm):
	print(counter)
	counter = counter + 1

	# Get all snps in window
	window_snp_indices = (snp_cm >= left_window_left) & (snp_cm <= right_window_right)
	if np.sum(window_snp_indices) == 0:
		# skip window
		# update window pos
		mid_window_left = mid_window_left + cm_window_size
		mid_window_right = mid_window_left + cm_window_size
		left_window_left = mid_window_left - cm_window_size
		right_window_right = mid_window_right + cm_window_size
		continue

	window_pos_tmp = snp_pos[window_snp_indices]
	window_rsids_tmp = ordered_rsids[window_snp_indices]
	window_cm_tmp = snp_cm[window_snp_indices]
	window_snp_integers_tmp = snp_integers[window_snp_indices]  # Used getting dosage

	# Further subset to hm3 snps
	valid_window_indices = []
	for ii, rsid in enumerate(window_rsids_tmp):
		if rsid in hm3_rsids:
			valid_window_indices.append(ii)
	valid_window_indices = np.asarray(valid_window_indices)

	if len(valid_window_indices) == 0:
		print('skipped window with 0 regression snps')
		mid_window_left = mid_window_left + cm_window_size
		mid_window_right = mid_window_left + cm_window_size
		left_window_left = mid_window_left - cm_window_size
		right_window_right = mid_window_right + cm_window_size
		continue

	window_pos = window_pos_tmp[valid_window_indices]
	window_rsids = window_rsids_tmp[valid_window_indices]
	window_cm = window_cm_tmp[valid_window_indices]
	window_snp_integers = window_snp_integers_tmp[valid_window_indices]

	if len(window_pos) <= 1:
		print('small window')
		pdb.set_trace()


	# Construct LD
	window_a0s = []
	window_a1s = []
	window_geno_mat = []
	for snp_integer in window_snp_integers:
		var = genotype_obj[snp_integer]
		dosage = var.alt_dosage # Note, not standardized but ok for ld
		window_geno_mat.append(dosage)
		window_a0s.append(var.alleles[0])
		window_a1s.append(var.alleles[1])  # This is alt allele!!
	window_geno_mat = np.asarray(window_geno_mat)
	geno_ss = window_geno_mat.shape[1]
	LD = np.corrcoef(window_geno_mat)
	squared_LD = np.square(LD)
	squared_adj_LD = squared_LD - ((1.0-squared_LD)/(geno_ss-2.0))
	window_a0s = np.asarray(window_a0s)
	window_a1s = np.asarray(window_a1s)

	# Loop through snps
	# If snp in middle and snp a HM3 snp: construct ld score (will involve subsetting to snps within 1 CM of focal snp)
	for window_snp_iter, window_snp_cm in enumerate(window_cm):
		window_snp_rsid = window_rsids[window_snp_iter]
		window_snp_pos = window_pos[window_snp_iter]
		window_snp_a0 = window_a0s[window_snp_iter]
		window_snp_a1 = window_a1s[window_snp_iter]

		# Disregard non hm3 snps
		if window_snp_rsid not in hm3_rsids:
			print('assumtpioneororr')
			continue
		# Disregard snps not in middle
		if window_snp_cm < mid_window_left:
			continue
		if window_snp_cm >= mid_window_right:
			continue
		# Extract indices of other snps in LD window that are within 1 CM of focal snp
		nearby_snp_indices = np.abs(window_cm - window_snp_cm) <= 1.0

		# Compute LD score
		snp_ldscore = np.sum(squared_adj_LD[window_snp_iter, nearby_snp_indices])
		#snp_anno_ldscore = np.dot(np.transpose(window_snp_anno[nearby_snp_indices,:]),squared_adj_LD[window_snp_iter, nearby_snp_indices])

		# print to output file
		t.write(str(args.chrom) + '\t' + window_snp_rsid + '\t' + str(window_snp_cm) + '\t' + str(window_snp_pos) + '\t' + window_snp_a0 + '\t' + window_snp_a1 + '\t' + str(snp_ldscore) + '\n')
		t.flush()

	# update window pos
	mid_window_left = mid_window_left + cm_window_size
	mid_window_right = mid_window_left + cm_window_size
	left_window_left = mid_window_left - cm_window_size
	right_window_right = mid_window_right + cm_window_size


t.close()









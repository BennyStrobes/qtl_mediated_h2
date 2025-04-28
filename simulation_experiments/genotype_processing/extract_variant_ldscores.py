import numpy as np
import os
import sys
import pdb
from bgen import BgenReader









def create_dictionary_list_of_hm3_rsids(hm3_rs_id_file):
	dicti = {}
	f = open(hm3_rs_id_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		dicti[line] = 1
	f.close()
	return dicti


def create_dictionary_mapping_from_rsid_to_cm_position(bim_file):
	dicti = {}
	f = open(bim_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		rsid =data[1]
		if rsid in dicti:
			print('assumption error')
			pdb.set_trace()
		cm_position = float(data[2])
		dicti[rsid] = cm_position
	f.close()
	return dicti

def extract_snp_cm_pos(ordered_rsids, rsid_to_cm):
	cm_arr = []
	for rsid in ordered_rsids:
		cm_arr.append(rsid_to_cm[rsid])
	return np.asarray(cm_arr)




#####################
# Command line args
#####################
chrom_num = sys.argv[1]
bgen_genotype_stem = sys.argv[2]
hm3_rs_id_file = sys.argv[3]
kg_genotype_dir = sys.argv[4]  # Used to get CM values for snps
variant_ld_score_file = sys.argv[5]  # Output file

# CM window size (impacts memory/speed trast)
cm_window_size = 1.0

# Create dictionary list of hm3rsids
hm3_rsids = create_dictionary_list_of_hm3_rsids(hm3_rs_id_file)
print(len(hm3_rsids))

# Create dictionary mapping from rsid to CM position
rsid_to_cm = create_dictionary_mapping_from_rsid_to_cm_position(kg_genotype_dir + '1000G.EUR.QC.' + str(chrom_num) + '.bim')

# Load in genotype object
genotype_obj = BgenReader(bgen_genotype_stem + '.bgen')
snp_pos = np.asarray(genotype_obj.positions())
ordered_rsids = np.asarray(genotype_obj.rsids())
snp_cm = extract_snp_cm_pos(ordered_rsids, rsid_to_cm)
snp_integers = np.arange(len(ordered_rsids))

# Open output file handle
t = open(variant_ld_score_file,'w')
t.write('chr\trsid\tcm\tpos\ta1\ta2\tldscore\n')

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

	window_pos = snp_pos[window_snp_indices]
	window_rsids = ordered_rsids[window_snp_indices]
	window_cm = snp_cm[window_snp_indices]
	window_snp_integers = snp_integers[window_snp_indices]  # Used getting dosage

	# Construct LD
	window_a0s = []
	window_a1s = []
	window_geno_mat = []
	for snp_integer in window_snp_integers:
		var = genotype_obj[snp_integer]
		dosage = var.minor_allele_dosage # Note, not standardized but ok for ld
		window_geno_mat.append(dosage)
		window_a0s.append(var.alleles[0])
		window_a1s.append(var.alleles[1])
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

		# print to output file
		t.write(str(chrom_num) + '\t' + window_snp_rsid + '\t' + str(window_snp_cm) + '\t' + str(window_snp_pos) + '\t' + window_snp_a0 + '\t' + window_snp_a1 + '\t' + str(snp_ldscore) + '\n')
		t.flush()

	# update window pos
	mid_window_left = mid_window_left + cm_window_size
	mid_window_right = mid_window_left + cm_window_size
	left_window_left = mid_window_left - cm_window_size
	right_window_right = mid_window_right + cm_window_size


t.close()











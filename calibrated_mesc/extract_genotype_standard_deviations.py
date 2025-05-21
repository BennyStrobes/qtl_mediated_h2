import numpy as np
import os
import sys
import pdb
from bgen import BgenReader
import argparse
import gzip



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

def load_in_alt_allele_sdevs(bfile, window_indices, ref_alt_alleles):
	sdevs = []
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

		# Compute standard deviation
		sdevs.append(np.std(dosage))

	return np.asarray(sdevs)

#####################
# Parse command line arguments
#####################
parser = argparse.ArgumentParser()
parser.add_argument('--chrom', default='None', type=str,
                    help='Chromosome number (e.g. 1)')
parser.add_argument('--bgen-file', default='None', type=str,
                    help='Absolute path to reference bgen genotype file')
parser.add_argument('--output-file', default='None', type=str,
                    help='Output file')
args = parser.parse_args()




ref_genotype_obj = BgenReader(args.bgen_file)
ref_G_obj_rsids = np.asarray(ref_genotype_obj.rsids())
ref_G_obj_positions = np.asarray(ref_genotype_obj.positions())
ref_geno_ref_alt_alleles = load_in_ref_alt_allele_arr(args.bgen_file.split('.bgen')[0] + '.pvar')
ref_geno_sdevs = load_in_alt_allele_sdevs(ref_genotype_obj, np.arange(len(ref_G_obj_rsids)), ref_geno_ref_alt_alleles)

t = open(args.output_file,'w')
t.write('variant_id\tchrom_num\tposition\ta1\ta2\tgenotype_dosage_standard_deviation\n')

for snp_iter in range(len(ref_G_obj_rsids)):

	t.write(ref_G_obj_rsids[snp_iter] + '\t' + str(args.chrom) + '\t' + str(ref_G_obj_positions[snp_iter]) + '\t')
	t.write(ref_geno_ref_alt_alleles[snp_iter][1] + '\t' + ref_geno_ref_alt_alleles[snp_iter][0] + '\t' + str(ref_geno_sdevs[snp_iter]) + '\n')

t.close()


print(args.output_file)




t.close()





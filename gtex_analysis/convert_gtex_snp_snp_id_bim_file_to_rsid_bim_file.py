import numpy as np
import os
import sys
import pdb



######################
# Command line args
######################
gtex_input_bim = sys.argv[1]
gwas_bim = sys.argv[2]
gtex_output_file = sys.argv[3]


gwas_snp_id_dicti = {}
f = open(gwas_bim)
repeat_snps = {}
for line in f:
	line = line.rstrip()
	data = line.split('\t')

	snpid1 = 'chr' + data[0] + '_' + data[3] + '_' + data[4] + '_' + data[5] + '_b38'
	snpid2 = 'chr' + data[0] + '_' + data[3] + '_' + data[5] + '_' + data[4] + '_b38'

	# Mismatch
	if snpid1 in gwas_snp_id_dicti or snpid2 in gwas_snp_id_dicti:
		repeat_snps[snpid1] = 1
		repeat_snps[snpid2] = 1
	gwas_snp_id_dicti[snpid1] = data[1]
	gwas_snp_id_dicti[snpid2] = data[1]
f.close()




f = open(gtex_input_bim)
t = open(gtex_output_file,'w')

for line in f:
	line = line.rstrip()
	data = line.split('\t')


	if data[1] not in gwas_snp_id_dicti:
		print('assumption eroror')
		pdb.set_trace()

	if data[1] in repeat_snps:
		print('assumption eroorro')
		pdb.set_trace()

	t.write(data[0] + '\t' + gwas_snp_id_dicti[data[1]] + '\t' + data[2] + '\t'+ data[3] + '\t' + data[4] + '\t' + data[5] + '\n')

f.close()
t.close()
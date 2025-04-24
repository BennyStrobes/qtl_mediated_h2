import numpy as np
import os
import sys
import pdb







######################
# Command line args
######################
gtex_bim = sys.argv[1]
gwas_bim = sys.argv[2]
output_file = sys.argv[3]


gwas_snp_id_dicti = {}
repeat_snps = {}
f = open(gwas_bim)
for line in f:
	line = line.rstrip()
	data = line.split('\t')

	snpid1 = 'chr' + data[0] + '_' + data[3] + '_' + data[4] + '_' + data[5] + '_b38'
	snpid2 = 'chr' + data[0] + '_' + data[3] + '_' + data[5] + '_' + data[4] + '_b38'

	if snpid1 in gwas_snp_id_dicti or snpid2 in gwas_snp_id_dicti:
		repeat_snps[snpid1] = 1
		repeat_snps[snpid2] = 1
		print('repeat')
	gwas_snp_id_dicti[snpid1] = 1
	gwas_snp_id_dicti[snpid2] = 1
f.close()


f = open(gtex_bim)
t = open(output_file,'w')

for line in f:
	line = line.rstrip()
	data = line.split('\t')


	if data[1] in gwas_snp_id_dicti and data[1] not in repeat_snps:
		t.write(data[1] + '\n')


f.close()
t.close()
import numpy as np
import os
import sys
import pdb








genotype_dir = sys.argv[1]

geno_file = genotype_dir + 'gwas_genotype_1.npy'
geno = np.load(geno_file)

n_chunks = 100
splits = np.array_split(np.arange(geno.shape[1]),n_chunks)

for chunk_iter in range(n_chunks):
	print(chunk_iter)
	splitter = splits[chunk_iter]

	new_geno_file = genotype_dir + 'gwas_genotype_1_' + str(chunk_iter) + '.npy'
	np.save(new_geno_file, geno[:,splitter])

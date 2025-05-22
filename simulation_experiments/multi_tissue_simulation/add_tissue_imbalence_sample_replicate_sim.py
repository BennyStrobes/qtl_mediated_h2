import sys





#####################
# Command line args
#####################
simulation_number = sys.argv[1]
chrom_string = sys.argv[2]
simulated_gene_expression_dir = sys.argv[3]
simulated_learned_gene_models_dir = sys.argv[4]
simulation_name_string = sys.argv[5]
processed_genotype_data_dir = sys.argv[6]
n_tissues = int(sys.argv[7])




global_output_file = simulated_learned_gene_models_dir + simulation_name_string + '_' + '100-1000' + '_replicate_eqtl_sumstats_xt_summary.txt'
t = open(global_output_file,'w')
t.write('tissue_name\tN_rep1\tN_rep2\teqtl_sumstat_rep1_file\teqtl_sumstat_rep1_file\n')
for tissue_iter in range(n_tissues):
	if tissue_iter == 0:
		eqtl_sample_size = 100
	else:
		eqtl_sample_size = 1000
	tissue_sumstat_output_file1 = simulated_learned_gene_models_dir + simulation_name_string + '_tissue' + str(tissue_iter) + '_' + str(eqtl_sample_size) + '_replicate1_eqtl_sumstats.txt'
	tissue_sumstat_output_file2 = simulated_learned_gene_models_dir + simulation_name_string + '_tissue' + str(tissue_iter) + '_' + str(eqtl_sample_size) + '_replicate2_eqtl_sumstats.txt'
	
	t.write('tissue' + str(tissue_iter) + '\t' + str(int(eqtl_sample_size/2)) + '\t' + str(int(eqtl_sample_size/2)) + '\t' + tissue_sumstat_output_file1 + '\t' + tissue_sumstat_output_file2 + '\n')
t.close()
print(global_output_file)
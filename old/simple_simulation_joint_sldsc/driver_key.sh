


output_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simple_simulation_joint_sldsc/"



gwas_ss="10000"
n_snps="800"
nm_h2=".2"
n_genes="60"
med_h2=".3"
snps_per_gene="5"
ge_h2=".2"
eqtl_ss="100"






sim_iter="1"
sh run_simple_simulation.sh $gwas_ss $n_snps $nm_h2 $n_genes $med_h2 $snps_per_gene $ge_h2 $eqtl_ss $sim_iter $output_dir
############################
# Input data
############################
# Real genotype data directories
simulation_genotype_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/genotype_processing/hm3_gwas_sample_size_100000/"

# Gencode hg19 gene annotation file
gencode_gene_annotation_file="/n/groups/price/ben/gene_annotation_files/gencode.v19.annotation.gtf.gz"

############################
# Output directories
############################
output_root_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/conditional_z_score_simulation_pca/"









simulation_number="1"

if false; then
sbatch run_simple_one_block_simulation.sh $simulation_genotype_dir $output_root_dir $simulation_number
fi










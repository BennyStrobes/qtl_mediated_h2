


####################
# Input data
####################

raw_ld_file="/n/scratch/users/b/bes710/qtl_mediated_h2/new_single_gene_simulation/input_data/ukb_b38_0.1_chr2.R_snp.98378738_101205867_ld.txt"

####################
# Output data
####################
output_root="/n/scratch/users/b/bes710/qtl_mediated_h2/new_single_gene_simulation/"

mediated_h2_results=${output_root}"mediated_h2_results/"

viz_dir=${output_root}"visualize_mediated_h2_results/"


##########################################################################
# Series of scripts devoted to understanding effects of correlated eqtls
##########################################################################

# Infinitesimal 
if false; then
python3 correlated_eqtl_sim_single_gene_infinitesimal_eqtl_shared_causal_snps.py
fi


if false; then
python3 correlated_eqtl_sim_single_gene_sparse_eqtl_shared_causal_snps.py
fi


if false; then
python3 correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps.py
fi

if false; then
python3 correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps_w_noise.py
fi


if false; then
python3 correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps_w_noise_windowed.py
fi


if false; then
python3 correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps_no_nm_var_in_gene.py
fi



if false ;then
python3 correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps_extra_nm_var_in_middle_of_gene.py
fi


if false; then
n_causal_eqtl_snps="2"
sbatch correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps_extra_nm_var_in_gene.sh $raw_ld_file $mediated_h2_results $n_causal_eqtl_snps
n_causal_eqtl_snps="5"
sbatch correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps_extra_nm_var_in_gene.sh $raw_ld_file $mediated_h2_results $n_causal_eqtl_snps
n_causal_eqtl_snps="25"
sbatch correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps_extra_nm_var_in_gene.sh $raw_ld_file $mediated_h2_results $n_causal_eqtl_snps


n_causal_eqtl_snps="2"
sbatch correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps_extra_nm_var_in_gene2.sh $raw_ld_file $mediated_h2_results $n_causal_eqtl_snps
n_causal_eqtl_snps="5"
sbatch correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps_extra_nm_var_in_gene2.sh $raw_ld_file $mediated_h2_results $n_causal_eqtl_snps
n_causal_eqtl_snps="25"
sbatch correlated_eqtl_sim_single_gene_sparse_eqtl_unshared_causal_snps_extra_nm_var_in_gene2.sh $raw_ld_file $mediated_h2_results $n_causal_eqtl_snps
fi


if false; then
module load R/3.5.1
fi
Rscript visualize_simulation_results.R $mediated_h2_results $viz_dir











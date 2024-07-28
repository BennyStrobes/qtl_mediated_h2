#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-40:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=12GB                         # Memory total in MiB (for all cores)






global_simulation_number="$1"
no_ld_simulation_res_dir="$2"

if false; then
module load gcc/6.2.0
module load python/3.6.0
source /n/groups/price/ben/environments/tensor_flow_cpu/bin/activate
fi

if false; then
module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate
fi


########################
# Simulation with no LD, no nm variants, and 3 causal tissue
# To estimate causal gene-trait effect sizes
########################
# Simulation parameters
n_sims="30"
gwas_ss="10000"
n_snps="7500"
eqtl_ss_2="500"
eqtl_ss_3="500"
med_h2_1="0.2"
med_h2_2="0.0"
med_h2_3="0.0"
ge_h2_1=".1"
ge_h2_2=".1"
ge_h2_3=".1"

eqtl_sample_size_arr=( "300" "1000") # 60 h
if false; then
for eqtl_ss_1 in "${eqtl_sample_size_arr[@]}"
do

	#####################
	# 1 causal tissue
	#####################
	date
	output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_no_ld_no_nm_variants_1_caus_tiss_eqtlss_"${eqtl_ss_1}
	python3 run_no_ld_simulation_with_no_ld_no_nm_var_and_1_causal_tissues_alt.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $eqtl_ss_2 $eqtl_ss_3 $med_h2_1 $med_h2_2 $med_h2_3 $ge_h2_1 $ge_h2_2 $ge_h2_3 $output_root $global_simulation_number
	date

	#####################
	# 3 causal un-correlated tissues
	#####################
	date
	sim_corr="0.0"
	output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_no_ld_no_nm_variants_3_caus_tiss_eqtlss_"${eqtl_ss_1}"_sim_corr_"${sim_corr}
	python3 run_no_ld_simulation_with_no_ld_no_nm_var_and_3_causal_tissues_alt.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $eqtl_ss_2 $eqtl_ss_3 $med_h2_1 $med_h2_2 $med_h2_3 $ge_h2_1 $ge_h2_2 $ge_h2_3 $output_root $global_simulation_number $sim_corr
	date

	#####################
	# 3 causal correlated tissues
	#####################
	date
	sim_corr="0.7"
	output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_no_ld_no_nm_variants_3_caus_tiss_eqtlss_"${eqtl_ss_1}"_sim_corr_"${sim_corr}
	python3 run_no_ld_simulation_with_no_ld_no_nm_var_and_3_causal_tissues_alt.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $eqtl_ss_2 $eqtl_ss_3 $med_h2_1 $med_h2_2 $med_h2_3 $ge_h2_1 $ge_h2_2 $ge_h2_3 $output_root $global_simulation_number $sim_corr
	date
done
fi

#####################
# 1 causal tissue
#####################
if false; then
eqtl_ss_1="300"
date
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_no_ld_no_nm_variants_1_caus_tiss_eqtlss_"${eqtl_ss_1}
python3 run_no_ld_simulation_with_no_ld_no_nm_var_and_1_causal_tissues_alt.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $eqtl_ss_2 $eqtl_ss_3 $med_h2_1 $med_h2_2 $med_h2_3 $ge_h2_1 $ge_h2_2 $ge_h2_3 $output_root $global_simulation_number
date
fi





########################
# Simulation with no LD, no nm variants, and 1 causal tissue
# To estimate squared causal gene trait effects
########################
# Simulation parameters
n_sims="30"
gwas_ss="10000"
n_snps="40000"
n_genes="100"
snps_per_gene="200"
med_h2_1="0.1"
nm_h2="0.3"
eqtl_ss_1="300"


#####################
# 1 causal tissue
#####################
if false; then
date
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_no_nm_variants_1_caus_tiss_eqtlss_"${eqtl_ss_1}
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_no_nm_var_and_1_causal_tissues.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $ge_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2
date
fi





#####################
# 1 causal tissue
#####################
date
eqtl_ss_1="300"
eqtl_architecture="polygenic"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_joint_reml"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_joint_reml.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture
date







if false; then
n_sims="1"
eqtl_ss_1="300"
eqtl_architecture="sparse"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_posterior_eqtl"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_posterior_eqtl.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture

eqtl_ss_1="300"
eqtl_architecture="polygenic"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_posterior_eqtl"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_posterior_eqtl.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture
fi



# Simulation parameters
n_sims="30"
gwas_ss="10000"
n_snps="50000"
n_genes="100"
snps_per_gene="200"
med_h2_1="0.1"
nm_h2="0.3"
eqtl_ss_1="300"

if false; then
module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/numpyro/bin/activate

n_sims="1"
eqtl_ss_1="300"
eqtl_architecture="sparse"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_posterior_eqtl_v2"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_posterior_eqtl_v2.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture

eqtl_ss_1="300"
eqtl_architecture="polygenic"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_posterior_eqtl_v2"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_posterior_eqtl_v2.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture
fi


if false; then
eqtl_ss_1="300"
eqtl_architecture="sparse"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_posterior_eqtl"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_posterior_eqtl_follow_up.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture
fi





if false; then
eqtl_ss_1="100"
eqtl_architecture="sparse"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_posterior_eqtl"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_posterior_eqtl.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture
eqtl_architecture="polygenic"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_posterior_eqtl"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_posterior_eqtl.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture

eqtl_ss_1="300"
eqtl_architecture="sparse"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_posterior_eqtl"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_posterior_eqtl.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture
eqtl_architecture="polygenic"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_posterior_eqtl"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_posterior_eqtl.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture

eqtl_ss_1="1000"
eqtl_architecture="sparse"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_posterior_eqtl"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_posterior_eqtl.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture
eqtl_architecture="polygenic"
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_eqtl_architecture_"$eqtl_architecture"_posterior_eqtl"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_and_1_causal_tissues_posterior_eqtl.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene $nm_h2 $eqtl_architecture
fi




#####################
# HE Regression stuff
#####################
if false; then
date
output_root=${no_ld_simulation_res_dir}"sim_"$global_simulation_number"_squared_gt_effects_no_ld_no_nm_variants_1_caus_tiss_eqtlss_"${eqtl_ss_1}"_he_regression"
python3 run_no_ld_simulation_for_squared_gene_trait_effects_with_no_ld_no_nm_var_and_1_causal_tissues_he_regression.py $global_simulation_number $n_sims $gwas_ss $n_snps $eqtl_ss_1 $med_h2_1 $ge_h2_1 $output_root $global_simulation_number $n_genes $snps_per_gene
date
fi















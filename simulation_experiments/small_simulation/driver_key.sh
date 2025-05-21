############################
# Input data
############################
# Real genotype data directories
simulation_genotype_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/genotype_processing/gwas_sample_size_100000/"




############################
# Output directories
############################
output_root_dir="/n/scratch/users/b/bes710/qtl_mediated_h2/simulation_experiments/small_simulation/"

ld_dir=${output_root_dir}"processed_LD/"

output_dir=${output_root_dir}"simulation_results/"


if false; then
module load gcc/9.2.0
module load python/3.9.14
module load cuda/12.1
source /n/groups/price/ben/environments/tf_new/bin/activate
fi

if false; then
python3 process_ld_data.py $simulation_genotype_dir $ld_dir
fi


# simulation parameters
gwas_ss="100000"
eqtl_ss="50"
nm_h2=".25"
med_h2=".1"
eqtl_h2=".05"
sim_sparsity="polygenic"

if false; then

eqtl_ss="50"
eqtl_h2=".05"
sh run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity

eqtl_ss="100"
eqtl_h2=".05"
sbatch run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity

eqtl_ss="300"
eqtl_h2=".05"
sbatch run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity

eqtl_ss="500"
eqtl_h2=".05"
sbatch run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity



eqtl_ss="50"
eqtl_h2=".025"
sbatch run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity

eqtl_ss="100"
eqtl_h2=".025"
sbatch run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity

eqtl_ss="300"
eqtl_h2=".025"
sbatch run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity

eqtl_ss="500"
eqtl_h2=".025"
sbatch run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity



eqtl_ss="50"
eqtl_h2=".1"
sbatch run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity

eqtl_ss="100"
eqtl_h2=".1"
sbatch run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity

eqtl_ss="300"
eqtl_h2=".1"
sbatch run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity

eqtl_ss="500"
eqtl_h2=".1"
sbatch run_simulation.sh $ld_dir $output_dir $gwas_ss $eqtl_ss $nm_h2 $med_h2 $eqtl_h2 $sim_sparsity
fi

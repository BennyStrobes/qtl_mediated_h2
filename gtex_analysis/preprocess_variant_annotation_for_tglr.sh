#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)

source ~/.bash_profile




ldsc_code_dir="$1"
hapmap3_rsid_file="$2"
baselineLD_anno_dir="$3"
ref_1kg_genotype_dir="$4"
tglr_variant_annotations_dir="$5"






num_chrom="22"


for chrom_num in $(seq 1 $(($num_chrom))); do 

	echo $chrom_num
	##################################
	# Create variant annotation files
	###################################
	# input file (baselineld annotation file)
	baseline_ld_annotation_stem=${baselineLD_anno_dir}baselineLD.${chrom_num}
	# output file (baselineld annotation file without any qtl maxcpp annotations)
	baseline_ld_no_qtl_annotation_stem=${tglr_variant_annotations_dir}baselineLD_no_qtl.${chrom_num}
	# Perform filtering
	source ~/.bash_profile
	python3 remove_qtl_annotations_from_baselineld_annotation_file.py $baseline_ld_annotation_stem $baseline_ld_no_qtl_annotation_stem

	# input file (baselineld annotation file)
	baseline_ld_annotation_stem=${baselineLD_anno_dir}baselineLD.${chrom_num}
	# output file (baselineld annotation file without any qtl maxcpp annotations)
	intercept_annotation_stem=${tglr_variant_annotations_dir}genotype_intercept.${chrom_num}
	# Perform filtering
	source ~/.bash_profile
	python3 get_genotype_intercept_annotation_file_from_baselineld_annotation_file.py $baseline_ld_annotation_stem $intercept_annotation_stem



	##################################
	# Create variant level LD-scores
	###################################
	# Load in LDSC module
	module load python/2.7.12
	# Run standard S-LDSC preprocessing on baselineLD annotations
	python ${ldsc_code_dir}ldsc.py\
		--l2\
		--bfile ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num}\
		--ld-wind-cm 1\
		--annot ${tglr_variant_annotations_dir}baselineLD_no_qtl.${chrom_num}.annot\
		--out ${tglr_variant_annotations_dir}baselineLD_no_qtl.${chrom_num}\
		--print-snps ${hapmap3_rsid_file}

	# Run standard S-LDSC preprocessing on genotype-intercept annotations
	python ${ldsc_code_dir}ldsc.py\
		--l2\
		--bfile ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num}\
		--ld-wind-cm 1\
		--annot ${tglr_variant_annotations_dir}genotype_intercept.${chrom_num}.annot\
		--out ${tglr_variant_annotations_dir}genotype_intercept.${chrom_num}\
		--print-snps ${hapmap3_rsid_file}




	##################################
	# Get weights for specific set of hapmap3 snps
	###################################
	source ~/.bash_profile
	# Filter genotype data to just regression snps in 1KG
	plink2 --bfile ${ref_1kg_genotype_dir}"1000G.EUR.hg38."${chrom_num} --extract ${hapmap3_rsid_file} --threads 1 --make-bed --keep-allele-order --out ${tglr_variant_annotations_dir}"_100G_regression_snps_only."${chrom_num}
	# Create regression snp weights
	source /n/groups/price/ben/environments/sldsc/bin/activate
	module load python/2.7.12
	python ${ldsc_code_dir}ldsc.py\
		--l2\
		--bfile ${tglr_variant_annotations_dir}"_100G_regression_snps_only."${chrom_num}\
		--ld-wind-cm 1\
		--out ${tglr_variant_annotations_dir}"regression_weights".${chrom_num}\
		--print-snps ${hapmap3_rsid_file}

	# Delete uncessary plink file
	rm ${tglr_variant_annotations_dir}"_100G_regression_snps_only."${chrom_num}*
done


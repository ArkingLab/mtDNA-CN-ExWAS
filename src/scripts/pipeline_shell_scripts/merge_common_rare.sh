#!/bin/sh
# Merge UK Biobank common variants (SNP array plink files)
# 	with rare variants (whole exome sequencing plink files)
# submit as:
#	qsub -cwd -t 1 -o merge-var.o -e merge-var.e ./merge_common_rare.sh

plink=/users/vpillala/bin/plink_linux_x86_64_20201019/plink
common_variants_dir="/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/hg38_lift/"
rare_variants_dir="/dcl01/arking/data/static/UKBiobank/WES/plink/"
output_dir="/dcl01/arking/data/static/UKBiobank/WES/plink/200k.merged.hg38lifted_common.rare.genotypes/"

# pare individuals in common fileset
# to only those with exome data (n=200,643)
# and merge variants such that 
 ${plink} \
 --bfile "${common_variants_dir}/ukb_chr${SGE_TASK_ID}.hg38_lift" \
 --bmerge "${rare_variants_dir}/UKBexomeOQFE_chr${SGE_TASK_ID}" \
 --make-bed \
 --out "${output_dir}/UKB_200k.all_genotypes_chr${SGE_TASK_ID}" \
 --keep /dcs01/arking/arkinglab/active/projects/mito_rare-variant/lmm/UKBexome.all_consent.IIDs

###${plink} \
#--bfile "${common_variants_dir}/ukb_chr24" \
#--bmerge "${rare_variants_dir}/UKBexomeOQFE_chrY" \
#--make-bed \
#--out "${output_dir}/UKB_200k.all_genotypes_chr24" \
#--keep UKBexome.all_consent.IIDs \
#--threads 8

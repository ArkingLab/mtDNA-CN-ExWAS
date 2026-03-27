#!/bin/sh
#$ -cwd
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=10G
plink2=/users/vpillala/bin/plink2
exome_plink_dir="/dcl01/arking/data/static/UKBiobank/WES/plink/"
exome_bgen_dir="/dcl01/arking/data/static/UKBiobank/WES/plink/convert_bgen/"
${plink2} --threads 8 --bfile "${exome_plink_dir}"/UKBexomeOQFE_chr${SGE_TASK_ID} --export bgen-1.2 bits=8 --out "${exome_bgen_dir}"/UKBexomeQQFE_chr${SGE_TASK_ID}
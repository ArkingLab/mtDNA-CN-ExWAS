#!/bin/sh

# Merge individual chromosome binary filesets -> single bfileset
plink=/users/vpillala/bin/plink_linux_x86_64_20201019/plink

${plink} --bfile UKBexomeOQFE_chr1 --merge-list merge-list.txt --make-bed --out UKBexomeQQFE_allChr --threads 8

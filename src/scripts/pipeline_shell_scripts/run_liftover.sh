#!/bin/sh
# Liftover script
# for chromosomes that DON'T have maps to alt chromosomes
# see https://github.com/sritchie73/liftOverPlink#usage
/users/vpillala/bin/liftover/liftOverPlink.py -e /users/vpillala/bin/liftover/liftOver \
-c /users/vpillala/bin/liftover/hg19ToHg38.over.chain.gz \
-p ukb_chr${SGE_TASK_ID}.ped \
-m ukb_chr${SGE_TASK_ID}.map \
-o ukb_chr${SGE_TASK_ID}.hg38_lift

/users/vpillala/bin/plink_linux_x86_64_20201019/plink \
--file ukb_chr${SGE_TASK_ID}.hg38_lift \
--make-bed \
--out /dcl01/arking/data/static/UKBiobank/GWAS_data/plink/hg38_lift/ukb_chr${SGE_TASK_ID}.hg38_lift
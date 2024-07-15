#!/bin/sh

bfile=ukb_chr${SGE_TASK_ID}

/users/vpillala/bin/plink_linux_x86_64_20201019/plink \
--bfile ${bfile} \
--recode tab \
--out ped/ukb_chr${SGE_TASK_ID}

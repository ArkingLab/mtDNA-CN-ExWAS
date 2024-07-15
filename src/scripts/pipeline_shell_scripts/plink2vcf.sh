#!/bin/sh

/users/vpillala/bin/plink_linux_x86_64_20201019/plink \
--bfile /dcl01/arking/data/static/UKBiobank/WES/plink/UKBexomeOQFE_chr${SGE_TASK_ID} \
--keep /dcs01/arking/arkinglab/active/projects/mito_rare-variant/lmm/UKBexome.all_consent.IIDs \
--recode vcf \
--out /dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/vcf/UKBexomeOQFE_chr${SGE_TASK_ID}

#!/bin/sh
module load conda_R/4.0.x
Rscript /users/vpillala/mito_rare-variant/scripts/ukb_gwas.plinkBED2gds.R ${SGE_TASK_ID} --num.cores 4

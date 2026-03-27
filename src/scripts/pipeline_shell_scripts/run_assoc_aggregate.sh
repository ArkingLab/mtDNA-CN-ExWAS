#!/bin/sh
# run as:
# qsub -cwd -t 2-21 -l mem_free=8G,h_vmem=8G ~/mito_rare-variant/scripts/run_assoc_aggregate.sh
chr=${SGE_TASK_ID}
#chr="3"

###module load conda_R/4.0.x
module load conda
source activate R_mkl

maf='0.01'
folder='cds_nonsyn_CADD18_0.01'
qvset='variants_in_cds_CADDphred_18'

time \
Rscript ~/mito_rare-variant/scripts/assoc_aggregate_chr2.R ${chr} \
--gds_dir /dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/ \
--output_dir /dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/run5/${folder}/ \
--AF.max ${maf} \
--qv ${qvset} \
--weight.beta_alpha 1 \
--weight.beta_beta 1 \
--invert FALSE \
--test Burden \
--nPCs 40 \
--min.mac 1 \
--score_type 0

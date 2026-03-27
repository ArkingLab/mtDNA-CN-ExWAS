## submit as:
## 	for i in {1..22}; do qsub -cwd -t 1-12 -F "${i}" /dcs01/arking/arkinglab/active/projects/mito_rare-variant/scripts/run_grm.chr.sh; done
#$ -l mem_free=2G,h_vmem=2G
#$ -pe local 4
bin=/users/vpillala/bin/
chr=$1
${bin}/gcta64 \
--bfile /dcl01/arking/data/static/UKBiobank/WES/plink/UKBexomeOQFE_chr${chr} \
--autosome \
--maf 0.01 \
--make-grm-part 12 ${SGE_TASK_ID} \
--out UKBexomeOQFE_chr${chr} \
--thread-num 4

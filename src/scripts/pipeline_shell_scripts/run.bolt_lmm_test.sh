#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe local 8
#$ -l mem_free=12G,h_vmem=12G
dir="/dcl01/arking/data/static/UKBiobank/WES/plink/convert_bgen/"
bolt=/users/vpillala/bin/BOLT-LMM_v2.3.4/bolt
# /dcs01/arking/arkinglab/software/src/BOLT-LMM/BOLT-LMM_v2.3.4/bolt

${bolt} \
--bfile=/dcl01/arking/data/static/UKBiobank/WES/plink/convert_bgen/ukb_chr7 \
--phenoFile=/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/bolt/mtDNA_ALLm2_pheno_BOLT.runfile.txt \
--phenoCol=mtDNA_CN \
--lmm \
--LDscoresFile=/dcs01/arking/arkinglab/software/src/BOLT-LMM/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz \
--statsFile=/dcl01/arking/data/static/UKBiobank/WES/plink/convert_bgen/stats.tab \
--verboseStats

# /users/vpillala/bin/BOLT-LMM_v2.3.4/bolt \
# --bed="${dir}/ukb_chr7.bed" \
# --bim="${dir}/ukb_chr7.bim" \
# --fam="${dir}/ukb_chr7.fam" \
# --remove=/dcl01/arking/data/static/UKBiobank/GWAS_data/pheno/ukb_nonwhites.txt \
# --exclude=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/snpQC/autosome_maf_lt_0.001.txt \
# --exclude=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/snpQC/autosome_missing_gt_0.1.txt \
# --phenoFile=/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/bolt/mtDNA_ALLm2_pheno_BOLT.runfile.txt \
# --phenoCol=mtDNA_CN \
# --covarFile=/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/bolt/mtDNA_ALLm2_covar_BOLT.runfile.txt \
# --covarCol=genotyping.array \
# --covarMaxLevels=30 \
# --LDscoresFile=/dcs01/arking/arkinglab/software/src/BOLT-LMM/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz \
# --geneticMapFile=/dcs01/arking/arkinglab/software/src/BOLT-LMM/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz \
# --bgenFile="${dir}/UKBexomeOQFE_chr7.bgen" \
# --sampleFile="${dir}/UKBexomeOQFE_chr7.sample" \
# --statsFile=ukb_exome_chr7_stats.gz \
# --statsFileBgenSnps=ukb_exome_chr7_bgen_stats.gz \
# --lmmForceNonInf \
# --numThreads=8 \
# --verboseStats

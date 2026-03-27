#!/bin/sh
#$ -cwd
#$ -pe local 8
#$ -l mem_free=15G,h_vmem=15G
bolt=/users/vpillala/bin/BOLT-LMM_v2.3.4/
dir1=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/
dir2=/dcl01/arking/data/static/UKBiobank/WES/plink/convert_bgen/


# ../bolt \
# ${bolt}/bolt \
bolt \
	--bed=${dir1}/ukb_chr{1:22}.bed \
	--bim=${dir1}/ukb_chr{1:22}.bim \
	--fam=/dcl01/arking/data/static/UKBiobank/WES/plink/convert_bgen/lmm/ukb_chr1.fam \
	--phenoFile=${dir2}/lmm/mtDNA_ALLm2_pheno_BOLT.runfile.txt \
	--phenoCol=mtDNA_CN \
	--lmm \
	--lmmForceNonInf \
	--remove=/dcl01/arking/data/static/UKBiobank/GWAS_data/pheno/ukb_nonwhites.txt \
	--remove=/dcl01/arking/data/static/UKBiobank/WES/plink/convert_bgen/lmm/bolt.in_plink_but_not_imputed.FID_IID.272472.txt \
	--exclude=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/snpQC/autosome_maf_lt_0.001.txt \
    --exclude=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/snpQC/autosome_missing_gt_0.1.txt \
    --covarFile=/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/bolt/mtDNA_ALLm2_covar_BOLT.runfile.txt \
    --covarCol=genotyping.array \
    --covarMaxLevels=30 \
	--LDscoresFile=${bolt}/tables/LDSCORE.1000G_EUR.tab.gz \
	--geneticMapFile=${bolt}/tables//genetic_map_hg38_withX.txt.gz \
	--statsFile=1.stats \
	--statsFileBgenSnps=2.stats \
	--bgenFile=${dir2}/UKBexomeQQFE_chr{1:22}.bgen \
	--sampleFile=${dir2}/UKBexomeQQFE_chr1.sample \
    --bgenMinMAF=1e-4 \
    --numThreads=24 \
	--verboseStats
	2>&1 | tee example.log # log output written to stdout and stderr


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

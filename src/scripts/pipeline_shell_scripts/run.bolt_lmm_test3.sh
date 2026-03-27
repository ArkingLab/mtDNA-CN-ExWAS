#!/bin/sh
#$ -cwd
#$ -pe local 16
#$ -l mem_free=5G,h_vmem=5G,h_fsize=10G
# run as: 
# 	qrsh -now y -pe local 16 -l h_vmem=5G,mem_free=5G,h_fsize=10G
# 	cd /dcs01/arking/arkinglab/software/src/BOLT-LMM/
# 	module load singularity/3.2
#	module load python/
# 	singularity shell -B /tmp,/dcl01,/dcs01 singularity-bolt-lmm.simg

bolt=/dcs01/arking/arkinglab/software/src/BOLT-LMM/BOLT-LMM_v2.3.4/
dir_geno=/dcl01/arking/data/static/UKBiobank/WES/plink/200k.merged.hg38lifted_common.rare.genotypes/
dir_lmm=/dcl01/arking/data/static/UKBiobank/WES/plink/200k.merged.hg38lifted_common.rare.genotypes/lmm/

# variant exclusions:
# model.snps:
# 	autosome maf <= 0.001
# 	autosome missing >= 10%
# sample exclusions:
# 	nonwhites (European-ancestry only, see LD)
# model covariates:
# 	wes batch (1=50k batch)

# ../bolt \
# ${bolt}/bolt \
#${bolt}/bolt \
${bolt}/bolt \
	--bed=${dir_geno}/UKB_200k.all_genotypes_chr{1:23}.bed \
	--bim=${dir_geno}/UKB_200k.all_genotypes_chr{1:23}.bim \
	--fam=${dir_geno}/UKB_200k.all_genotypes_chr1.fam \
	--modelSnps=${dir_lmm}/modelsnps.rsids \
	--phenoFile=${dir_lmm}/200k_wes.lmm.pheno.txt \
	--phenoCol=mtDNA_CN \
	--lmm \
	--lmmForceNonInf \
	--remove=/dcl01/arking/data/static/UKBiobank/GWAS_data/pheno/ukb_nonwhites.txt \
	--exclude=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/snpQC/autosome_maf_lt_0.001.txt \
    --exclude=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/snpQC/autosome_missing_gt_0.1.txt \
    --covarFile=${dir_lmm}/200k_wes.lmm.covar.txt \
    --covarCol=wes.batch \
	--LDscoresFile=${bolt}/tables/LDSCORE.1000G_EUR.tab.gz \
	--geneticMapFile=/users/vpillala/bin/BOLT-LMM_v2.3.4/tables/genetic_map_hg38_withX.txt.gz \
	--statsFile=${dir_lmm}/chr1-22_XY.stats.gz \
    --numThreads=16 \
	--verboseStats \
	2>&1 | tee ${dir_lmm}/bolt_lmm.run.log # log output written to stdout and stderr


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

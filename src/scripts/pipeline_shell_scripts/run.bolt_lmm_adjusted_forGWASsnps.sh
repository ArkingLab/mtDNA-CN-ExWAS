#!/bin/sh
##$ -cwd
##$ -pe local 16
##$ -l mem_free=5G,h_vmem=5G,h_fsize=10G
# run as: 
# 	qrsh -now y -pe local 16 -l h_vmem=5G,mem_free=5G,h_fsize=10G
# OPTIONAL:
# 	cd /dcs01/arking/arkinglab/software/src/BOLT-LMM/
# 	module load singularity/3.2
#	module load python/
# 	singularity shell -B /tmp,/dcl01,/dcs01 singularity-bolt-lmm.simg

bolt=/dcs01/arking/arkinglab/software/src/BOLT-LMM/BOLT-LMM_v2.3.4/
dir_geno=/dcl01/arking/data/static/UKBiobank/WES/plink/200k.merged.hg38lifted_common.rare.genotypes/
##dir_lmm=/dcl01/arking/data/static/UKBiobank/WES/plink/200k.merged.hg38lifted_common.rare.genotypes/lmm_adjusted/
dir_lmm=/dcs04/arking/data/active/projects/mito_rare-variant/lmm_adjusted/

# variant exclusions:
# model.snps:
# 	autosome maf <= 0.001
# 	autosome missing >= 10%
# sample exclusions:
# 	nonwhites (European-ancestry only, see LD)
# model covariates:
# 	wes batch (1=50k batch)
# 	133 GWAS SNPs (129 autosomal, 4 Xchr)

# ../bolt \
# ${bolt}/bolt \
#${bolt}/bolt \
${bolt}/bolt \
	--bed=${dir_geno}/UKB_200k.all_genotypes_chr{1:23}.bed \
	--bim=${dir_geno}/UKB_200k.all_genotypes_chr{1:23}.bim \
	--fam=${dir_geno}/UKB_200k.all_genotypes_chr1.fam \
	--modelSnps=${dir_lmm}/modelsnps.rsids \
	--phenoFile=${dir_lmm}/200k_wes.lmm.pheno.wPerm.txt \
	--phenoCol=mtDNA_CN.perm \
	--lmm \
	--lmmForceNonInf \
	--remove=/dcl01/arking/data/static/UKBiobank/GWAS_data/pheno/ukb_nonwhites.txt \
	--exclude=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/snpQC/autosome_maf_lt_0.001.txt \
	--exclude=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/snpQC/autosome_missing_gt_0.1.txt \
    --covarFile=${dir_lmm}/200k_wes.lmm.covar.txt \
    --covarCol=wes.batch \
    --qCovarCol=SNP_{1:133} \
	--LDscoresFile=${bolt}/tables/LDSCORE.1000G_EUR.tab.gz \
	--geneticMapFile=/users/vpillala/bin/BOLT-LMM_v2.3.4/tables/genetic_map_hg38_withX.txt.gz \
	--statsFile=${dir_lmm}/run2.alt/chr1-22_X.stats.adjustedforGWA_snps.gz \
    --numThreads=16 \
	--verboseStats \
	2>&1 | tee ${dir_lmm}/run3.perm/bolt_lmm.run.log # log output written to stdout and stderr


cp ${dir_lmm}/run3.perm/chr1-22_X.stats.adjustedforGWA_snps.gz ${dir_lmm}/run3.perm/chr1-22_X.stats.adjustedforGWA_snps.backup.gz

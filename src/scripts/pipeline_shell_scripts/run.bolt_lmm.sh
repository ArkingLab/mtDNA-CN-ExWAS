#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe local 8
#$ -l mem_free=12G,h_vmem=12G
/dcs01/arking/arkinglab/software/src/BOLT-LMM/BOLT-LMM_v2.3.4/bolt \
    --bed=/dcl01/arking/data/static/UKBiobank/WES/plink/200k.merged.common_rare_genotypes/UKB_200k.all_genotypes_chr{1:24}.bed \
    --bim=/dcl01/arking/data/static/UKBiobank/WES/plink/200k.merged.common_rare_genotypes/UKB_200k.all_genotypes_chr{1:24}.bim \
    --fam=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/ukb_chr22_fixCol6.fam \
    --remove=/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/bolt/bolt.in_plink_but_not_imputed.FID_IID.805.txt \
    --remove=/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/bolt/bolt.in_plink_but_not_imputed.FID_IID.102.txt \
    --remove=/dcl01/arking/data/static/UKBiobank/GWAS_data/pheno/ukb_nonwhites.txt \

    --exclude=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/snpQC/autosome_maf_lt_0.001.txt \
    --exclude=/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/snpQC/autosome_missing_gt_0.1.txt \
    
    --phenoFile=/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/bolt/mtDNA_ALLm2_pheno_BOLT.runfile.txt \
    --phenoCol=mtDNA_CN \
    --covarFile=/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/bolt/mtDNA_ALLm2_covar_BOLT.runfile.txt \
    --covarCol=genotyping.array \
    --covarMaxLevels=30 \
    --LDscoresFile=/dcs01/arking/arkinglab/software/src/BOLT-LMM/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=/dcs01/arking/arkinglab/software/src/BOLT-LMM/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz \
    --lmmForceNonInf \
    --numThreads=22 \
    --statsFile=/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/bolt/mtDNA_CN.ALLm2.stats.gz \


    # --bgenFile=/dcl01/arking/data/static/UKBiobank/GWAS_data/EGAD00010001474/ukb_imp_chr{1:22}_v3.bgen \
    # --bgenMinMAF=1e-3 \
    # --bgenMinINFO=0.3 \
    # --sampleFile=/dcl01/arking/data/static/UKBiobank/GWAS_data/pheno/ukb1773_imp_chr17_v2_s487406.sample \
    # --statsFileBgenSnps=/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/bolt/mtDNA_CN.ALLm2.bgen.stats.gz \

    --verboseStats


# /dcl01/arking/data/static/UKBiobank/WES/plink/UKBexomeOQFE_chr{1:22}.bed
# /dcl01/arking/data/static/UKBiobank/WES/plink/UKBexomeOQFE_chr{1:22}.bim
# /dcl01/arking/data/static/UKBiobank/WES/plink/UKBexomeOQFE_chr{1:22}.fam

# /dcl01/arking/data/static/UKBiobank/WES/plink/UKBexomeOQFE_chrX.bed
# /dcl01/arking/data/static/UKBiobank/WES/plink/UKBexomeOQFE_chrX.bim
# /dcl01/arking/data/static/UKBiobank/WES/plink/UKBexomeOQFE_chrX.fam

# /dcl01/arking/data/static/UKBiobank/WES/plink/UKBexomeOQFE_chrY.bed
# /dcl01/arking/data/static/UKBiobank/WES/plink/UKBexomeOQFE_chrY.bim
# /dcl01/arking/data/static/UKBiobank/WES/plink/UKBexomeOQFE_chrY.fam

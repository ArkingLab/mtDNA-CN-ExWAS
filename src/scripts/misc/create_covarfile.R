# load
library(dplyr)
library(readr)

# Read in data ----
# GWAS autosomal SNPs
ukbPheno <- # takes 2-3 mins
  as_tibble(readRDS("/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/ukbPheno_12022020_genotypes.rds"))
# GWAS XChr SNPs
xchr <-
  read_tsv(
    '/dcs04/arking/data/active/projects/mito_rare-variant/lmm_adjusted/xchr_only_variants.raw',
    col_names = T
  )
# UKB exome seq batch (50k first-pass = 1, rest 2)
wes.batch <- 
  read_tsv(
    '/dcs04/arking/data/active/projects/mito_rare-variant/lmm_adjusted/200k_wes.lmm.covar.txt',
    col_names = T
  )
wes.batch <- wes.batch %>% select(IID, wes.batch)

# Combine data to final COVARIATE file for exomes data ----
gwas <- ukbPheno[,c(1,193:321)] %>% inner_join(xchr, by=c('IID')) %>%
  mutate(FID=IID) %>% select(FID, IID, everything()) %>% select(-PAT, -MAT, -SEX, -PHENOTYPE)
snps.list <- colnames(gwas)[3:ncol(gwas)]
# write actual snp names to file (for reference)
write_lines(x=snps.list, 
  file='/dcs04/arking/data/active/projects/mito_rare-variant/lmm_adjusted/n199834.covar.n133gwasSNPIDs.list.txt', 
  sep='\n')
# reformat colnames for easier input into BOLT_LMM ("SNP_#" format)
colnames(gwas) <- c('FID','IID',c(paste0('SNP_',c(1:133))))

# construct covariate table
covar <-
  gwas %>%
  inner_join(wes.batch, by = 'IID')
  # %>%
  #rename(FID = FID.x) %>%
  #select(-FID.y)

# write to file
write_tsv(x = covar,
          file = '/dcs04/arking/data/active/projects/mito_rare-variant/lmm_adjusted/run2/200k_wes.lmm.covar.txt',
          col_names = T)

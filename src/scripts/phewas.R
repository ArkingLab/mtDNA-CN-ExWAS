# PheWAS Analysis
# Author: Vamsee Pillalamarri
# For an individual, pull out the ICD codes, 
# then convert to the format PheWAS::createPhenotypes() accepts

# Pre Reqs:
# alldata<-readRDS(file="/dcl01/arking/data/static/UKBiobank/GWAS_data/pheno/ukbCombined_031521.rds")
# alldata <- alldata %>% as_tibble()
# Filter to IIDs in study (n=180256)
# k <- read.table('~/mtrv/burden/genesis/gds/n180256.keep.IIDs')
# k <- k$V2
# alldata <- alldata %>% filter(f.eid %in% k)
# alldata <- alldata %>% select(f.eid, colnames(alldata)[stringr::str_detect(colnames(alldata), "f.41270")])

source('src/phewas_functions.R')

# run as:
icd10_codes <- list()
for(i in 1:nrow(alldata)) {
  message(i)
  icd10_codes[[i]] <- get.pheno_matrix(codes = alldata[i, ])
}
icd10_codes <- as_tibble(bind_rows(icd10_codes))

# Get UKB data
ukb.age_sex <- alldata %>% select(id=f.eid,age=f.21022.0.0,sex=f.31.0.0) %>% 
  mutate(sex2=ifelse(sex=="Male",'M','F')) %>% select(id, age, sex=sex2) %>%
  as_tibble()
ukb.pcs <- alldata %>% select(id=f.eid, starts_with("PC")) %>% as_tibble()
ukb.mtDNA <- alldata %>% select(id=f.eid, contains('mtDNA'), starts_with('arrayCN')) %>% as_tibble()
ukb.misc <- alldata %>% select(id=f.eid, genotyping.array, in.white.British.ancestry.subset, used.in.pca.calculation) %>% as_tibble()

# id.sex info
#load('~/mtrv/burden/genesis/gds/n180256.keep_samples.AnnotatedDataFrame.RData')
#id.sex <- as_tibble(annot@data)
#id.sex <- id.sex %>% mutate(sex2=ifelse(sex=="Male","M","F")) %>% select(id=sample.id, sex=sex2)


# convert icd codes to phecodes
# use a custom map which removes the '.' within the ICD codes because the UKBB input
# actually has no '.' in its representation of ICD10 codes
my_icd_map = PheWAS::phecode_map_icd10 %>% mutate(code=gsub("[.]","",code))
icd10_phecodes <- createPhenotypes(id.vocab.code.index=icd10_codes, id.sex=id.sex, 
                                   full.population.ids=id.sex$id, min.code.count=1, 
                                   vocabulary.map=my_icd_map) %>% as_tibble()

# run phewas
# phecodes
i <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n180256.icd10_codes_phecodes.rds')
i <- i$icd10_phecodes
# covariates
k <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n180256.keep_samples.AnnotatedDataFrame.RData')
k <- as_tibble(k@data)
k <- k %>% mutate(sex=as.factor(sex)) %>%
  mutate(wes.batch=as.factor(wes.batch)) %>%
  select(id=sample.id, mtDNA_CN, age, sex, PC1, PC2, PC3, PC4, PC5, PC6,PC7,PC8,PC9,PC10,PC11,
         PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,PC24,PC25,
         PC26,PC27,PC28,PC29,PC30,PC31,PC32,PC33,PC34,PC35,PC36,PC37,PC38,PC39,
         PC40,wes.batch)

# SAMHD1 carriers
samhd1.carrier <-
  read.csv(
    '~/mtrv/phesant/samhd1.carriers/SAMHD1.carriers.mac_gteq_6',
    col.names = c('id', 'SAMHD1.carrier')
  ) %>% as_tibble()
data <- i %>% inner_join(k, by='id') %>% inner_join(samhd1.carrier, by='id')
results_samhd1.carrier_mac_gteq6=phewas_ext(phenotypes=colnames(i)[-1],genotypes='SAMHD1.carrier',
                                          data=data, cores=20, 
                                          covariates=colnames(k)[-c(1,2)], 
                                          additive.genotypes=F)

samhd1.carrier <-
  read.csv(
    '~/mtrv/phesant/samhd1.carriers/SAMHD1.carriers.mac_lt_6',
    col.names = c('id', 'SAMHD1.carrier')
  ) %>% as_tibble()
data <- i %>% inner_join(k, by='id') %>% inner_join(samhd1.carrier, by='id')
results_samhd1.carrier_mac_lt6=phewas_ext(phenotypes=colnames(i)[-1],genotypes='SAMHD1.carrier',
                                          data=data, cores=20, 
                                          covariates=colnames(k)[-c(1,2)], 
                                          additive.genotypes=F)

# TWNK carriers
twnk.carrier <-
  read.csv(
    '~/mtrv/phesant/twnk.carriers/TWNK.carriers.mac_gteq_6',
    col.names = c('id', 'TWNK.carrier')
  ) %>% as_tibble()
data <- i %>% inner_join(k, by='id') %>% inner_join(twnk.carrier, by='id')
results_twnk.carrier_mac_gteq6=phewas_ext(phenotypes=colnames(i)[-1],genotypes='TWNK.carrier',
                                          data=data, cores=15, 
                                          covariates=colnames(k)[-c(1,2)], 
                                          additive.genotypes=F)

twnk.carrier <-
  read.csv(
    '~/mtrv/phesant/twnk.carriers/TWNK.carriers.mac_lt_6',
    col.names = c('id', 'TWNK.carrier')
  ) %>% as_tibble()
data <- i %>% inner_join(k, by='id') %>% inner_join(twnk.carrier, by='id')
results_twnk.carrier_mac_lt6=phewas_ext(phenotypes=colnames(i)[-1],genotypes='TWNK.carrier',
                                          data=data, cores=15, 
                                          covariates=colnames(k)[-c(1,2)], 
                                          additive.genotypes=F)


# TFAM carriers
tfam.carrier <-
  read.csv(
    '~/mtrv/phesant/tfam.carriers/TFAM.carriers.mac_gteq_6',
    col.names = c('id', 'TFAM.carrier')
  ) %>% as_tibble()
data <- i %>% inner_join(k, by='id') %>% inner_join(tfam.carrier, by='id')
results_twnk.carrier_mac_lt6=phewas_ext(phenotypes=colnames(i)[-1],genotypes='TFAM.carrier',
                                        data=data, cores=15, 
                                        covariates=colnames(k)[-c(1,2)], 
                                        additive.genotypes=F)

tfam.carrier <-
  read.csv(
    '~/mtrv/phesant/tfam.carriers/TFAM.carriers.mac_lt_6',
    col.names = c('id', 'TFAM.carrier')
  ) %>% as_tibble()
data <- i %>% inner_join(k, by='id') %>% inner_join(tfam.carrier, by='id')
results_tfam.carrier_mac_lt6=phewas_ext(phenotypes=colnames(i)[-1],genotypes='TFAM.carrier',
                                        data=data, cores=15, 
                                        covariates=colnames(k)[-c(1,2)], 
                                        additive.genotypes=F)

# fdr <- qvalue(results_tfam$p)
# results_tfam$q <- fdr$qvalues
# sig_results_tfam <- results_tfam[results_tfam$bonferroni&!is.na(results_tfam$p),]


results_mtDNA_CN=phewas_ext(phenotypes=colnames(i)[-1],genotypes='mtDNA_CN',
                            data=data,cores=15,covariates=colnames(k)[-c(1,2)],
                            additive.genotypes=F)


# UKB mtDNA-CN metric PHEWAS ----
i <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n502485.icd10_phecodes.rds')
ukb.data <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n502485.ukb.data.rds')
ukb.data <- ukb.data$ukb.data
icd10_phecodes=i

data_m1 <- ukb.data %>% filter(in.white.British.ancestry.subset==1 & used.in.pca.calculation==1) %>% 
  select(id, age, sex, arrayCN_PCA_m1, paste0('PC',1:40,separate=''), genotyping.array) %>% 
  inner_join(icd10_phecodes, by='id')
data_m2 <- ukb.data %>% filter(in.white.British.ancestry.subset==1 & used.in.pca.calculation==1) %>% 
  select(id, age, sex, arrayCN_PCA_m2, paste0('PC',1:40,separate=''), genotyping.array) %>% 
  inner_join(icd10_phecodes, by='id')
data_m3 <- ukb.data %>% filter(in.white.British.ancestry.subset==1 & used.in.pca.calculation==1) %>% 
  select(id, age, sex, arrayCN_PCA_m3, paste0('PC',1:40,separate=''), genotyping.array) %>% 
  inner_join(icd10_phecodes, by='id')

results_mtDNA_CN_m1=phewas_ext(phenotypes=colnames(i)[-1],genotypes='arrayCN_PCA_m1',
                            data=data_m1,cores=10,covariates=c('age','sex',
                                                            paste0('PC',1:40,separate=''),
                                                            'genotyping.array'),
                            additive.genotypes=F)
# saveRDS(object=results_mtDNA_CN_m1, file='~/mtrv/phesant/phewas_results_mtDNA_CN_m1.rds')

results_mtDNA_CN_m2=phewas_ext(phenotypes=colnames(i)[-1],genotypes='arrayCN_PCA_m2',
                            data=data_m2,cores=15,covariates=c('age','sex',
                                                               paste0('PC',1:40,separate=''),
                                                               'genotyping.array'),
                            additive.genotypes=F)

results_mtDNA_CN_m3=phewas_ext(phenotypes=colnames(i)[-1],genotypes='arrayCN_PCA_m3',
                            data=data_m3,cores=10,covariates=c('age','sex',
                                                               paste0('PC',1:40,separate=''),
                                                               'genotyping.array'),
                            additive.genotypes=F)
results_mtDNA_CN_m3 %>% addPhecodeInfo() %>% as_tibble() %>% arrange(p) %>% 
  select(phenotype, beta, SE, OR, p, n_total, n_cases, n_controls, description, group)

# view results
phewas_m1 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_mtDNA_CN_m1.rds') %>%
  filter(!is.na(p)) %>% as_tibble() %>% addPhecodeInfo() %>% mutate(case.ctrl_ratio=(n_cases/n_controls)) %>% arrange(p)
phewas_m1 %>% 
  filter(p <= (0.05/nrow(phewas_m1))) %>%
  select(phenotype, OR, description, p, n_cases, n_controls,case.ctrl_ratio) %>% arrange(p) %>% gt() %>%
  fmt_number(columns = 'case.ctrl_ratio', n_sigfig = 3) %>% 
  fmt_scientific(columns='p', decimals = 2) %>% 
  fmt_number(columns='OR', n_sigfig = 3)

phewas_m2 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_mtDNA_CN_m2.rds') %>%
  filter(!is.na(p)) %>% as_tibble() %>% addPhecodeInfo() %>% mutate(case.ctrl_ratio=(n_cases/n_controls)) %>% arrange(p)
phewas_m2 %>% 
  filter(p <= (0.05/nrow(phewas_m2))) %>%
  select(phenotype, OR, description, p, n_cases, n_controls,case.ctrl_ratio) %>% arrange(p) %>% gt() %>%
  fmt_number(columns = 'case.ctrl_ratio', n_sigfig = 3) %>% 
  fmt_scientific(columns='p', decimals = 2) %>% 
  fmt_number(columns='OR', n_sigfig = 3)

phewas_m3 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_mtDNA_CN_m3.rds') %>%
  filter(!is.na(p)) %>% as_tibble() %>% addPhecodeInfo() %>% mutate(case.ctrl_ratio=(n_cases/n_controls)) %>% arrange(p)
phewas_m3 %>% 
  filter(p <= (0.05/nrow(phewas_m3))) %>%
  select(phenotype, OR, description, p, n_cases, n_controls,case.ctrl_ratio) %>% arrange(p) %>% gt() %>%
  fmt_number(columns = 'case.ctrl_ratio', n_sigfig = 3) %>% 
  fmt_scientific(columns='p', decimals = 2) %>% 
  fmt_number(columns='OR', n_sigfig = 3)

# case-control ratio (all of m1,m2,m3 are same, so choose m1)
phewas_m1 %>% ggplot(aes(x=case.ctrl_ratio)) + geom_histogram() + ggtitle('Case Control Ratio for 1536 Phecodes')


# GWAS SNPs ----
phewas_gwas.snps1 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results.gwas_snp_phewas1.unrel.rds')
phewas_gwas.snps2 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results.gwas_snp_phewas2.unrel.rds')
phewas_gwas.snps <- rbind(phewas_gwas.snps1, phewas_gwas.snps2) %>% as_tibble() # %>% addPhecodeInfo()
phewas_gwas.snps %>% 
  addPhecodeInfo() %>% 
  arrange(p) %>% 
  mutate(case.ctrl_ratio=(n_cases/n_controls)) %>%
  filter(p <= (0.05 / 238908)) %>%
  dplyr::select(phenotype, snp, beta, SE, OR, p, n_cases, n_controls, 
                case.ctrl_ratio, description, group) %>% 
  gt() %>% 
  fmt_number(columns = 'case.ctrl_ratio', n_sigfig = 3) %>% 
  fmt_scientific(columns='p', decimals = 2) %>% 
  fmt_number(columns='OR', n_sigfig = 3)

phewasManhattan(phewas_gwas.snps, OR.size=T, OR.direction=T, annotate.size=3) + 
  theme(legend.position = 'none')




# MTDNA rv carrier (any SAMHD1, TWNK, TFAM carrier (each MAC all)) ----
samhd1.carrier.mac_gteq6 <-
  read.csv(
    '~/mtrv/phesant/samhd1.carriers/SAMHD1.carriers.mac_gteq_6',
    col.names = c('id', 'SAMHD1.carrier.mac_gteq6')
  ) %>% as_tibble()
samhd1.carrier.mac_lt6 <-
  read.csv(
    '~/mtrv/phesant/samhd1.carriers/SAMHD1.carriers.mac_lt_6',
    col.names = c('id', 'SAMHD1.carrier.mac_lt6')
  ) %>% as_tibble()
samhd1.carrier <- samhd1.carrier.mac_gteq6 %>% 
  inner_join(samhd1.carrier.mac_lt6, by='id') %>% 
  mutate(SAMHD1.carrier.mac_all=ifelse(SAMHD1.carrier.mac_gteq6==1 | SAMHD1.carrier.mac_lt6==1, 1,0))

twnk.carrier.mac_gteq6 <-
  read.csv(
    '~/mtrv/phesant/twnk.carriers/TWNK.carriers.mac_gteq_6',
    col.names = c('id', 'TWNK.carrier.mac_gteq6')
  ) %>% as_tibble()
twnk.carrier.mac_lt6 <-
  read.csv(
    '~/mtrv/phesant/twnk.carriers/TWNK.carriers.mac_lt_6',
    col.names = c('id', 'TWNK.carrier.mac_lt6')
  ) %>% as_tibble()
twnk.carrier <- twnk.carrier.mac_gteq6 %>% inner_join(twnk.carrier.mac_lt6, by='id') %>% 
  mutate(TWNK.carrier.mac_all=ifelse(TWNK.carrier.mac_gteq6==1 | TWNK.carrier.mac_lt6==1, 1,0))

tfam.carrier.mac_gteq6 <-
  read.csv(
    '~/mtrv/phesant/tfam.carriers/TFAM.carriers.mac_gteq_6',
    col.names = c('id', 'TFAM.carrier.mac_gteq6')
  ) %>% as_tibble()
tfam.carrier.mac_lt6 <-
  read.csv(
    '~/mtrv/phesant/tfam.carriers/TFAM.carriers.mac_lt_6',
    col.names = c('id', 'TFAM.carrier.mac_lt6')
  ) %>% as_tibble()
tfam.carrier <- tfam.carrier.mac_gteq6 %>% inner_join(tfam.carrier.mac_lt6, by='id') %>% 
  mutate(TFAM.carrier.mac_all=ifelse(TFAM.carrier.mac_gteq6==1 | TFAM.carrier.mac_lt6==1, 1,0))

mtDNA.rv.carrier <- samhd1.carrier %>% inner_join(twnk.carrier, by='id') %>% inner_join(tfam.carrier, by='id')
# mtDNA.rv.carrier %>% count(SAMHD1.carrier.mac_gteq6,TWNK.carrier.mac_gteq6,TFAM.carrier.mac_gteq6)
mtDNA.rv.carrier <- mtDNA.rv.carrier %>%
  mutate(mtDNA.rv.carrier.mac_gteq6=ifelse(SAMHD1.carrier.mac_gteq6==1 | TWNK.carrier.mac_gteq6==1 | TFAM.carrier.mac_gteq6==1, 1, 0))
mtDNA.rv.carrier <- mtDNA.rv.carrier %>% inner_join(k, by=c('id'='sample.id'))

# 1. using all MAC variants
data <- mtDNA.rv.carrier %>%  
  select(id, mtDNA.rv.carrier, age, sex, mtDNA_CN, paste0('PC',1:40,separate=''), wes.batch) %>% 
  inner_join(icd10_phecodes, by='id')


results_mtDNA_rv.carrier=phewas_ext(phenotypes=colnames(i)[-1],genotypes='mtDNA.rv.carrier',
                               data=data,cores=15,covariates=c('age','sex',
                                                                  paste0('PC',1:40,separate=''),
                                                                  'wes.batch'))

# 1. using all MAC variants for SAMHD1,TWNK,TFAM mac all, separately
data <- mtDNA.rv.carrier %>%  
  select(id, SAMHD1.carrier.mac_all, age, sex, paste0('PC',1:40,separate=''), wes.batch) %>% 
  inner_join(icd10_phecodes, by='id')
data <- data %>% inner_join(twnk.carrier, by='id') %>% inner_join(tfam.carrier)

results_mtDNA_rv.carrier_SAMHD1=phewas_ext(phenotypes=colnames(icd10_phecodes)[-1],genotypes='SAMHD1.carrier.mac_all',
                                    data=data,cores=15,covariates=c('age','sex',
                                                                    paste0('PC',1:40,separate=''),
                                                                    'wes.batch'))
results_mtDNA_rv.carrier_TWNK=phewas_ext(phenotypes=colnames(icd10_phecodes)[-1],genotypes='TWNK.carrier.mac_all',
                                           data=data,cores=15,covariates=c('age','sex',
                                                                           paste0('PC',1:40,separate=''),
                                                                           'wes.batch'))
results_mtDNA_rv.carrier_TFAM=phewas_ext(phenotypes=colnames(icd10_phecodes)[-1],genotypes='TFAM.carrier.mac_all',
                                           data=data,cores=15,covariates=c('age','sex',
                                                                           paste0('PC',1:40,separate=''),
                                                                           'wes.batch'))


# 2. using only MAC >= 6 variants
data <- mtDNA.rv.carrier %>%  
  select(id, mtDNA.rv.carrier.mac_gteq6, age, sex, mtDNA_CN, paste0('PC',1:40,separate=''), wes.batch) %>% 
  inner_join(icd10_phecodes, by='id')


results_mtDNA_rv.carrier=phewas_ext(phenotypes=colnames(i)[-1],genotypes='mtDNA.rv.carrier.mac_gteq6',
                                    data=data,cores=15,covariates=c('age','sex',
                                                                    paste0('PC',1:40,separate=''),
                                                                    'wes.batch'))


# SAMHD1 meta analysis (mac < and >= 6)
samhd1.mac_lt_6 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_samhd1_mac_lt6.rds') %>% 
  addPhecodeInfo() %>% as_tibble()
samhd1.mac_gteq_6 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_samhd1_mac_gteq6.rds') %>%
  as_tibble()
samhd1.mac_lt_6$study <- 'mac_lt6'
samhd1.mac_gteq_6$study <- 'mac_gteq6'
samhd1.rv.phewas_res <- rbind(samhd1.mac_lt_6, samhd1.mac_gteq_6)
samhd1.rv.phewas_res$adjustment <- NA
samhd1.rv.phewas_res <- phewasMeta(samhd1.rv.phewas_res, fixed = F, keep.both = T, cores=1)
samhd1.rv.phewas_res %>% addPhecodeInfo() %>% select(p, description, group) %>% arrange(p) %>% filter(p <= (0.05 / 1436))
phewasManhattan(samhd1.rv.phewas_res, annotate.phenotype.description = T, OR.direction=T, OR.size=T) + theme(legend.position = 'bottom')

# TWNK meta analysis (mac < and >= 6)
twnk.mac_lt_6 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_twnk_mac_lt6.rds') %>% 
  addPhecodeInfo() %>% as_tibble()
twnk.mac_gteq_6 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_twnk_mac_gteq6.rds') %>%
  addPhecodeInfo() %>% as_tibble()
twnk.mac_lt_6$study <- 'mac_lt6'
twnk.mac_gteq_6$study <- 'mac_gteq6'
twnk.rv.phewas_res <- rbind(twnk.mac_lt_6, twnk.mac_gteq_6)
twnk.rv.phewas_res$adjustment <- NA
twnk.rv.phewas_res <- phewasMeta(twnk.rv.phewas_res, fixed = F, keep.both = T, cores=1)
twnk.rv.phewas_res %>% addPhecodeInfo() %>% select(p, description, group) %>% arrange(p) # %>% filter(p <= (0.05 / 1436))
phewasManhattan(twnk.rv.phewas_res, annotate.phenotype.description = T, OR.direction=T, OR.size=T) +
  theme(legend.position = 'bottom')

# TFAM meta analysis (mac < and >= 6)
tfam.mac_lt_6 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_tfam_mac_lt6.rds') %>% 
  addPhecodeInfo() %>% as_tibble()
tfam.mac_gteq_6 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_tfam_mac_gteq6.rds') %>%
  addPhecodeInfo() %>% as_tibble()
tfam.mac_lt_6$study <- 'mac_lt6'
tfam.mac_gteq_6$study <- 'mac_gteq6'
tfam.rv.phewas_res <- rbind(tfam.mac_lt_6, tfam.mac_gteq_6)
tfam.rv.phewas_res$adjustment <- NA
tfam.rv.phewas_res <- phewasMeta(tfam.rv.phewas_res, fixed = F, keep.both = T, cores=1)
tfam.rv.phewas_res %>% addPhecodeInfo() %>% select(p, description, group) %>% arrange(p) # %>% filter(p <= (0.05 / 1436))
phewasManhattan(tfam.rv.phewas_res, annotate.phenotype.description = T, OR.direction=T, OR.size=T) +
  theme(legend.position = 'bottom')





# meta analysis of MAC >= 6 variants across SAMHD1, TWNK, and TFAM phewases
samhd1.mac_gteq_6 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_samhd1_mac_gteq6.rds') %>%
  as_tibble()
twnk.mac_gteq_6 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_twnk_mac_gteq6.rds') %>%
  addPhecodeInfo() %>% as_tibble()
tfam.mac_gteq_6 <- TopmedPipeline::getobj('~/mito_rare-variant/resources/phewas_results_tfam_mac_gteq6.rds') %>%
  addPhecodeInfo() %>% as_tibble()

samhd1.mac_gteq_6$snp <- 'carrier'
twnk.mac_gteq_6$snp <- 'carrier'
tfam.mac_gteq_6$snp <- 'carrier'

samhd1.mac_gteq_6$study <- 'samhd1'
twnk.mac_gteq_6$study <- 'twnk'
tfam.mac_gteq_6$study <- 'tfam'

samhd1.mac_gteq_6$adjustment <- NA
twnk.mac_gteq_6$adjustment <- NA
tfam.mac_gteq_6$adjustment <- NA


mtDNA.rv.carrier.meta <- rbind(samhd1.mac_gteq_6,
                               twnk.mac_gteq_6,
                               tfam.mac_gteq_6)
mtDNA.rv.carrier.meta <- phewasMeta(mtDNA.rv.carrier.meta, fixed = F, 
                                    keep.both = T, cores=1)
mtDNA.rv.carrier.meta %>% addPhecodeInfo() %>% select(p, description, group) %>% 
  arrange(p) # %>% filter(p <= (0.05 / 1436))
phewasManhattan(mtDNA.rv.carrier.meta, annotate.phenotype.description = T, OR.direction=T, OR.size=T) +
  theme(legend.position = 'bottom')




# testing LONP1
gwas_rs11085147_T=phewas_ext(phenotypes=colnames(i)[-1],genotypes='rs11085147_T',
                                    data=data2,cores=15,covariates=c('age','sex',
                                                                    paste0('PC',1:40,separate=''),
                                                                    'genotyping.array'))






# testing SAMHD1,TWNK,TFAM unrelated
unrel <- read_tsv('~/mtrv/burden/genesis/gds/n176064.unrelated.plink.to_keep.IIDs', 
                  col_names=c('id'))
k <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n180256.keep_samples.AnnotatedDataFrame.RData')
k <- as_tibble(k@data)
k <- k %>% select('sample.id', 'wes.batch')
i <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n180256.icd10_codes_phecodes.rds')
i <- i$icd10_phecodes
ukb.dat <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n502485.ukb.data.rds')
samhd1.carrier <-
  read.csv(
    '~/mtrv/phesant/samhd1.carriers/SAMHD1.carriers.mac_gteq_6',
    col.names = c('id', 'SAMHD1.carrier')
  ) %>% as_tibble()
samhd1.carrier <-
  read.csv(
    '~/mtrv/phesant/samhd1.carriers/SAMHD1.carriers.mac_lt_6',
    col.names = c('id', 'SAMHD1.carrier')
  ) %>% as_tibble()

twnk.carrier <-
  read.csv(
    '~/mtrv/phesant/twnk.carriers/TWNK.carriers.mac_gteq_6',
    col.names = c('id', 'TWNK.carrier.mac_gteq6')
  ) %>% as_tibble()
twnk.carrier <-
  read.csv(
    '~/mtrv/phesant/twnk.carriers/TWNK.carriers.mac_lt_6',
    col.names = c('id', 'TWNK.carrier')
  ) %>% as_tibble()

tfam.carrier <-
  read.csv(
    '~/mtrv/phesant/tfam.carriers/TFAM.carriers.mac_gteq_6',
    col.names = c('id', 'TFAM.carrier.mac_gteq6')
  ) %>% as_tibble()
tfam.carrier <-
  read.csv(
    '~/mtrv/phesant/tfam.carriers/TFAM.carriers.mac_lt_6',
    col.names = c('id', 'TFAM.carrier')
  ) %>% as_tibble()


# samhd1
data <- ukb.dat$ukb.data %>% 
  inner_join(unrel, by='id') %>%
  inner_join(i, by='id') %>% 
  inner_join(k, by=c('id'='sample.id'), ) %>% 
  inner_join(samhd1.carrier, by='id')

samhd1.mac_gteq_6.unrel=phewas_ext(phenotypes=colnames(i)[-1],genotypes='SAMHD1.carrier',
                             data=data,cores=1,covariates=c('age','sex',
                                                              paste0('PC',1:40,separate=''),
                                                              'wes.batch'))
samhd1.mac_lt_6.unrel=phewas_ext(phenotypes=colnames(i)[-1],genotypes='SAMHD1.carrier',
                                   data=data,cores=20,covariates=c('age','sex',
                                                                  paste0('PC',1:40,separate=''),
                                                                  'wes.batch'))
# getting C.I. for specific significant phe-codes
tmp=phewas_ext(phenotypes='442.2',genotypes='SAMHD1.carrier',
           data=data,cores=1,covariates=c('age','sex',
                                          paste0('PC',1:40,separate=''),
                                          'wes.batch'),
           MASS.confint.level = 0.95, additive.genotypes = F)

# twnk
data <- ukb.dat$ukb.data %>% 
  inner_join(unrel, by='id') %>%
  inner_join(i, by='id') %>% 
  inner_join(k, by=c('id'='sample.id'), ) %>% 
  inner_join(twnk.carrier, by='id')

twnk.mac_gteq_6.unrel=phewas_ext(phenotypes=colnames(i)[-1],genotypes='TWNK.carrier.mac_gteq6',
                                   data=data,cores=1,covariates=c('age','sex',
                                                                  paste0('PC',1:40,separate=''),
                                                                  'wes.batch'))
twnk.mac_lt_6.unrel=phewas_ext(phenotypes=colnames(i)[-1],genotypes='TWNK.carrier',
                                 data=data,cores=20,covariates=c('age','sex',
                                                                paste0('PC',1:40,separate=''),
                                                                'wes.batch'))
# tfam
data <- ukb.dat$ukb.data %>% 
  inner_join(unrel, by='id') %>%
  inner_join(i, by='id') %>% 
  inner_join(k, by=c('id'='sample.id'), ) %>% 
  inner_join(tfam.carrier, by='id')

tfam.mac_gteq_6.unrel=phewas_ext(phenotypes=colnames(i)[-1],genotypes='TFAM.carrier',
                                   data=data,cores=20,covariates=c('age','sex',
                                                                  paste0('PC',1:40,separate=''),
                                                                  'wes.batch'))
tfam.mac_lt_6.unrel=phewas_ext(phenotypes=colnames(i)[-1],genotypes='TFAM.carrier',
                                 data=data,cores=20,covariates=c('age','sex',
                                                                 paste0('PC',1:40,separate=''),
                                                                 'wes.batch'))










# meta analysis of SAMHD1 variants after firth logistic regression correction
# MAC >= 6
samhd1.carrier <-
  read.csv(
    '~/mtrv/phesant/samhd1.carriers/SAMHD1.carriers.mac_gteq_6',
    col.names = c('id', 'SAMHD1.carrier')
  ) %>% as_tibble()
unrel <- read_tsv('~/mtrv/burden/genesis/gds/n176064.unrelated.plink.to_keep.IIDs', 
                  col_names=c('id'))
k <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n180256.keep_samples.AnnotatedDataFrame.RData')
k <- as_tibble(k@data)
k <- k %>% select('sample.id', 'wes.batch')
i <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n180256.icd10_codes_phecodes.rds')
i <- i$icd10_phecodes
ukb.dat <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n502485.ukb.data.rds')
data <- ukb.dat$ukb.data %>% 
  inner_join(unrel, by='id') %>%
  inner_join(i, by='id') %>% 
  inner_join(k, by=c('id'='sample.id'), ) %>% 
  inner_join(samhd1.carrier, by='id')

s1 <- TopmedPipeline::getobj('~/mtrv/phesant/samhd1.carriers/phewas_results_samhd1_mac_gteq6.unrel.rds')
s1 <- s1 %>% as_tibble() %>% filter(!is.na(p)) %>% filter(p <= 0.05) %>% addPhecodeInfo() %>% arrange(p)
x1 <- process.firth(s=s1, ncores=15)
# if errors in x1, then:
# x1.tmp <- process.firth(s=s1[which(sapply(x1, is.null, simplify='array')),], ncores=15)
# x1[which(sapply(x1, is.null, simplify='array'))] <- x1.tmp
s1$study <- 's1'
s1$beta <- sapply(x1, function(x) x$coefficients[2], simplify='array')
s1$SE <- sapply(x1, function(x) sqrt(diag(vcov(x))[2]), simplify='array')
s1$p <- sapply(x1, function(x) x$prob[2], simplify='array')
  
# MAC < 6
samhd1.carrier <-
  read.csv(
    '~/mtrv/phesant/samhd1.carriers/SAMHD1.carriers.mac_lt_6',
    col.names = c('id', 'SAMHD1.carrier')
  ) %>% as_tibble()
unrel <- read_tsv('~/mtrv/burden/genesis/gds/n176064.unrelated.plink.to_keep.IIDs', 
                  col_names=c('id'))
k <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n180256.keep_samples.AnnotatedDataFrame.RData')
k <- as_tibble(k@data)
k <- k %>% select('sample.id', 'wes.batch')
i <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n180256.icd10_codes_phecodes.rds')
i <- i$icd10_phecodes
ukb.dat <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n502485.ukb.data.rds')
data <- ukb.dat$ukb.data %>% 
  inner_join(unrel, by='id') %>%
  inner_join(i, by='id') %>% 
  inner_join(k, by=c('id'='sample.id'), ) %>% 
  inner_join(samhd1.carrier, by='id')
rm(list=c('k','i','ukb.dat','unrel'))
s2 <- TopmedPipeline::getobj('~/mtrv/phesant/samhd1.carriers/phewas_results_samhd1_mac_lt6.unrel.rds')
s2 <- s2 %>% as_tibble() %>% filter(!is.na(p)) %>% filter(p <= 0.05) %>% addPhecodeInfo() %>% arrange(p)
x2 <- process.firth(s=s2, ncores=15)
s2$study <- 's2'
s2$beta <- sapply(x2, function(x) x$coefficients[2], simplify='array')
s2$SE <- sapply(x2, function(x) sqrt(diag(vcov(x))[2]), simplify='array')
s2$p <- sapply(x2, function(x) x$prob[2], simplify='array')

s <- rbind(s1,s2)
s$adjustment <- 'NA'

s.meta <- phewasMeta(s, fixed=F, keep.both = T)





















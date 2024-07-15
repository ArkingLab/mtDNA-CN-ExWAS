## Read in genotypes from merged gds format
# Genotypes
dir="/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/"
gds <- seqOpen(paste0(dir,'UKBexomeQQFE_chr1-24.gds'))

## Read in covariates, create AnnotatedDataFrame, attach to merged gds
# ldl
library(TopmedPipeline)
library(dplyr)
library(Biobase)
library(RNOmni)
ukb <- tibble(getobj('/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/ukbPheno_12022020_genotypes.rds'))
# ldl <- ukb %>% select(IID, LDL)
crp <- ukb %>% dplyr::select(IID, crp)


# mtDNA CN phenotype data
mt.cn <-
  read_tsv(
    '/dcs01/arking/arkinglab/active/projects/mito_rare-variant/resources/200k_wes.lmm.pheno.txt',
    col_names = T
  )
# covariates (wes.batch)
covar_wes.batch <-
  read_tsv(
    '/dcs01/arking/arkinglab/active/projects/mito_rare-variant/resources/200k_wes.lmm.covar.txt',
    col_names = T
  ) %>% 
  dplyr::select(IID, wes.batch)
# other covariates (age, sex)
covar_age_sex <-
  ukb %>% dplyr::select(IID, age, sex)
    # read_tsv(
  #   '/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/covar/n180256.age_sex.covar.txt',
  #   col_names=T
  # )
# pcs 1-10
covar_pcs <-
  ukb %>% dplyr::select(IID, starts_with('PC'))
  # read_tsv(
  #   '/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/covar/ukb_sample_qc.pcs1-10.n180256.txt',
  #   col_names=T
  # )
# final covariates
covar <- covar_age_sex %>%
  inner_join(covar_pcs, by='IID') %>%
  inner_join(covar_wes.batch, by='IID') 
  #%>%
  # dplyr::select(-FID.x,-FID.y) %>%
  # dplyr::select(-starts_with('SNP_'))
# nonwhites
nonwhites <-
  read_tsv(
    '/dcs01/arking/arkinglab/active/projects/mito_rare-variant/resources/ukb_nonwhites.txt',
    col_names = T
  ) %>% dplyr::select(IID)

# Create AnnotatedDataFrame with metadata
# NOTE: `SeqVarTools::SeqVarData()` does a very strict identity check 
# using `identical()` - hence, annot$sample.id and gds$sample.id vectors
# need to be same class -- specifically, since sample.ids are represented 
# as.integer in the gds, annot$sample.id should also be as.integer
annot <- tibble(IID=seqGetData(gds,"sample.id")) %>%
  # left_join(mt.cn, by='IID') %>%
  left_join(crp, by='IID') %>% 
  left_join(covar, by='IID') %>%
  dplyr::select('sample.id'='IID','crp','age','sex',paste0('PC',1:40,sep=''),'wes.batch')
annot$sample.id <- as.integer(annot$sample.id) 

# metadata to describe columns
metadata <- data.frame(labelDescription=c("IID (as.integer)",
                                          'crp',
                                          'Age',
                                          'Sex',
                                          paste0('PC',1:40,sep=''),
                                          "WES Batch Number (1=50k, 2=Remaining Samples)"),
                       row.names=names(annot),
                       stringsAsFactors = F)
annot <- AnnotatedDataFrame(as.data.frame(annot), metadata)
save(list='annot', file='/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/n200643.all_samples.CRP.AnnotatedDataFrame.RData')

# add AnnotatedDataFrame to gds
stopifnot(identical(annot$sample.id, seqGetData(gds, "sample.id")))
seqData <- SeqVarData(gds, sampleData=annot)

# 180256
samples <- seqGetData(gds, 'sample.id')
keep.samples <-
  tibble(IID=samples) %>%
  # inner_join(mt.cn, by='IID') %>%
  # select(-FID) %>%
  anti_join(nonwhites, by='IID') %>%
  # filter(!is.na(mtDNA_CN)) %>%
  inner_join(covar, by='IID') %>%
  #select(FID, IID, mtDNA_CN, wes.batch)
  inner_join(crp, by='IID') %>% 
  # dplyr::select(IID, crp) %>%
  # inner_join(covar, by='IID') %>% # attach covariates for subset
  filter(!is.na(crp))

keep.samples$crp <- RNOmni::RankNorm(keep.samples$crp)

# write to file (for PLINK)
runonce=TRUE
if(!runonce){
  keep.samples %>% mutate(FID=IID) %>% dplyr::select(FID, IID) %>%
    write_tsv(paste0('/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/n', nrow(keep.samples), '.crp.keep.IIDs'), col_names = F)
}


# Create new AnnotatedDataFrame for filtered samples
keep.samples_annot <-
  #keep.samples %>% select('sample.id' = 'IID', 'mtDNA_CN', 'wes.batch') %>%
  keep.samples %>% dplyr::select('sample.id'='IID','crp','age','sex',paste0('PC',1:40,sep=''),'wes.batch') %>%
  mutate(Population=rep('EA', nrow(keep.samples))) %>%
  mutate(Population.Description=rep('European Ancestry', nrow(keep.samples)))
keep.samples_annot$sample.id <- as.integer(keep.samples_annot$sample.id)
metadata <- data.frame(labelDescription=c("IID (as.integer)",
                                          'CRP',
                                          'Age',
                                          'Sex',
                                          paste0('PC',1:40,sep=''),
                                          "WES Batch Number (1=50k, 2=Remaining Samples)",
                                          "Population",
                                          "Population Description"),
                       row.names=names(keep.samples_annot),
                       stringsAsFactors = F)

# Create new sample subset's AnnotatedDataFrame
annot <- keep.samples_annot <- Biobase::AnnotatedDataFrame(as.data.frame(keep.samples_annot), metadata)
dir='/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/'
save(list='annot', file=paste0(dir,'/n',nrow(keep.samples),'.keep_samples.wCRP.untransformed.AnnotatedDataFrame.RData'))
## Read in genotypes from merged gds format
# Genotypes
dir="/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/"
gds <- seqOpen(paste0(dir,'UKBexomeQQFE_chr1-24.gds'))

## Read in covariates, create AnnotatedDataFrame, attach to merged gds
# ldl
library(TopmedPipeline)
library(dplyr)
ukb <- tibble(getobj('/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/ukbPheno_12022020_genotypes.rds'))
ldl <- ukb %>% select(IID, LDL)

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
  )
# other covariates (age, sex)
covar_age_sex <-
  read_tsv(
    '/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/covar/n180256.age_sex.covar.txt',
    col_names=T
  )
# pcs 1-10
covar_pcs <-
  read_tsv(
    '/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/covar/ukb_sample_qc.pcs1-10.n180256.txt',
    col_names=T
  )
# final covariates
covar <- covar_age_sex %>%
  inner_join(covar_pcs, by='IID') %>%
  inner_join(covar_wes.batch, by='IID') %>%
  select(-FID.x,-FID.y) %>%
  select(-starts_with('SNP_'))
# nonwhites
nonwhites <-
  read_tsv(
    '/dcs01/arking/arkinglab/active/projects/mito_rare-variant/resources/ukb_nonwhites.txt',
    col_names = T
  )

# Create AnnotatedDataFrame with metadata
# NOTE: `SeqVarTools::SeqVarData()` does a very strict identity check 
# using `identical()` - hence, annot$sample.id and gds$sample.id vectors
# need to be same class -- specifically, since sample.ids are represented 
# as.integer in the gds, annot$sample.id should also be as.integer
annot <- tibble(IID=seqGetData(gds,"sample.id")) %>%
  left_join(mt.cn, by='IID') %>%
  left_join(ldl, by='IID') %>% 
  left_join(covar, by='IID') %>%
  dplyr::select('sample.id'='IID','mtDNA_CN','LDL','age','sex',paste0('PC',1:10,sep=''),'wes.batch')
annot$sample.id <- as.integer(annot$sample.id) 

# metadata to describe columns
metadata <- data.frame(labelDescription=c("IID (as.integer)",
                                          "mtDNA_CN (Units are standard deviation)",
                                          'LDL (mg/dL)',
                                          'Age',
                                          'Sex',
                                          paste0('PC',1:10,sep=''),
                                          "WES Batch Number (1=50k, 2=Remaining Samples)"),
                       row.names=names(annot),
                       stringsAsFactors = F)
annot <- AnnotatedDataFrame(as.data.frame(annot), metadata)
save(list='annot', file='/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/n200643.all_samples.wLDL.AnnotatedDataFrame.RData')

# add AnnotatedDataFrame to gds
stopifnot(identical(annot$sample.id, seqGetData(gds, "sample.id")))
seqData <- SeqVarData(gds, sampleData=annot)

# 180256
samples <- seqGetData(gds, 'sample.id')
keep.samples <-
  tibble(IID=samples) %>%
  inner_join(mt.cn, by='IID') %>%
  # select(-FID) %>%
  anti_join(nonwhites, by='IID') %>%
  filter(!is.na(mtDNA_CN)) %>%
  inner_join(covar, by='IID') %>%
  #select(FID, IID, mtDNA_CN, wes.batch)
  select(IID, mtDNA_CN) %>%
  inner_join(covar, by='IID') %>% # attach covariates for subset
  inner_join(ldl, by='IID') %>% 
  filter(!is.na(LDL))


# write to file (for PLINK)
runonce=TRUE
if(!runonce){
  keep.samples %>% select(FID, IID) %>%
    write_tsv('/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/n180256.keep.IIDs', col_names = F)
}

# Create new AnnotatedDataFrame for filtered samples
keep.samples_annot <-
  #keep.samples %>% select('sample.id' = 'IID', 'mtDNA_CN', 'wes.batch') %>%
  keep.samples %>% select('sample.id'='IID','mtDNA_CN','LDL','age','sex',paste0('PC',1:10,sep=''),'wes.batch') %>%
  mutate(Population=rep('EA', nrow(keep.samples))) %>%
  mutate(Population.Description=rep('European Ancestry', nrow(keep.samples)))
keep.samples_annot$sample.id <- as.integer(keep.samples_annot$sample.id)
metadata <- data.frame(labelDescription=c("IID (as.integer)",
                                          "mtDNA_CN (Units are standard deviation)",
                                          'LDL',
                                          'Age',
                                          'Sex',
                                          paste0('PC',1:10,sep=''),
                                          "WES Batch Number (1=50k, 2=Remaining Samples)",
                                          "Population",
                                          "Population Description"),
                       row.names=names(keep.samples_annot),
                       stringsAsFactors = F)

# Create new sample subset's AnnotatedDataFrame
annot <- keep.samples_annot <- Biobase::AnnotatedDataFrame(as.data.frame(keep.samples_annot), metadata)
dir='/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/'
save(list='annot', file=paste0(dir,'/n172158.keep_samples.wLDL.AnnotatedDataFrame.RData'))
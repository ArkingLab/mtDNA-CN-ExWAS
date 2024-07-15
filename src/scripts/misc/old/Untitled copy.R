# get covar of people with specific variant
# Get subset AnnotatedDataFrame
source('/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/analysis_pipeline/TopmedPipeline/R/utils.R')
gds.dir='/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/'
k <- getobj(paste0(gds.dir,
                   'n180256.keep_samples.AnnotatedDataFrame.RData'))

ukb.pheno <- readRDS('/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/covar/ukbPheno_12022020_genotypes.rds')

variant <- 
  as.integer(readLines('/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/single_snp/indel_20:33663664:D:6.IIDs'))

variant.pheno <- 
  tibble(k@data) %>% inner_join(tibble(variant), by=c('sample.id'='variant')) %>%
  left_join(ukb.pheno, by=c('sample.id'='IID'))
  
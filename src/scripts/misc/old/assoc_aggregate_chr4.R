# setup
suppressPackageStartupMessages(library(argparser, quietly = T))
p <- argparser::arg_parser("Aggregate Association testing using GENESIS")
p <- argparser::add_argument(p, "chr",
                             help="Chromosome, values: 1-22,23=X,24=Y")
p <- argparser::add_argument(p, "--test",
                             default='SMMAT',
                             help="Pass through argument to `test` argument of `assocTestAggregate`")
p <- argparser::add_argument(p, "--AF.max",
                             default=0.01,
                             help="Pass through argument to `AF.max` argument of `assocTestAggregate`")
p <- argparser::add_argument(p, "--qv",
                             default="nonsyn.impact_mod_high",
                             help="Variant selection for aggregate testing. 
                             For options, see folder names in 
                             ~/mtrv/resources/qual.var/
                             Must match one of these folder names.")
p <- argparser::add_argument(p, "--genes",
                             default=NA,
                             help="Specific genes to test. List as GeneA,GeneB,etc. Must be on selected chromosome.")
p <- argparser::add_argument(p, "--genes.tested",
                             default=NA,
                             help="Path to file for genes previously tested/not-tested.")
p <- argparser::add_argument(p, "--max.nvar",
                             default=1000,
                             help="Upper threshold of # variants in a gene for it to be tested. Any genes with nvar > max.nvar will not be tested. Set this to speed up runtime/memory req. Use `Inf` to set no-limit.")
p <- argparser::add_argument(p, "--ncores",
                             default=1,
                             help="Number of cores to use for parallel processing")
p <- argparser::add_argument(p, "--res.dir",
                             help="Path for outputting results.")
argv <- argparser::parse_args(p)

message('-----------------------------------------------------')
# Process input args
chr=argv$chr
AF.max=argv$AF.max
test=argv$test
qv=argv$qv
ncores=argv$ncores
genes=argv$genes
res.dir=argv$res.dir; 

  # chr='19'
  # AF.max=0.0001
  # test='SMMAT'
  # qv='nonsyn.impact_mod_high'
  # ncores=1

message('chr: ', chr)
message('AF.max: ', AF.max)
message('test: ', test)
message('qualifying variant set: ', qv)
message('ncores: ', ncores)
if(is.na(res.dir)){
  # no res.dir provided, create default:
  res.dir = paste0('~/mtrv/resources/res/', 
                   qv,'.','AF.max.',AF.max,'_chr',chr, '/')
}
message('Results directory: \n', res.dir)

message('Loading packages..')
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(CompQuadForm))
suppressPackageStartupMessages(library(GenomeInfoDb))

# Helper functions
getobj <- function (Rdata) {
  require(tools)
  if (tolower(file_ext(Rdata)) == "rds") {
    return(readRDS(Rdata))
  }
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata, 
                  "\nReturning only the first object"))
  }
  return(get(objname))
}
write_tsv_custom <- function(x, file, append = TRUE, ...) {
  # SOURCE: https://github.com/tidyverse/readr/issues/832#issuecomment-384031454
  col_names <- !append
  if (append) {
    if (!file.exists(file)) {
      append <- FALSE
      col_names <- TRUE
    }
  }
  readr::write_tsv(x = x, file = file, append = append, col_names = col_names, ...)
}
setFilter <- function(seqData, sample.selection, chr, var, AF.max, reset=T){
# sets sample / variant filters on gds
  if(reset){SeqArray::seqResetFilter(seqData)}
  message('Setting sample & variant filters ',
          '[qualifying variants, then MAF filter of ',
          '[singletons=1/(2*N) .. to .. MAF < ',AF.max*100,'%]]..')
  # sample filter
  samples <- seqGetData(seqData, 'sample.id')
  seqSetFilter(seqData,
               sample.sel = which(samples %in% sample.selection),
               verbose = T)
  
  # variant filters
  seqSetFilterChrom(seqData, include = chr, verbose = T)
  pos <- seqGetData(seqData, 'position') # set qualifying variants
  variant.id <- seqGetData(seqData, 'variant.id')
  seqSetFilter(seqData, variant.id = variant.id[pos %in% var])
  
  # Filter variants to min MAF = MAF @ singleton MAC = (1 / (2*N)) where N=415422 and max-MAF = 0.01 (1%)
  # NOTE: MAF = the prp. of minor allele in pool of all alleles at any given site,
  #             and `all alleles` defined for diploid indv. as 2 * N=#individuals
  #  thus, MAF = MAC / (2*N),
  #  and   MAC = MAF * N * 2
  #  and   N = MAC / (2*MAF)
  # NOTE: can use MAF = [1/(2*N), 0.01] or MAC=[1, 8308.44=MAC at 1% MAF]
  # NOTE: either filter also (correctly) removes MAF/MAC=0 variants, 
  #       those that are unique to non-white populations (i.e. not in 415422 set)
  seqSetFilterCond(seqData, maf = c((1/(2*length(sample.selection))), AF.max), 
                   .progress = FALSE, verbose=T, parallel=(parallel::detectCores()-2))
  # maf <- seqAlleleFreq(seqData, minor=T, .progress=T, parallel=(parallel::detectCores()-2))
  # mac <- seqAlleleCount(seqData, minor=T, .progress=T, parallel=(parallel::detectCores()-2))
  return(seqData)
}
get.gtf <- function(chr){
  require(GenomeInfoDb)
# loads gtf
  # if(chr=='2'){
  #   # load a gtf that doesn't have TTN in it - TTN is too big to do analysis on
  #   # using this approach of single-gene testing (need to implement TopmedPipeline's
  #   # approach of testing segments then combining)
  #   gtf <- getobj(paste0('~/mtrv/resources/gtf/gencode_v38.chr',chr,'_protein_coding_genes.gtf2.rds'))
  # } else if (chr=='23'){
  #   gtf <- getobj(paste0('~/mtrv/resources/gtf/gencode_v38.chr','X','_protein_coding_genes.gtf.rds'))
  # } else {
  #   gtf <- getobj(paste0('~/mtrv/resources/gtf/gencode_v38.chr',chr,'_protein_coding_genes.gtf.rds'))
  # }
  message('Loading gtf')
  
  if(chr=='23'){
    gtf <- getobj(paste0('~/mtrv/resources/gtf/gencode_v38.chr','X','_protein_coding_genes.gtf.rds'))
    GenomeInfoDb::seqlevels(gtf) <- c('X'="23") # to match style with UKB data
  } else if(chr=='24'){
    gtf <- getobj(paste0('~/mtrv/resources/gtf/gencode_v38.chr','Y','_protein_coding_genes.gtf.rds'))
    seqlevels(gtf) <- c('Y'="24") # to match style with UKB data
  } else {
    gtf <- getobj(paste0('~/mtrv/resources/gtf/gencode_v38.chr',chr,'_protein_coding_genes.gtf.rds'))
  }
  
  # chromosome specific exclusions for some genes that are too big / dense variants
  if(chr==2){
    gtf <- gtf[-which(names(gtf) == 'TTN')]
  } else if(chr==5){
    gtf <- gtf[!stringr::str_detect(names(gtf), 'PCDH')]
  }
  return(gtf)
}
get.qual_var <- function(qv.dir, qv, chr){
  # qualifying variants
  qv.dir <- '~/mtrv/resources/qual.var/'
  if(chr=='23') {
    chr_='X'
    # var <- as.integer(readLines(paste0(qv.dir,'/', qv,'/',qv,'.chr','X','.uniq_pos')))  
    if (qv == 'nonsyn.impact_mod_high') {
      var <- as.integer(readLines(paste0(qv.dir,'/', qv,'/ukb23149_c', chr_,
                                         '_b0_v1.pvar.vep.out.impact_mod_high.uniqpos')))
    } else if (qv == 'nonsyn.impact_mod_high_CADD18') {
      var <- as.integer(readLines(paste0(qv.dir,'/', qv,'/ukb23149_c', 
                                         chr_,'_b0_v1.pvar.vep.out.impact_mod_high_CADD18.uniqpos')))
    } else if (qv == 'nonsyn.impact_high_pLOF') {
      var <- as.integer(readLines(paste0(qv.dir,'/', qv,'/ukb23149_c', 
                                         chr_,'_b0_v1.pvar.vep.out.impact_high_pLOF.uniqpos')))
    }
  } else {
    if (qv == 'nonsyn.impact_mod_high') {
      var <- as.integer(readLines(paste0(qv.dir,'/', qv,'/ukb23149_c', chr,
                                         '_b0_v1.pvar.vep.out.impact_mod_high.uniqpos')))
    } else if (qv == 'nonsyn.impact_mod_high_CADD18') {
      var <- as.integer(readLines(paste0(qv.dir,'/', qv,'/ukb23149_c', 
                                         chr,'_b0_v1.pvar.vep.out.impact_mod_high_CADD18.uniqpos')))
    } else if (qv == 'nonsyn.impact_high_pLOF') {
      var <- as.integer(readLines(paste0(qv.dir,'/', qv,'/ukb23149_c', 
                                         chr,'_b0_v1.pvar.vep.out.impact_high_pLOF.uniqpos')))
    }
    # var <- as.integer(readLines(paste0(qv.dir,'/', qv,'/',qv,'.chr',chr,'.uniq_pos')))
  }
}
precompute <- function() {
  # Loads precomputed data
  # `a` = all samples AnnotatedDF
  # `sample.selection` = selected samples after filtering to test
  #  (note: GENESIS might further filter depending on results of `complete.case` covariate filter)
  # `nullmod` = Random Effecsts Fitted null model
  # `seqData` gds file of genotype data
  message('Loading precomputed data..')
  sample.selection <- readLines('~/mtrv/resources/n428k.exomes.nonwhites_nonconsent.removed.sample.IIDs.list')
  
  seqData <- seqOpen(paste0('~/mtrv/resources/gds/UKBexomeQQFE_chr',chr,'.gds'), allow.duplicate=T)
  smpls <- seqGetData(seqData, 'sample.id')
  
  a <- getobj(paste0('~/mtrv/resources/n450k.Full_unfiltered.exomes_covar_AnnotatedDataFrame.rds'))
  a_ <- a@data %>% as_tibble()
  a_ <- a_ %>% arrange(match(sample.id, smpls)) # rearrange 'sample.id' column
  a <- Biobase::AnnotatedDataFrame(data = as.data.frame(a_), varMetadata = a@varMetadata) # recreate AnnotatedDataFrame
  
  seqData <- SeqVarData(seqData, sampleData=a)
  
  # pare sample.selection to those samples with !is.na(mtDNA_CN), n=415422
  a_ <- a@data %>% as_tibble() %>% filter(sample.id %in% sample.selection) %>% filter(!is.na(arrayCN_PCAodd_m2))
  sample.selection <- a_$sample.id # N = 415,422
  nullmod <- getobj(paste0('~/mtrv/resources/n450k.exomes.w_nonwhites_removed.nullmod.re.rds'))
  return(list(seqData, nullmod, sample.selection))
}

p <- precompute()
seqData <- p[[1]]
nullmod <- p[[2]]
sample.selection <- p[[3]]

gtf <- get.gtf(chr)

# if genes is not NA, select genes listed in `genes`
if(!is.na(genes)) {
  if (genes == 'non-tested') {
    # genes tested will be those not previously tested due to mem limits, etc.
    # these are expected to be in a file, which, when read with `read_tsv()`, 
    # should produce a tibble of cols: <gene, tested>
    genes.tested <- read_tsv(argv$genes.tested, show_col_types = FALSE)
    genes.tested <- genes.tested %>% filter(tested==0)
    genes <- genes.tested$gene
    message('Genes selected will be non-tested genes for chr', chr, '.')
    message('Genes selected: ', paste0(genes, separate='', collapse=', '))
    
    gtf <- gtf[which(names(gtf) %in% genes)]
    genes.tested <-
      tibble(gene = names(gtf), tested = rep(1, length(names(gtf))))
  } else{
    # assume a list of comma-separated genes was provided
    # as a single string argument --genes <geneA,geneB,..>
    message('Genes selected: ', genes)
    gtf <- gtf[which(names(gtf) %in% strsplit(genes, split = ',')[[1]])]
    genes.tested <-
      tibble(gene = names(gtf), tested = rep(1, length(names(gtf))))
  }
}

# get qual var
var <- get.qual_var(qv.dir, qv, chr)

# set sample & variant filters
seqData <- setFilter(seqData, sample.selection=sample.selection, chr=chr, var=var, 
                     AF.max=AF.max, reset=T)

# create iterator
iterator1 <-
  SeqVarRangeIterator(seqData,
                      variantRanges = gtf,
                      verbose = F)
nvar <- sapply(iterator1@variantFilter, length, simplify='array')
message('summary of variant block sizes:')
print(summary(nvar))

# OLD NOTE: Only test genes with number of variants < 99th percentile
# because > 99th percentile nvar generally are too big and take too long to run
# additionally, memory requirements scale with nvar, s.t. [N x nvar] becomes
# inconsistent, leading to some iterations where memory requirement is low
# and some have very high requirement (i.e. those with 1000s of variants)
# such as TTN which has >20k! and on chr15, a series of overlapping genes encoding
# protocadherins (PCHD*) that lead to many thousands per gene.
# Until I come up with a better method for testing these (perhaps after initial run,
# I can test them separately using a single run), I will choose to not test these
# using a quantile(nvar, 0.99) cutoff. Only do this if genes vector == NA or if 
# any nvar > 1000 (arbritrary threshold)
if(!is.na(genes)){
  iterator <- iterator1
  gtf_ <- gtf
  # print nvar sizes for genes selected
  nvar <- sapply(iterator1@variantFilter, length, simplify='array')
  message('summary of variant block sizes:')
  print(summary(nvar))
} else {
  if((argv$max.nvar !="Inf") & any(nvar > argv$max.nvar)){
    message('\n',
            '>> IMPT NOTE: removing any genes from assoc testing with nvar > max.nvar threshold = ',
            argv$max.nvar,
            ' to speed up runtime / mem req. <<',
            '\n'
    )
    message('Variant Block sizes > max nvar: ', 
            paste0(nvar[which(nvar > argv$max.nvar)], separate=' ', collapse=' '))
    gtf_ <- gtf[which(nvar < argv$max.nvar)]
    
    # set genes not tested to 0 in `genes.tested` tbl
    diff.genes <- setdiff(names(gtf), names(gtf_))
    genes.tested <-
      tibble(gene = names(gtf), tested = rep(1, length(names(gtf))))
    genes.tested$tested[which(genes.tested$gene %in% diff.genes)] <- 0
    
    # For some reason, SeqVarTools::SeqVarRangeIterator() needs 
    # to reset seqData and re-filter for variants/samples before creating new iterator..
    # so, resetting...
    message('RESETTING sample/variant filters:')
    seqData <- setFilter(seqData, sample.selection, chr, var, AF.max, reset=T)
    
    rm(iterator1)
    iterator2 <- # recreate iterator
      SeqVarRangeIterator(seqData,
                          variantRanges = gtf_,
                          verbose = F)
    # SeqVarTools::resetIterator(iterator2)
    message('summary of new variant block sizes:')
    nvar <- sapply(iterator2@variantFilter, length, simplify='array')
    print(summary(nvar))
    iterator <- iterator2
  } else {
    iterator <- iterator1
    genes.tested <-
      tibble(gene = names(gtf), tested = rep(1, length(names(gtf))))
    gtf_ <- gtf
  }
}

# writeLines(
#   con = paste0(res.dir, '/genes.tested.txt'),
#   text = names(iterator@variantRanges)
# )

# run association..
message('---Chromosome: ', chr, ' | AF.max: ', AF.max, ' | QV: ', qv, '---')
message('Number of genes to test: ', length(iterator@variantRanges))
if(ncores==1) {
  assoc <-
    assocTestAggregate(
      iterator,
      nullmod,
      test = test,
      verbose = T,
      BPPARAM = BiocParallel::SerialParam(progressbar = TRUE,
                                          stop.on.error = FALSE)
    )
} else {
  assoc <-
    assocTestAggregate(
      iterator,
      nullmod,
      test = test,
      verbose = T,
      BPPARAM = BiocParallel::MulticoreParam(
        workers = ncores,
        progressbar = TRUE,
        tasks = length(gtf_),
        stop.on.error = FALSE
      )
    )
}

# harvest results
dir.create(file.path(res.dir), showWarnings = FALSE)
# reformat test statistics table
assoc[[1]] <- 
  tibble::tibble(assoc[[1]]) %>%
  dplyr::mutate(unit=seq(1,nrow(assoc[[1]]))) %>%
  dplyr::mutate(cfreq=sapply(assoc[[2]], 
                             FUN=function(x)sum(x$freq), 
                             simplify='array')) %>%
  dplyr::mutate(gene=names(gtf_))

# print assoc summary stats to file
assoc[[1]] %>% write_tsv_custom(file =  paste0(res.dir,'/rv_assoc.aggr.',qv,'.', 
                                               test,'.','AF.max.', AF.max,'.',
                                               'chr',chr,'.txt'), append = T)
# save full assoc results as RData
if(!is.na(genes)){
  if(genes=='non-tested'){
    save(list = c('assoc'),
         file = paste0(res.dir,'/rv_assoc.aggr.',qv,'.', test,'.','AF.max.',
                       AF.max,'.','chr',chr,'.second_pass_genes.RData'))    
  } else {
    save(list = c('assoc'),
         file = paste0(res.dir,'/rv_assoc.aggr.',qv,'.', test,'.','AF.max.',
                       AF.max,'.','chr',chr,'.individual_genes.RData'))
  }
} else {
  save(list = c('assoc'),
       file = paste0(res.dir,'/rv_assoc.aggr.',qv,'.', test,'.','AF.max.',
                     AF.max,'.','chr',chr,'.RData'))
}

# write genes.tested to file
if(!is.na(genes) & genes=='non-tested'){
  genes.tested_old <- read_tsv(argv$genes.tested)
  genes.tested_old$tested[which(genes.tested_old$gene %in% names(gtf_))] <- 1
  genes.tested <- genes.tested_old
}
genes.tested %>% 
  write_tsv_custom(file = paste0(res.dir,'/rv_assoc.aggr.',qv,'.', test,'.','AF.max.', 
                                 AF.max,'.','chr',chr,'.genes_tested.txt'), append = T)

message('-----------------------------------------------------')














# GENESIS aggregate association testing script for single chromosome
# for mitochondrial copy number
sessionInfo()

suppressPackageStartupMessages(library(argparser, quietly = T))

p <- argparser::arg_parser("Aggregate Association testing using GENESIS")
p <- argparser::add_argument(p, "chr",
                             help="Chromosome, values: 1-22,23=X,24=Y")
p <- argparser::add_argument(p, "--gds_dir",
                             default='/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/', 
                             help='Top level of gds directory.')
p <- argparser::add_argument(p, "--output_dir",
                             default="/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/",
                             help="Output directory for assoc test results")
p <- argparser::add_argument(p, "--test",
                             default='Burden',
                             help="Pass through argument to `test` argument of `assocTestAggregate`")
p <- argparser::add_argument(p, "--AF.max",
                             default=0.05,
                             help="Pass through argument to `AF.max` argument of `assocTestAggregate`")
p <- argparser::add_argument(p, "--singleton",
                             default=NULL,
                             help="minor allele count filter (1=singleton)")
p <- argparser::add_argument(p, "--qv",
                             default="variants_in_cds",
                             help="Variant selection for aggregate testing. 
                             For options, see folder names in 
                             /dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/qual.var/ 
                             Must match one of these folder names.")
p <- argparser::add_argument(p, "--invert",
                             default=FALSE,
                             help="For variant_in_cds, invert positions,
                             such that variants_not_in_cds are selected instead.
                             Only run with argument qv=variants_in_cds")
p <- argparser::add_argument(p, "--weight.beta_alpha",
                             default=1,
                             help="Pass through argument to `weight.beta` argument of `assocTestAggregate` - Alpha in dbeta(alpha,beta)")
p <- argparser::add_argument(p, "--weight.beta_beta",
                             default=25,
                             help="Pass through argument to `weight.beta` argument of `assocTestAggregate` - Beta in dbeta(alpha,beta)")
p <- argparser::add_argument(p, "--nPCs",
                             default=10,
                             help="Number of PCs to select for adjusting within the null model.")
argv <- argparser::parse_args(p)

# setup
message('Loading packages')
# library(dplyr)
# library(readr)
# library(stringr)
# library(SNPRelate)
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(TopmedPipeline))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(dplyr))
# suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

dir=argv$gds_dir
chr=argv$chr
qv = argv$qv
AF.max = argv$AF.max
test = argv$test
invert = argv$invert
output_dir = argv$output_dir
wAlpha = argv$weight.beta_alpha
wBeta = argv$weight.beta_beta
weight.beta = c(wAlpha,wBeta)
nPCs = argv$nPCs

message("Output directory: ", output_dir)
message("AF.max: ", AF.max)

# dir='/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/'
# output_dir = '/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/'
# chr=20
# # qv = 'variants_in_cds_CADDphred_90thPercentile'
# qv = 'variants_in_cds'
# AF.max = 0.001
# MAC.max = 1
# test = 'SMMAT'
# invert = FALSE
# wAlpha = 1
# wBeta = 1
# weight.beta = c(wAlpha,wBeta)

# Setup ----
# Get genotypes -> gds
message('Loading gds')
gds <- seqOpen(paste0(dir,'UKBexomeQQFE_chr1-24.gds'))

# load annotated sample & variant data (cadd scores) for each sample.id / variant.id
# Convert cadd scores to AnnotatedDataFrame and attach both to gds -> `seqData`
a <- getobj(paste0(dir, '/n200643.all_samples.AnnotatedDataFrame.RData'))
t <- TopmedPipeline::getobj(paste0(dir,'UKBexomeQQFE_chr1-24.cadd_scores_for_variant_ids.rds'))
cadd.weights <- Biobase::AnnotatedDataFrame(data.frame(variant.id=t$variant.id, weight=t$CADD_PHRED))
seqData <- SeqVarData(gds, sampleData=a, variantData=cadd.weights)

# Get subsetted sample.id AnnotatedDataFrame
k <- getobj(paste0(dir,'n180256.keep_samples.AnnotatedDataFrame.RData'))
                   # 'n176064.keep_samples.AnnotatedDataFrame.RData'))

# Fit NULL MODEL -
# using subset covar/pheno/outcome data
# NOTE: Null Model was pre-computed using GRM via kingToMatrix
# NOTE: grm was created by `kingToMatrix` as sparse matrix using 
# data of pairwise kinship coefficients from UKB `gfetch rel` relatedness dat
outcome='mtDNA_CN'
# covars=c('age','sex',paste0('PC',1:nPCs,sep=''),'wes.batch')
# covars=c('age','sex','wes.batch')
# cov.mat <- GENESIS::kingToMatrix(king=paste0(dir,'/ukb_rel_a17731_s15854.kin'), estimator='Kinship', sample.include=k@data$sample.id, verbose=T)
# nullmod.noPCs <- fitNullModel(x=k, outcome=outcome, covars=covars, cov.mat=cov.mat, verbose=T)
# saveRDS(object = nullmod.noPCs, file = paste0(dir, '/n180256.nullmod.w_grm.noPCs.rds'))
# nullmod <- TopmedPipeline::getobj(paste0(dir,'/n180256.nullmod.w_grm.rds'))
nullmod <- TopmedPipeline::getobj(paste0(dir,'/n180256.nullmod.w_grm.noPCs.rds'))
message('Null Model Covariates:\n', paste0(colnames(nullmod$model.matrix), collapse = ','))

# Set SAMPLES FILTER on `seqData` (n=180256)
# samples <- seqGetData(gds, 'sample.id')
# seqSetFilter(gds, sample.sel = which(samples %in% k@data$sample.id))
samples <- seqGetData(gds, 'sample.id')
seqSetFilter(seqData, sample.sel = which(samples %in% k@data$sample.id),verbose=T)
# seqSetFilterAnnotID(seqData, id = as.character(k@data$sample.id), verbose=T)

# Set CHROMOSOME FILTER (1-22,'X','Y') (for parallelization)
seqSetFilterChrom(seqData, include=chr, verbose=T)
message('Selecting chromosome: ', chr)

# Set QUALIFYING VARIANTS FILTER
message('Selecting qualifying variants: ', qv)
qv.dir <- '/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/qual.var/'
var <- as.integer(readLines(paste0(qv.dir,'/',qv,'/',qv,'.chr',chr,'.uniq_pos')))
# if(invert){
#   # note: only should be run using `variants_in_cds` set
#   # gets the set diff of variants -- i.e. variants NOT in coding sequence ("cds")
#   pos <- seqGetData(seqData, 'position')
#   invert.var = setdiff(pos, var)
#   message('inverted variant selection.')
#   seqSetFilterPos(seqData, chr=chr, pos = invert.var, 
#                   verbose=T, multi.pos=F, intersect=F)
#   message('Selecting variants with MAF.max ', AF.max)
#   seqSetFilterCond(seqData, maf = c(0, AF.max), .progress = TRUE, verbose=T)
# } else {
# set filter positions
pos <- seqGetData(seqData, 'position')
variant.id <- seqGetData(seqData, 'variant.id')
seqSetFilter(seqData, variant.id=variant.id[pos %in% var])
## TODO - make sure we're not double counting multi allelic variants ## 

# CADD weighting - OLD code
# variant.id=variant.id[pos %in% var]
# pos=pos[pos %in% var]
# cadd <- readr::read_tsv(paste0('~/mtrv/annotations/cadd_scored/ukb23156_c',chr,'_cadd.scores.txt'))
# cadd <- cadd %>% right_join(tibble(pos), by='pos')
# weights <- data.frame(variant.id, weight=cadd$CADD_PHRED)
# seqData@variantData <- AnnotatedDataFrame(weights)

# seqSetFilterPos(
#   seqData,
#   chr = chr,
#   pos = var,
#   verbose = T,
#   multi.pos = F,
#   intersect = F
# )

# if(!is.null(argv$singleton)){
#   # filter to singletons
#   # mac=c(1,2) sets singleton filter (all remaining variants are singletons)
#   seqSetFilterCond(seqData, mac = c(1, 2), .progress = TRUE, verbose=T)
# } else {
#   # filter on minor allele frequency
message('Selecting variants with MAF.max ', AF.max)
seqSetFilterCond(seqData, maf = c(0, AF.max), .progress = TRUE, verbose=T)
# }
# }
# seqSetFilter(seqData, variant.sel = which(pos %in% var), verbose=T)
# note: interset=F and multi.pos=F are defaults
# meaning... duplicate positions indicate multiple alleles for a variant
#   So, matches of pos in var can have duplicates.
#   We're choosing to NOT allow this i.e. forcing unique matches
##   since all duplicate matches will match against a variant with the prescribed
##   property (e.g. a coding seq variant or CADD.threshold, etc)

# Association testing ----
# get gene definitions
# NEW - using GENCODE v36 GTF
# filters:
# select `protein_coding` `genes` on chosen chromosome
gtf <- rtracklayer::import('~/mtrv/burden/genesis/resources/gencode.v38.annotation.gtf.gz')
gtf <- gtf[seqnames(gtf)==paste0('chr',chr) & 
             gtf$type=='gene' & 
             gtf$gene_type=='protein_coding']
gene_symbols <- mcols(gtf)$gene_name
mcols(gtf) <- NULL
mcols(gtf) <- DataFrame(gene_name=gene_symbols)
seqlevelsStyle(gtf) <- 'NCBI' # change 'chr#' -> #
names(gtf) <- make.names(gene_symbols, unique=TRUE)

# Create an iterator where each successive unit is a different gene
# create iterator for `assoc_aggregate`
iterator <- SeqVarRangeIterator(seqData, variantRanges=gtf, verbose=F)

# Check that variants x samples won't be greater than mem limit of R
# Set cutoff for max variants / max array size issue
# Remove any gtf where NxM > 2^31 - 1
# where N = variants in gene to test, M = samples (180256); thus N*M <= (2^31)-1
# Therefore, cutoff for M=180256 is N=floor(((2^31)-1) / 180256) = 11913
max_var=11913
nvar <- sapply(iterator@variantFilter, length, simplify='array')
# nvar <- countOverlaps(gtf, pos)
if(any(nvar > max_var)){
  # if any blocks are too big, prune them; then rebuild iterator..
  message('any(nvar > max_var)=TRUE; rebuilding iterator..')
  seqResetFilter(seqData)
  rm(iterator)
  gtf <- gtf[nvar < max_var] # remove genes w/ num variants > max_var
  seqSetFilter(seqData, sample.sel = which(samples %in% k@data$sample.id),verbose=T)
  seqSetFilterChrom(seqData, include=chr, verbose=T)
  seqSetFilter(seqData, variant.id=variant.id[pos %in% var])
  seqSetFilterCond(seqData, maf = c(0, AF.max), .progress = TRUE, verbose=T)
  iterator <- SeqVarRangeIterator(seqData, variantRanges=gtf, verbose=F)
}

# do an aggregate association test on the rare variants in each gene
message('running association test: ', test)
# message('max MAF: ', AF.max)
# message('variant weights: dBeta(MAF, alpha,beta) = dBeta(MAF,', 
#         paste(weight.beta,collapse=','),')')
message('variant weights: user-specified CADD scores')
assoc <- assocTestAggregate(iterator, nullmod,test=test,weight.user='weight',genome.build='hg38',verbose=T) #weight.beta = weight.beta,
assoc[[1]] <- 
  tibble::tibble(assoc[[1]]) %>%
  dplyr::mutate(unit=seq(1,nrow(assoc[[1]]))) %>%
  dplyr::mutate(cfreq=sapply(assoc[[2]], FUN=function(x)sum(x$freq), simplify='array')) %>%
  dplyr::mutate(gene=names(gtf)) %>%
  dplyr::arrange(pval_SMMAT)

# assoc[[1]] %>% filter(cfreq >= 0.0001) # filter to CAF of 0.0001

# save results
save(list = c('assoc'),
     file = paste0(output_dir,'/assoc_aggregate.', test,'.','AF.max.',
                   AF.max,'.',chr,'.results.RData'))

message('Done.')






# OLD CODE----
# OLD - using txdb
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# seqlevels(txdb) <- seqlevels0(txdb) # reset seqlevels
# seqlevels(txdb) <- paste0('chr', chr) # select specific chromosome
# # txdb.chr <- keepSeqlevels(txdb, paste0('chr', chr))
# genes <- sort(genes(txdb))
# genes <- renameSeqlevels(x=genes, value=as.character(chr))

# remove large genes > 1Mb
# genes <- genes[-which((end(genes) - start(genes))> 1e6)]
# iterator <- SeqVarRangeIterator(seqData, variantRanges=genes, verbose=F)

# gr <- granges(gds)
#gr <- granges(seqData)

# define aggregate units:
# find variants that overlap with each gene in hg38
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#gr <- renameSeqlevels(gr, paste0("chr", seqlevels(gr)))
# ts <- transcriptsByOverlaps(txdb, gr, columns="GENEID")

# define genes as overlapping transcripts 
# (possibly over simplistic definition of aggregate unit)
# genes <- reduce(ts)
# genes <- renameSeqlevels(genes, sub("chr", "", seqlevels(genes)))






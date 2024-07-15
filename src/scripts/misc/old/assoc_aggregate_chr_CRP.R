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
argv <- argparser::parse_args(p)

# setup
message('Loading packages')
# library(dplyr)
# library(readr)
# library(stringr)
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(TopmedPipeline))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(SeqVarTools))
# library(SNPRelate)
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

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

# dir='/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/'
# output_dir = '/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/'
# chr=21
# qv = 'variants_in_cds_CADDphred'
# AF.max = 0.01
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

# Attach AnnotatedDataFrame to gds -> `seqData`
a <- getobj(paste0(dir,
                   '/n200643.all_samples.AnnotatedDataFrame.RData'))
seqData <- SeqVarData(gds, sampleData=a)

# Get subset AnnotatedDataFrame
k <- getobj(paste0(dir,
                   'n178330.keep_samples.wCRP.untransformed.AnnotatedDataFrame.RData'))

# Fit NULL MODEL - fixed effects (note: grm not provided)
# using subset covar/pheno/outcome data
outcome='crp'
covars=c('age','sex',paste0('PC',1:40,sep=''),'wes.batch')
nullmod.fe <- 
  fitNullModel(x=k, outcome=outcome, covars=covars, verbose=T)
# message('Null Model Formula:\n', nullmod.fe$model$formula)

# Set SAMPLES FILTER on `seqData` (n=180256)
# samples <- seqGetData(gds, 'sample.id')
# seqSetFilter(gds, sample.sel = which(samples %in% k@data$sample.id))
samples <- seqGetData(gds, 'sample.id')
seqSetFilter(seqData, sample.sel = which(samples %in% k@data$sample.id), 
             verbose=T)
# seqSetFilterAnnotID(seqData, id = as.character(k@data$sample.id), verbose=T)

# Set CHROMOSOME FILTER (1-22,'X','Y') (for parallelization)
#seqSetFilterChrom(seqData, include=chr, verbose=T)
message('Selecting chromosome: ', chr)

# Set QUALIFYING VARIANTS FILTER
message('Selecting qualifying variants: ', qv)
qv.dir <- '/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/qual.var/'
var <- as.integer(readLines(paste0(qv.dir,'/',qv,'/',qv,'.chr',chr,'.uniq_pos')))
if(invert){
  # note: only should be run using `variants_in_cds` set
  # gets the set diff of variants -- i.e. variants NOT in coding sequence ("cds")
  pos <- seqGetData(seqData, 'position')
  invert.var = setdiff(pos, var)
  message('inverted variant selection.')
  seqSetFilterPos(seqData, chr=chr, pos = invert.var, 
                  verbose=T, multi.pos=F, intersect=F)
  message('Selecting variants with MAF.max ', AF.max)
  seqSetFilterCond(seqData, maf = c(0, AF.max), .progress = TRUE, verbose=T)
} else {
  # set filter positions
  seqSetFilterPos(
    seqData,
    chr = chr,
    pos = var,
    verbose = T,
    multi.pos = F,
    intersect = F
  )
  # if(!is.null(argv$singleton)){
  #   # filter to singletons
  #   # mac=c(1,2) sets singleton filter (all remaining variants are singletons)
  #   seqSetFilterCond(seqData, mac = c(1, 2), .progress = TRUE, verbose=T)
  # } else {
  #   # filter on minor allele frequency
    message('Selecting variants with MAF.max ', AF.max)
    seqSetFilterCond(seqData, maf = c(0, AF.max), .progress = TRUE, verbose=T)
  # }
}
# seqSetFilter(seqData, variant.sel = which(pos %in% var), verbose=T)
# note: interset=F and multi.pos=F are defaults
# meaning... duplicate positions indicate multiple alleles for a variant
#   So, matches of pos in var can have duplicates.
#   We're choosing to NOT allow this i.e. forcing unique matches
##   since all duplicate matches will match against a variant with the prescribed
##   property (e.g. a coding seq variant or CADD.threshold, etc)

# Association testing ----
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(txdb) <- seqlevels0(txdb) # reset seqlevels
seqlevels(txdb) <- paste0('chr', chr) # select specific chromosome
# txdb.chr <- keepSeqlevels(txdb, paste0('chr', chr))
genes <- sort(genes(txdb))
genes <- renameSeqlevels(x=genes, value=as.character(chr))

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

# create an iterator where each successive unit is a different gene
iterator <- SeqVarRangeIterator(seqData, variantRanges=genes, verbose=F)

# do a aggregate association test on the rare variants in each gene
message('running association test: ', test)
message('max MAF: ', AF.max)
message('variant weights: dBeta(MAF, alpha,beta) = dBeta(MAF,', 
        paste(weight.beta,collapse=','),')')
assoc <- assocTestAggregate(iterator, nullmod.fe,
                            test=test,
                            genome.build='hg38',
                            weight.beta = weight.beta,
                            verbose=T)

# save results
save(list = c('assoc'),
     file = paste0(output_dir,'/assoc_aggregate.', test,'.','AF.max.',
                   AF.max,'.',chr,'.results.RData'))

message('Done.')















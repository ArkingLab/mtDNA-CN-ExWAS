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
p <- argparser::add_argument(p, "--min.mac",
                             default=1,
                             help="Min minor allele count of a variant to be included in assocTestAggregate")
p <- argparser::add_argument(p, "--outcome",
                             default="mtDNA_CN",
                             help="Outcome for association testing")
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
p <- argparser::add_argument(p, "--genes",
                             default=NA,
                             help="Specific genes to test. List as GeneA,GeneB,etc. Must be on selected chromosome.")
p <- argparser::add_argument(p, "--score_type",
                             default='c',
                             help="Type of score to weight variants by.")
argv <- argparser::parse_args(p)

# setup
message('Loading packages')
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(CompQuadForm))

# read in arguments
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
min.mac=argv$min.mac
nPCs = argv$nPCs
outcome = argv$outcome
genes=argv$genes
score_type=argv$score_type

# input args (hard-coded)
chr='11'
AF.max=0.01
test='SMMAT'
qv='variants_in_cds'

# Helper function
getobj <- function (Rdata)
{
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

message("Output directory: ", output_dir)
message('Variant weights: user-specified CADD scores')
message('Selecting chromosome: ', chr)
message('Selecting qualifying variants: ', qv)
message('Selecting variants with MAF.max ', AF.max)
message('Selecting variants with minimum MAC ', min.mac)
message('Running association test: ', test)
if(!is.na(genes)){message('Specific genes selected: ', genes)}

# Load in precomputed data
# NOTE: gtf splits were precomuted using `~/mito_rare-variant/scripts/split_gencode_v38_gtf.R`
# gtf = gene info on this chromosome
#`a` = all samples covar AnnotatedDF
# `sample.selection` = sample subset to run analysis on
# `nullmod` = # Null model, wGRM + PCs
# `t` optional CADD scores for variants
message('Loading precomputed data')
gtf <- getobj(paste0('~/mtrv/resources/gtf/gencode_v38.chr',chr,'_protein_coding_genes.gtf.rds'))
a <- getobj(paste0('~/mtrv/resources/n450k.Full_unfiltered.exomes_covar_AnnotatedDataFrame.rds')) # all samples AnnotatedDF
sample.selection <- readLines('~/mtrv/resources/n428k.exomes.nonwhites_nonconsent.removed.sample.IIDs.list')
nullmod <- getobj(paste0('~/mtrv/resources/n450k.exomes.w_nonwhites_removed.nullmod.re.rds'))
# t <- getobj(paste0('~/mtrv/resources/UKBexomeQQFE_chr1-24.cadd_scores_for_variant_ids.rds'))

# Compute dynamic data:
message('Computing dynamic data (chromosome specific)')
if(!is.na(genes)){
  # if genes is not NA, select genes listed in `genes`
  gtf <- gtf[strsplit(genes, split = ',')[[1]]]
}
qv.dir <- '~/mtrv/resources/qual.var/' # qualifying variants
var <- as.integer(readLines(paste0(qv.dir,'/',qv,'/',qv,'.chr',chr,'.uniq_pos')))

# Pick variant weight score type
weights = switch(
  score_type,
  "0"=NULL,
  "a"= Biobase::AnnotatedDataFrame(data.frame(variant.id=t$variant.id, weight=t$CADD_PHRED)),  
  "b"= Biobase::AnnotatedDataFrame(data.frame(variant.id=t$variant.id, weight=t$cadd_maf.combined_weight)),  
  "c"= Biobase::AnnotatedDataFrame(data.frame(variant.id=t$variant.id, weight=t$CADD_RAW)),
  "d"= Biobase::AnnotatedDataFrame(data.frame(variant.id=t$variant.id, weight=t$CADD_RAW.scaled)),
  "e"= Biobase::AnnotatedDataFrame(data.frame(variant.id=t$variant.id, weight=t$cadd.raw_maf_combined_weight)),
  "f"= Biobase::AnnotatedDataFrame(data.frame(variant.id=t$variant.id, weight=t$cadd.raw_neg_log10maf_combined_weight))
)

# load gds with annotated samples & supply weights if available
gds <- seqOpen(paste0('~/mtrv/resources/gds/UKBexomeQQFE_chr1-22.gds'))
if(!is.null(weights)){
  seqData <- SeqVarData(gds, sampleData=a, variantData=weights)
} else {
  seqData <- SeqVarData(gds, sampleData=a)
}

# set sample & variant filters
message('Setting sample & variant filters')
samples <- seqGetData(gds, 'sample.id')
seqSetFilter(seqData, sample.sel = which(samples %in% sample.selection),verbose=T) # 429002
seqSetFilterChrom(seqData, include=chr, verbose=T)
pos <- seqGetData(seqData, 'position') # set qualifying variants
variant.id <- seqGetData(seqData, 'variant.id')
seqSetFilter(seqData, variant.id=variant.id[pos %in% var])
seqSetFilterCond(seqData, maf = c(0, AF.max), .progress = TRUE, verbose=T,  parallel=30)

if(min.mac > 1){
  seqSetFilterCond(seqData, mac = min.mac, .progress = TRUE, verbose=T, parallel=30)
}

# create iterator, test if any iterations will exceed memory limit
iterator <- SeqVarRangeIterator(seqData, variantRanges=gtf, verbose=F)
max_var=11913
nvar <- sapply(iterator@variantFilter, length, simplify='array')
if(any(nvar > max_var)){
  # if any blocks are too big, prune them; then, rebuild iterator..
  message('any(nvar > max_var)=TRUE; rebuilding iterator..')
  seqResetFilter(seqData)
  rm(iterator)
  gtf <- gtf[nvar < max_var] # remove genes w/ num variants > max_var
  seqSetFilter(seqData, sample.sel = which(samples %in% k@data$sample.id),verbose=T)
  seqSetFilterChrom(seqData, include=chr, verbose=T)
  seqSetFilter(seqData, variant.id=variant.id[pos %in% var])
  seqSetFilterCond(seqData, maf = c(0, AF.max), .progress = TRUE, verbose=T)
  seqSetFilterCond(seqData, mac = min.mac, .progress = TRUE, verbose=T)
  iterator <- SeqVarRangeIterator(seqData, variantRanges=gtf, verbose=F)
}

# run association testing
# multicoreParam <- MulticoreParam(workers = multicoreWorkers(), progressbar=T)
# BPPARAM = multicoreParam
res.dir = paste0('~/mtrv/resources/res/', qv,'_','AF.max.',AF.max,'_chr',chr, '/')
dir.create(file.path(res.dir), showWarnings = TRUE)
if(!is.null(weights)){
  assoc <-
    assocTestAggregate(
      iterator,
      nullmod,
      test = test,
      weight.user = 'weight',
      genome.build = 'hg38',
      verbose = T
    ) #weight.beta = weight.beta,
} else {
  assoc <-
    assocTestAggregate(
      iterator,
      nullmod,
      test = test,
      genome.build = 'hg38',
      verbose = T,
      BPPARAM = BiocParallel::MulticoreParam(
        workers = 5,
        tasks = length(gtf),
        progressbar = T,
        stop.on.error = FALSE,
        resultdir = res.dir
      )
    )
}
assoc[[1]] <- 
  tibble::tibble(assoc[[1]]) %>%
  dplyr::mutate(unit=seq(1,nrow(assoc[[1]]))) %>%
  dplyr::mutate(cfreq=sapply(assoc[[2]], FUN=function(x)sum(x$freq), simplify='array')) %>%
  dplyr::mutate(gene=names(gtf)) # %>%
  #dplyr::arrange(pval_SMMAT)
  #dplyr::arrange(pval_burden)

# save results
save(list = c('assoc'),
     file = paste0(output_dir,'/assoc_aggregate.', test,'.','AF.max.',
                   AF.max,'.',chr,'.results.RData'))
message('Done.')

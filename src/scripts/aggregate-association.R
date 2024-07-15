#!/usr/bin/env Rscript

# Aggregate Association Testing using GENESIS
# Author: Vamsee Pillalamarri

# Load required packages
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(glue)
  library(SeqArray)
  library(SeqVarTools)
  library(GENESIS)
  library(Biobase)
  library(survey)
  library(CompQuadForm)
  library(GenomeInfoDb)
  library(future)
  library(furrr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--test"), default = "SMMAT", help = "Test type for assocTestAggregate"),
  make_option(c("--AF.max"), type = "double", default = 0.01, help = "Maximum allele frequency"),
  make_option(c("--qv"), default = "nonsyn.impact_mod_high", help = "Variant selection for aggregate testing"),
  make_option(c("--genes"), default = NA, help = "Specific genes to test (comma-separated)"),
  make_option(c("--genes.tested"), default = NA, help = "Path to file for genes previously tested/not-tested"),
  make_option(c("--max.nvar"), type = "integer", default = 1000, help = "Upper threshold of # variants in a gene"),
  make_option(c("--ncores"), type = "integer", default = 1, help = "Number of cores for parallel processing"),
  make_option(c("--res.dir"), help = "Path for outputting results")
)

opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = 1)
chr <- opt$args[1]

# Helper functions
get_gtf <- function(chr) {
  gtf_file <- glue("~/mtrv/resources/gtf/gencode_v38.chr{chr}_protein_coding_genes.gtf.rds")
  if (chr == "23") {
    gtf_file <- glue("~/mtrv/resources/gtf/gencode_v38.chrX_protein_coding_genes.gtf.rds")
    chr_name <- "X"
  } else if (chr == "24") {
    gtf_file <- glue("~/mtrv/resources/gtf/gencode_v38.chrY_protein_coding_genes.gtf.rds")
    chr_name <- "Y"
  } else {
    chr_name <- chr
  }
  
  gtf <- readRDS(gtf_file)
  seqlevels(gtf) <- setNames(chr, chr_name)
  
  # Chromosome-specific exclusions
  if (chr == "2") {
    gtf <- gtf[names(gtf) != "TTN"]
  } else if (chr == "5") {
    gtf <- gtf[!stringr::str_detect(names(gtf), "PCDH")]
  }
  
  return(gtf)
}

get_qual_var <- function(qv, chr) {
  qv_dir <- "~/mtrv/resources/qual.var/"
  chr_name <- ifelse(chr == "23", "X", chr)
  
  file_pattern <- glue("ukb23149_c{chr_name}_b0_v1.pvar.vep.out.{qv}.uniqpos")
  var_file <- list.files(path = file.path(qv_dir, qv), pattern = file_pattern, full.names = TRUE)
  
  if (length(var_file) == 0) {
    stop(glue("Qualifying variant file not found for chr{chr} and qv={qv}"))
  }
  
  return(as.integer(readLines(var_file)))
}

setFilter <- function(seqData, sample.selection, chr, var, AF.max, reset = TRUE) {
  if (reset) SeqArray::seqResetFilter(seqData)
  message("Setting sample & variant filters...")
  
  data <- seqGetData(seqData, c("sample.id", "position", "variant.id"))
  seqSetFilter(seqData,
               sample.sel = which(data$sample.id %in% sample.selection),
               variant.sel = which(data$position %in% var),
               verbose = TRUE)
  
  seqSetFilterChrom(seqData, include = chr, verbose = TRUE)
  
  seqSetFilterCond(seqData, 
                   maf = c(1 / (2 * length(sample.selection)), AF.max),
                   .progress = FALSE, 
                   verbose = TRUE, 
                   parallel = parallel::detectCores() - 2)
  
  return(seqData)
}

precompute <- function(chr) {
  message("Loading precomputed data...")
  sample.selection <- readLines("~/mtrv/resources/n428k.exomes.nonwhites_nonconsent.removed.sample.IIDs.list")
  
  seqData <- seqOpen(glue("~/mtrv/resources/gds/UKBexomeQQFE_chr{chr}.gds"), allow.duplicate = TRUE)
  smpls <- seqGetData(seqData, "sample.id")
  
  a <- readRDS("~/mtrv/resources/n450k.Full_unfiltered.exomes_covar_AnnotatedDataFrame.rds")
  a_data <- as_tibble(a@data) %>%
    arrange(match(sample.id, smpls)) %>%
    filter(sample.id %in% sample.selection, !is.na(arrayCN_PCAodd_m2))
  
  sample.selection <- a_data$sample.id
  
  a <- AnnotatedDataFrame(data = as.data.frame(a_data), varMetadata = a@varMetadata)
  seqData <- SeqVarData(seqData, sampleData = a)
  
  nullmod <- readRDS("~/mtrv/resources/n450k.exomes.w_nonwhites_removed.nullmod.re.rds")
  
  return(list(seqData = seqData, nullmod = nullmod, sample.selection = sample.selection))
}

# Main execution
main <- function() {
  message(glue("-----------------------------------------------------"))
  message(glue("chr: {chr}"))
  message(glue("AF.max: {opt$AF.max}"))
  message(glue("test: {opt$test}"))
  message(glue("qualifying variant set: {opt$qv}"))
  message(glue("ncores: {opt$ncores}"))
  
  if (is.null(opt$res.dir)) {
    opt$res.dir <- glue("~/mtrv/resources/res/{opt$qv}.AF.max.{opt$AF.max}_chr{chr}/")
  }
  message(glue("Results directory: {opt$res.dir}"))
  
  # Load precomputed data
  pre_data <- precompute(chr)
  seqData <- pre_data$seqData
  nullmod <- pre_data$nullmod
  sample.selection <- pre_data$sample.selection
  
  # Get GTF data
  gtf <- get_gtf(chr)
  
  # Process gene selection
  if (!is.na(opt$genes)) {
    if (opt$genes == "non-tested") {
      genes.tested <- read_tsv(opt$genes.tested, show_col_types = FALSE) %>%
        filter(tested == 0)
      genes <- genes.tested$gene
      message(glue("Genes selected will be non-tested genes for chr{chr}."))
      message(glue("Genes selected: {paste(genes, collapse = ', ')}"))
    } else {
      genes <- strsplit(opt$genes, split = ",")[[1]]
      message(glue("Genes selected: {paste(genes, collapse = ', ')}"))
    }
    gtf <- gtf[names(gtf) %in% genes]
    genes.tested <- tibble(gene = names(gtf), tested = 1)
  }
  
  # Get qualifying variants
  var <- get_qual_var(opt$qv, chr)
  
  # Set sample & variant filters
  seqData <- setFilter(seqData, sample.selection, chr, var, opt$AF.max)
  
  # Create iterator
  iterator <- SeqVarRangeIterator(seqData, variantRanges = gtf, verbose = FALSE)
  nvar <- sapply(iterator@variantFilter, length)
  message("Summary of variant block sizes:")
  print(summary(nvar))
  
  # Filter genes based on max.nvar
  if (!is.infinite(opt$max.nvar) && any(nvar > opt$max.nvar)) {
    message(glue("\n>> IMPT NOTE: removing genes with nvar > max.nvar threshold = {opt$max.nvar} <<\n"))
    message(glue("Variant Block sizes > max nvar: {paste(nvar[nvar > opt$max.nvar], collapse = ' ')}"))
    gtf <- gtf[nvar <= opt$max.nvar]
    genes.tested <- tibble(gene = names(gtf), tested = 1)
    
    # Reset seqData and create new iterator
    seqData <- setFilter(seqData, sample.selection, chr, var, opt$AF.max)
    iterator <- SeqVarRangeIterator(seqData, variantRanges = gtf, verbose = FALSE)
    nvar <- sapply(iterator@variantFilter, length)
    message("Summary of new variant block sizes:")
    print(summary(nvar))
  }
  
  # Run association test
  message(glue("---Chromosome: {chr} | AF.max: {opt$AF.max} | QV: {opt$qv}---"))
  message(glue("Number of genes to test: {length(iterator@variantRanges)}"))
  
  plan(multiprocess, workers = opt$ncores)
  
  assoc <- future_map(iterator@variantRanges, ~{
    assocTestAggregate(
      iterator,
      nullmod,
      test = opt$test,
      verbose = TRUE
    )
  }, .progress = TRUE)
  
  # Process results
  assoc[[1]] <- tibble(assoc[[1]]) %>%
    mutate(
      unit = seq(1, nrow(assoc[[1]])),
      cfreq = sapply(assoc[[2]], function(x) sum(x$freq)),
      gene = names(gtf)
    )
  
  # Save results
  dir.create(opt$res.dir, showWarnings = FALSE, recursive = TRUE)
  
  result_file <- glue("{opt$res.dir}/rv_assoc.aggr.{opt$qv}.{opt$test}.AF.max.{opt$AF.max}.chr{chr}.txt")
  write_tsv(assoc[[1]], file = result_file, append = TRUE)
  
  rdata_file <- glue("{opt$res.dir}/rv_assoc.aggr.{opt$qv}.{opt$test}.AF.max.{opt$AF.max}.chr{chr}.RData")
  save(assoc, file = rdata_file)
  
  genes_tested_file <- glue("{opt$res.dir}/rv_assoc.aggr.{opt$qv}.{opt$test}.AF.max.{opt$AF.max}.chr{chr}.genes_tested.txt")
  write_tsv(genes.tested, file = genes_tested_file, append = TRUE)
  
  message("-----------------------------------------------------")
}

# Run the main function
main()

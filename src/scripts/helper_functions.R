# helper_functions.R
# Optimized helper functions for genetic analysis

library(tools)
library(GenomeInfoDb)
library(dplyr)
library(SeqArray)
library(Biobase)
library(SeqVarTools)

#' Load R objects from RData or RDS files
#'
#' @param Rdata Path to the R data file
#' @return Loaded R object
#' @export
getobj <- function(Rdata) {
  if (tolower(tools::file_ext(Rdata)) == "rds") {
    return(readRDS(Rdata))
  }
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning("Multiple objects stored in file ", Rdata, ". Returning only the first object.")
  }
  get(objname[1])
}

#' Set filters on GDS file
#'
#' @param seqData SeqVarData object
#' @param sample_selection Vector of sample IDs
#' @param chr Chromosome number
#' @param var Vector of variant positions
#' @param AF_max Maximum allele frequency
#' @param gene Gene name (optional)
#' @param reset Boolean to reset filters
#' @return Filtered SeqVarData object
#' @export
setFilter <- function(seqData, sample_selection, chr, var, AF_max, gene = NULL, reset = TRUE) {
  if (reset) SeqArray::seqResetFilter(seqData)
  
  message(sprintf("Setting sample & variant filters [qualifying variants, then MAF filter of [singletons=1/(2*N) .. to .. MAF < %.2f%%]]", AF_max * 100))
  
  if (!is.null(gene)) {
    message("Setting filter to ", gene, " gene boundaries.")
    gtf <- get_gtf(chr)
    gtf <- gtf[strsplit(gene, split = ',')[[1]]]
    seqSetFilterChrom(seqData, include = chr, from.bp = start(gtf), to.bp = end(gtf))
  } else {
    seqSetFilterChrom(seqData, include = chr)
  }
  
  samples <- seqGetData(seqData, 'sample.id')
  seqSetFilter(seqData, sample.sel = which(samples %in% sample_selection), verbose = TRUE)
  
  pos <- seqGetData(seqData, 'position')
  variant_id <- seqGetData(seqData, 'variant.id')
  seqSetFilter(seqData, variant.id = variant_id[pos %in% var])
  
  seqSetFilterCond(seqData, 
                   maf = c(1 / (2 * length(sample_selection)), AF_max), 
                   .progress = FALSE, 
                   verbose = TRUE, 
                   parallel = max(1, parallel::detectCores() - 2))
  
  seqData
}

#' Get GTF data for a chromosome
#'
#' @param chr Chromosome number
#' @return GTF data for the chromosome
#' @export
get_gtf <- function(chr) {
  message('Loading gtf')
  
  gtf_path <- "resources/gtf/gencode_v38.chr%s_protein_coding_genes.gtf.rds"
  
  if (chr == '23') {
    gtf <- getobj(sprintf(gtf_path, 'X'))
    GenomeInfoDb::seqlevels(gtf) <- c('X' = "23")
  } else if (chr == '24') {
    gtf <- getobj(sprintf(gtf_path, 'Y'))
    GenomeInfoDb::seqlevels(gtf) <- c('Y' = "24")
  } else {
    gtf <- getobj(sprintf(gtf_path, chr))
  }
  
  # Chromosome-specific exclusions
  if (chr == 2) {
    gtf <- gtf[names(gtf) != 'TTN']
  } else if (chr == 5) {
    gtf <- gtf[!stringr::str_detect(names(gtf), 'PCDH')]
  }
  
  return(gtf)
}

#' Get qualifying variants
#'
#' @param qv Qualifying variant type
#' @param chr Chromosome number
#' @return Vector of qualifying variant positions
#' @export
get_qual_var <- function(qv, chr) {
  qv_dir <- 'resources/qual.var/'
  chr_ <- ifelse(chr == '23', 'X', chr)
  
  file_pattern <- paste0(qv_dir, '/', qv, '/ukb23149_c', chr_, 
                         '_b0_v1.pvar.vep.out.%s.uniqpos')
  
  file_path <- switch(qv,
                      "nonsyn.impact_mod_high" = sprintf(file_pattern, "impact_mod_high"),
                      "nonsyn.impact_mod_high_CADD18" = sprintf(file_pattern, "impact_mod_high_CADD18"),
                      "nonsyn.impact_high_pLOF" = sprintf(file_pattern, "impact_high_pLOF"),
                      stop("Invalid qv parameter"))
  
  as.integer(readLines(file_path))
}

#' Precompute data for analysis
#'
#' @param chr Chromosome number
#' @return List containing seqData, nullmod, and sample_selection
#' @export
precompute <- function(chr) {
  message('Loading precomputed data...')
  
  sample_selection <- readLines('resources/n428k.exomes.nonwhites_nonconsent.removed.sample.IIDs.list')
  
  seqData <- seqOpen(paste0('resources/gds/UKBexomeQQFE_chr', chr, '.gds'), allow.duplicate = TRUE)
  samples <- seqGetData(seqData, 'sample.id')
  
  a <- getobj('resources/n450k.Full_unfiltered.exomes_covar_AnnotatedDataFrame.rds')
  a_data <- a@data %>% 
    as_tibble() %>% 
    arrange(match(sample.id, samples))
  
  a <- Biobase::AnnotatedDataFrame(data = as.data.frame(a_data), varMetadata = a@varMetadata)
  
  seqData <- SeqVarData(seqData, sampleData = a)
  
  sample_selection <- a_data %>% 
    filter(sample.id %in% sample_selection, !is.na(arrayCN_PCAodd_m2)) %>% 
    pull(sample.id)
  
  nullmod <- getobj('resources/n450k.exomes.w_nonwhites_removed.nullmod.re.rds')
  
  list(seqData = seqData, nullmod = nullmod, sample_selection = sample_selection)
}

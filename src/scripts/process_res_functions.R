# Process association aggregate results
suppressPackageStartupMessages({
  library(TopmedPipeline)
  library(dplyr)
  library(ggplot2)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GeneAnswers)
})

source('src/scripts/plots.R')

#' Process single chromosome results
#'
#' @param assoc.dir Directory containing association results
#' @param test Test type
#' @param chr Chromosome number
#' @param AF.max Maximum allele frequency
#' @param txdb TxDb object (optional)
#' @return List containing results and variants
#' @export
process_single_chrom <- function(assoc.dir, test, chr, AF.max, txdb = NULL) {
  file_path <- file.path(assoc.dir, sprintf("assoc_aggregate.%s.AF.max.%s.%s.results.RData", test, AF.max, chr))
  r <- getobj(file_path)
  res <- tibble(r[[1]]) %>% mutate(unit = row_number())
  var <- r[[2]]
  list(res = res, var = var)
}

#' Process and arrange results for a single chromosome
#'
#' @param assoc.dir Directory containing association results
#' @param test Test type
#' @param chr Chromosome number
#' @param AF.max Maximum allele frequency
#' @param txdb TxDb object (optional)
#' @param order_by Column to order results by
#' @param print.var Whether to print variant information
#' @return List containing processed results and variants
#' @export
process_and_arrange_results <- function(assoc.dir, test, chr, AF.max = '0.01', txdb = NULL, order_by = NULL, print.var = FALSE) {
  assoc <- process_single_chrom(assoc.dir, test, chr, AF.max, txdb)
  res <- assoc$res
  var <- assoc$var
  
  # Arrange results based on test type or specified column
  if (!is.null(order_by)) {
    res <- res %>% arrange(!!sym(order_by))
  } else {
    res <- switch(test,
                  'Burden' = res %>% arrange(Score.pval),
                  'SKAT' = res %>% arrange(pval),
                  'fastSKAT' = res %>% arrange(pval),
                  'SMMAT' = res %>% arrange(pval_SMMAT),
                  stop("Invalid test type"))
  }
  
  unit.of.interest <- res$unit[1]
  
  if (test == 'SMMAT') {
    print_min_pvals(res)
  } else {
    message('Unit corresponding to lowest p-val: ', unit.of.interest)
  }
  
  if (print.var) {
    print(tibble(var[[unit.of.interest]]) %>% arrange(desc(MAC)), n = 5)
  }
  
  list(res = res, var = var)
}

#' Print minimum p-values for SMMAT test
#'
#' @param res Results data frame
#' @export
print_min_pvals <- function(res) {
  min_pvals <- c(burden = min(res$pval_burden, na.rm = TRUE),
                 theta = min(res$pval_theta, na.rm = TRUE),
                 SMMAT = min(res$pval_SMMAT, na.rm = TRUE))
  
  for (test in names(min_pvals)) {
    cat(sprintf('Lowest p-val_%s: %g | Unit: %d\n', 
                test, 
                min_pvals[test], 
                res$unit[which.min(res[[paste0("pval_", test)]])]))
  }
}

#' Process results for all chromosomes
#'
#' @param assoc.dir Directory containing association results
#' @param test Test type
#' @param AF.max Maximum allele frequency
#' @param txdb TxDb object (optional)
#' @param cfreq.cutoff Cutoff for cumulative frequency
#' @param plot Whether to generate QQ plots
#' @return List containing all results and variants
#' @export
process_all_chromosomes <- function(assoc.dir = './', test = 'SMMAT', AF.max = 1e-02, txdb = NULL, cfreq.cutoff = NULL, plot = FALSE) {
  r.all <- list()
  v.all <- list()
  
  if (plot) {
    pdf(paste0(basename(assoc.dir), '_qqplots.pdf'))
    on.exit(dev.off())
  }
  
  for (chr in 1:22) {
    message('\n------chr: ', chr, '-------')
    res_file <- file.path(assoc.dir, sprintf("assoc_aggregate.%s.AF.max.%s.%s.results.RData", test, AF.max, chr))
    
    if (file.exists(res_file)) {
      message('res file: ', res_file)
      r <- process_and_arrange_results(assoc.dir, test, chr, AF.max, txdb)
      res <- r$res
      
      if (!is.null(txdb)) {
        res <- annotate_with_txdb(res, txdb, chr)
      }
      
      r.all[[chr]] <- res
      v.all[[chr]] <- r$var
      
      if (plot) {
        qqPlot(pval = res$pval_SMMAT,
               title = paste0('p-value SMMAT for chr: ', chr),
               subtitle = basename(assoc.dir))
      }
    } else {
      message('!!Chr res file not found. Skipping..')
    }
  }
  
  r.all <- bind_rows(r.all)
  
  if (!is.null(cfreq.cutoff)) {
    r.all <- r.all %>% 
      filter(cfreq >= cfreq.cutoff) %>% 
      arrange(pval_SMMAT)
  } else {
    r.all <- r.all %>% arrange(pval_SMMAT)
  }
  
  list(results = r.all, variants = v.all)
}

#' Annotate results with TxDb information
#'
#' @param res Results data frame
#' @param txdb TxDb object
#' @param chr Chromosome number
#' @return Annotated results data frame
#' @export
annotate_with_txdb <- function(res, txdb, chr) {
  seqlevels(txdb) <- paste0('chr', chr)
  genes <- sort(genes(txdb))
  genes <- renameSeqlevels(x = genes, value = as.character(chr))
  
  res %>%
    mutate(entrez.id = names(genes),
           chr = chr,
           start = start(genes),
           end = end(genes))
}

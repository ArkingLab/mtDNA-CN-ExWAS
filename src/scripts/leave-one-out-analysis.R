# Leave-One-Out Analysis for Rare Variant Burden Testing
# Author: Vamsee Pillalamarri

library(optparse)
library(SeqArray)
library(SeqVarTools)
library(GenomicRanges)
library(dplyr)
library(GENESIS)
library(CompQuadForm)
library(survey)
library(ggplot2)
library(hrbrthemes)
library(ggpubr)

# Helper functions
source("src/helper_functions.R")

#' Main function to conduct leave-one-out analysis
#'
#' @param gene Name of gene
#' @param chr Chromosome number
#' @param AF.max Maximum allele frequency
#' @param qv Qualifying variant type
#' @return Variant-level results of leave-one-out analysis
#' @export
loo <- function(gene, chr, AF.max, qv = "nonsyn.impact_mod_high", cores = 1) {
  p <- precompute(chr)
  seqData <- p$seqData
  nullmod <- p$nullmod
  sample_selection <- p$sample_selection
  
  qvar <- get_qual_var(qv = qv, chr = chr)
  
  gds <- setFilter(seqData, sample_selection = sample_selection, chr = chr, 
                   var = qvar, AF.max = AF.max, reset = TRUE)
  
  gtf <- get_gtf(chr)
  gtf <- gtf[strsplit(gene, split = ",")[[1]]]  # set gene
  
  pos <- seqGetData(gds, "position")
  var <- seqGetData(gds, "variant.id")
  mac <- SeqVarTools::alleleCount(gds, parallel = cores, n = 1)
  
  res <- pbapply::pblapply(seq_along(var), function(i) {
    message(sprintf("i: %d / %d -- Leaving out: chr%s:%d | MAC: %d", 
                    i, length(var), chr, pos[i], mac[i]))
    var_test <- rep(TRUE, length(var))
    var_test[i] <- FALSE  # Leave one out
    
    seqSetFilter(gds, action = "push+intersect", variant.sel = var_test)
    
    iterator <- SeqVarRangeIterator(gds, variantRanges = gtf, verbose = FALSE)
    assoc <- assocTestAggregate(
      iterator,
      nullmod,
      test = "SMMAT",
      verbose = TRUE,
      BPPARAM = BiocParallel::SerialParam(progressbar = TRUE, stop.on.error = FALSE)
    )
    
    seqSetFilter(gds, action = "pop")
    seqSetFilter(gds, action = "pop")
    
    c(assoc[[1]], list(variant = i, var = var[i], pos = pos[i], mac = mac[i]))
  }, cl = cores)
  
  res <- tibble::as_tibble(dplyr::bind_rows(res))
  return(res)
}

#' Function to plot LOO results
#'
#' @param gene Name of gene
#' @param res_loo Results from loo()
#' @param plot_type Plot type: SMMAT, Burden, or SKAT
#' @param save_dir Directory to save plot
#' @export
plot_loo <- function(gene, res_loo, plot_type = NULL, save_dir = NULL) {
  g1 <- ggplot(res_loo, aes(x = seq_len(nrow(res_loo)), y = -log10(pval_SMMAT))) +
    geom_line(alpha = 2/3) +
    labs(x = "Variant Left Out in Aggregate Test", y = "-log10 P-Value SMMAT") +
    theme_minimal() +
    theme(text = element_text(size = 8))
  
  g2 <- ggplot(res_loo, aes(x = seq_len(nrow(res_loo)), y = -log10(pval_burden))) +
    geom_line(alpha = 2/3) +
    labs(y = "-log10 P-Value Burden") +
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), text = element_text(size = 8))
  
  g3 <- ggplot(res_loo, aes(x = seq_len(nrow(res_loo)), y = -log10(pval_theta))) +
    geom_line(alpha = 2/3) +
    labs(x = "Variant Left Out in Aggregate Test", y = "-log10 P-Value SKAT", 
         caption = gene) +
    theme_minimal() +
    theme(text = element_text(size = 8))
  
  if (!is.null(plot_type)) {
    g <- switch(plot_type,
                SMMAT = g1,
                Burden = g2,
                SKAT = g3,
                stop("Invalid plot type"))
  } else {
    g <- ggarrange(g1, g2, g3, ncol = 1, nrow = 3, align = "v", labels = "AUTO")
  }
  
  if (!is.null(save_dir)) {
    ggsave(plot = g, filename = file.path(save_dir, paste0("loo_plot.", gene, ".pdf")), 
           device = "pdf", scale = 1)
  }
  
  return(g)
}

# Main function
main <- function(opts) {
  res_loo <- loo(gene = opts$gene, chr = opts$chr, AF.max = opts$af_max, 
                 qv = opts$qv, cores = opts$cores)
  
  if (!is.null(opts$output)) {
    saveRDS(res_loo, file = opts$output)
  }
  
  if (opts$plot) {
    plot <- plot_loo(gene = opts$gene, res_loo = res_loo, 
                     plot_type = opts$plot_type, save_dir = opts$save_dir)
    print(plot)
  }
}

# Parse command line options
option_list <- list(
  make_option(c("-g", "--gene"), type = "character", help = "Gene name"),
  make_option(c("-c", "--chr"), type = "character", help = "Chromosome number"),
  make_option(c("-a", "--af_max"), type = "numeric", default = 0.01, help = "Maximum allele frequency"),
  make_option(c("-q", "--qv"), type = "character", default = "nonsyn.impact_mod_high", 
              help = "Qualifying variant type"),
  make_option(c("-o", "--output"), type = "character", help = "Output file for results"),
  make_option(c("-p", "--plot"), action = "store_true", default = FALSE, help = "Generate plot"),
  make_option(c("-t", "--plot_type"), type = "character", help = "Plot type (SMMAT, Burden, or SKAT)"),
  make_option(c("-s", "--save_dir"), type = "character", help = "Directory to save plot"),
  make_option(c("-n", "--cores"), type = "integer", default = 1, help = "Number of cores to use")
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Run main function
main(opts)

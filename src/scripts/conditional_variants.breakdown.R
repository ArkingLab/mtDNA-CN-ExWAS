# Conditional Variant-level Breakdown of Rare variant Burden tests
# Author: Vamsee Pillalamarri

library(dplyr)
library(SeqArray)
library(SeqVarTools)
library(qvalue)
library(tools)

# Helper function to load R objects
getobj <- function(Rdata) {
  if (tolower(file_ext(Rdata)) == "rds") {
    return(readRDS(Rdata))
  }
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata, "\nReturning only the first object"))
  }
  return(get(objname))
}

# Function to perform conditional variant breakdown of rare variant burden signal
# Returns conditional variant-level effect sizes for variants within a given gene
# Collapse rare vars w/ MAC < 6 into single variable for testing
var.breakdown2 <- function(gds, t, gene) {
  seqResetFilter(gds)
  var.id <- seqGetData(gds, 'variant.id')
  samp.id <- seqGetData(gds, 'sample.id')

  k <- getobj('~/mtrv/n428k.exomes.nonwhites_nonconsent.removed.AnnotatedDataFrame.rds')

  x <- t[[2]][gene][[1]]
  mac_lt_6_idx <- which(x$MAC < 6)
  mac_gt_6_idx <- which(x$MAC >= 6)
  
  message('# Variants with MAC >= 6: ', length(mac_gt_6_idx))
  message('# Variants with MAC < 6: ', length(mac_lt_6_idx))
  message('Table of variants with MAC < 6: ')
  print(table(x$MAC[mac_lt_6_idx]))

  # Get genotype matrix MAC >= 6
  seqSetFilter(gds, sample.sel = samp.id %in% k@data$sample.id, variant.sel = var.id %in% x[mac_gt_6_idx,]$variant.id)
  geno_gt_6 <- altDosage(gds)
  colnames(geno_gt_6) <- paste0('var.', x[mac_gt_6_idx,]$variant.id,
                                '.pos.', paste0(x[mac_gt_6_idx,]$chr, '_', x[mac_gt_6_idx,]$pos),
                                '_mac_', x[mac_gt_6_idx,]$MAC)
  geno_gt_6 <- tibble::rownames_to_column(as.data.frame(geno_gt_6), "rowname") %>% as_tibble()
  geno_gt_6$rowname <- as.numeric(geno_gt_6$rowname)

  # Get genotype matrix MAC < 6
  seqResetFilter(gds)
  seqSetFilter(gds, sample.sel = samp.id %in% k@data$sample.id, variant.sel = var.id %in% x[mac_lt_6_idx,]$variant.id)
  geno_lt_6 <- altDosage(gds)
  colnames(geno_lt_6) <- paste0('var.', x[mac_lt_6_idx,]$variant.id,
                                '.pos.', paste0(x[mac_lt_6_idx,]$chr, '_', x[mac_lt_6_idx,]$pos),
                                '_mac_', x[mac_lt_6_idx,]$MAC)
  geno_lt_6 <- tibble::rownames_to_column(as.data.frame(geno_lt_6), "rowname") %>% as_tibble()
  geno_lt_6$rowname <- as.numeric(geno_lt_6$rowname)

  # Collapse rare vars w/ MAC < 6 into single variable
  mac_lt_6 <- tibble(rowname = geno_lt_6$rowname, mac_lt_6 = rowSums(geno_lt_6[,-1], na.rm = TRUE))

  # Combine phenotype, covariates & genotype data in single tibble
  k <- tibble(k@data) %>% 
    inner_join(geno_gt_6, by = c('sample.id' = 'rowname')) %>% 
    inner_join(mac_lt_6, by = c('sample.id' = 'rowname')) %>%
    select(-arrayCN_PCAodd, -arrayCN_PCAodd_m1, -arrayCN_PCAodd_m3)

  k$sex <- as.factor(k$sex)
  k$wes.batch <- as.factor(k$wes.batch)

  # Fit null & full models
  k <- k %>% filter(!is.na(arrayCN_PCAodd_m2))
  k_1 <- k[,1:45] # for null
  k_2 <- k[,c(1, 46:ncol(k))] %>% select(-haplo) # for full
  
  message('n.samples: ', nrow(k))
  fit.null <- lm('arrayCN_PCAodd_m2 ~ . -sample.id', data = k_1)

  k_2$mtDNA_CN <- k$arrayCN_PCAodd_m2
  k_2$mtDNA_CN.resid <- resid(fit.null)

  fit.full <- lm('mtDNA_CN.resid ~ . -mtDNA_CN -sample.id', data = k_2)

  singularities <- summary(fit.full)$aliased
  if(any(singularities)) {
    message('WARNING: ', sum(singularities), 
            ' singularities exist, check whether some variants are multi-colinear. 
            Could be due to multi-allelic variants. NOTE: THESE ARE DROPPED from results!')
  }

  coef_summary <- summary(fit.full)$coef[-1,]
  var.names <- rownames(coef_summary)
  beta <- coef_summary[,1]
  se <- coef_summary[,2]
  p <- coef_summary[,4]
  q <- qvalue(p, pi0 = 1, fdr.level = 0.1)
  
  res <- tibble(
    var = var.names,
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p), 
    q = as.numeric(q$qvalues),
    fdr_sig = q$significant
  ) %>% mutate(bonf.sig = p <= (0.05 / (length(var.names))))
  
  seqClose(gds)

  return(list('var.list' = x, 'geno_gteq6' = k_2, 'fit.full' = fit.full, 'var.bkdwn' = res, 'geno_lt6' = geno_lt_6))
}

# Function to get synthetic allele
get.synthetic_allele <- function(x, level = 'all', include_lt6 = FALSE) {
  v <- x$var.bkdwn %>% arrange(p)
  g <- x$geno_gteq6 %>% select(sample.id, contains('var'), mac_lt_6)
  
  include_lt6 <- level == 'mac_lt_6' || include_lt6
  if (!include_lt6) g <- g %>% select(-mac_lt_6)
  
  process_variants <- function(vars, betas = NULL) {
    if (length(vars) == 0) return(NULL)
    
    message('Variants: ', paste(vars, collapse = ', '))
    vars_geno <- g %>% select(sample.id, all_of(vars))
    
    syn_allele <- apply(vars_geno[,-1], 1, function(x) any(x == 1, na.rm = TRUE)) %>% 
      as.integer()
    
    result <- tibble(IID = vars_geno$sample.id, syn.allele = syn_allele)
    
    if (!is.null(betas)) {
      scored_geno <- vars_geno
      for (i in seq_along(betas)) {
        scored_geno[,i+1] <- scored_geno[,i+1] * betas[i]
      }
      list(syn.allele = result, alleles = vars_geno, scored.alleles = scored_geno)
    } else {
      list(syn.allele = result, alleles = vars_geno, scored.alleles = NULL)
    }
  }

  result <- switch(
    level,
    bonf = process_variants(v$var[v$bonf.sig == 1], v$beta[v$bonf.sig == 1]),
    fdr = process_variants(v$var[v$fdr_sig], v$beta[v$fdr_sig]),
    mac_lt_6 = list(syn.allele = tibble(IID = g$sample.id, syn.allele = g$mac_lt_6), 
                    alleles = NULL, scored.alleles = NULL),
    all = process_variants(v$var, v$beta)
  )
  
  if (is.null(result)) {
    message("No variants at level: ", level)
    return(list(syn.allele = NA, alleles = NA, scored.alleles = NA))
  }
  
  result
}

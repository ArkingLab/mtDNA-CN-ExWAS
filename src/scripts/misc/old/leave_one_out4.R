# Leave one out analysis functions for a single gene
loo <- function(gene, chr, AF.max){
  require(SeqArray)
  require(SeqVarTools)
  require(GenomicRanges)
  require(dplyr)
  require(GENESIS)
  require(CompQuadForm)
  require(survey)
  
  # Helper function s
  getobj <- function (Rdata){
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
  setFilter <- function(seqData, sample.selection, chr, var, AF.max, gene=NULL, reset=T){
    # sets sample / variant filters on gds
    if(reset){SeqArray::seqResetFilter(seqData)}
    message('Setting sample & variant filters ',
            '[qualifying variants, then MAF filter of ',
            '[singletons=1/(2*N) .. to .. MAF < ',AF.max*100,'%]]..')
    
    # gene filter (optional)
    if(!is.null(gene)){
      message('Setting filter to ', gene, ' gene boundaries.')
      gtf <- get.gtf(chr)
      gtf <- gtf[strsplit(gene, split = ',')[[1]]] # set gene
      seqSetFilterChrom(seqData, include=chr,from.bp=start(gtf),to.bp=end(gtf))
    } else {
      seqSetFilterChrom(seqData, include=chr)
    }
    
    # sample filter
    samples <- seqGetData(seqData, 'sample.id')
    seqSetFilter(seqData,
                 sample.sel = which(samples %in% sample.selection),
                 verbose = T)
    
    # variant filters
    # seqSetFilterChrom(seqData, include = chr, verbose = T)
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
  
  #tmp <- seqAlleleCount(seqData, minor=T, parallel=14)
  
  p <- precompute()
  seqData <- p[[1]]
  nullmod <- p[[2]]
  sample.selection <- p[[3]]
  
  # get qual var
  qvar <- get.qual_var(qv.dir='~/mtrv/resources/qual.var/', qv='nonsyn.impact_mod_high', chr)
  
  # set sample & variant filters
  gds <- setFilter(seqData, sample.selection=sample.selection, chr=chr, var=qvar, 
                       AF.max=AF.max, reset=T)
  
  # gtf, select gene
  gtf <- get.gtf(chr)
  gtf <- gtf[strsplit(gene, split = ',')[[1]]] # set gene
  
  # Association testing using leave-one-out
  pos <- seqGetData(gds, 'position')
  var <- seqGetData(gds, 'variant.id')
  mac <- SeqVarTools::alleleCount(gds, parallel=(parallel::detectCores()-2), n=1)
  res <- list()
  nvar <- length(var)
  for (i in 1:length(var)) {
    message('i: ', i, ' / ', nvar,' -- Leaving out: chr', chr, ':', pos[i], ' | MAC: ', mac[i])
    var.test <- rep(TRUE, length(var))
    var.test[i] <- FALSE # Leave one out
    
    seqSetFilter(gds, action = 'push+intersect', variant.sel = var.test)
    
    # do a aggregate association test on the rare variants in each gene
    # after leaving one variant out
    iterator <-
      SeqVarRangeIterator(gds, variantRanges = gtf, verbose = F)
    assoc <-
      assocTestAggregate(
        iterator,
        nullmod,
        test = 'SMMAT',
        verbose = T,
        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE,
                                            stop.on.error = FALSE)
      )
    # assoc <- assocTestAggregate(
    #   iterator,
    #   nullmod,
    #   test = 'SMMAT',
    #   genome.build = 'hg38',
    #   verbose = F
    # )
    res[[i]] <- assoc
    message('p-values (burden, adj-SKAT, SMMAT): ', paste0(
      c(assoc[[1]]$pval_burden, assoc[[1]]$pval_theta, assoc[[1]]$pval_SMMAT),
      collapse = ' , '
    ))
    seqSetFilter(gds, action = 'pop')
    seqSetFilter(gds, action = 'pop') # Magnitude.
  }
  
  res <- tibble::as_tibble(dplyr::bind_rows(lapply(res, function(x) x[[1]])) %>%
                             mutate(variant=1:length(pos), var=var, pos=pos, mac=mac))
  # 
  # # - save to file
  # savefile = paste0(save.dir,
  #                   '/',
  #                   gene,
  #                   '.',
  #                   qv,
  #                   '.',
  #                   sprintf(as.character(AF.max)),
  #                   '.LOO.rds')
  # saveRDS(object = res, file = savefile)
  return(res)  
}
plot.loo <- function(gene, res.loo, plot.type=NULL, save.dir=NULL){
  require(ggplot2)
  require(dplyr)
  require(hrbrthemes)
  require(ggpubr)
  
  g1 <- res.loo %>%
    ggplot(aes(x=seq(1,nrow(res.loo)), y=-log10(pval_SMMAT))) +
    # ggplot(aes(x=seq(1,nrow(res.loo)), y=-log10(pval_burden))) +
    # ggplot(aes(x = seq(1, nrow(res.loo)), y = -log10(pval_theta))) +
    geom_line(alpha = 2 / 3) +
    # geom_point(alpha = 1 / 2, size = 0.75) +
    xlab('Variant Left Out in Aggregate Test') +
    ylab('-log10 P-Value SMMAT') +
    # theme_ipsum_rc() +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size=8))
  
  g2 <- res.loo %>%
    # ggplot(aes(x=seq(1,nrow(res.loo)), y=-log10(pval_SMMAT))) +
    ggplot(aes(x=seq(1,nrow(res.loo)), y=-log10(pval_burden))) +
    # ggplot(aes(x = seq(1, nrow(res.loo)), y = -log10(pval_theta))) +
    geom_line(alpha = 2 / 3) +
    ylab('-log10 P-Value Burden') +
    # geom_point(alpha = 1 / 2, size = 0.75) +
    # xlab('Variant Left Out in Aggregate Test') +
    # ggtitle(gene) +
    # theme_ipsum_rc() +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size=8))
  
  g3 <- res.loo %>%
    # ggplot(aes(x=seq(1,nrow(res.loo)), y=-log10(pval_SMMAT))) +
    # ggplot(aes(x=seq(1,nrow(res.loo)), y=-log10(pval_burden))) +
    ggplot(aes(x = seq(1, nrow(res.loo)), y = -log10(pval_theta))) +
    geom_line(alpha = 2 / 3) +
    # geom_point(alpha = 1 / 2, size = 0.75) +
    theme_minimal() +
    theme(text = element_text(size=8)) +
    labs(caption=gene) +
    xlab('Variant Left Out in Aggregate Test') +
    ylab('-log10 P-Value SKAT')
  # ggtitle(gene) +
  # theme_ipsum_rc()
  
  if(!is.null(plot.type)){
    if(plot.type=='SMMAT'){
      g <- g1 + xlab('Variant Left Out in Aggregate Test')
    } else if(plot.type=='Burden'){
      g <- g2 + xlab('Variant Left Out in Aggregate Test')
    } else if(plot.type=='SKAT'){
      g <- g3 + xlab('Variant Left Out in Aggregate Test')
    }
  } else {
    g <- ggarrange(g1,g2,g3, ncol=1, nrow=3, align='v', labels='AUTO')
  }
  
  if(!is.null(save.dir)){
    ggsave(plot = g, filename=paste0(save.dir,'/loo_plot.', gene,'.pdf'), device = 'pdf', scale = 1)
  }
  g
  return(g)
}

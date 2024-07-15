# what does this code do? Run phewas?
# genes <- c('SAMHD1','TWNK','TFAM','JAK2','CLPX','LY75.CD302','AK2','MGME1','CHEK2','CAVIN2')
# chrs <- c(20, 10, 10, 9, 15, 2, 1, 20, 22, 2)
genes <- c('AK2','MGME1','CHEK2','CAVIN2')
chrs <- c(1, 20, 22, 2)
for(i in 1:length(genes)){
  # get gene variants
  x.gene <- var.breakdown2(gds=seqOpen(paste0('~/mtrv/resources/gds/UKBexomeQQFE_chr',chrs[i],'.gds'), allow.duplicate=T), 
                             t= getobj(paste0('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr',chrs[i],'.RData')),
                             gene=genes[i])
  pass_vars <- x.gene[[4]]$var[which((abs(x.gene[[4]]$beta) >= 0.3) & (x.gene[[4]]$fdr_sig==TRUE))]
  
  if(length(pass_vars) > 0){
    v.gene <- get.synthetic_allele(x.gene, level='fdr', include_lt6=T)
    vars <- v.gene$alleles %>% select(all_of(pass_vars))
    IID <- v.gene$alleles$sample.id
    vars$syn.allele <- apply(vars, 1, function(x) ifelse(any(x!=0, na.rm=T), 1, 0))
    vars$IID <- IID
    vars <- vars %>% select(IID, everything())
    
    # create phewas data input
    var.carriers <- vars %>% filter(syn.allele==1)
    data <- create_phewas_data.input(var.carriers, subset.to.unrelated=TRUE, subset.to.white.British.ancestry=FALSE)
    phecodes = data[[3]]
    sex.check=data[[2]]
    message('number of phecodes: ', length(phecodes))
    
    # Create a bigmatrix object from `data`
    require(bigmemory)
    require(biganalytics)
    data <- data[[1]]
    k2 <- bigmemory::as.big.matrix(as.data.frame(data))
    data_desc <- bigmemory::describe(k2)
    
    # cleanup
    rm(x.gene)
    
    # run pheWAS for gene
    phewas <- process.glm(phecodes=sample(phecodes), sex_check=sex.check, 
                          data=data_desc, ncores=5,var='var',
                          file.conn=paste0('./phewas.corr_cutoff_genes.fdr10pct_include_lt6.ESgt0.3.',genes[i],'.inFullExomesSet.txt'))
    }
  closeAllConnections()
}


# get list of fdr10% vars that are with abs(betas) > 0.3
genes <- c('SAMHD1','TWNK','TFAM','JAK2','CLPX','LY75.CD302','AK2','MGME1','CHEK2','CAVIN2')
chrs <- c(20, 10, 10, 9, 15, 2, 1, 20, 22, 2)
pass_vars <- list()
for(i in 1:length(genes)){
  message('Gene: ', genes[[i]])
  # get gene variants
  gds=seqOpen(paste0('~/mtrv/resources/gds/UKBexomeQQFE_chr',chrs[i],'.gds'), allow.duplicate=T)
  t=getobj(paste0('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr',chrs[i],'.RData'))
  x.gene <- var.breakdown2(gds, t, gene=genes[i])
  vars <- x.gene[[4]] %>% filter((abs(beta) >= 0.3) & (fdr_sig==TRUE))
  vars$gene <- rep(genes[[i]], nrow(vars))
  pass_vars[[i]] <- vars
}
pass_vars <- bind_rows(pass_vars)












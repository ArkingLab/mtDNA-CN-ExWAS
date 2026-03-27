


# gather Burden results for CADD18 MAF 0.01 run
for(i in 1:22){res.all_burden[[i]] <-
  TopmedPipeline::getobj(paste0('assoc_aggregate.Burden.AF.max.0.01.',i,'.results.RData'))[[1]]; message(i)}
res.all_burden <- bind_rows(res.all_burden)
res.all_burden <- res.all_burden %>% filter(!is.na(Score))

mito.carta <- read_tsv('/users/vpillala/mito_rare-variant/resources/gene_lists/Human_MitoPathways3.0.txt', col_names=T)

# Tests mean t-vaules between in set and out set of genes
test.enrich.mitocarta <- function(line, res){
  require(dplyr)
  require(stringr)
  
  # line, format: [id,link,gene1,gene2, .., geneN]
  line.split = str_split(line, ', ')[[1]]
  #line.split = line.split[3:length(line.split)]
  
  # in.stat <- res %>% filter(gene %in% line.split) %>% select(Score.Stat)
  # out.stat <- res %>% filter(gene %notin% line.split) %>% select(Score.Stat)
  
  # in.stat <- res %>% filter(gene %in% line.split) %>% select(z.burden)
  # out.stat <- res %>% filter(gene %notin% line.split) %>% select(z.burden)
  
  # in.stat <- res %>% filter(gene %in% line.split) %>% select(z.skat)
  # out.stat <- res %>% filter(gene %notin% line.split) %>% select(z.skat)
  
  in.stat <- res %>% filter(gene %in% line.split) %>% select(z.smmat)
  out.stat <- res %>% filter(gene %notin% line.split) %>% select(z.smmat)
  
  if(nrow(in.stat) > 1){
    # test <- t.test(abs(in.stat$z.burden), abs(out.stat$z.burden), var.equal=T)
    # test <- t.test(abs(in.stat$z.skat), abs(out.stat$z.skat), var.equal=T)
    test <- t.test(abs(in.stat$z.smmat), abs(out.stat$z.smmat), var.equal=T)
    return(c(nrow(in.stat), nrow(out.stat), test$estimate, test$p.value))
  } else {
    return(NA)
  }
}

k_1 <- k %>% filter(varset=='cds_nonsyn_cadd18_0.01') %>% filter(!is.na(err)) %>% filter(err==0)
# k_1 <- TopmedPipeline::getobj('k_1.rds')
# k_1$z.smmat <- abs(qnorm(k_1$pval_SMMAT, mean=0, sd=1))
# k_1$z.skat <- abs(qnorm(k_1$pval_theta, mean=0, sd=1))
# k_1$z.burden <- abs(qnorm(k_1$pval_burden, mean=0, sd=1))
k_1$z.smmat <- gCMAP::zScores(pval=k_1$pval_SMMAT)
k_1$z.burden <- gCMAP::zScores(pval=k_1$pval_burden)
k_1$z.skat <- gCMAP::zScores(pval=k_1$pval_theta)

k_2 <- k_1 %>% filter(gene != 'SAMHD1')
k_3 <- k_1
k_3$z.smmat[which(k_3$gene=='SAMHD1')] = max(k_2$z.smmat, na.rm=T)
k_3$z.burden[which(k_3$gene=='SAMHD1')] = max(k_2$z.burden, na.rm=T)
k_3$z.skat[which(k_3$gene=='SAMHD1')] = max(k_2$z.skat, na.rm=T)

res.all_burden_noSAMHD1 <- res.all_burden %>% filter(gene!='SAMHD1')

enrich.pval <- matrix(nrow = nrow(mito.carta), ncol = 5)
# [n1, n2, est1, est2, p.value]
for(i in 1:nrow(mito.carta)){
  message('i: ', i)
  # enrich.pval[i,] <- test.enrich.msigdb(line=msig[i], res=res.all_burden)
  enrich.pval[i,] <- test.enrich.mitocarta(line=mito.carta[i,3], res=res)
}

x <- tibble(cbind(mito.carta$MitoPathway, as.data.frame(enrich.pval)))
colnames(x) <- c('id', 'n.in','n.out','est.in','est.out','pvalue')
#x <- x %>% filter(!is.na(n.in))
x.burden_1 <- x
x.smmat <- x
x.burden_2 <- x
x.skat <- x




test.enrich <- function(gene.set, gene.set_name='msig',res){
  require(dplyr)
  require(stringr)
  
  `%notin%` <- Negate(`%in%`)
    
  enrich.pval <- matrix(nrow = nrow(gene.set), ncol = 5)
  
  for(i in 1:nrow(gene.set)){
    message('i: ', i)
    
    if(gene.set_name=='msig'){
      # line, format: [id,link,gene1,gene2, .., geneN]
      line=gene.set$value[i]
      line.split = str_split(line, '\t')[[1]]
      line.split = line.split[3:length(line.split)]
    } else {
      # gene.sets = data.frame(id, genes)
      # where `id` is name of gene set
      # `genes` is list of genes that are comma separated
      line=gene.set[i,2]
      line.split = str_split(line, ', ')[[1]]
    }
    
    in.stat <- res %>% filter(gene %in% line.split) %>% select(z.smmat)
    out.stat <- res %>% filter(gene %notin% line.split) %>% select(z.smmat)
    
    if(nrow(in.stat) > 1){
      test <- t.test(abs(in.stat$z.smmat), abs(out.stat$z.smmat), var.equal=T)
      enrich.pval[i,] <- c(nrow(in.stat), nrow(out.stat), test$estimate, test$p.value)
    } else {
      enrich.pval[i,] <- rep(NA, 5)
    }
  }
  return(enrich.pval)
}












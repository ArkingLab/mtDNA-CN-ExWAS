
# gather Burden results for CADD18 MAF 0.01 run
for(i in 1:22){res.all_burden[[i]] <-
  TopmedPipeline::getobj(paste0('assoc_aggregate.Burden.AF.max.0.01.',i,'.results.RData'))[[1]]; message(i)}
res.all_burden <- bind_rows(res.all_burden)
res.all_burden <- res.all_burden %>% filter(!is.na(Score))

msig <- readLines('~/mtrv/resources/msigdb.v7.4.symbols.gmt.txt')
msig.IDs <- readLines('~/mito_rare-variant/resources/gene_lists/msig.IDs')
disgenenet.ids <- readLines('~/mito_rare-variant/resources/gene_lists/disgenet.curated.v7.symbols.IDs')


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


res <- k %>% filter(err==0 & !is.na(err) & cfreq >= 1e-4) %>% arrange(pval_SMMAT) %>% distinct(gene, .keep_all=T)
res <- k %>% filter(err==0 & !is.na(err) & cfreq >= 1e-4) %>% arrange(pval_SMMAT) %>% filter(varset=='cds_nonsyn_cadd18_0.0001')
res$z.smmat <- qnorm(res$pval_SMMAT/2, 0, 1, lower.tail=FALSE)
res$z.skat <- qnorm(res$pval_theta/2, 0, 1, lower.tail=FALSE)
res$z.burden <- qnorm(res$pval_burden/2, 0, 1, lower.tail=FALSE)


# Tests mean t-vaules between in set and out set of genes
test.enrich <- function(gene.set, gene.set_name='msig',gene.set_ids, res){
  require(dplyr)
  require(stringr)
  
  `%notin%` <- Negate(`%in%`)
  
  # enrich.pval <- matrix(nrow = nrow(gene.set), ncol = 6)
  enrich.pval <- matrix(nrow = nrow(gene.set), ncol = 5)
  
  for(i in 1:nrow(gene.set)){
    message('i: ', i)
    
    if(gene.set_name=='msig'){
      # line, format: [id,link,gene1,gene2, .., geneN]
      line=gene.set$value[i]
      line.split = str_split(line, '\t')[[1]]
      line.split = line.split[3:length(line.split)]
    } else if (gene.set_name=='disgenenet'){
      line=gene.set$value[i]
      line.split = str_split(line, '\t')[[1]]
      line.split = line.split[3:length(line.split)]
    } else if (gene.set_name=='mitocarta') {
      line=gene.set$Genes[i]
      line.split=str_split(line, ', ')[[1]]
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
      tt <- t.test(abs(in.stat$z.smmat), abs(out.stat$z.smmat), var.equal=T)
      # wt <- wilcox.test(abs(in.stat$z.smmat), abs(out.stat$z.smmat), var.equal=T)
      # enrich.pval[i,] <- c(nrow(in.stat), nrow(out.stat), tt$estimate, tt$p.value, wt$p.value)
      enrich.pval[i,] <- c(nrow(in.stat), nrow(out.stat), tt$estimate, tt$p.value)
    } else {
      enrich.pval[i,] <- rep(NA, 5)
    }
  }
  
  x <- tibble(cbind(gene.set_ids, as.data.frame(enrich.pval)))
  # colnames(x) <- c('id', 'n.in','n.out','tt.est.in','tt.est.out', 'tt.pvalue','wt.pvalue')
  colnames(x) <- c('id', 'n.in','n.out','tt.est.in','tt.est.out', 'tt.pvalue')
  
  return(x)
}

a <- test.enrich(gene.set=msig, gene.set_name='msig', res)
a <- test.enrich(gene.set=r, gene.set_name='disgenenet', res)

#x <- x %>% filter(!is.na(n.in))
x.burden_1 <- x
x.smmat <- x
x.burden_2 <- x
x.skat <- x

# just the GO terms
msig.go
enrich.pval.go <- matrix(nrow = length(msig.go), ncol = 5)
# [n1, n2, est1, est2, p.value]
for(i in 1:length(msig)){
  message('i: ', i)
  # enrich.pval[i,] <- test.enrich.msigdb(line=msig[i], res=res.all_burden)
  enrich.pval[i,] <- test.enrich.msigdb(line=msig[i], res=k_3)
}

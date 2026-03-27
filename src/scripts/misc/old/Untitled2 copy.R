gtf <- list()
for(chr in 1:22){
  gtf[[chr]] <- getobj(paste0('~/mtrv/resources/gtf/gencode_v38.chr',chr,'_protein_coding_genes.gtf.rds'))
}
gtf <- unlist(as(gtf, 'GRangesList'))

genes <- readLines('genes.list')
genes <- genes[!(genes %in% res$gene)]
gtf <- gtf[genes] %>% as_tibble()

res2 <- list()
for(i in 1:nrow(gtf)){
  skip=F
  chr=gtf$seqnames[i]
  gene=gtf$gene_name[i]
  
  message('gene: ', gene)
  t <- getobj(paste0('~/mtrv/resources/res/rv_assoc.aggr.variants_in_cds.SMMAT.AF.max.0.01.chr',chr,'.RData'))
  gds <- seqOpen(paste0('resources/gds/UKBexomeQQFE_chr', chr, '.gds'), allow.duplicate=T)
  tryCatch(x <- var.breakdown2(gds, t, gene)[[4]], error = function(e){skip=T})
  if(skip){next} else{
    res2[[i]] <- 
      x %>% arrange(p) %>% slice(1) %>% mutate(gene=gene)
  }
}
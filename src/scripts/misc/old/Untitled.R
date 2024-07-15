res.chrs <- list()
for(i in seq(1,23)){
  message('chr: ', i)
  gtf <- TopmedPipeline::getobj(paste0('~/mtrv/burden/genesis/gds/gencode_v38.chr',i,'_protein_coding_genes.gtf.rds'))
  chr_res <- rep(NA, length(gtf))
  for (k in 1:length(gtf)) {
    gene=names(gtf[k])
    message('gene: ', gene)
    pvals = t.run4$P_BOLT_LMM[t.run4$P_BOLT_LMM < 1 & t.run4$CHR==i & t.run4$BP >= start(gtf[k]) & t.run4$BP <= end(gtf[k])]
    # weights = t.run4$A1FREQ[t.run4$P_BOLT_LMM < 1 & t.run4$CHR==i & t.run4$BP >= start(gtf[k]) & t.run4$BP <= end(gtf[k])]
    # weights = -log10(weights)
    # chr_res[k] <- CCT(pvals, weights)
    chr_res[k] <- CCT(pvals)
  }
  res.chrs[[i]] <- tibble(gene=names(gtf),pval_CCT=chr_res)
}
# uses SeqArray to find MAC < 6 variants
# assumes .gds files exist for each chromosome
require(SeqArray)
for(i in 1:22){
  fname=paste0('UKBexomeQQFE_chr',i,'.gds', separate='')
  gds <- seqOpen(fname)
  seqSetFilterCond(gds, mac=c(0,5), parallel=10) # ~2.1 min for chr1
  a <- seqGetData(gds, c(chr="chromosome", pos="position", allele="allele"))
  
  message('Chr: ', i, ': variants with MAC < 6: ', length(a[[1]]))
  as_tibble(a) %>%
    separate(allele, c('A1','A2'), sep=',') %>% 
    mutate(var=paste0(chr, ':', pos, ':', A1, ':',A2)) %>% 
    select(var) %>% 
    write_tsv('UKBexomeQQFE_chr1-22.mac_lt6.variant.ids', col_names=F, append=T)
  
  seqClose(gds)
}
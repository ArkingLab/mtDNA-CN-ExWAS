library(SeqArray)
gene='SLC25A37'
t <- TopmedPipeline::getobj('~/mtrv/burden/genesis/run5/cds_nonsyn_0.01/assoc_aggregate.SMMAT.AF.max.0.01.10.results.RData')
x <- var.breakdown(gds=gds, t=t, gene=gene)
seqSetFilter(gds, variant.sel=x[[1]]$variant.id)
x[[1]]$allele <- seqGetData(gds, 'allele')

x[[1]] %>% 
  filter(MAC >= 6) %>%
  # filter(MAC < 6) %>%
  mutate(reformat=paste0(chr,':',pos,':', stringr::str_replace(allele,',',':'))) %>% 
  arrange(MAC) %>% 
  select(reformat, MAC) %>% 
  write_tsv(paste0('~/mtrv/phesant/',tolower(gene),'.carriers/',gene,'_mtDNA_CN_n180k_MAC_gteq6_variants.withMAC.txt'))
  # write_tsv(paste0('~/mtrv/phesant/',tolower(gene),'.carriers/',gene,'_mtDNA_CN_n180k_MAC_lt6_variants.withMAC.txt'))
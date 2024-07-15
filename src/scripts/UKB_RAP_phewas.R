# PheWAS

# library(devtools)
# devtools::install_github("PheWAS/PheWAS")
# library(PheWAS)
library(tools)
library(dplyr)
library(readr)
library(parallel)
library(splines)
library(doParallel)
library(logistf)
library(bigmemory)
library(biganalytics)
library(doSNOW)
# install.packages(c('bigmemory','biganalytics','doSNOW'))

# Main functions
# see phewas_functions.R
# source('phewas_functions.R')

# gene-level variant breakdown ----
x.samhd1 <- var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr20.gds'), 
               # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr20/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
               t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
               gene='SAMHD1')
x.twnk <- var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr10.gds'), 
                         # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr10/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                         t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                         gene='TWNK')
x.tfam <- var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr10.gds'), 
                         # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr10/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                         t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                         gene='TFAM')
x.jak2 <- var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr9.gds'), 
                         # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr9/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr9.RData'),
                         t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr9.RData'),
                         gene='JAK2')
x.clpx <- var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr15.gds'), 
                         # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr15/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr15.RData'),
                         t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr15.RData'),
                         gene='CLPX')
x.ly75.cd302 <- var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr2.gds'), 
                               # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr2/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr2.RData'),
                               t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr2.RData'),
                               gene='LY75.CD302')
x.ak2 <- var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr1.gds'), 
                        # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr1/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr1.RData'),
                        t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr1.RData'),
                        gene='AK2')
x.mgme1 <- var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr20.gds'), 
                          # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr20/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                          t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                          gene='MGME1')
x.polrmt <- var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr19.gds'), 
                          # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr19/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr19.RData'),
                          t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr19.RData'),
                          gene='POLRMT')
x.chek2 <- var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr22.gds'), 
                           # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr19/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr19.RData'),
                           t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr22.RData'),
                           gene='CHEK2')
x.cavin2 <- var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr2.gds'), 
                           # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr19/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr19.RData'),
                           t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr2.RData'),
                           gene='CAVIN2')

# list of variants (bonf and fdr10pct)
bonf.vars <-
  rbind(x.samhd1$var.bkdwn %>% filter(bonf.sig==1),
      x.twnk$var.bkdwn %>% filter(bonf.sig==1),
      x.tfam$var.bkdwn %>% filter(bonf.sig==1),
      x.jak2$var.bkdwn %>% filter(bonf.sig==1),
      x.clpx$var.bkdwn %>% filter(bonf.sig==1),
      x.ly75.cd302$var.bkdwn %>% filter(bonf.sig==1),
      x.ak2$var.bkdwn %>% filter(bonf.sig==1),
      x.mgme1$var.bkdwn %>% filter(bonf.sig==1),
      x.polrmt$var.bkdwn %>% filter(bonf.sig==1),
      x.chek2$var.bkdwn %>% filter(bonf.sig==1),
      x.cavin2$var.bkdwn %>% filter(bonf.sig==1))
bonf.vars$gene <- c(rep('SAMHD1',nrow(x.samhd1[[4]] %>% filter(bonf.sig==TRUE))),
                        rep('TWNK',nrow(x.twnk[[4]] %>% filter(bonf.sig==TRUE))),
                        rep('TFAM',nrow(x.tfam[[4]] %>% filter(bonf.sig==TRUE))),
                        rep('JAK2',nrow(x.jak2[[4]] %>% filter(bonf.sig==TRUE))),
                        rep('CLPX',nrow(x.clpx[[4]] %>% filter(bonf.sig==TRUE))),
                        rep('LY75.CD302',nrow(x.ly75.cd302[[4]] %>% filter(bonf.sig==TRUE))),
                        rep('AK2',nrow(x.ak2[[4]] %>% filter(bonf.sig==TRUE))),
                        rep('MGME1',nrow(x.mgme1[[4]] %>% filter(bonf.sig==TRUE))),
                        rep('POLRMT',nrow(x.polrmt[[4]] %>% filter(bonf.sig==TRUE))),
                        rep('CHEK2',nrow(x.chek2[[4]] %>% filter(bonf.sig==TRUE))),
                        rep('CAVIN2',nrow(x.cavin2[[4]] %>% filter(bonf.sig==TRUE))))
# bonf.vars$gene <- c(rep('SAMHD1',4),rep('TWNK',4),rep('TFAM',3),'JAK2',rep('CLPX',2), rep('LY75.CD302',2), 'AK2','MGME1','CAVIN2')

fdr10pct.vars <-
  rbind(x.samhd1$var.bkdwn %>% filter(fdr_sig==1),
        x.twnk$var.bkdwn %>% filter(fdr_sig==1),
        x.tfam$var.bkdwn %>% filter(fdr_sig==1),
        x.jak2$var.bkdwn %>% filter(fdr_sig==1),
        x.clpx$var.bkdwn %>% filter(fdr_sig==1),
        x.ly75.cd302$var.bkdwn %>% filter(fdr_sig==1),
        x.ak2$var.bkdwn %>% filter(fdr_sig==1),
        x.mgme1$var.bkdwn %>% filter(fdr_sig==1),
        x.polrmt$var.bkdwn %>% filter(fdr_sig==1),
        x.chek2$var.bkdwn %>% filter(fdr_sig==1),
        x.cavin2$var.bkdwn %>% filter(fdr_sig==1))
fdr10pct.vars$gene <- c(rep('SAMHD1',nrow(x.samhd1[[4]] %>% filter(fdr_sig==TRUE))),
                        rep('TWNK',nrow(x.twnk[[4]] %>% filter(fdr_sig==TRUE))),
                        rep('TFAM',nrow(x.tfam[[4]] %>% filter(fdr_sig==TRUE))),
                        rep('JAK2',nrow(x.jak2[[4]] %>% filter(fdr_sig==TRUE))),
                        rep('CLPX',nrow(x.clpx[[4]] %>% filter(fdr_sig==TRUE))),
                        rep('LY75.CD302',nrow(x.ly75.cd302[[4]] %>% filter(fdr_sig==TRUE))),
                        rep('AK2',nrow(x.ak2[[4]] %>% filter(fdr_sig==TRUE))),
                        rep('MGME1',nrow(x.mgme1[[4]] %>% filter(fdr_sig==TRUE))),
                        rep('POLRMT',nrow(x.polrmt[[4]] %>% filter(fdr_sig==TRUE))),
                        rep('CHEK2',nrow(x.chek2[[4]] %>% filter(fdr_sig==TRUE))),
                        rep('CAVIN2',nrow(x.cavin2[[4]] %>% filter(fdr_sig==TRUE))))

# gene-level synthetic alleles from bonf-sig variants in each gene ----
v.samhd1 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr20.gds'), 
                                                # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr20/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                                t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                                gene='SAMHD1'), 
                                 level='bonf', include_lt6=T)
v.twnk <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr10.gds'), 
                                                # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr10/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                                              t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                                                gene='TWNK'), 
                                 level='bonf', include_lt6=T)
v.tfam <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr10.gds'), 
                                                # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr10/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                                              t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                                                gene='TFAM'), 
                                 level='bonf', include_lt6=T)
v.jak2 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr9.gds'), 
                                                # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr9/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr9.RData'),
                                              t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr9.RData'),
                                                gene='JAK2'),
                                 level='bonf', include_lt6=T)
v.ly75.cd302 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr2.gds'), 
                                                # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr2/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr2.RData'),
                                                t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr2.RData'),
                                                gene='LY75.CD302'), 
                                 level='bonf', include_lt6=T)
v.clpx <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr15.gds'), 
                                                # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr15/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr15.RData'),
                                              t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr15.RData'),
                                                gene='CLPX'), 
                                 level='bonf', include_lt6=T)
v.ak2 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr1.gds'), 
                                                # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr1/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr1.RData'),
                                             t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr1.RData'),
                                                gene='AK2'), 
                                 level='bonf', include_lt6=T)
v.mgme1 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr20.gds'), 
                                                # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr20/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                               t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                                gene='MGME1'), 
                                 level='bonf', include_lt6=T)

v.CHEK2 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr22.gds'), 
                                               # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr20/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                               t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr22.RData'),
                                               gene='CHEK2'), 
                                level='bonf', include_lt6=T)

v.CAVIN2 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr2.gds'), 
                                               # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr20/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                               t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr2.RData'),
                                               gene='CAVIN2'), 
                                level='bonf', include_lt6=T)

# gene-level synthetic alleles from FDR10pct-sig variants in each gene ----
v.samhd1 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr20.gds', allow.duplicate=T), 
                                                # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr20/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                                t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                                gene='SAMHD1'), 
                                 level='fdr', include_lt6=T)
v.twnk <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr10.gds'), 
                                              # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr10/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                                              t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                                              gene='TWNK'), 
                               level='fdr', include_lt6=T)
v.tfam <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr10.gds'), 
                                              # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr10/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                                              t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr10.RData'),
                                              gene='TFAM'), 
                               level='fdr', include_lt6=T)
v.jak2 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr9.gds'), 
                                              # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr9/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr9.RData'),
                                              t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr9.RData'),
                                              gene='JAK2'),
                               level='fdr', include_lt6=T)
v.clpx <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr15.gds'), 
                                              # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr15/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr15.RData'),
                                              t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr15.RData'),
                                              gene='CLPX'), 
                               level='fdr', include_lt6=T)
v.ly75.cd302 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr2.gds'), 
                                                    # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr2/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr2.RData'),
                                                    t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr2.RData'),
                                                    gene='LY75.CD302'), 
                                     level='fdr', include_lt6=T)
v.ak2 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr1.gds'), 
                                             # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr1/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr1.RData'),
                                             t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr1.RData'),
                                             gene='AK2'), 
                              level='fdr', include_lt6=T)
v.mgme1 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr20.gds'), 
                                               # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr20/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                               t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                               gene='MGME1'), 
                                level='fdr', include_lt6=T)
v.chek2 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr22.gds'), 
                                               # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr20/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                               t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr22.RData'),
                                               gene='CHEK2'), 
                                level='fdr', include_lt6=T)
v.cavin2 <- get.synthetic_allele(var.breakdown2(gds=seqOpen('~/mtrv/resources/gds/UKBexomeQQFE_chr2.gds'), 
                                               # t= getobj('~/mtrv/resources/res/nonsyn.impact_mod_high.AF.max.0.01_chr20/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr20.RData'),
                                               t= getobj('~/mtrv/resources/res/rv_assoc.aggr.nonsyn.impact_mod_high.SMMAT.AF.max.0.01.chr2.RData'),
                                               gene='CAVIN2'), 
                                level='fdr', include_lt6=T)
v.samhd1 <- get.synthetic_allele(x.samhd1, level='fdr', include_lt6=T)
v.twnk <- get.synthetic_allele(x.twnk, level='fdr', include_lt6=T)
v.tfam <- get.synthetic_allele(x.tfam, level='fdr', include_lt6=T)
v.jak2 <- get.synthetic_allele(x.jak2, level='fdr', include_lt6=T)
v.clpx <- get.synthetic_allele(x.clpx, level='fdr', include_lt6=T)
v.ly75.cd302 <- get.synthetic_allele(x.ly75.cd302, level='fdr', include_lt6=T)
v.ak2 <- get.synthetic_allele(x.ak2, level='fdr', include_lt6=T)
v.mgme1 <- get.synthetic_allele(x.mgme1, level='fdr', include_lt6=T)
# v.polrmt <- get.synthetic_allele(x.polrmt, level='fdr', include_lt6=T)
v.chek2 <- get.synthetic_allele(x.chek2, level='fdr', include_lt6=T)
v.cavin2 <- get.synthetic_allele(x.cavin2, level='fdr', include_lt6=T)


# combine synthetic alleles ----
v <- v.samhd1$syn.allele %>%
  left_join(v.twnk$syn.allele, by='IID') %>% 
  left_join(v.tfam$syn.allele, by='IID') %>% 
  left_join(v.jak2$syn.allele, by='IID') %>%
  left_join(v.clpx$syn.allele, by='IID') %>%
  left_join(v.ly75.cd302$syn.allele, by='IID') %>%
  left_join(v.ak2$syn.allele, by='IID') %>%
  left_join(v.mgme1$syn.allele, by='IID') %>% 
  # left_join(v.chek2$syn.allele, by='IID') %>% 
  left_join(v.cavin2$syn.allele, by='IID')
  # colnames(v) <- c('IID','SAMHD1','TWNK','TFAM','JAK2','CLPX','LY75.CD302','AK2','MGME1','CHEK2','CAVIN2')
colnames(v) <- c('IID','SAMHD1','TWNK','TFAM','JAK2','CLPX','LY75.CD302','AK2','MGME1','CAVIN2')

v$syn.allele <- apply(v, 1, function(x) ifelse(any(x==1, na.rm=T), 1, 0))
var.carriers <- v %>% filter(syn.allele==1)
n.carriers <- apply(var.carriers[,2:(ncol(var.carriers)-1)], 1, sum)
table(n.carriers) # 446 IIDs with 2 hits, 3 with 3 hits (bonf)
                  # 684 IIDs with 2 hits, 8 with 3 hits (fdr 10%)

# positive betas
v.tmp <- v %>% select(IID, SAMHD1, JAK2, CLPX,LY75.CD302,AK2,MGME1)
v.tmp$syn.allele <- apply(v.tmp, 1, function(x) ifelse(any(x==1, na.rm=T), 1, 0))
var.carriers <- v.tmp %>% filter(syn.allele==1)
data <- create_phewas_data.input(var.carriers, subset.to.unrelated=TRUE, subset.to.white.British.ancestry=FALSE)
phecodes = data[[3]]
sex.check=data[[2]]
length(phecodes)
data <- data[[1]]
k2 <- bigmemory::as.big.matrix(as.data.frame(data))
data_desc <- bigmemory::describe(k2)
phewas <- process.glm(phecodes=sample(phecodes), sex_check=sex.check, data_desc=data_desc, ncores=10, var='var', file.conn='./phewas.corr_cutoff_genes.mac_gteq6_bonf_include_lt6.positive_betas.inFullExomesSetk.txt')


# negative betas
v.tmp <- v %>% select(IID, TWNK, TFAM)
v.tmp$syn.allele <- apply(v.tmp, 1, function(x) ifelse(any(x==1, na.rm=T), 1, 0))
var.carriers <- v.tmp %>% filter(syn.allele==1)
data <- create_phewas_data.input(var.carriers, subset.to.unrelated=TRUE, subset.to.white.British.ancestry=FALSE)
phecodes = data[[3]]
sex.check=data[[2]]
length(phecodes)
data <- data[[1]]
k2 <- bigmemory::as.big.matrix(as.data.frame(data))
data_desc <- bigmemory::describe(k2)
phewas <- process.glm(phecodes=sample(phecodes), sex_check=sex.check, data_desc=data_desc, ncores=10, var='var', file.conn='./phewas.corr_cutoff_genes.mac_gteq6_bonf_include_lt6.negative_betas.inFullExomesSetk.txt')

# no JAK2 (recompute syn.allele)
v1 <- v %>% select(-JAK2, -syn.allele)
v1$syn.allele <- apply(v1, 1, function(x) ifelse(any(x==1, na.rm=T), 1, 0))
var.carriers <- v1 %>% filter(syn.allele==1)
n.carriers <- apply(var.carriers[,2:(ncol(var.carriers)-1)], 1, sum)
table(n.carriers)

# combine mac>=6 bonf for TWNK, TFAM, MGME1 (depletion syndrome genes)
v_depletion <- v.twnk %>% 
  left_join(v.tfam, by='IID') %>% 
  left_join(v.mgme1, by='IID')
colnames(v_depletion) <- c('IID','TWNK','TFAM','MGME1')
v_depletion$syn.allele <- apply(v_depletion, 1, function(x) ifelse(any(x==1, na.rm=T), 1, 0))
var.carriers <- v_depletion %>% filter(syn.allele==1)
n.carriers <- apply(var.carriers[,2:(ncol(var.carriers)-1)], 1, sum)
table(n.carriers) # 14 with 2 hits

#### rv score from synthetic alleles ----####
v <- v.samhd1$scored.alleles %>%
  left_join(v.twnk$scored.alleles, by='sample.id') %>% 
  left_join(v.tfam$scored.alleles, by='sample.id') %>% 
  left_join(v.jak2$scored.alleles, by='sample.id') %>%
  left_join(v.ly75.cd302$scored.alleles, by='sample.id') # %>%
  # left_join(v.clpx$scored.alleles, by='sample.id') %>%
  # left_join(v.ak2$scored.alleles, by='sample.id') %>%
  # left_join(v.mgme1$scored.alleles, by='sample.id') %>% 
  # left_join(v.CAVIN2$scored.alleles, by='sample.id')

# downselect to those alleles with effect size > +/- 0.3
es_0.3_vars <- c('var.330053.pos.20_36893060_mac_304',
                 'var.330419.pos.20_36905431_mac_104',
                 'var.330977.pos.20_36930783_mac_197',
                 'var.768340.pos.10_100990873_mac_11',
                 'var.768362.pos.10_100990943_mac_6',
                 'var.768474.pos.10_100993257_mac_48',
                 'var.361598.pos.10_58388676_mac_320',
                 'var.361871.pos.10_58394946_mac_81',
                 'var.44367.pos.9_5073770_mac_284', # Jak2
                 'var.1141606.pos.2_159854509_mac_8',
                 'mac_lt_6.y') # should this be included?
# create weighted score
v <- v %>% select(sample.id, all_of(es_0.3_vars))
v$scored.allele <- apply(v[,2:ncol(v)], 1, sum, na.rm=T) # non absolute-value
# v$scored.allele <- abs(apply(v[,2:ncol(v)], 1, sum, na.rm=T)) # absolute-value
var.carriers <- v %>% select(IID=sample.id, scored.allele) %>% filter(scored.allele!=0)


# Only variants with effect sizes > 0.1339 (largest abs ES from common var GWAS) ----
# fdr10pct_ES <- fdr10pct.vars %>% filter(abs(beta) > 0.1339) %>% arrange(beta)
fdr10pct_ES <- fdr10pct.vars %>% filter(abs(beta) >= 0.30) %>% arrange(beta)
vars <- list()
IID <- v.samhd1$alleles$sample.id
vars[[1]] <- v.samhd1$alleles %>% select(all_of(fdr10pct_ES$var[which(fdr10pct_ES$gene=='SAMHD1')]))
vars[[2]] <- v.twnk$alleles %>% select(all_of(fdr10pct_ES$var[which(fdr10pct_ES$gene=='TWNK')]))
vars[[2]] <- rename(vars[[2]], mac_lt_6_TWNK=mac_lt_6)
vars[[3]] <- v.tfam$alleles %>% select(all_of(fdr10pct_ES$var[which(fdr10pct_ES$gene=='TFAM')]))
vars[[3]] <- rename(vars[[3]], mac_lt_6_TFAM=mac_lt_6)
# vars[[4]] <- v.jak2$alleles %>% select(all_of(fdr10pct_ES$var[which(fdr10pct_ES$gene=='JAK2')]))
vars[[5]] <- v.clpx$alleles %>% select(all_of(fdr10pct_ES$var[which(fdr10pct_ES$gene=='CLPX')]))
vars[[6]] <- v.ly75.cd302$alleles %>% select(all_of(fdr10pct_ES$var[which(fdr10pct_ES$gene=='LY75.CD302')]))
vars[[7]] <- v.ak2$alleles %>% select(all_of(fdr10pct_ES$var[which(fdr10pct_ES$gene=='AK2')]))
vars[[8]] <- v.mgme1$alleles %>% select(all_of(fdr10pct_ES$var[which(fdr10pct_ES$gene=='MGME1')]))
vars[[9]] <- v.chek2$alleles %>% select(all_of(fdr10pct_ES$var[which(fdr10pct_ES$gene=='CHEK2')]))
vars[[10]] <- v.cavin2$alleles %>% select(all_of(fdr10pct_ES$var[which(fdr10pct_ES$gene=='CAVIN2')]))
vars <- bind_cols(vars)
vars$syn.allele <- apply(vars, 1, function(x) ifelse(any(x!=0, na.rm=T), 1, 0))
vars$IID <- IID
vars <- vars %>% select(IID, everything())
var.carriers <- vars %>% filter(syn.allele==1)
n.carriers <- apply(var.carriers[,2:(ncol(var.carriers)-1)], 1, sum)
table(n.carriers)



# create data set ----
# data <- create_phewas_data.input(var.carriers, subset.to.unrelated=TRUE, subset.to.white.British.ancestry=FALSE)
data <- create_phewas_data.input2(var.carriers)

# data <- create_phewas_data.input2(var.carriers)
phecodes = data[[3]]
sex.check=data[[2]]
length(phecodes)

## Create a bigmatrix object from `data`
require(bigmemory)
require(biganalytics)
data <- data[[1]]
k2 <- bigmemory::as.big.matrix(as.data.frame(data))
data_desc <- bigmemory::describe(k2)

# Run PheWAS (examples)
# phewas <- process.glm(phecodes=sample(phecodes), sex_check=sex.check, 
phewas <- process.glm3.1(phecodes=sample(phecodes), both_sexes=sex.check, 
                       data_=data, ncores=1,var='var',
                       file.conn='./phewas.rv_scored_allele.es_0.3.noJAK2.inFullExomesSet.txt')

# RV SCORE (using brute-force GLM logistic regressions)
phewas <- process.glm3.1_rvscore(phecodes=sample(phecodes), data=data, var='var', both_sexes=sex.check, ncores=1, file.conn='~/mtrv/phewas_results.rvscore_es0.3.glm.txt')
phewas <- process.glm3.1_rvscore(phecodes=sample(phecodes), data=data, var='var', both_sexes=sex.check, ncores=1, file.conn='~/mtrv/phewas_results.rvscore_es0.3.glm2_wTFAM_mac_lt_6.txt')
phewas <- process.glm3.1_rvscore(phecodes=sample(phecodes), data=data, var='var', both_sexes=sex.check, ncores=1, file.conn='~/mtrv/phewas_results.rvscore_es0.3.glm2_wTFAM_mac_lt_6_wJAK2.txt')

phewas <- # parallel rv score GLM
  process.glm2_rvscore(
    phecodes = sample(phecodes),
    both_sexes = sex.check,
    data_desc = data_desc,
    ncores = 15,
    var = 'var',
    back_correct_FirthSE = TRUE,
    file.conn = './phewas_results.rvscore_es0.3.glm2_wTFAM_mac_lt_6_wJAK2.txt'
  )

phewas <- process.glm2(phecodes=sample(phecodes), sex_check=sex.check, 
                      data=data_desc, ncores=7,var='var',
                      file.conn='./phewas.rv_scored_allele.es_0.3.noJAK2.inFullExomesSet.txt')
                      # file.conn='./phewas.rv_scored_allele.noJAK2.inFullExomesSet.txt')
                      # file.conn='./phewas.corr_cutoff_genes.fdr10pct_include_lt6.ESgt0.3.noJAK2.inFullExomesSet.txt')
                      # file.conn='./phewas.depletion_genes.mac_gteq6_bonf_include_lt6_TWNK_TFAM_MGME1.inFullExomesSetk.txt2')
                      # file.conn='./phewas.corr_cutoff_genes.mac_gteq6_fdr10pct_include_lt6.inFullExomesSetk2.txt')
                      # file.conn='./phewas.corr_cutoff_genes.mac_gteq6_bonf_include_lt6.inFullExomesSetk2.txt')
                      # file.conn='./phewas.corr_cutoff_genes.mac_gteq6_fdr10pct_include_lt6.noJAK2.inFullExomesSetk.txt')
                      # file.conn='./phewas.corr_cutoff_genes.mac_gteq6_fdr10pct_include_lt6.positive_betas.inFullExomesSetk.txt')
                      # file.conn='./phewas.corr_cutoff_genes.mac_gteq6_fdr10pct_include_lt6.negative_betas.inFullExomesSetk.txt')
                      # file.conn='./phewas.corr_cutoff_genes.mac_gteq6_fd`r10pct_include_lt6.noJAK2_SAMHD1_TWNK_TFAM.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_gteq6_bonf_include_lt6_TFAM.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_gteq6_bonf_MGME1.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_gteq6_bonf_AK2.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_gteq6_bonf_include_lt6_TWNK.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_gteq6_bonf_LY75_CD302.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_gteq6_bonf_JAK2.inFullExomesSetk.txt3')
                      # file.conn='./phewas.all_genes.mac_gteq6_bonf_CLPX.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_gteq6_bonf_TFAM.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_gteq6_bonf_SAMHD1.inFullExomesSetk.txt')
                      # file.conn='./phewas3.all_genes.mac_all_gteq6_FDR10pct_SAMHD1_SLC25A37.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_all_gteq6_FDR10pct_TWNK_TFAM.inFullExomesSetk.txt')
                      # file.conn='./phewas.SLC25A37.mac_lt6.inFullExomesSetk.txt')
                      # file.conn='./phewas.TWNK.mac_lt6.inFullExomesSetk.txt')
                      # file.conn='./phewas.SAMHD1.mac_lt6.inFullExomesSetk.txt')
                      # file.conn='./phewas.TFAM.mac_gteq6.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_gteq6_FDR10pctVars.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_gteq6.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_lt6.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_all_gteq6.inFullExomesSetk.txt')
                      # file.conn='./phewas.all_genes.mac_all_gteq6_FDR10pct.inFullExomesSetk.txt')
phewas2 <- process.glm(phecodes=sample(phecodes[which(!(phecodes %in% phewas$phecode))]), sex_check=sex.check, 
                      data=data_desc, ncores=10,var='var', 
                      file.conn='./phewas.corr_cutoff_genes.fdr10pct_include_lt6.ESgt0.3.noJAK2.inFullExomesSet.txt')

phewas <- process.glm2(phecodes=sample(phecodes), sex_check=sex.check, 
                      data=data_desc, ncores=12,var='var',
                      file.conn='./phewas.corr_cutoff_genes.mac_gteq6_bonf_include_lt6.scored_allele_absolute_value.inFullExomesSetk.txt')
closeAllConnections()



# gene by gene
v.tmp <- v %>% select(IID, JAK2)
var.carriers <- v.tmp %>% filter(JAK2==1)
data <- create_phewas_data.input(var.carriers, subset.to.unrelated=TRUE, subset.to.white.British.ancestry=FALSE)
phecodes = data[[3]]
sex.check=data[[2]]
length(phecodes)
data <- data[[1]]
k2 <- bigmemory::as.big.matrix(as.data.frame(data))
data_desc <- bigmemory::describe(k2)
phewas <- process.glm(phecodes=sample(phecodes), sex_check=sex.check, data_desc=data_desc, ncores=10, var='var', file.conn='./phewas.corr_cutoff_genes.mac_gteq6_fdr10pct_include_lt6.JAK2.inFullExomesSetk.txt')




## Test phegroups (rough code) ----
i <- getobj('~/ref/UKB -mtDNA_2022/data_tables//n502485.icd10_phecodes.rds') # phecodes
phecode_def <- read_csv('~/mtrv/phecode_definitions1.2.csv')
phegroups <- as_tibble(getobj('~/ref/UKB -mtDNA_2022/vkp/n502485.phecode_cases_by_phecode_category.rds'), rownames='id')
phegroups$id <- as.numeric(phegroups$id)

var.carriers = var.carriers_gteq6_FDR10pct_SAMHD1_SLC25A37
data <- create_phewas_data.input(var.carriers)
# data <- data[[1]] %>% select(-any_of(colnames(i)[-1]))
data <- data %>% inner_join(phegroups, by='id')
var='var'
# phecodes <- colnames(phegroups)[2:ncol(phegroups)]
phecodes <- phecode_def$phecode[which(phecode_def$category_number==0)]# & phecode_def$leaf==1)]
data <- data[[1]] %>% select(any_of(c('var','age','sex','Center',paste0('PC',1:40,separate=''), phecodes)))
k2 <- bigmemory::as.big.matrix(as.data.frame(data))
data_desc <- bigmemory::describe(k2)
# for(i in 1:length(phecodes)){
#   if(!(phecodes[i] %in% colnames(data))) { next }
#   f=paste0('`', phecodes[i],'`',' ~ ', var, ' + ns(age, df=2) + sex + Center + ', paste0('PC',1:40,separate='',collapse='+'));
#   tmp <- glm(formula(f), data=data, family='binomial'); 
#   message(phecodes[i], 
#           ' beta= ', summary(tmp)$coef[2,1], 
#           ' SE= ', summary(tmp)$coef[2,2], 
#           ' p= ', summary(tmp)$coef[2,4])
# }
phewas <-
  process.glm(
    phecodes = sample(phecodes),
    data_desc = data_desc,
    var = 'var',
    sex_check = NULL,
    file.conn = './tmp4.txt',
    ncores = 8
  )

# post-processing code ----
## meta analysis by-gene
require(dplyr)
require(readr)
require(meta)
require(poolr)
source('~/mtrv/cct.R')

options(width=200)

gene='SAMHD1'
setwd(paste0('~/mtrv/',gene,'/'))
phecode_def <- 
  read_csv('~/mito_rare-variant/resources/phecode_definitions1.2.csv')

# dx cd vkp/${gene}/
# dx download phewas\*
s1 <- read_tsv(paste0('phewas.',gene,'.mac_lt6.inFullExomesSetk.txt'))
s2 <- read_tsv(paste0('phewas.',gene,'.mac_gteq6_FDR10pct.inFullExomesSetk.txt'))
# s2 <- read_tsv(paste0('phewas.',gene,'.mac_gteq6.inFullExomesSetk.txt'))

super.set.to_test <- intersect(s1$phecode, s2$phecode)
s_1 <- s1 %>% filter(phecode %in% super.set.to_test) %>% arrange(phecode)
s_2 <- s2 %>% filter(phecode %in% super.set.to_test) %>% arrange(phecode)

# combine each mac gteq6 and lt6 results into single gene dfs
s <- s_1 %>% inner_join(s_2, by='phecode') %>% arrange(phecode)

x <- s
meta.res.df <- data.frame(phecode=rep('', nrow(x)),
                          beta=rep(NA, nrow(x)),
                          SE=rep(NA, nrow(x)),
                          # p.mac_lt6=rep(NA, nrow(x)),
                          # p.mac_gteq6=rep(NA, nrow(x)),
                          pval.fixed=rep(NA, nrow(x)),
                          pval.rand=rep(NA, nrow(x)),
                          p.fisher=rep(NA, nrow(x)),
                          p.cct=rep(NA,nrow(x))) %>% as_tibble()
for(i in 1:nrow(x)) {
  mres <- metagen(TE=c(x$beta.x[i], x$beta.y[i]),
                  seTE=c(x$SE.x[i], x$SE.y[i]))
  meta.res.df$phecode[i] <- x$phecode[i]
  meta.res.df$beta[i] <- mres$TE.random
  meta.res.df$SE[i] <- mres$seTE.random
  meta.res.df$pval.fixed[i] <- mres$pval.fixed
  meta.res.df$pval.rand[i] <- mres$pval.random
  meta.res.df$p.fisher[i] <- poolr::fisher(p = c(x$p.x[i], x$p.y[i]))$p
  meta.res.df$p.cct[i] <- CCT(pvals = c(x$p.x[i], x$p.y[i]))
}
meta.res.df %>% 
  arrange(p.cct) %>% ## arrange(pval.fixed)
  inner_join(phecode_def, by='phecode')

## meta analysis across-genes
require(dplyr)
require(readr)
require(meta)
require(poolr)
source('~/mtrv/cct.R')
options(width=200)
setwd(paste0('~/mtrv//'))
phecode_def <- 
  read_csv('~/mito_rare-variant/resources/phecode_definitions1.2.csv')

# gene='TFAM'
# dx cd vkp/${gene}/
# dx download phewas\*
genes=c('SAMHD1','SLC25A37','TWNK','TFAM')
# gt6 fdr10pct
s1 <- read_tsv(paste0('~/mtrv/',genes[1],'/phewas.',genes[1],'.mac_gteq6_FDR10pct.inFullExomesSetk.txt'))
s2 <- read_tsv(paste0('~/mtrv/',genes[2],'/phewas.',genes[2],'.mac_gteq6_FDR10pct.inFullExomesSetk.txt'))
s3 <- read_tsv(paste0('~/mtrv/',genes[3],'/phewas.',genes[3],'.mac_gteq6_FDR10pct.inFullExomesSetk.txt'))
s4 <- read_tsv(paste0('~/mtrv/',genes[4],'/phewas.',genes[4],'.mac_gteq6_FDR10pct.inFullExomesSetk.txt'))
# gt6 (any)
s1 <- read_tsv(paste0('~/mtrv/',genes[1],'/phewas.',genes[1],'.mac_gteq6.inFullExomesSetk.txt'))
s2 <- read_tsv(paste0('~/mtrv/',genes[2],'/phewas.',genes[2],'.mac_gteq6.inFullExomesSetk.txt'))
s3 <- read_tsv(paste0('~/mtrv/',genes[3],'/phewas.',genes[3],'.mac_gteq6.inFullExomesSetk.txt'))
s4 <- read_tsv(paste0('~/mtrv/',genes[4],'/phewas.',genes[4],'.mac_gteq6.inFullExomesSetk.txt'))
# lt6
s1 <- read_tsv(paste0('~/mtrv/',genes[1],'/phewas.',genes[1],'.mac_lt6.inFullExomesSetk.txt'))
s2 <- read_tsv(paste0('~/mtrv/',genes[2],'/phewas.',genes[2],'.mac_lt6.inFullExomesSetk.txt'))
s3 <- read_tsv(paste0('~/mtrv/',genes[3],'/phewas.',genes[3],'.mac_lt6.inFullExomesSetk.txt'))
s4 <- read_tsv(paste0('~/mtrv/',genes[4],'/phewas.',genes[4],'.mac_lt6.inFullExomesSetk.txt'))
# mac_all (w/ gteq6 FDR10pct)
s1 <- read_table(paste0('~/mtrv/',genes[1],'/phewas.',genes[1],'.mac_all_wFDR10pctGteq6_.inFullExomesSetk.txt'), col_types=cols('phecode'=col_character()))
s2 <- read_table(paste0('~/mtrv/',genes[2],'/phewas.',genes[2],'.mac_all_gteq6_FDR10pct.inFullExomesSetk.txt'), col_types=cols('phecode'=col_character()))
s3 <- read_table(paste0('~/mtrv/',genes[3],'/phewas.',genes[3],'.mac_all_gteq6_FDR10pct.inFullExomesSetk.txt'), col_types=cols('phecode'=col_character()))
s4 <- read_table(paste0('~/mtrv/',genes[4],'/phewas.',genes[4],'.mac_all_gteq6_FDR10pct.inFullExomesSetk.txt'), col_types=cols('phecode'=col_character()))
# mac_all (w/ any gteq6)
s1 <- read_table(paste0('~/mtrv/',genes[1],'/phewas.',genes[1],'.mac_all.inFullExomesSetk.txt'), col_types=cols('phecode'=col_character()))
s2 <- read_table(paste0('~/mtrv/',genes[2],'/phewas.',genes[2],'.mac_all_gteq6.inFullExomesSetk.txt'), col_types=cols('phecode'=col_character()))
s3 <- read_table(paste0('~/mtrv/',genes[3],'/phewas.',genes[3],'.mac_all_gteq6.inFullExomesSetk.txt'), col_types=cols('phecode'=col_character()))
s4 <- read_table(paste0('~/mtrv/',genes[4],'/phewas.',genes[4],'.mac_all.inFullExomesSetk.txt'), col_types=cols('phecode'=col_character()))

super.set.to_test_all4genes <- Reduce(intersect,list(s1$phecode, s2$phecode, s3$phecode, s4$phecode))
# super.set.to_test_SAMHD1_TWNK_TFAM <- Reduce(intersect,list(s1$phecode, s3$phecode, s4$phecode))
# super.set.to_test_SAMHD1_SLC25A37 <- Reduce(intersect,list(s1$phecode, s2$phecode))
# super.set.to_test_TWNK_TFAM <- Reduce(intersect,list(s3$phecode, s4$phecode))
# super.set.to_test_SAMHD1_TWNK <- Reduce(intersect,list(s1$phecode, s3$phecode))
# super.set.to_test_SAMHD1_TFAM <- Reduce(intersect,list(s1$phecode, s4$phecode))

super.set.to_test = super.set.to_test_all4genes
# super.set.to_test = super.set.to_test_SAMHD1_SLC25A37
# super.set.to_test = super.set.to_test_TWNK_TFAM
# super.set.to_test = super.set.to_test_SAMHD1_TWNK
# super.set.to_test = super.set.to_test_SAMHD1_TFAM
# super.set.to_test = super.set.to_test_SAMHD1_TWNK_TFAM

s_1 <- s1 %>% filter(phecode %in% super.set.to_test) %>% arrange(phecode)
s_2 <- s2 %>% filter(phecode %in% super.set.to_test) %>% arrange(phecode)
s_3 <- s3 %>% filter(phecode %in% super.set.to_test) %>% arrange(phecode)
s_4 <- s4 %>% filter(phecode %in% super.set.to_test) %>% arrange(phecode)

# combine each mac gteq6 and lt6 results
s_ <- s_1 %>%
  inner_join(s_2, by='phecode') %>%
  inner_join(s_3, by='phecode') %>%
  inner_join(s_4, by='phecode') %>%
  arrange(phecode)
# s_ <- s_1 %>%
#   inner_join(s_2, by='phecode') %>%
#   arrange(phecode)
# s_ <- s_3 %>%
#   inner_join(s_4, by='phecode') %>%
#   arrange(phecode)
# s_ <- s_1 %>%
#   inner_join(s_3, by='phecode') %>%
#   arrange(phecode)
# s_ <- s_1 %>%
#   inner_join(s_4, by='phecode') %>%
#   arrange(phecode)
# s_ <- s_1 %>%
#   inner_join(s_3, by='phecode') %>%
#   inner_join(s_4, by='phecode') %>%
#   arrange(phecode)

x <- s_
meta.res.df <- data.frame(phecode=rep('', nrow(x)),
                          beta=rep(NA, nrow(x)),
                          SE=rep(NA, nrow(x)),
                          pval.fixed=rep(NA, nrow(x)),
                          pval.rand=rep(NA, nrow(x)),
                          p.fisher=rep(NA, nrow(x)),
                          p.cct=rep(NA,nrow(x))) %>% as_tibble()
for(i in 1:nrow(x)) {
  # inv variance meta analysis
  mres <- metagen(TE=c(x$beta.x[i], x$beta.y[i], x$beta.x.x[i], x$beta.y.y[i]),
                  seTE=c(x$SE.x[i], x$SE.y[i], x$SE.x.x[i], x$SE.y.y[i]))
  # mres <- metagen(TE=c(x$beta.x[i], x$beta.y[i]),
  #                 seTE=c(x$SE.x[i], x$SE.y[i]))
  # mres <- metagen(TE=c(x$beta.x[i], x$beta.y[i], x$beta[i]),
  #                 seTE=c(x$SE.x[i], x$SE.y[i], x$SE[i]))
  
  meta.res.df$phecode[i] <- x$phecode[i]
  meta.res.df$beta[i] <- mres$TE.random
  meta.res.df$SE[i] <- mres$seTE.random
  meta.res.df$pval.fixed[i] <- mres$pval.fixed
  meta.res.df$pval.rand[i] <- mres$pval.random

  # meta-analysis using p-value combinations
  meta.res.df$p.fisher[i] <- poolr::fisher(p = c(x$p.x[i], x$p.y[i], x$p.x.x[i], x$p.y.y[i]))$p
  meta.res.df$p.cct[i] <- CCT(pvals = c(x$p.x[i], x$p.y[i], x$p.x.x[i], x$p.y.y[i]))
  # meta.res.df$p.fisher[i] <- poolr::fisher(p = c(x$p.x[i], x$p.y[i]))$p
  # meta.res.df$p.cct[i] <- CCT(pvals = c(x$p.x[i], x$p.y[i]))
  # meta.res.df$p.fisher[i] <- poolr::fisher(p = c(x$p.x[i], x$p.y[i], x$p[i]))$p
  # meta.res.df$p.cct[i] <- CCT(pvals = c(x$p.x[i], x$p.y[i], x$p[i]))
}
meta.res.df %>% arrange(pval.fixed) %>% inner_join(phecode_def, by='phecode')
meta.res.df %>% arrange(pval.rand) %>% inner_join(phecode_def, by='phecode')
meta.res.df %>% arrange(p.fisher) %>% inner_join(phecode_def, by='phecode')
meta.res.df %>% arrange(p.cct) %>% inner_join(phecode_def, by='phecode')






















# Load variant carrier lists (OLD) ----
# Individual gene carrier lists
# lt carriers
var.carriers_lt6_SAMHD1=tibble(IID=readLines('~/mtrv/SAMHD1/SAMHD1_mac_lt6_IIDs.in_450kExomes.txt'))
var.carriers_lt6_TWNK=tibble(IID=readLines('~/mtrv/TWNK/TWNK_mac_lt6_IIDs.in_450kExomes.txt'))
var.carriers_lt6_TFAM=tibble(IID=readLines('~/mtrv/TFAM/TFAM_mac_lt6_IIDs.in_450kExomes.txt'))
var.carriers_lt6_SLC25A37=tibble(IID=readLines('~/mtrv/SLC25A37/SLC25A37_mac_lt6_IIDs.in_450kExomes.txt'))

# gteq carriers
var.carriers_gteq6_SAMHD1=tibble(IID=readLines('~/mtrv/SAMHD1/SAMHD1_mac_gteq6_IIDs.in_450kExomes.txt'))
var.carriers_gteq6_TWNK=tibble(IID=readLines('~/mtrv/TWNK/TWNK_mac_gteq6_IIDs.in_450kExomes.txt'))
var.carriers_gteq6_TFAM=tibble(IID=readLines('~/mtrv/TFAM/TFAM_mac_gteq6_IIDs.in_450kExomes.txt'))
var.carriers_gteq6_SLC25A37=tibble(IID=readLines('~/mtrv/SLC25A37/SLC25A37_mac_gteq6_IIDs.in_450kExomes.txt'))

# gteq6 (FDR10pct only) carriers
var.carriers_gteq6_FDR10pct_SAMHD1=tibble(IID=readLines('~/mtrv/SAMHD1/SAMHD1_mac_gteq6_FDR10pct_IIDs.in_450kExomes.txt'))
var.carriers_gteq6_FDR10pct_TWNK=tibble(IID=readLines('~/mtrv/TWNK/TWNK_mac_gteq6_FDR10pct_IIDs.in_450kExomes.txt'))
var.carriers_gteq6_FDR10pct_TFAM=tibble(IID=readLines('~/mtrv/TFAM/TFAM_mac_gteq6_FDR10pct_IIDs.in_450kExomes.txt'))
var.carriers_gteq6_FDR10pct_SLC25A37=tibble(IID=readLines('~/mtrv/SLC25A37/SLC25A37_mac_gteq6_FDR10pct_IIDs.in_450kExomes.txt'))

# mac_all (with FDR10pct gteq6) carriers
var.carriers_mac_all_gteq6_FDR10pct_SAMHD1=tibble(IID=readLines('~/mtrv/SAMHD1/SAMHD1_mac_all_FDR10pct_IIDs.in_450kExomes.txt'))
var.carriers_mac_all_gteq6_FDR10pct_TWNK=tibble(IID=readLines('~/mtrv/TWNK/TWNK_mac_all_FDR10pct_IIDs.in_450kExomes.txt'))
var.carriers_mac_all_gteq6_FDR10pct_TFAM=tibble(IID=readLines('~/mtrv/TFAM/TFAM_mac_all_FDR10pct_IIDs.in_450kExomes.txt'))
var.carriers_mac_all_gteq6_FDR10pct_SLC25A37=tibble(IID=readLines('~/mtrv/SLC25A37/SLC25A37_mac_all_FDR10pct_IIDs.in_450kExomes.txt'))

# combination alternates
var.carriers_lt6 <-
  rbind(var.carriers_lt6_SAMHD1,
        var.carriers_lt6_TWNK,
        var.carriers_lt6_TFAM,
        var.carriers_lt6_SLC25A37)
var.carriers_gteq6 <-
  rbind(var.carriers_gteq6_SAMHD1,
        var.carriers_gteq6_TWNK,
        var.carriers_gteq6_TFAM,
        var.carriers_gteq6_SLC25A37)
var.carriers_gteq6_FDR10pct <-
  rbind(var.carriers_gteq6_FDR10pct_SAMHD1,
        var.carriers_gteq6_FDR10pct_TWNK,
        var.carriers_gteq6_FDR10pct_TFAM,
        var.carriers_gteq6_FDR10pct_SLC25A37)
var.carriers_gteq6_FDR10pct_SAMHD1_SLC25A37 <-
  rbind(var.carriers_gteq6_FDR10pct_SAMHD1,
        var.carriers_gteq6_FDR10pct_SLC25A37)
var.carriers_gteq6_FDR10pct_SAMHD1_TWNK <-
  rbind(var.carriers_gteq6_FDR10pct_SAMHD1,
        var.carriers_gteq6_FDR10pct_TWNK)
var.carriers_gteq6_FDR10pct_SAMHD1_TFAM <-
  rbind(var.carriers_gteq6_FDR10pct_SAMHD1,
        var.carriers_gteq6_FDR10pct_TFAM)
var.carriers_gteq6_FDR10pct_TWNK_TFAM <-
  rbind(var.carriers_gteq6_FDR10pct_TWNK,
        var.carriers_gteq6_FDR10pct_TFAM)
var.carriers_mac.all_any_gteq6 <-
  rbind(rbind(var.carriers_lt6_SAMHD1, var.carriers_gteq6_SAMHD1),
        rbind(var.carriers_lt6_TWNK, var.carriers_gteq6_TWNK),
        rbind(var.carriers_lt6_TFAM, var.carriers_gteq6_TFAM),
        rbind(var.carriers_lt6_SLC25A37, var.carriers_gteq6_SLC25A37))
var.carriers_mac.all_any_gteq6_FDR10pct <-
  rbind(rbind(var.carriers_lt6_SAMHD1, var.carriers_gteq6_FDR10pct_SAMHD1),
        rbind(var.carriers_lt6_TWNK, var.carriers_gteq6_FDR10pct_TWNK),
        rbind(var.carriers_lt6_TFAM, var.carriers_gteq6_FDR10pct_TFAM),
        rbind(var.carriers_lt6_SLC25A37, var.carriers_gteq6_FDR10pct_SLC25A37))
var.carriers_mac.all_any_gteq6_FDR10pct_SAMHD1_SLC25A37 <-
  rbind(rbind(var.carriers_lt6_SAMHD1, var.carriers_gteq6_FDR10pct_SAMHD1),
        rbind(var.carriers_lt6_SLC25A37, var.carriers_gteq6_FDR10pct_SLC25A37))
var.carriers_mac.all_any_gteq6_FDR10pct_TWNK_TFAM <-
  rbind(rbind(var.carriers_lt6_TWNK, var.carriers_gteq6_FDR10pct_TWNK),
        rbind(var.carriers_lt6_TFAM, var.carriers_gteq6_FDR10pct_TFAM))

# var.carriers = var.carriers_mac.all_any_gteq6_FDR10pct_SAMHD1_SLC25A37
var.carriers = var.carriers_lt6_SAMHD1

# mtDNA-CN GWAS (common variants)
gwas <- read_csv('~/mito_rare-variant/resources/mtDNA_CN_GWAS_results.csv')
gwas$A1FREQ.2 <-
  ifelse(gwas$A1FREQ > 0.5, 1 - gwas$A1FREQ, gwas$A1FREQ)
gwas$`EFFECT SIZE ESTIMATE` <-
  ifelse(gwas$A1FREQ > 0.5, -1*gwas$`EFFECT SIZE ESTIMATE`, gwas$`EFFECT SIZE ESTIMATE`)


# rv results
library(vroom)
MAC.cutoff=6
res.dir='/Users/vkp/projects/mito_rare-variant/resources/results/single_var/'
N=414665

# Run2 = non inverse rank-normalized phenotype
res <- vroom(paste0(res.dir,'bolt.lmm_n450k_exomes.chr1-22.adjustedforGWA_SNPs.stats.gz'), progress=F)
res <- res %>% mutate(MAC=(A1FREQ*N)*2)
res <- res %>% filter(MAC >= MAC.cutoff)
res <- res %>% mutate(bonf.sig=ifelse(P_BOLT_LMM <= (0.05 / nrow(res)),TRUE,FALSE))
res <- res %>% mutate(BETA.dir=ifelse(BETA>0,'+','-'))
qobj <- qvalue(res$P_BOLT_LMM)
res$qvalue <- qobj$qvalues
res$fdr.sig <- qobj$significant
rv <- res %>% filter(A1FREQ < 0.002 & qvalue <= 0.1)
#rv %>% write_tsv('~/mito_rare-variant/resources/mtrv.bolt.fdr0.1_maf_lt0.002.hits.txt') # manually remove closely positioned variants, keeping top P-value in each locus
rv <- read_tsv('~/mito_rare-variant/resources/mtrv.bolt.fdr0.1_maf_lt0.002.hits.txt') # post-manual paring of rv

gwas_ <- gwas %>% select(SNP, A1FREQ=A1FREQ.2, BETA=`EFFECT SIZE ESTIMATE`)
rv_ <- rv %>% select(SNP, A1FREQ, BETA)
gwas_rv <- rbind(gwas_, rv_)
gwas_rv <- gwas_rv %>% mutate(BETA.dir=ifelse(BETA>0,'+','-'))
gwas_rv <- gwas_rv %>% mutate(BETA_abs=abs(BETA))

# Plot both
gwas_rv %>%
  ggplot(aes(x=A1FREQ, y=BETA)) +
  geom_point(size=0.8) +
  hrbrthemes::theme_ipsum_rc() +
  scale_x_log10() +
  geom_vline(xintercept = min(gwas$A1FREQ.2), color='grey', linetype='dashed') +
  labs(title='Genetic variants influencing mtDNA-CN')
  
# titrate effect size, count beta directions











# IN OTHER PAPERS
w <- read_tsv('/Users/vkp/projects/mito_rare-variant/analyses/investigation_bolt_beta.signs/Barton_etal_2021_ukb_exome_res.txt') %>% 
  tidyr::separate(TRAIT, into=c('TRAIT_CAT','TRAIT'), sep="_", extra = 'merge')

w %>% count(TRAIT_CAT)
# # A tibble: 6 × 2
# TRAIT_CAT        n
# < chr >        < int >
#   1 biochemistry   646
# 2 blood          442
# 3 bmd             11
# 4 body            65
# 5 bp              18
# 6 lung             7

chisq.betas <- function(x){
  # trait <- x$TRAIT[1]
  trait_cat <- x$TRAIT_CAT[1]
  betas <- x$BETA_Sign
  n = length(betas)
  npos = length(which(betas=="+"))
  nneg = length(which(betas=="-"))
  test <- chisq.test(c(npos,nneg))
  val <- tibble(trait_cat = trait_cat,
                chisq = test$statistic,
                pval = test$p.value)
  return(val)
}

chisq.betas1 <- function(x){
  betas <- x$BETA_Sign
  n = length(betas)
  npos = length(which(betas=="+"))
  nneg = length(which(betas=="-"))
  test <- chisq.test(c(npos,nneg))
  return(test$p.value)
}

trait.chisq <- w %>% group_by(TRAIT) %>% group_map(~ chisq.betas1(.x)) %>% simplify2array()
names(trait.chisq) <- unique(w$TRAIT)

trait_cat.chisq <- w %>% group_by(TRAIT_CAT) %>% group_map(~ chisq.betas1(.x))
names(trait_cat.chisq) <- unique(w$TRAIT_CAT)




# Individual Traits
traits <- unique(w$TRAIT)
# test.res <- 
#   tibble(
#   trait = traits,
#   npos = rep(NA, length(traits)),
#   nneg = rep(NA, length(traits)),
#   chisq = rep(NA, length(traits)),
#   p.value = rep(NA, length(traits))
# )
test.res <- list()
for(i in 1:length(traits)){
  trait <- traits[i]
  tmp <- w %>% filter(TRAIT==trait)
  trait_cat <- tmp$TRAIT_CAT[1]
  ntot = nrow(tmp)
  npos = length(which(tmp$BETA_Sign=="+"))
  nneg = length(which(tmp$BETA_Sign=="-"))
  npos_nneg_ratio <- npos/nneg
  test <- chisq.test(c(npos,nneg))
  test.res[[i]] <- tibble(trait=trait, ntot=ntot,npos=npos, nneg=nneg, npos_nneg_ratio=npos_nneg_ratio, chisq=test$statistic, pvalue=test$p.value, trait_cat=trait_cat)
  message('trait: ', trait, ' | npos: ', npos, ' | nneg: ', nneg, ' | p.value = ', test$p.value)
}
test.res <- bind_rows(test.res)


# Trait categories
trait_cats <- unique(w$TRAIT_CAT)
test.res <- list()
for(i in 1:length(trait_cats)){
  trait_cat <- trait_cats[i]
  tmp <- w %>% filter(TRAIT_CAT==trait_cat)
  ntot = nrow(tmp)
  npos = length(which(tmp$BETA_Sign=="+"))
  nneg = length(which(tmp$BETA_Sign=="-"))
  test <- chisq.test(c(npos,nneg))
  test.res[[i]] <- tibble(trait_cat=trait_cat, ntot=ntot,npos=npos, nneg=nneg, chisq=test$statistic, pvalue=test$p.value)
  message('trait_cat: ', trait_cat, ' | npos: ', npos, ' | nneg: ', nneg, ' | p.value = ', test$p.value)
}
test.res <- bind_rows(test.res)





# Vuckovic et al. 2020 blood cell traits
# Individual Traits
t <- read_tsv('~/mito_rare-variant/analyses/investigation_bolt_beta.signs/Vuckovic_etal_2020_blood_cell_traits_ukb.txt')
x <- t
x <- t %>% filter(`Minor Allele Frequency` <= 0.01)
traits <- unique(x$`Associated Blood Index`)
test.res <- list()
for(i in 1:length(traits)){
  trait <- traits[i]
  tmp <- x %>% filter(`Associated Blood Index`==trait)
  # trait_cat <- tmp$TRAIT_CAT[1]
  ntot = nrow(tmp)
  # npos = length(which(tmp$`(UNIVAR) Estimate of Additive Allelic Effect (REF=Baseline, ALT=Effect)`>0))
  # nneg = length(which(tmp$`(UNIVAR) Estimate of Additive Allelic Effect (REF=Baseline, ALT=Effect)`<0))
  npos = length(which(tmp$`(MULTI) Estimate of Additive Allelic Effect (REF=Baseline, ALT=Effect)`>0))
  nneg = length(which(tmp$`(MULTI) Estimate of Additive Allelic Effect (REF=Baseline, ALT=Effect)`<0))
  npos_nneg_ratio <- npos/nneg
  test <- chisq.test(c(npos,nneg))
  test.res[[i]] <- tibble(trait=trait, ntot=ntot,npos=npos, nneg=nneg, npos_nneg_ratio=npos_nneg_ratio, chisq=test$statistic, pvalue=test$p.value)
  message('trait: ', trait, ' | npos: ', npos, ' | nneg: ', nneg, ' | p.value = ', test$p.value)
}
test.res <- bind_rows(test.res)
qobj <- qvalue(test.res$pvalue)
test.res$qvalue <- qobj$qvalues

# table
require(gt)
require(gtExtras)
test.res %>% 
  arrange(pvalue) %>%
  gt() %>%
  fmt_scientific(columns = c('pvalue','qvalue'), decimals = 2) %>%
  tab_style(
    style = list(cell_text(color = "red")),
    locations = cells_body(columns = c('pvalue'),
                           rows = pvalue < (0.05/29))) %>% 
  tab_header(
    title = md("BETA directionality in Vuckovic et al. 2020 blood cell traits"),
    subtitle = md("MAF < 0.01 (1%)")
  )
  gt_theme_nytimes()


trait='HLSR#'
trait='RET#'
trait='MRV'
tmp <- x %>% filter(`Associated Blood Index`==trait)
tmp %>% ggplot(aes(x=`Minor Allele Frequency`,y=`(MULTI) Estimate of Additive Allelic Effect (REF=Baseline, ALT=Effect)`)) +
  geom_point(size=1) + 
  scale_x_log10() +
  hrbrthemes::theme_ipsum_rc(axis_title_just = 'r') +
  labs(title=paste0(trait))
  





















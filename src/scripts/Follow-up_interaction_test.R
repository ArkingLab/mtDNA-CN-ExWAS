library(dplyr)
library(ggplot2)
library(readr)
library(hrbrthemes)
library(ggtext)
library(ggsignif)

# Helper function
getobj <- function (Rdata)
{
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

# Load info
k <- getobj('~/mito_rare-variant/resources/n428k.Subset_removed.nonwhites.nonconsent.exomes_covar_AnnotatedDataFrame.rds')
k <- as_tibble(k@data) %>% filter(!is.na(arrayCN_PCAodd_m2))


# SNP 1 - SPAG4 (20:35619007:C:G) (upstream)
# Read in SPAG4 variants + diff metrics of mtDNA CN + covars
snp1 <- read_delim('~/mito_rare-variant/analyses/lmm_adjusted/single-variant_followups/spag4.samhd1/snp.20:35619007:C:G.raw', delim=' ')
snp1 <- snp1 %>% dplyr::select(IID,het=`20:35619007:C:G_G`)

# SNP 2 - SAMHD1 (20_36893060_C_T) (downstream)
# Read in SAMHD1 variants + diff metrics of mtDNA CN + covars
snp2 <- read_delim('~/mito_rare-variant/analyses/lmm_adjusted/single-variant_followups/spag4.samhd1/snp.20:36893060:C:T.raw', delim=' ')
snp2 <- snp2 %>% dplyr::select(IID,het=`20:36893060:C:T_T`)


# Compare individuals in snp1 to snp2
# snp1.iids <- snp1 %>% filter(het==1) %>% dplyr::select(IID)
# snp2.iids <- snp2 %>% filter(het==1) %>% dplyr::select(IID)
snp1.iids <- snp1 %>% filter(het!=0) %>% dplyr::select(IID)
snp2.iids <- snp2 %>% filter(het!=0) %>% dplyr::select(IID)

message('# IIDs overlap between SNP1 (SPAG4, upstream) and SNP2 (SAMHD1, downstream): ', 
        length(which(snp1.iids$IID %in% snp2.iids$IID)), 
        '\n out of SNP1 HETs: ', length(snp1.iids$IID),
        '\n and SNP2 HETs: ', length(snp2.iids$IID))

# join snp1 and snp2
j <- k %>%
  left_join(snp1, by=c('sample.id'='IID')) %>%
  left_join(snp2, by=c('sample.id'='IID')) # %>% 
# left_join(ukb, by=c('sample.id'='IID'))
j$het.x <- as.factor(j$het.x)
j$het.y <- as.factor(j$het.y)
j <- j %>% mutate(het.xy=as.factor(ifelse(het.x==1 & het.y==1, 1, 0)))
# j <- j %>% mutate(het.xmiss.y=as.factor(ifelse(is.na(het.x) & het.y==1, 1, 0)))
# j <- j %>% mutate(het.x.ymiss=as.factor(ifelse(het.x==1 & is.na(het.y), 1, 0)))
# j <- j %>% mutate(het.x_only=as.factor(ifelse(het.x==1 & het.y==0, 1, 0)))
# j <- j %>% mutate(het.y_only=as.factor(ifelse(het.x==0 & het.y==1, 1, 0)))

# Haplotype Group INDICES
hap.group_idx <- j %>% group_by(het.x, het.y) %>% group_indices()
j <- j %>% mutate(haplotype=hap.group_idx)

# Adjust mtDNA CN for age,sex,PC1:40,wes.batch
j$mtDNA_CN.adj <- resid(lm(formula=paste0('arrayCN_PCAodd_m2 ~ age+sex+',paste('PC',1:40,sep='',collapse='+'),'+wes.batch'), data=j))





# pc plots of snp1,2 vs rest of individuals
# g <- j %>% ggplot() + geom_point(data=subset(j, het.x==0 & het.y==0), aes(x=PC1,y=PC2), size=0.5, alpha=0.5, color='blue') + geom_point(data=subset(j, het.x==1), aes(x=PC1,y=PC2), size=0.8, alpha=0.8) + geom_point(data=subset(j, het.y==1), aes(x=PC1,y=PC2), size=0.8, alpha=0.8,) #+ 
#   #theme_ipsum_rc() + scale_color_discrete(name = "SNP", labels = c("20:36893060_C_T", "20:35619007_C_G")) + theme(legend.position = 'bottom') + xlim(c(-25,25)) + ylim(c(-25,25))
# g

j %>%
  # inner_join(ukb, by=c('sample.id'='IID')) %>%
  ggplot(aes(x=PC1,y=PC2)) + 
  geom_point(alpha=0.5,size=0.2) + 
  # geom_point(data=subset(j,het.xy==1),aes(x=PC1,y=PC2, color='red'),size=0.3) +
  geom_point(data=subset(j,haplotype==5 | haplotype==6),aes(x=PC1,y=PC2, color='red'),size=0.3) +
  # scale_color_continuous(type = 'viridis') +
  # theme_modern_rc() +
  theme_bw() +
  labs(title='Principal Components Analysis of n=415,422 IIDs', subtitle = 'Red dots: n=262 IIDs with SPAG4 .. SAMHD1 haplotype')+
  theme(legend.position = 'none') +
  xlim(-25,25) + 
  ylim(-25,25)

# j %>% filter(het.xy==1) %>%
#   # inner_join(ukb, by=c('sample.id'='IID')) %>%
#   ggplot(aes(x=PC1,y=PC2,color=mtDNA_CN)) + 
#   geom_point() + 
#   scale_color_continuous(type = 'viridis') +
#   theme_modern_rc()






# mtDNA create haplotype of SNP1 -> SNP2
# Summarise mtDNA by group

j %>% group_by(het.x) %>% count()
j %>% group_by(het.y) %>% count()

# unadjused cn
j %>% 
  group_by(het.x,het.y) %>%
  summarise(n(), round(mean(mtDNA_CN),1)) %>%
  dplyr::select(SNP1=het.x, SNP2=het.y, n=`n()`, 
                `Mean mtDNA_CN`=`round(mean(mtDNA_CN), 1)`)

# adjusted cn
j %>% 
  group_by(het.x,het.y) %>%
  dplyr::summarise(n(), round(mean(mtDNA_CN.adj),3)) %>%
  dplyr::select(SNP1=het.x, SNP2=het.y, n=`n()`, 
                `Mean Adj mtDNA_CN`=`round(mean(mtDNA_CN.adj), 3)`)


# t.test adjusted mtDNA_CN
t.test(j$mtDNA_CN.adj[j$haplotype==1], j$mtDNA_CN.adj[j$haplotype==2])# comparing 0-0 and 0-1
t.test(j$mtDNA_CN.adj[j$haplotype==1], j$mtDNA_CN.adj[j$haplotype==4])
t.test(j$mtDNA_CN.adj[j$haplotype==1], j$mtDNA_CN.adj[j$haplotype==5]) # comparing 0-0 and 1-1
# t.test(j$mtDNA_CN.adj[j$haplotype==1], j$mtDNA_CN.adj[j$haplotype==6]) 

# Combine 1-1 and 2-2 together (adding in the homozygote) and test CN association
j$haplotype2 <- rep(0,nrow(j))
j$haplotype2[which(j$haplotype==5)] <- 1
j$haplotype2[which(j$haplotype==6)] <- 1
t.test.res <- t.test(j$mtDNA_CN.adj[j$haplotype==1], j$mtDNA_CN.adj[j$haplotype2==1]) # comparing 0-0 and [1-1 + 2-2]
# The association to mtDNA-CN was only significant among carriers of both variants
t.test.res$p.value



# Anova
t <- j %>% dplyr::select(SNV1=het.x, SNV2=het.y, both_var=haplotype2, mtDNA_CN.adj)
t <- t %>% filter(!is.na(SNV1) & !is.na(SNV2))

# base_model <- lm(formula=paste0('arrayCN_PCAodd_m2 ~ age+sex+',paste('PC',1:10,sep='',collapse='+'),'+wes.batch','+het.y'), data=tmp)
# full_model <- lm(formula=paste0('arrayCN_PCAodd_m2 ~ age+sex+',paste('PC',1:10,sep='',collapse='+'),'+wes.batch','+het.y','+het.x'), data=tmp)

base_model <- lm(mtDNA_CN.adj ~ SNV2, data=t)
full_model <- lm(mtDNA_CN.adj ~ SNV2 + SNV1, data=t)
full_model2.1 <- lm(mtDNA_CN.adj ~ SNV2 + SNV1 + both_var, data=t)
full_model2.2 <- lm(mtDNA_CN.adj ~ SNV2 + both_var, data=t)
interact_model <- lm(mtDNA_CN.adj ~ SNV2 * SNV1, data=t)
a1 <- anova(base_model, full_model)
a2.1 <- anova(base_model, full_model2.1)
a2.2 <- anova(base_model, full_model2.2)

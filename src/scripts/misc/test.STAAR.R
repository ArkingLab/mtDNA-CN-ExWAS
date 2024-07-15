# Fit null model for STAAR
# nullmodel.glmmkin <- fit_null_glmmkin(fixed=as.formula(paste('mtDNA_CN ~ age + sex + ', paste0('PC', 1:40, collapse='+'), ' + wes.batch')),
#                                       data = k, family=gaussian(link='identity'), 
#                                       id='sample.id', kins=grm, verbose=T)
nullmodel.glmmkin <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n180256.nullmod.STAAR.w_grm.w_PCs.rds')
k <- TopmedPipeline::getobj('~/mtrv/burden/genesis/gds/n180256.keep_samples.AnnotatedDataFrame.RData')
t <- TopmedPipeline::getobj(paste0(dir,'UKBexomeQQFE_chr1-24.cadd_scores_for_variant_ids.rds')) # CADD scores
gds <- seqOpen('~/mtrv/burden/genesis/gds/UKBexomeQQFE_chr1-24.gds')
seqSetFilter(gds, sample.sel = which(samples %in% k@data$sample.id),verbose=T, action='push+set')
var <- seqGetData(gds, 'variant.id')

# assoc <- load('~/mtrv/burden/genesis/run4/cds_nonsyn_0.01/assoc_aggregate.SMMAT.AF.max.0.01.20.results.RData')
assoc <- TopmedPipeline::getobj('~/mtrv/burden/genesis/run4/cds_nonsyn_0.01/assoc_aggregate.SMMAT.AF.max.0.01.10.results.RData')
pval_staar_o <- rep(NA, nrow(assoc[[1]]))
for(i in 1:nrow(assoc[[1]])){
  skip_to_next <- FALSE
  gene=assoc[[1]]$gene[i]
  message(gene)
  set.var_ids <- tibble(assoc[[2]][gene][[1]])$variant.id
  cadd_phred =t$CADD_PHRED[set.var_ids]
  cadd_phred[which(is.na(cadd_phred))] = 0
  pvalues = tryCatch(test.STAAR(nullmodel.glmmkin=nullmodel.glmmkin, seqData=gds, var=var, set.var_ids=set.var_ids, min.maf=0.00001, weights = cadd_phred, weight.cutoff=NULL), 
           error = function(e) { message('gene ', gene, ' failed.'); skip_to_next <<- TRUE; pvalues=data.frame(results_STAAR_O=1)})
  pval_staar_o[i] <- pvalues$results_STAAR_O
}
res <- tibble(gene=assoc[[1]]$gene, pval_staar_o)

# Run STAAR
test.STAAR <- function(nullmodel.glmmkin, seqData, var, set.var_ids, min.maf=0.01, weights=NULL, weight.cutoff=NULL){
  require(STAAR)
  require(SeqVarTools)
  if(!is.null(weight.cutoff)){
    set.var_ids=set.var_ids[-which(weights < weight.cutoff)]
    weights=weights[-which(weights < weight.cutoff)]
  }
  seqSetFilter(gds, variant.sel=which(var %in% set.var_ids), action='push+set')
  geno <- altDosage(gds)
  #message('# samples: ', length(seqGetData(seqData, 'sample.id')))
  # message('Dimensions of geno matrix: ', dim(geno)[1], ' x ', dim(geno)[2])
  seqSetFilter(gds, action='pop')
  pvalues <- STAAR(genotype=geno, obj_nullmodel=nullmodel.glmmkin, rare_maf_cutoff=min.maf, annotation_phred=weights)
  return(pvalues)
}
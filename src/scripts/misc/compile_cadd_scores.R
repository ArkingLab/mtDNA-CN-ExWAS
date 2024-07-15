library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
files <- list.files(path='/dcs04/arking/data/active/projects/mito_rare-variant/annotations/cadd_scored/', pattern='ukb23156_c([0-9]+|X)_cadd_annotated.txt.gz.col_1_14.txt')

for(i in 2:length(files)){
	message(files[i])
	t <- as_tibble(fread(files[i], header=F, data.table=F, showProgress=T))
	
	#t <- t %>% dplyr::mutate(CADD_PHRED=stringr::str_extract(V2, "CADD_PHRED=[0-9]\\d*(\\.\\d+)?")) %>% dplyr::mutate(CADD_PHRED2=stringr::str_extract(CADD_PHRED, "[0-9]\\d*(\\.\\d+)?"))
	t <- t %>% dplyr::mutate(CADD_RAW=stringr::str_extract(V2, "CADD_RAW=[-]?[0-9]\\d*(\\.\\d+)?")) %>% dplyr::mutate(CADD_RAW2=stringr::str_extract(CADD_RAW, "[-]?[0-9]\\d*(\\.\\d+)?"))

	#t <- t %>% dplyr::select(-V2, -CADD_PHRED) %>% dplyr::mutate(CADD_PHRED=as.numeric(CADD_PHRED2)) %>% dplyr::select(-CADD_PHRED2)
	t <- t %>% dplyr::select(-V2, -CADD_RAW) %>% dplyr::mutate(CADD_RAW=as.numeric(CADD_RAW2)) %>% dplyr::select(-CADD_RAW2)

	t <- t %>% tidyr::separate(col=V1, sep='_', into=c('chr','pos','ref','alt'))
	
	#t <- t %>% dplyr::filter(!is.na(CADD_PHRED)) %>% distinct(pos, ref, alt, .keep_all=T)
	t <- t %>% dplyr::filter(!is.na(CADD_RAW)) %>% distinct(pos, ref, alt, .keep_all=T)
	
	fname <- paste0(stringr::str_remove(files[i], '_annotated.txt.gz.col_1_14.txt'),'_scores.txt')
	t %>% readr::write_tsv(fname)
}



# later..
# compile list of cadd scores
cadd <- list()
cadd.scores <- list.files(path='~/mtrv/annotations/cadd_scored/',pattern='*_cadd_scores.txt', full.names=T)
for(i in 1:length(cadd.scores)){message(cadd.scores[i]); cadd[[i]]<-read_tsv(cadd.scores[i], col_names=T, progress=T)}
cadd <- bind_rows(cadd) %>% arrange(chr,pos)

# line it up with gds-specific variant.ids
x <- tibble(variant.id=seqGetData(seqData,'variant.id'), chr=seqGetData(seqData,'chromosome'), pos=seqGetData(seqData,'position'))
t <- x %>% left_join(cadd) %>% distinct(variant.id, .keep_all=T)
# save to file
saveRDS(object=t, file='../gds/UKBexomeQQFE_chr1-24.cadd_scores_for_variant_ids.rds')

# load cadd weights for each variant.id
t <- TopmedPipeline::getobj(paste0(dir,'UKBexomeQQFE_chr1-24.cadd_scores_for_variant_ids.rds'))

# attach cadd weights as AnnotatedDataFrame to variantData of gds
cadd.weights <-  Biobase::AnnotatedDataFrame(data.frame(variant.id=t$variant.id, weight=t$CADD_PHRED))
seqData <- SeqVarData(gds, sampleData=a, variantData=cadd.weights)

x <- tibble(variant.id=seqGetData(seqData,'variant.id'), chr=seqGetData(seqData,'chromosome'), pos=seqGetData(seqData,'position'))
#!/usr/bin/env Rscript
# processes VEP output to filter on CADD cutoffs
# expected:
# VEP output filtered for CADD18

suppressPackageStartupMessages(library(argparser, quietly = T))

p <- argparser::arg_parser("Filter VEP output for variants that are >= CADD score.")
p <- argparser::add_argument(p, "chr",
                             help="Chromosome, values: 1-22,23=X,24=Y")
p <- argparser::add_argument(p, "--vep_dir",
                             default='/dcs04/arking/data/active/projects/mito_rare-variant/annotations/cadd_scored/', 
                             help='Top level of vep directory.')
p <- argparser::add_argument(p, "--out_dir",
                             default='/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/qual.var/variants_in_cds_CADDphred_90thPercentiile/', 
                             help='Top level of out directory.')
p <- argparser::add_argument(p, "--quantile",
                             default=0.9, 
                             help='Quantile of CADD score to select (0 - 1)')
argv <- argparser::parse_args(p)

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

# dir = '~/projects/mito_rare-variant/'
dir=argv$vep_dir
chr=argv$chr
t <- read_tsv(paste0(dir,'/ukb23156_c',chr,'_cadd_annotated.txt.gz'), 
              col_names = F,
              progress = T)

# extract CADD score
cadd.exp = 'CADD_PHRED=([0-9]*.[0-9]*|[0-9])'
t$cadd <- as.numeric(str_remove(str_extract(t$X14, cadd.exp),'CADD_PHRED='))
cadd.cutoff <- quantile(t$cadd, c(argv$quantile), na.rm = T) # 90th percentile
message('CADD cutoff corresponding to quantile ',
        argv$quantile,
        ' selected: ',
        cadd.cutoff)

# filter on CADD and nonsynonymous ("/" in column 11 indicates amino acid change)
t <- t %>%
  filter(cadd >= cadd.cutoff) %>%
  filter(str_detect(X11, '/')) 

pos <- t %>%
  separate(col=X2, into = c('chr','pos'), sep = ':') %>%
  dplyr::select(pos) %>%
  distinct(pos)
# split pos with '-' (which indicate ranges) into 
pos <- str_replace(pos$pos, '-[0-9]*', '')

message('Selected ', length(pos), ' unique positions that match filter.')

out_dir=argv$out_dir
write_tsv(x = as.data.frame(pos),
          file = paste0(out_dir, '/variants_in_cds_CADDphred_18.chr',chr,'.uniq_pos'),
          col_names=F)
message('File written.')

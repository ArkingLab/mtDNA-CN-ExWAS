# Leave one out analysis for a single gene
suppressPackageStartupMessages(library(argparser, quietly = T))
p <- argparser::arg_parser("Leave-one-out association testing using GENESIS")
p <- argparser::add_argument(p, "gene",
                             help="Gene name")
p <- argparser::add_argument(p, "chr",
                             help="Chromosome, values: 1-22,23=X,24=Y")
p <- argparser::add_argument(p, "--AF.max",
                             default=0.01,
                             help="Upper threshold of MAF for variants to include")
p <- argparser::add_argument(p, "--res.dir",
                             default='./',
                             help="Directory for outputting results.")
argv <- argparser::parse_args(p)

source('leave_one_out3.R')

message('-----------------------------------------------------')
# process input args
gene=argv$gene
chr=argv$chr
AF.max=argv$AF.max
res.dir=argv$res.dir

message('gene: ', gene)
message('chr: ', chr)
message('AF.max: ', AF.max)
message('Results directory: \n', res.dir)

res <- loo(gene=gene, chr=chr, AF.max=AF.max)

# save
save(list = c('res'),
     file = paste0(res.dir,'/loo_assoc.aggr.AF.max.', AF.max,'.','chr',chr,'.gene.',gene,'.RData'))
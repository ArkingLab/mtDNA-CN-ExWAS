#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparser)
  library(SeqArray)
})

# Set up argument parser
p <- arg_parser("UKBiobank GWAS plink binary fileset to GDS for a single chromosome.")
p <- add_argument(p, "chr", help="Chromosome, values: 1-22,X,Y [or 23,24]")
p <- add_argument(p, "--bed_dir", 
                  default='/dcl01/arking/data/static/UKBiobank/GWAS_data/plink/', 
                  help='Top level of PLINK binary fileset directory.')
p <- add_argument(p, "--output_dir",
                  default="/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/grm.gwas.2/",
                  help="Output directory for gds")
p <- add_argument(p, "--num_cores", 
                  default = 10,
                  help = "Number of cores for running SeqArray::seqBED2GDS")
p <- add_argument(p, "--compression", 
                  default = "LZMA_RA",
                  help = "Compression method for GDS file")

argv <- parse_args(p)

# Construct file paths
chr <- argv$chr
file_prefix <- file.path(argv$bed_dir, paste0("ukb_chr", chr))
bed_file <- paste0(file_prefix, ".bed")
bim_file <- paste0(file_prefix, ".bim")
fam_file <- paste0(file_prefix, ".fam")
output_file <- file.path(argv$output_dir, paste0("ukb_chr", chr, ".gds"))

# Ensure output directory exists
dir.create(argv$output_dir, showWarnings = FALSE, recursive = TRUE)

# Convert PLINK files to GDS
tryCatch({
  seqBED2GDS(
    bed.fn = bed_file,
    bim.fn = bim_file,
    fam.fn = fam_file,
    out.gdsfn = output_file,
    storage.option = argv$compression,
    parallel = as.integer(argv$num_cores),
    verbose = TRUE
  )
  cat(sprintf("Successfully converted chromosome %s to GDS format.\n", chr))
}, error = function(e) {
  cat(sprintf("Error converting chromosome %s: %s\n", chr, e$message))
})

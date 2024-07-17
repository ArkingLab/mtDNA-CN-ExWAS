#!/usr/bin/env Rscript
# Formats a VCF to a GDS
# Uses `SeqArray::seqVCF2GDS`
# Combines UK Biobank WES pVCF blocks (b0 -- bN) for a single chromosome
# Creates merged chromosome-specific pVCF --> GDS
# Author: Vamsee Pillalamarri

suppressPackageStartupMessages({
  library(argparser)
  library(SeqArray)
})

# Set up argument parser
p <- arg_parser("UKBiobank WES pVCF to GDS")
p <- add_argument(p, "chr", help="Chromosome, values: 1-22,X,Y")
p <- add_argument(p, "--pVCF_dir", 
                  default='/dcl01/arking/data/static/UKBiobank/WES/pVCF/', 
                  help='Top level of pVCF directory.')
p <- add_argument(p, "--output_dir",
                  default="/dcs04/arking/data/active/projects/mito_rare-variant/burden/genesis/gds/",
                  help="Output directory for gds")
p <- add_argument(p, "--num_cores", 
                  default = 10,
                  help = "Number of cores for running SeqArray::seqVCF2GDS")
p <- add_argument(p, "--compression", 
                  default = "LZ4_RA",
                  help = "Compression method for GDS file")

argv <- parse_args(p)

# Construct file paths and ensure output directory exists
chr <- argv$chr
input_dir <- argv$pVCF_dir
output_dir <- argv$output_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Find all relevant pVCF files
files <- list.files(
  path = input_dir,
  pattern = sprintf('ukb23156_c%s_b[0-9]+_v1.vcf.gz$', chr),
  full.names = TRUE
)

if (length(files) == 0) {
  stop(sprintf("No matching pVCF files found for chromosome %s", chr))
}

# Construct output file path
output_file <- file.path(output_dir, sprintf('ukb23156_c%s_blockmerged.gds', chr))

# Convert pVCF to GDS
tryCatch({
  seqVCF2GDS(
    files,
    output_file,
    storage.option = argv$compression,
    parallel = as.integer(argv$num_cores),
    verbose = TRUE
  )
  cat(sprintf("Successfully converted chromosome %s pVCF blocks to GDS format.\n", chr))
}, error = function(e) {
  cat(sprintf("Error converting chromosome %s: %s\n", chr, e$message))
})

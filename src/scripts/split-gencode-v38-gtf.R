suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(parallel)
})

# Define constants
DIR <- '~/mtrv/burden/genesis/gds/'
GTF_PATH <- '~/mtrv/burden/genesis/resources/gencode.v38.annotation.gtf.gz'

# Import GTF file
gtf <- rtracklayer::import(GTF_PATH)

# Filter GTF for protein-coding genes
protein_coding_genes <- gtf[gtf$type == 'gene' & gtf$gene_type == 'protein_coding']

# Get unique chromosomes
chrs <- unique(seqnames(protein_coding_genes))

# Function to process each chromosome
process_chromosome <- function(chr) {
  message(chr)
  
  x <- protein_coding_genes[seqnames(protein_coding_genes) == chr]
  gene_symbols <- mcols(x)$gene_name
  
  # Simplify metadata
  mcols(x) <- DataFrame(gene_name = gene_symbols)
  
  # Set sequence style and names
  seqlevelsStyle(x) <- 'NCBI'
  names(x) <- make.names(gene_symbols, unique = TRUE)
  
  # Save as RDS file
  output_file <- file.path(DIR, paste0('gencode_v38.', chr, '_protein_coding_genes.gtf.rds'))
  saveRDS(object = x, file = output_file)
  
  return(NULL)
}

# Process chromosomes in parallel
mclapply(chrs, process_chromosome, mc.cores = parallel::detectCores() - 1)

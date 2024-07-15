# SPAG4 and SAMHD1 SNP Analysis
# Author: Vamsee Pillalamarri

# Load required libraries
library(dplyr)
library(readr)
library(SeqArray)
library(Rfast)
library(ggplot2)
library(RNOmni)

# Helper functions
load_data <- function(file_path) {
  TopmedPipeline::getobj(file_path)
}

process_snp_data <- function(file_path, snp_id) {
  snp_data <- read_delim(file_path, delim=' ')
  snp_data <- snp_data %>% dplyr::select(IID, het = !!snp_id)
  snp_iids <- snp_data %>% filter(het == 1) %>% dplyr::select(IID)
  list(data = snp_data, iids = snp_iids)
}

# 1. Data Loading and Preparation
# Load sample data for 180256 IIDs
sample_data <- load_data('~/mtrv/burden/genesis/gds/n180256.keep_samples.AnnotatedDataFrame.RData')
bases_mapped <- read_tsv('/dcs04/arking/data/ukb_wes/qcstats/n200k.bases_mapped.txt', 
                         col_names=c('IID','bases.mapped'))

# Open GDS file and filter samples
gds_file <- '/dcs04/arking/data/active/projects/mito_rare-variant/lmm_adjusted/single-variant_followups/spag4.samhd1/ukb23156_c20_b11_v1.vcf.gz.gds'
gds <- seqOpen(gds_file)
samples <- seqGetData(gds, 'sample.id')
seqSetFilter(gds, sample.sel = which(samples %in% sample_data@data$sample.id), verbose=TRUE)
positions <- seqGetData(gds, 'position')

# Load SNP data
snp1 <- process_snp_data('~/mtrv/lmm_adjusted/single-variant_followups/spag4.samhd1/snp_20:35619007:C:G.raw', 
                         '20:35619007:C:G_HET')
snp2 <- process_snp_data('~/mtrv/lmm_adjusted/single-variant_followups/spag4.samhd1/20:36893060:C:T.raw', 
                         '20:36893060:C:T_HET')

# 2. SNP Overlap Analysis
message('# IIDs overlap between SNP1 (SPAG4, upstream) and SNP2 (SAMHD1, downstream): ', 
        length(intersect(snp1$iids$IID, snp2$iids$IID)), 
        '\n out of SNP1 HETs: ', nrow(snp1$iids),
        '\n and SNP2 HETs: ', nrow(snp2$iids))

# 3. Data Integration
k <- as_tibble(sample_data@data) %>% mutate(wes.batch = as.factor(wes.batch))
j <- k %>%
  left_join(snp1$data, by = c('sample.id' = 'IID')) %>%
  left_join(snp2$data, by = c('sample.id' = 'IID')) %>%
  mutate(
    het.x = as.factor(het.x),
    het.y = as.factor(het.y),
    het.xy = as.factor(ifelse(het.x == 1 & het.y == 1, 1, 0)),
    het.y_only = as.factor(ifelse(het.x == 0 & het.y == 1, 1, 0)),
    het.x_only = as.factor(ifelse(het.x == 1 & het.y == 0, 1, 0)),
    het.x_miss = as.factor(ifelse(is.na(het.x), 1, 0))
  ) %>%
  left_join(bases_mapped, by = c('sample.id' = 'IID')) %>%
  mutate(cov = (bases.mapped / (38997831 * 2)))

# 4. Coverage Analysis
g.cov <- j %>% 
  ggplot(aes(x = cov)) + 
  geom_histogram(bins = 50) +
  geom_vline(xintercept = median(j$cov, na.rm = TRUE), col = 'red', linetype = 'dashed') +
  labs(title = "Coverage Distribution", x = "Coverage", y = "Count")
ggsave(plot = g.cov, filename = 'n180256_wes.cov.pdf', device = 'pdf')

# 5. Read Depth Analysis
snp1_pos <- 35619007
snp2_pos <- 36893060

# Get read depth data
depth_data <- seqGetData(gds, 'annotation/format/DP')
colnames(depth_data) <- seqGetData(gds, 'position')
rownames(depth_data) <- seqGetData(gds, 'sample.id')

# Normalize read depth by coverage
depth_data_norm <- depth_data / j$cov
# Remove samples with NA coverage
na_cov_idx <- which(is.na(j$cov))
depth_data_norm <- depth_data_norm[-na_cov_idx, ]
j <- j[-na_cov_idx, ]

# 6. Regression Analysis Functions
run_regression <- function(depth_data, variant, covariates, positions) {
  results <- apply(depth_data, 2, function(x) {
    model <- lm(x ~ variant + covariates)
    summary(model)$coefficients[2, c(1, 4)]  # beta and p-value for variant
  })
  
  tibble(
    pos = positions,
    beta = results[1, ],
    pvalue = results[2, ]
  ) %>%
    mutate(sig = factor(ifelse(pvalue <= (0.05 / ncol(depth_data)), 
                               'p <= Bonferroni', 'p > Bonferroni')))
}

# 7. Run Regressions
covariates <- model.matrix(~ wes.batch + age + sex, data = j)[, -1]
het_xy_results <- run_regression(depth_data_norm, j$het.xy, covariates, positions)
het_x_miss_results <- run_regression(depth_data_norm, j$het.x_miss, covariates, positions)
het_y_only_results <- run_regression(depth_data_norm, j$het.y_only, covariates, positions)

# 8. Plotting Function
plot_regression_results <- function(results, title, x_range = NULL) {
  p <- ggplot(results, aes(x = 1:nrow(results), y = beta)) +
    geom_vline(xintercept = which(results$pos == snp1_pos), linetype = 'dashed', col = 'red', size = 0.5) +
    geom_vline(xintercept = which(results$pos == snp2_pos), linetype = 'dashed', col = 'red', size = 0.5) +
    geom_point(aes(y = 0, color = sig), size = 0.2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_step(alpha = 1, size = 0.3) +
    theme_bw() +
    labs(title = title, x = "Position", y = "Beta")
  
  if (!is.null(x_range)) {
    p <- p + xlim(x_range)
  }
  
  return(p)
}

# 9. Generate Plots
plot_het_xy <- plot_regression_results(het_xy_results, "Het XY Results", c(15450, 15550))
plot_het_x_miss <- plot_regression_results(het_x_miss_results, "Het X Miss Results", c(15450, 15550))
plot_het_y_only <- plot_regression_results(het_y_only_results, "Het Y Only Results", c(15450, 15550))

# Save plots
ggsave("het_xy_plot.pdf", plot_het_xy)
ggsave("het_x_miss_plot.pdf", plot_het_x_miss)
ggsave("het_y_only_plot.pdf", plot_het_y_only)

# 10. Additional Analysis: Read Depth at Specific SNPs
snp1_depth <- depth_data[, which(positions == snp1_pos)]
snp2_depth <- depth_data[, which(positions == snp2_pos)]

g_snp1 <- ggplot(j, aes(x = het.xy, y = snp1_depth)) +
  geom_boxplot() +
  labs(title = "Read Depth at SNP1", x = "Het XY", y = "Read Depth") +
  theme_bw()

g_snp2 <- ggplot(j, aes(x = het.xy, y = snp2_depth)) +
  geom_boxplot() +
  labs(title = "Read Depth at SNP2", x = "Het XY", y = "Read Depth") +
  theme_bw()

ggsave("snp1_depth.pdf", g_snp1)
ggsave("snp2_depth.pdf", g_snp2)

# 11. Clean up
seqClose(gds)

# 12. Session Info
sessionInfo()

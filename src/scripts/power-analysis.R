# Power Analysis for Quantitative Trait
# This script calculates and visualizes power for a quantitative trait association study.
# Author: Vamsee Pillalamarri, adapted from various sources
# Sources: 
# https://web.pdx.edu/~newsomj/uvclass/ho_power.pdf
# https://www.sfu.ca/~lockhart/richard/350/08_2/lectures/PowerSampleSize/web.pdf
# etc.

library(ggplot2)

# Constants and parameters
GENOME_WIDE_P_VALUE <- 5.0e-8
SAMPLE_SIZE <- 10500
MAF_RANGE <- seq(0.05, 0.5, 0.05)
EFFECT_SIZE_RANGE <- seq(0, 0.3, 0.005)

#' Calculate non-centrality parameter
#'
#' @param n Sample size
#' @param maf Minor allele frequency
#' @param effect_size Effect size in standard deviation units
#' @return Non-centrality parameter
calculate_ncp <- function(n, maf, effect_size) {
  2 * maf * (1 - maf) * n * effect_size^2
}

#' Calculate power
#'
#' @param critical_chisq Critical chi-square value
#' @param ncp Non-centrality parameter
#' @return Power
calculate_power <- function(critical_chisq, ncp) {
  pchisq(critical_chisq, df = 1, ncp = ncp, lower.tail = FALSE)
}

#' Generate power matrix
#'
#' @param maf_range Range of minor allele frequencies
#' @param effect_size_range Range of effect sizes
#' @param n Sample size
#' @param critical_chisq Critical chi-square value
#' @return Matrix of power values
generate_power_matrix <- function(maf_range, effect_size_range, n, critical_chisq) {
  power_matrix <- matrix(NA, nrow = length(maf_range), ncol = length(effect_size_range))
  
  for (i in seq_along(maf_range)) {
    for (j in seq_along(effect_size_range)) {
      ncp <- calculate_ncp(n, maf_range[i], effect_size_range[j])
      power_matrix[i, j] <- calculate_power(critical_chisq, ncp)
    }
  }
  
  return(power_matrix)
}

#' Create power heatmap
#'
#' @param maf_range Range of minor allele frequencies
#' @param effect_size_range Range of effect sizes
#' @param power_matrix Matrix of power values
#' @return ggplot object
create_power_heatmap <- function(maf_range, effect_size_range, power_matrix) {
  data <- expand.grid(MAF = maf_range, EffectSize = effect_size_range)
  data$Power <- as.vector(power_matrix)
  
  ggplot(data, aes(x = MAF, y = EffectSize, fill = Power)) +
    geom_tile() +
    scale_fill_gradientn(colors = gray.colors(10, start = 0.9, end = 0.1)) + # start = 0.7, end = 0
    labs(x = "Cumulative Minor Allele Frequency",
         y = "Effect Size (sd)",
         fill = "Power") +
    theme_minimal() +
    theme(legend.position = "right")
}

# Main execution
main <- function() {
  critical_chisq <- qchisq(GENOME_WIDE_P_VALUE, df = 1, lower.tail = FALSE)
  power_matrix <- generate_power_matrix(MAF_RANGE, EFFECT_SIZE_RANGE, SAMPLE_SIZE, critical_chisq)
  
  plot <- create_power_heatmap(MAF_RANGE, EFFECT_SIZE_RANGE, power_matrix)
  
  # Display the plot
  print(plot)
  
  # Optionally save the plot
  # ggsave("quantitative_power_analysis.pdf", plot, width = 10, height = 8)
}

# Run the analysis
main()

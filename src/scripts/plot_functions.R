# Plotting functions for genetic analysis
# Author: Vamsee Pillalamarri

library(ggplot2)

#' Create a QQ plot
#'
#' @param pval Vector of p-values
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param conf_level Confidence level for the envelope (default: 0.95)
#' @return A ggplot object
#' @export
qqPlot <- function(pval, title = 'QQ Plot', subtitle = 'subtitle', conf_level = 0.95) {
  if (!is.numeric(pval) || any(pval < 0 | pval > 1, na.rm = TRUE)) {
    stop("pval must be a numeric vector with values between 0 and 1")
  }
  
  pval <- pval[!is.na(pval)]
  n <- length(pval)
  
  if (n == 0) {
    stop("No valid p-values provided")
  }
  
  alpha <- (1 - conf_level) / 2
  
  dat <- data.frame(
    obs = sort(pval),
    exp = seq_len(n) / n,
    upper = qbeta(alpha, seq_len(n), rev(seq_len(n))),
    lower = qbeta(1 - alpha, seq_len(n), rev(seq_len(n)))
  )
  
  ggplot(dat, aes(x = -log10(exp), y = -log10(obs))) +
    geom_ribbon(aes(ymin = -log10(lower), ymax = -log10(upper)), fill = "grey80", alpha = 0.5) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    scale_x_continuous(name = expression(paste(-log[10], "(expected P)")),
                       expand = c(0, 0)) +
    scale_y_continuous(name = expression(paste(-log[10], "(observed P)")),
                       expand = c(0, 0)) +
    labs(title = title, subtitle = subtitle) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}

#' Plot Manhattan plot for association results
#'
#' @param assoc Association results
#' @param signif Significance threshold
#' @param out_file Output file path
#' @return A ggplot object
#' @export
plot_assoc_manhattan <- function(assoc, signif = 5e-8, out_file = NULL) {
  # Placeholder for future implementation
  message("plot_assoc_manhattan function is not fully implemented yet.")
  
  # Example of how the function might look when implemented:
  # 
  # require(dplyr)
  # 
  # reduce_assoc <- function(chr_assoc, var) {
  #   data.frame(
  #     chr = chr_assoc$chr,
  #     pos = var$pos,
  #     pval = chr_assoc$pval
  #   )
  # }
  # 
  # plot_data <- bind_rows(mapply(reduce_assoc, assoc[[1]], assoc[[2]], SIMPLIFY = FALSE))
  # 
  # cmap <- setNames(rainbow(length(unique(plot_data$chr))), unique(plot_data$chr))
  # 
  # p <- ggplot(plot_data, aes(x = chr, y = -log10(pval), color = factor(chr))) +
  #   geom_point(position = position_dodge(0.8)) +
  #   scale_color_manual(values = cmap, guide = "none") +
  #   geom_hline(yintercept = -log10(signif), linetype = 'dashed') +
  #   theme_minimal() +
  #   theme(panel.grid.major.x = element_blank(),
  #         panel.grid.minor.x = element_blank(),
  #         axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   labs(x = "Chromosome", y = expression(-log[10](p)))
  # 
  # if (!is.null(out_file)) {
  #   ggsave(out_file, plot = p, width = 10, height = 5)
  # }
  # 
  # return(p)
}

#' Wrapper function to create and save a QQ plot
#'
#' @param pval Vector of p-values
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param out_file Output file path (optional)
#' @param width Plot width (default: 7)
#' @param height Plot height (default: 7)
#' @export
create_and_save_qqplot <- function(pval, title = 'QQ Plot', subtitle = 'subtitle', 
                                   out_file = NULL, width = 7, height = 7) {
  plot <- qqPlot(pval, title, subtitle)
  
  if (!is.null(out_file)) {
    ggsave(out_file, plot = plot, width = width, height = height)
    message("QQ plot saved to: ", out_file)
  }
  
  return(plot)
}

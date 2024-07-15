# Plotting functions
qqPlot <- function(pval, title='QQPLOT', subtitle='subtitle') {
  require(ggplot2)
  pval <- pval[!is.na(pval)]
  n <- length(pval)
  x <- 1:n
  dat <- data.frame(obs=sort(pval),
                    exp=x/n,
                    upper=qbeta(0.025, x, rev(x)),
                    lower=qbeta(0.975, x, rev(x)))
  
  x <- ggplot(dat, aes(-log10(exp), -log10(obs))) +
    geom_line(aes(-log10(exp), -log10(upper)), color="gray") +
    geom_line(aes(-log10(exp), -log10(lower)), color="gray") +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red") +
    xlab(expression(paste(-log[10], "(expected P)"))) +
    ylab(expression(paste(-log[10], "(observed P)"))) +
    labs(title = title, subtitle = subtitle) +
    theme_bw()
  print(x)
  message('plotted.')
}

# DOESN't WORK YET
plot_assoc_manhattan <- function(assoc) {
  # code borrowed from:
  # https://github.com/UW-GAC/analysis_pipeline/blob/master/R/assoc_plots.R
  require(ggplot2)
  # assoc <- assoc %>% filter(pval <= plot_max_p)
  reduce_assoc <- function(chr.assoc, var){
    # reduces single chr assoc results from running GENESIS::assocAggregateTest()
    chr=chr.assoc$chr
    pos = var[]
    data.frame(chr=chr.assoc$chr)
  }
  sapply(assoc[[1]], reduce_assoc, var)
  p <- ggplot(assoc, aes(chr, -log10(pval), group=interaction(chr, pos), color=chr)) +
    geom_point(position=position_dodge(0.8)) +
    scale_color_manual(values=cmap, breaks=names(cmap)) +
    geom_hline(yintercept=-log10(signif), linetype='dashed') +
    theme_bw() +
    theme(legend.position="none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("chromosome") + ylab(expression(-log[10](p)))
  ggsave(config["out_file_manh"], plot=p, width=10, height=5)  
}

# qqPlot(assoc$Score.pval)

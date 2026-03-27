# Supplementary tables

# Gene-specific variant breakdowns ----

# TWNK
x <- readRDS('~/mito_rare-variant/analyses/burden/SAMHD1/TWNK.fit.full.rds')
gt1.TWNK <- 
  tbl_regression(x, pvalue_fun = function(x) prettyNum(x, digits=2)) %>%
  bold_p(t=(0.05 / (length(coef(x))-1))) %>%
  add_glance_source_note() %>%
  add_q() %>%
  add_significance_stars(hide_ci = F, hide_p = F) %>%
  modify_caption("Rare variants in **_TWNK_**") %>%
  modify_header(update = list(label ~ "**Variant**"))

# TFAM
x <- readRDS('~/mito_rare-variant/analyses/burden/SAMHD1/TFAM.fit.full.rds')
gt1.TFAM <- 
  tbl_regression(x, pvalue_fun = function(x) prettyNum(x, digits=2)) %>%
  bold_p(t=(0.05 / (length(coef(x))-1))) %>%
  add_glance_source_note() %>%
  add_q() %>%
  add_significance_stars(hide_ci = F, hide_p = F) %>%
  modify_caption("Rare variants in **_TFAM_**") %>%
  modify_header(update = list(label ~ "**Variant**"))
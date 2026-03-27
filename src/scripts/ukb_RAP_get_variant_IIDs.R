## bash pre-req code (run once)
# dx cd project-G7KB5zQJz55qG0158ZXb2J5p:/vkp/SAMHD1
# dx download SAMHD1_mtDNA_CN_n180k_MAC_lt6_variants.withMAC.txt
# dx download SAMHD1_mtDNA_CN_n180k_MAC_gteq6_variants.withMAC.txt
# dx download SAMHD1_mtDNA_CN_n180k_MAC_gteq6_variants.FDR10pct.txt
# 
# plink to get variants
# less SAMHD1_mtDNA_CN_n180k_MAC_lt6_variants.withMAC.txt | cut -f1 | tail -n +2 > SAMHD1_mtDNA_CN_n180k_MAC_lt6_variants.txt
# plink --bfile ~/ref/UKB\ -mtDNA_2022/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ interim\ 450k\ release/ukb23149_c8_b0_v1 --extract SAMHD1_mtDNA_CN_n180k_MAC_lt6_variants.txt --recode AD --out SAMHD1_mac_lt6_variants
# 
# less SAMHD1_mtDNA_CN_n180k_MAC_gteq6_variants.withMAC.txt | cut -f1 | tail -n +2 > SAMHD1_mtDNA_CN_n180k_MAC_gteq6_variants.txt
# plink --bfile ~/ref/UKB\ -mtDNA_2022/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ interim\ 450k\ release/ukb23149_c8_b0_v1 --extract SAMHD1_mtDNA_CN_n180k_MAC_gteq6_variants.txt --recode AD --out SAMHD1_mac_gteq6_variants
# 
# rm *.nosex *.log

# R
library(dplyr)
library(readr)

# Constants
RECODED_PLINK_LT6 <- 'SAMHD1_mac_lt6_variants.raw'
RECODED_PLINK_GTEQ6 <- 'SAMHD1_mac_gteq6_variants.raw'

VAR_LIST_LT6 <- 'SAMHD1_mtDNA_CN_n180k_MAC_lt6_variants.txt'
VAR_LIST_GTEQ6 <- 'SAMHD1_mtDNA_CN_n180k_MAC_gteq6_variants.txt'
VAR_LIST_GTEQ6_FDR10PCT <- 'SAMHD1_mtDNA_CN_n180k_MAC_gteq6_variants.FDR10pct.txt'

# Function to get variant IIDs
get_var_iids <- function(recoded_plink, var_list) {
  t <- read_table(recoded_plink)
  var <- readLines(var_list)
  iids <- t$IID
  
  var_iids <- lapply(var, function(v) {
    message(v)
    t_ <- t %>% select(all_of(paste0(v, '_HET')))
    iids[which(t_[,1] != 0)]
  })
  
  names(var_iids) <- var
  return(var_iids)
}

# Function to process variant data
process_variants <- function(recoded_plink, var_list, output_prefix) {
  var_iids <- get_var_iids(recoded_plink, var_list)
  var_iids_counts <- sapply(var_iids, length)
  var_iids_tbl <- tibble(variant = names(var_iids_counts), MAC_450k = var_iids_counts)
  
  var_iids_flat <- unlist(var_iids)
  names(var_iids_flat) <- NULL
  
  writeLines(as.character(var_iids_flat), paste0(output_prefix, '_IIDs.in_450kExomes.txt'))
  
  return(list(iids = var_iids_flat, counts = var_iids_tbl))
}

# Process variants
lt6_results <- process_variants(RECODED_PLINK_LT6, VAR_LIST_LT6, 'SAMHD1_mac_lt6')
gteq6_results <- process_variants(RECODED_PLINK_GTEQ6, VAR_LIST_GTEQ6, 'SAMHD1_mac_gteq6')
gteq6_fdr10pct_results <- process_variants(RECODED_PLINK_GTEQ6, VAR_LIST_GTEQ6_FDR10PCT, 'SAMHD1_mac_gteq6_FDR10pct')

# Combine results
var_iids_mac_all <- c(lt6_results$iids, gteq6_results$iids)
var_iids_mac_all_FDR10pct <- c(lt6_results$iids, gteq6_fdr10pct_results$iids)

writeLines(as.character(var_iids_mac_all), 'SAMHD1_mac_all_IIDs.in_450kExomes.txt')
writeLines(as.character(var_iids_mac_all_FDR10pct), 'SAMHD1_mac_all_FDR10pct_IIDs.in_450kExomes.txt')

# Uncomment the following lines if you want to generate plots
# library(ggplot2)
# 
# generate_plot <- function(data, filename) {
#   g <- data %>% 
#     ggplot(aes(x = MAC, y = MAC_450k)) + 
#     geom_point() + 
#     geom_smooth(method = lm) + 
#     labs(x = "MAC in 180k cohort", y = "MAC in 450k cohort")
#   ggsave(plot = g, filename = filename, device = 'pdf')
# }
# 
# var_mac_lt6 <- read_tsv('SAMHD1_mtDNA_CN_n180k_MAC_lt6_variants.withMAC.txt') %>%
#   inner_join(lt6_results$counts, by = c('reformat' = 'variant'))
# 
# var_mac_gteq6 <- read_tsv('SAMHD1_mtDNA_CN_n180k_MAC_gteq6_variants.withMAC.txt') %>%
#   inner_join(gteq6_results$counts, by = c('reformat' = 'variant'))
# 
# var_mac_all <- bind_rows(var_mac_lt6, var_mac_gteq6)
# 
# generate_plot(var_mac_all, 'SAMHD1_mac_all_cohort_comparison.pdf')

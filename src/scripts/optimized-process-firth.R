# Rewritten function to process PheWAS of rare-variant carriers
# Author: Vamsee Pillalamarri

library(logistf)
library(dplyr)
library(furrr)
library(purrr)
library(glue)

# Rewritten function to process PheWAS of rare-variant carriers
# Using `logistif::fit.firth()` and `furrr` parallelization mapping
process.firth <- function(s, data, gene, ncores = 1) {
  # Input validation
  required_cols <- c("phenotype", "note")
  if (!all(required_cols %in% colnames(s))) {
    stop(glue("Missing required columns in 's': {paste(setdiff(required_cols, colnames(s)), collapse = ', ')}"))
  }
  
  if (!all(c(gene, "age", "sex", "wes.batch", paste0("PC", 1:40)) %in% colnames(data))) {
    stop(glue("Missing required columns in 'data'"))
  }
  
  # Set up parallelization
  plan(multisession, workers = ncores)
  
  # Define the function to process each row
  process_row <- function(row) {
    tryCatch({
      phecode <- row$phenotype
      
      # Create formula
      covariates <- if (grepl("sex", row$note)) {
        c("age", "wes.batch", paste0("PC", 1:40))
      } else {
        c("age", "sex", "wes.batch", paste0("PC", 1:40))
      }
      
      formula <- glue("{phecode} ~ {gene}.carrier + {paste(covariates, collapse = ' + ')}")
      
      # Fit Firth's logistic regression
      fit.firth <- logistf(
        formula = as.formula(formula),
        data = data,
        control = logistf.control(maxit = 100, maxstep = 50),
        plcontrol = logistpl.control(maxit = 400, maxstep = 20),
        plconf = c(2),
        model = FALSE
      )
      
      # Extract and return results
      tibble(
        phenotype = phecode,
        beta = fit.firth$coefficients[2],
        SE = sqrt(diag(vcov(fit.firth))[2]),
        p = fit.firth$prob[2]
      )
    }, error = function(e) {
      tibble(
        phenotype = row$phenotype,
        beta = NA_real_,
        SE = NA_real_,
        p = NA_real_,
        error = as.character(e)
      )
    })
  }
  
  # Process all rows in parallel
  results <- future_map_dfr(1:nrow(s), ~process_row(s[.x, ]), .progress = TRUE)
  
  # Clean up parallel backend
  plan(sequential)
  
  return(results)
}

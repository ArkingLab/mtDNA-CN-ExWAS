# Functions for PheWAS analysis
# Author: Vamsee Pillalamarri

library(PheWAS)
library(dplyr)
library(logistf)
library(splines)
library(bigmemory)
library(biganalytics)
library(doParallel)
library(snow)
library(doSNOW)

# Helper function to load data
load_data <- function(path) {
  if (tools::file_ext(path) == "rds") {
    return(readRDS(path))
  } else {
    loaded_obj <- load(path)
    if (length(loaded_obj) > 1) {
      warning("Multiple objects in ", path, ". Returning the first.")
    }
    return(get(loaded_obj[1]))
  }
}

# Function to get phenotype matrix
get_pheno_matrix <- function(codes) {
  IID <- codes$f.eid
  codes <- t(codes[, -1])
  codes <- codes[!is.na(codes)]
  data.frame(
    id = rep(IID, length(codes)),
    vocabulary_id = rep('ICD10', length(codes)),
    code = codes,
    index = rep(1, length(codes))
  )
}

# Function to prepare UKB data
prepare_ukb_data <- function(alldata) {
  list(
    age_sex = alldata %>% 
      select(id = f.eid, age = f.21022.0.0, sex = f.31.0.0) %>% 
      mutate(sex = ifelse(sex == "Male", 'M', 'F')) %>% 
      select(id, age, sex),
    pcs = alldata %>% select(id = f.eid, starts_with("PC")),
    mtDNA = alldata %>% select(id = f.eid, contains('mtDNA'), starts_with('arrayCN')),
    misc = alldata %>% select(id = f.eid, genotyping.array, in.white.British.ancestry.subset, used.in.pca.calculation)
  )
}

# Function to convert ICD codes to phecodes
convert_to_phecodes <- function(icd10_codes, id_sex) {
  my_icd_map <- PheWAS::phecode_map_icd10 %>% mutate(code = gsub("[.]", "", code))
  createPhenotypes(
    id.vocab.code.index = icd10_codes,
    id.sex = id_sex,
    full.population.ids = id_sex$id,
    min.code.count = 1,
    vocabulary.map = my_icd_map
  ) %>% as_tibble()
}

# Function to run PheWAS
run_phewas <- function(phenotypes, genotypes, data, covariates, cores = 1) {
  phewas_ext(
    phenotypes = phenotypes,
    genotypes = genotypes,
    data = data,
    cores = cores,
    covariates = covariates,
    additive.genotypes = FALSE
  )
}

# Function to create analysis data frame for passing into process.glm()
create_phewas_data.input <- function(var.carriers,
                                     subset.to.unrelated = TRUE,
                                     subset.to.white.British.ancestry = TRUE,
                                     filter.to.exomes = TRUE,
                                     ukb.center.path = '/dcl01/arking/data/active/ukb_mitoscore_pheWAS/ukb.center.rds',
                                     ukb.icd10_phecodes.path = '/dcl01/arking/data/active/ukb_mitoscore_pheWAS/n502485.icd10_phecodes.rds',
                                     ukb.dat.path = '/dcl01/arking/data/active/ukb_mitoscore_pheWAS/n502485.ukb.data.rds') {
  # Helper function to load objects from RData or RDS files
  getobj <- function(Rdata) {
    if (tolower(tools::file_ext(Rdata)) == "rds") {
      return(readRDS(Rdata))
    }
    objname <- load(Rdata)
    if (length(objname) > 1) {
      warning("Multiple objects stored in file ", Rdata, "\nReturning only the first object")
    }
    return(get(objname[1]))
  }
  
  # Read pre-req data
  message('Reading pre-req data')
  ukb.center <- getobj(ukb.center.path)
  i <- getobj(ukb.icd10_phecodes.path) # phecodes
  dat <- getobj(ukb.dat.path)$ukb.data # ukb data
  
  # Create synthetic allele of rv carriers
  message('Creating synthetic allele from var.carriers')
  var <- tibble(id = dat$id, var = as.factor(ifelse(dat$id %in% var.carriers$IID, 1, 0)))
  
  # Create dataset
  message('Creating dataset')
  data <- dat %>%
    inner_join(i, by = 'id') %>%
    inner_join(ukb.center, by = c('id' = 'IID')) %>%
    inner_join(var, by = 'id') %>%
    dplyr::select(
      id, age, sex, Center, var, genotyping.array,
      in.white.British.ancestry.subset, used.in.pca.calculation,
      starts_with('PC'), everything()
    ) %>%
    mutate(
      across(c(Center, sex, genotyping.array, in.white.British.ancestry.subset, used.in.pca.calculation), as.factor),
      age2 = age^2
    )
  
  # Filter dataset to exome IIDs only
  if (filter.to.exomes) {
    message('Filtering dataset')
    exomes.iids <- readLines('~/mtrv/450k_exomes.IIDs.list')
    data <- data %>% filter(id %in% exomes.iids)
  }
  
  message('Number of var: ')
  print(data %>% count(var))
  
  # Subset to unrelated IIDs only (optional)
  if (subset.to.unrelated) {
    data <- data %>% filter(used.in.pca.calculation == 1)
    message('Number of var after subsetting to unrelated: ')
    print(data %>% count(var))
  }
  
  # Subset to white IIDs only (optional)
  if (subset.to.white.British.ancestry) {
    data <- data %>% filter(in.white.British.ancestry.subset == 1)
    message('Number of var after subsetting to white British ancestry: ')
    print(data %>% count(var))
  }
  
  # Compute additional data and filter phecodes
  message('Computing additional data and filtering phecodes')
  sex.check <- sapply(select(data, -id, -age, -sex, -Center, -var, -genotyping.array, -in.white.British.ancestry.subset, -used.in.pca.calculation, -starts_with('PC'), -age2), 
                      function(x) length(unique(data$sex[which(x)])) == 2)
  
  n.var.cases <- colSums(data %>% select(-id, -age, -sex, -Center, -var, -genotyping.array, -in.white.British.ancestry.subset, -used.in.pca.calculation, -starts_with('PC'), -age2) & data$var == 1)
  
  phecodes <- names(n.var.cases)
  phecodes <- phecodes[n.var.cases > 0]
  
  return(list(data, sex.check, phecodes))
}

# Function to create analysis data frame for passing into process.glm()
# Input is a tibble with a column `IID` and a corresponding scored synthetic allele column
create_phewas_data.input2 <- function(var.carriers, 
                                      ukb_center_path = '~/ukb/UKB -mtDNA_2022/data_tables/ukb.center.rds',
                                      ukb_icd10_phecodes_path = '~/ukb/UKB -mtDNA_2022/data_tables/n502485.icd10_phecodes.rds',
                                      ukb_data_path = '~/ukb/UKB -mtDNA_2022/data_tables/n502485.ukb.data.rds',
                                      exomes_iids_path = '~/mtrv/450k_exomes.IIDs.list') {
  
  library(dplyr)
  library(tidyr)
  library(readr)
  
  # Helper function to load data
  load_data <- function(path) {
    if (tools::file_ext(path) == "rds") {
      return(readRDS(path))
    } else {
      loaded_obj <- load(path)
      if (length(loaded_obj) > 1) {
        warning("Multiple objects in ", path, ". Returning the first.")
      }
      return(get(loaded_obj[1]))
    }
  }
  
  # Load required data
  message('Reading pre-requisite data')
  ukb.center <- load_data(ukb_center_path)
  i <- load_data(ukb_icd10_phecodes_path)
  dat <- load_data(ukb_data_path)$ukb.data
  
  # Create synthetic allele of rv carriers
  message('Attaching scored allele from var.carriers')
  var <- tibble(id = dat$id) %>%
    left_join(var.carriers, by = c("id" = "IID")) %>%
    mutate(var = coalesce(scored.allele, 0))
  
  # Create and process dataset
  message('Creating and processing dataset')
  data <- dat %>%
    inner_join(i, by = 'id') %>%
    inner_join(ukb.center, by = c('id' = 'IID')) %>%
    inner_join(var, by = 'id') %>%
    select(id, age, sex, Center, var, genotyping.array,
           in.white.British.ancestry.subset, used.in.pca.calculation,
           starts_with('PC'), everything()) %>%
    mutate(across(c(Center, sex, genotyping.array, 
                    in.white.British.ancestry.subset, used.in.pca.calculation), as.factor),
           age2 = age ^ 2)
  
  # Filter dataset to exome IIDs only
  message('Filtering dataset')
  exomes.iids <- read_lines(exomes_iids_path)
  data <- data %>% 
    filter(id %in% exomes.iids) %>%
    filter(in.white.British.ancestry.subset == 1 & used.in.pca.calculation == 1)
  
  # Compute sex check and n.var.cases
  message('Computing additional data and filtering phecodes')
  phecode_cols <- colnames(i)[-1]
  sex_check <- sapply(data[phecode_cols], function(x) length(unique(data$sex[x])) == 2)
  
  n_var_cases <- sapply(data[phecode_cols], function(x) sum(x & data$var != 0))
  
  # Filter phecodes
  phecodes <- phecode_cols[n_var_cases > 0]
  
  return(list(data = data, sex_check = sex_check, phecodes = phecodes))
}

# Function to create analysis data frame for passing into process.glm()
# Input is a tibble with a column `IID` and a corresponding 'qvar' column
# that contains a quantitative phenotype of interest
# This is the independent variable in a model of phecode ~ qvar + cov
create_phewas_data.input3 <- function(var.carriers,
                                      filter.to.exomes = FALSE,
                                      subset.to.unrelated = FALSE,
                                      subset.to.white.British.ancestry = FALSE,
                                      ukb.center.path = '~/ukb/UKB -mtDNA_2022/data_tables/ukb.center.rds',
                                      ukb.icd10_phecodes.path = '~/ukb/UKB -mtDNA_2022/data_tables/n502485.icd10_phecodes.rds',
                                      ukb.dat.path = '~/ukb/UKB -mtDNA_2022/data_tables/n502485.ukb.data.rds',
                                      exomes.iids.path = '~/mtrv/450k_exomes.IIDs.list') {
  
  library(dplyr)
  library(tidyr)
  library(readr)
  
  # Helper function to load data
  load_data <- function(path) {
    if (tools::file_ext(path) == "rds") {
      return(readRDS(path))
    } else {
      loaded_obj <- load(path)
      if (length(loaded_obj) > 1) {
        warning("Multiple objects in ", path, ". Returning the first.")
      }
      return(get(loaded_obj[1]))
    }
  }
  
  # Load required data
  message('Reading pre-requisite data')
  ukb.center <- load_data(ukb.center.path)
  i <- load_data(ukb.icd10_phecodes.path)
  dat <- load_data(ukb.dat.path)$ukb.data
  
  # Create synthetic allele of rv carriers
  message('Attaching qvar from var.carriers')
  var <- tibble(id = dat$id) %>%
    left_join(var.carriers, by = c('id' = 'IID')) %>%
    select(IID = id, var = qvar)
  
  # Create dataset
  message('Creating dataset')
  data <- dat %>%
    inner_join(i, by = 'id') %>%
    inner_join(ukb.center, by = c('id' = 'IID')) %>%
    inner_join(var, by = c('id' = 'IID')) %>%
    select(
      id, age, sex, Center, var, genotyping.array,
      in.white.British.ancestry.subset, used.in.pca.calculation,
      starts_with('PC'), everything()
    ) %>%
    mutate(
      across(c(Center, sex, genotyping.array, 
               in.white.British.ancestry.subset, used.in.pca.calculation), as.factor),
      age2 = age ^ 2
    )
  
  # Apply filters
  if (filter.to.exomes) {
    message('Filtering dataset to exomes')
    exomes.iids <- read_lines(exomes.iids.path)
    data <- data %>% filter(id %in% exomes.iids)
  }
  
  if (subset.to.unrelated) {
    message('Subsetting to unrelated IIDs')
    data <- data %>% filter(used.in.pca.calculation == 1)
    message('Number of var != 0:')
    print(data %>% count(var != 0))
  }
  
  if (subset.to.white.British.ancestry) {
    message('Subsetting to white British ancestry')
    data <- data %>% filter(in.white.British.ancestry.subset == 1)
    message('Number of var != 0:')
    print(data %>% count(var != 0))
  }
  
  # Compute sex.check and n.var.cases
  message('Computing additional data and filtering phecodes')
  phecode_cols <- colnames(i)[-1]
  sex_check <- sapply(data[phecode_cols], function(x) length(unique(data$sex[which(x)])) == 2)
  
  n_var_cases <- sapply(data[phecode_cols], function(x) sum(x & data$var != 0))
  
  # Filter phecodes
  phecodes <- phecode_cols[n_var_cases > 0]
  
  return(list(data = data, sex_check = sex_check, phecodes = phecodes))
}

# Function to run PheWAS on synthetic allele
process.glm <- function(phecodes,
                        both_sexes = NULL,
                        data_desc,
                        ncores = 1,
                        var = 'var',
                        file.conn = NULL,
                        back_correct_FirthSE = TRUE,
                        debug_log = FALSE) {
  
  # Write header to file if streaming results
  if (!is.null(file.conn)) {
    header <- c('phecode', 'beta', 'SE', 'chisq', 'OR', 'p', 'n.cases', 'n.var.cases', 'n.controls', 'n.var.controls', 'case.ctrl_ratio', 'firth_correct')
    readr::write_lines(paste(header, collapse = '\t'), file = file.conn, append = TRUE)
  }
  
  cl <- parallel::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)
  
  results <- foreach(phecode = phecodes,
                     .packages = c('splines', 'logistf', 'biganalytics'),
                     .errorhandling = 'pass') %dopar% {
    
    data_ <- bigmemory::attach.big.matrix(data_desc)
    
    # Construct formula
    f <- paste0('`', phecode, '`', ' ~ ', var, ' + ns(age, df=2) + sex + Center + ', paste0('PC', 1:40, collapse = '+'))
    if (!is.null(both_sexes) && !both_sexes[as.character(phecode)]) {
      f <- gsub(' \\+ sex', '', f)
    }
    message(f)
    
    # Run GLM
    fit <- biganalytics::bigglm.big.matrix(
      formula(f),
      data_,
      family = binomial(link = 'logit'),
      fc = if (!is.null(both_sexes) && !both_sexes[as.character(phecode)]) c('var', 'Center') else c('var', 'sex', 'Center'),
      quiet = TRUE
    )
    
    # Firth Correct GLM p-value if it is <= 0.05
    firth_correct <- FALSE
    if (summary(fit)$mat[2, 'p'] <= 0.05) {
      message(phecode, ' Firth Correcting...')
      firth_correct <- TRUE
      fit.firth <- logistf(
        formula(f),
        data = as.data.frame(data_[, c(phecode, 'var', 'age', 'sex', 'Center', paste0('PC', 1:40))]),
        control = logistf.control(maxit = 100, maxstep = 50),
        plcontrol = logistpl.control(maxit = 400, maxstep = 20),
        plconf = c(2),
        model = FALSE
      )
    }
    
    # Extract values from GLM or Firth model
    val <- setNames(rep(NA, 12), c('phecode', 'beta', 'SE', 'chisq', 'OR', 'p', 'n.cases', 'n.var.cases', 'n.controls', 'n.var.controls', 'case.ctrl_ratio', 'firth_correct'))
    val['phecode'] <- phecode
    
    if (firth_correct) {
      val['beta'] <- fit.firth$coefficients[2]
      val['OR'] <- exp(fit.firth$coefficients[2])
      val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
      val['chisq'] <- qchisq(1 - fit.firth$prob, df = 1)[2]
      val['p'] <- fit.firth$prob[2]
      if (back_correct_FirthSE) {
        val['SE'] <- abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
      }
    } else {
      val['beta'] <- summary(fit)$mat[2, 'Coef']
      val['OR'] <- exp(summary(fit)$mat[2, 'Coef'])
      val['SE'] <- summary(fit)$mat[2, 'SE']
      val['chisq'] <- qchisq(1 - summary(fit)$mat[2, 'p'], df = 1)
      val['p'] <- summary(fit)$mat[2, 'p']
    }
    
    val['n.cases'] <- sum(data_[, phecode] == 1)
    val['n.var.cases'] <- sum(data_[, 'var'] == 2 & data_[, phecode] == 1)
    val['n.controls'] <- sum(data_[, phecode] == 0)
    val['n.var.controls'] <- sum(data_[, 'var'] == 2 & data_[, phecode] == 0)
    val['case.ctrl_ratio'] <- val['n.cases'] / val['n.controls']
    val['firth_correct'] <- firth_correct
    
    # Write to file if streaming results
    if (!is.null(file.conn)) {
      readr::write_lines(paste(val, collapse = '\t'), file = file.conn, append = TRUE)
    }
    
    val
  }
  
  parallel::stopCluster(cl)
  bind_rows(results) %>% as_tibble()
}

# Function to run PheWAS on scored synthetic allele
process.glm2 <- function(phecodes,
                         both_sexes = NULL,
                         data_desc,
                         ncores = 1,
                         var = 'var',
                         file.conn = NULL,
                         back_correct_FirthSE = TRUE,
                         debug_log = F) {
  require(dplyr)
  require(logistf)
  require(doParallel)
  require(parallel)
  require(splines)
  require(logistf)
  require(bigmemory)
  require(biganalytics)
  require(snow)
  require(doSNOW)
  
  # pre-compute some info about dataset
  # message('Precomputing data: n.var in cases, controls, etc.')
  # sex.check <-
  #   apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) {
  #     ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
  #   })
  # var <- data %>% select(var)
  # n.cases <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x)))
  # n.var.cases <- apply(data %>% dplyr::select(var,all_of(phecodes)), 2, function(x) length(which(x==TRUE & var==1)))
  # n.var.cases <- n.var.cases[-1]
  # n.controls <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x==FALSE)))
  # case.ctrl_ratio <- n.cases / n.controls
  
  # print(str(data[,1:10]))
  # print(str(n.var.cases))
  # print(str(n.var.controls))
  
  # filter phecodes to test
  # filter: number of cases  with variant (`var`) > 0
  # phecodes <- phecodes[-which(n.var.cases == 0)]
  # message('Number of phecodes: ', length(phecodes))
  
  # if streaming results to file, write a file with a header
  if (!is.null(file.conn)) {
    header = c(
      'phecode',
      'beta',
      'SE',
      'chisq',
      'OR',
      'p',
      'n.cases',
      'n.var.cases',
      'n.controls',
      'n.var.controls',
      'case.ctrl_ratio',
      'firth_correct'
    )
    readr::write_lines(
      x = paste(header, sep = '\t', collapse = '\t'),
      file = file.conn,
      append = TRUE
    )
  }
  
  # ## Create a bigmatrix object from `data`
  # data2 <- bigmemory::as.big.matrix(as.data.frame(data))
  # data_desc <- bigmemory::describe(data2)
  # rm(data)
  # doParallel::registerDoParallel(ncores)
  # ifelse(debug_log,
  #        cl = parallel::makeCluster(ncores, outfile = ''),
  #        cl = parallel::makeCluster(ncores))
  cl = parallel::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)
  
  x <- foreach(
    i = 1:length(phecodes),
    .errorhandling = 'pass',
    .noexport = c('sex.check'),
    # .export = c('both_sexes','phecodes','data_desc','back_correct_FirthSE','file.conn','var'),
    # .combine = 'rbind',
    .packages = c('splines', 'logistf', 'biganalytics')
  ) %dopar% {
    data_ <- bigmemory::attach.big.matrix(data_desc)
    
    phecode <- phecodes[i]
    if (is.null(both_sexes)) {
      f = paste0(
        '`',
        phecode,
        '`',
        ' ~ ',
        var,
        ' + ns(age, df=2) + sex + Center + ',
        paste0('PC', 1:40, separate = '', collapse =
                 '+')
      )
    } else if (!both_sexes[as.character(phecode)]) {
      f = paste0(
        '`',
        phecode,
        '`',
        ' ~ ',
        var,
        ' + ns(age,df=2) + Center + ',
        paste0('PC', 1:40, separate = '', collapse =
                 '+')
      )
    } else {
      f = paste0(
        '`',
        phecode,
        '`',
        ' ~ ',
        var,
        ' + ns(age, df=2) + sex + Center + ',
        paste0('PC', 1:40, separate = '', collapse =
                 '+')
      )
    }
    message(f)
    
    # Run GLM
    fit <-
      biganalytics::bigglm.big.matrix(
        formula(f),
        data_,
        family = gaussian(link = 'identity'),
        fc = ifelse(
          !is.null(both_sexes),
          yes = ifelse(!both_sexes[as.character(phecode)], # <--- DIFF
                       c('Center'), # <--- DIFF
                       c('sex', 'Center')),
          # <--- DIFF
          no = c('sex', 'Center')
        ),
        # <--- DIFF
        quiet = T
      )
    # Firth Correct GLM p-value if it is <= 0.05
    firth_correct = FALSE
    # if(summary(fit)$coef[2,4] <= 0.05){
    if (summary(fit)$mat[2, 'p'] <= 0.05) {
      message(phecode, ' Firth Correcting...')
      firth_correct = TRUE
      fit.firth <-
        logistf(
          formula(f),
          as.data.frame(data_[, c(phecode,
                                  'var',
                                  'age',
                                  'sex',
                                  'Center',
                                  paste0('PC', 1:40, separate = ''))]),
          control = logistf.control(maxit = 100, maxstep = 50),
          plcontrol = logistpl.control(maxit = 400, maxstep = 20),
          plconf = c(2),
          model = F
        )
    }
    
    # extract values from GLM or Firth model -> write to file, and/or return
    val <- rep(NA, 12)
    names(val) <- c(
      'phecode',
      'beta',
      'SE',
      'chisq',
      'OR',
      'p',
      'n.cases',
      'n.var.cases',
      'n.controls',
      'n.var.controls',
      'case.ctrl_ratio',
      'firth_correct'
    )
    val['phecode'] <- phecode
    if (firth_correct) {
      val['beta'] <- fit.firth$coefficients[2]
      val['OR'] <- exp(fit.firth$coefficients[2])
      val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
      val['chisq'] <- qchisq(1 - fit.firth$prob, df = 1)[2]
      val['p'] <- fit.firth$prob[2]
      if (back_correct_FirthSE) {
        val['SE'] <- # recompute SE using beta and Chisq value
          abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
      }
    } else {
      # val['beta'] <- coef(fit)[2]
      # val['OR'] <- exp(coef(fit)[2])
      # val['SE'] <- summary(fit)$coef[2, 2]
      # val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
      # val['p'] <- summary(fit)$coef[2, 4]
      val['beta'] <- summary(fit)$mat[2, 'Coef']
      val['OR'] <- exp(summary(fit)$mat[2, 'Coef'])
      val['SE'] <- summary(fit)$mat[2, 'SE']
      val['chisq'] <- qchisq(1 - summary(fit)$mat[2, 'p'], df =
                               1)
      val['p'] <- summary(fit)$mat[2, 'p']
    }
    # val['n.cases'] <- n.cases[phecode]
    # data <- data %>% select(var, all_of(phecode))
    # val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
    # val['n.controls'] <- n.controls[phecode]
    # val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
    # val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
    # val['firth_correct'] <- firth_correct
    # message('table var: ', unique(data_[,'var']))
    val['n.cases'] <- length(which(data_[, phecode] == 1)) # n.cases[phecode]
    # val['n.var.cases'] <- length(which(data_[,'var']==1 & data_[, phecode]==1)) # note: `var`==1 indicates presence of variant
    # val['n.var.cases'] <- length(which(data_[,'var']==2 & data_[, phecode]==1)) # note: `var`==2 indicates presence of variant
    # val['n.var.cases'] <- length(which(data_[,'var']!=0 & data_[, phecode]==1)) # note: `var`==2 indicates presence of variant <--- DIFF from process.glm()
    val['n.var.cases'] <- mean(data_$var[which(data_[, 'var'] !=
                                                 0 &
                                                 data_[, phecode] == 1)], na.rm = T) # <--- DIFF from process.glm()
    
    val['n.controls'] <- length(which(data_[, phecode] == 0))
    # val['n.var.controls'] <- length(which(data_[,'var']==1 & data_[,phecode]==0)) # note: `var`==1 indicates presence of variant
    # val['n.var.controls'] <- length(which(data_[,'var']==2 & data_[,phecode]==0)) # note: `var`==2 indicates presence of variant
    val['n.var.controls'] <- length(which(data_[, 'var'] !=
                                            0 &
                                            data_[, phecode] == 0)) # note: `var`==2 indicates presence of variant <--- DIFF from process.glm()
    
    val['case.ctrl_ratio'] <-
      (length(which(data_[, phecode] == 1))) / (length(which(data_[, phecode] ==
                                                               0)))
    val['firth_correct'] <- firth_correct
    
    # optionally, stream to file
    if (!is.null(file.conn)) {
      val_ <- paste(c(val), sep = '', collapse = '\t')
      readr::write_lines(x = val_,
                         file = file.conn,
                         append = TRUE)
    }
    return(val)
  }
  
  # doParallel::stopImplicitCluster()
  # stopCluster(cl)
  closeAllConnections()
  x <- bind_rows(x) %>% as_tibble()
  return(x)
}

# Function to run PheWAS on scored synthetic allele # NEWEST (as of 6/25/22)
process.glm2_rvscore <- function(phecodes,
                                 both_sexes = NULL,
                                 data_desc,
                                 ncores = 1,
                                 var = 'var',
                                 file.conn = NULL,
                                 back_correct_FirthSE = TRUE,
                                 var.is_factor = FALSE,
                                 debug_log = F) {
  require(dplyr)
  require(logistf)
  # require(doParallel)
  require(parallel)
  require(splines)
  require(logistf)
  require(bigmemory)
  require(biganalytics)
  require(snow)
  require(doSNOW)
  
  # pre-compute some info about dataset
  # message('Precomputing data: n.var in cases, controls, etc.')
  # sex.check <-
  #   apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) {
  #     ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
  #   })
  # var <- data %>% select(var)
  # n.cases <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x)))
  # n.var.cases <- apply(data %>% dplyr::select(var,all_of(phecodes)), 2, function(x) length(which(x==TRUE & var==1)))
  # n.var.cases <- n.var.cases[-1]
  # n.controls <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x==FALSE)))
  # case.ctrl_ratio <- n.cases / n.controls
  
  # print(str(data[,1:10]))
  # print(str(n.var.cases))
  # print(str(n.var.controls))
  
  # filter phecodes to test
  # filter: number of cases  with variant (`var`) > 0
  # phecodes <- phecodes[-which(n.var.cases == 0)]
  # message('Number of phecodes: ', length(phecodes))
  
  # if streaming results to file, write a file with a header
  if (!is.null(file.conn)) {
    header = c(
      'phecode',
      'beta',
      'SE',
      'chisq',
      'OR',
      'p',
      'n.cases',
      'n.var.cases',
      'n.controls',
      'n.var.controls',
      'case.ctrl_ratio',
      'firth_correct'
    )
    readr::write_lines(
      x = paste(header, sep = '\t', collapse = '\t'),
      file = file.conn,
      append = TRUE
    )
  }
  
  # ## Create a bigmatrix object from `data`
  # data2 <- bigmemory::as.big.matrix(as.data.frame(data))
  # data_desc <- bigmemory::describe(data2)
  # rm(data)
  # doParallel::registerDoParallel(ncores)
  # ifelse(debug_log,
  #        cl = parallel::makeCluster(ncores, outfile = ''),
  #        cl = parallel::makeCluster(ncores))
  cl = parallel::makeCluster(ncores, outfile = "")
  doSNOW::registerDoSNOW(cl)
  # registerDoParallel(cl)
  message('ncores: ', ncores)
  
  x <- foreach(
    i = 1:length(phecodes),
    .errorhandling = 'pass',
    .noexport = c('sex.check'),
    # .export = c('both_sexes','phecodes','data_desc','back_correct_FirthSE','file.conn','var'),
    # .combine = 'rbind',
    .packages = c('splines', 'logistf', 'biganalytics')
  ) %dopar% {
    data_ <- bigmemory::attach.big.matrix(data_desc)
    
    phecode <- phecodes[i]
    if (is.null(both_sexes)) {
      f = paste0(
        '`',
        phecode,
        '`',
        ' ~ ',
        var,
        ' + ns(age, df=2) + sex + Center + ',
        paste0('PC', 1:40, separate = '', collapse =
                 '+')
      )
    } else if (!both_sexes[as.character(phecode)]) {
      f = paste0(
        '`',
        phecode,
        '`',
        ' ~ ',
        var,
        ' + ns(age,df=2) + Center + ',
        paste0('PC', 1:40, separate = '', collapse =
                 '+')
      )
    } else {
      f = paste0(
        '`',
        phecode,
        '`',
        ' ~ ',
        var,
        ' + ns(age, df=2) + sex + Center + ',
        paste0('PC', 1:40, separate = '', collapse =
                 '+')
      )
    }
    message(f)
    # cat(f)
    
    # vector of factor variable names for bigglm.big.matrix
    var.is_factor = T
    get.fc <- function(both_sexes, phecode, var.is_factor) {
      fc = NULL
      if (!both_sexes[as.character(phecode)]) {
        fc = c('Center')
      } else {
        fc = c('sex', 'Center')
      }
      if (var.is_factor) {
        fc <- c(fc, 'var')
      }
      return(fc)
    }
    
    # Run GLM
    fit <-
      biganalytics::bigglm.big.matrix(
        formula(f),
        data_,
        # family=gaussian(link='identity'),
        family = binomial(link = 'logit'),
        # fc=ifelse(!is.null(both_sexes),
        #           yes = ifelse(!both_sexes[as.character(phecode)], # <--- DIFF
        #                        c('Center'),  # <--- DIFF
        #                        c('sex','Center')), # <--- DIFF
        #           no = c('sex','Center')), # <--- DIFF
        fc = get.fc(both_sexes, phecode, var.is_factor),
        quiet = F
      )
    # Firth Correct GLM p-value if it is <= 0.05
    firth_correct = FALSE
    # if(summary(fit)$coef[2,4] <= 0.05){
    if (summary(fit)$mat[2, 'p'] <= 0.0005) {
      message('FIRTH Correcting...', f)
      firth_correct = TRUE
      fit.firth <-
        logistf(
          formula(f),
          as.data.frame(data_[, c(phecode,
                                  'var',
                                  'age',
                                  'sex',
                                  'Center',
                                  paste0('PC', 1:40, separate = ''))]),
          control = logistf.control(maxit = 100, maxstep = 50),
          plcontrol = logistpl.control(maxit = 400, maxstep = 20),
          plconf = c(2),
          model = F
        )
    }
    
    # extract values from GLM or Firth model -> write to file, and/or return
    val <- rep(NA, 12)
    names(val) <- c(
      'phecode',
      'beta',
      'SE',
      'chisq',
      'OR',
      'p',
      'n.cases',
      'n.var.cases',
      'n.controls',
      'n.var.controls',
      'case.ctrl_ratio',
      'firth_correct'
    )
    val['phecode'] <- phecode
    if (firth_correct) {
      val['beta'] <- fit.firth$coefficients[2]
      val['OR'] <- exp(fit.firth$coefficients[2])
      val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
      val['chisq'] <- qchisq(1 - fit.firth$prob, df = 1)[2]
      val['p'] <- fit.firth$prob[2]
      if (back_correct_FirthSE) {
        val['SE'] <- # recompute SE using beta and Chisq value
          abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
      }
    } else {
      # val['beta'] <- coef(fit)[2]
      # val['OR'] <- exp(coef(fit)[2])
      # val['SE'] <- summary(fit)$coef[2, 2]
      # val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
      # val['p'] <- summary(fit)$coef[2, 4]
      val['beta'] <- summary(fit)$mat[2, 'Coef']
      val['OR'] <- exp(summary(fit)$mat[2, 'Coef'])
      val['SE'] <- summary(fit)$mat[2, 'SE']
      val['chisq'] <- qchisq(1 - summary(fit)$mat[2, 'p'], df =
                               1)
      val['p'] <- summary(fit)$mat[2, 'p']
    }
    # val['n.cases'] <- n.cases[phecode]
    # data <- data %>% select(var, all_of(phecode))
    # val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
    # val['n.controls'] <- n.controls[phecode]
    # val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
    # val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
    # val['firth_correct'] <- firth_correct
    # message('table var: ', unique(data_[,'var']))
    val['n.cases'] <- length(which(data_[, phecode] == 1)) # n.cases[phecode]
    # val['n.var.cases'] <- length(which(data_[,'var']==1 & data_[, phecode]==1)) # note: `var`==1 indicates presence of variant
    # val['n.var.cases'] <- length(which(data_[,'var']==2 & data_[, phecode]==1)) # note: `var`==2 indicates presence of variant
    # val['n.var.cases'] <- length(which(data_[,'var']!=0 & data_[, phecode]==1)) # note: `var`==2 indicates presence of variant <--- DIFF from process.glm()
    # val['n.var.cases'] <- mean(data_[,'var'][which(data_[,'var']!=0 & data_[, phecode]==1)], na.rm=T) # <--- DIFF from process.glm()
    
    val['n.controls'] <- length(which(data_[, phecode] == 0))
    # val['n.var.controls'] <- length(which(data_[,'var']==1 & data_[,phecode]==0)) # note: `var`==1 indicates presence of variant
    # val['n.var.controls'] <- length(which(data_[,'var']==2 & data_[,phecode]==0)) # note: `var`==2 indicates presence of variant
    # val['n.var.controls'] <- length(which(data_[,'var']!=0 & data_[,phecode]==0)) # note: `var`==2 indicates presence of variant <--- DIFF from process.glm()
    
    val['case.ctrl_ratio'] <-
      (length(which(data_[, phecode] == 1))) / (length(which(data_[, phecode] ==
                                                               0)))
    val['firth_correct'] <- firth_correct
    
    # optionally, stream to file
    if (!is.null(file.conn)) {
      val_ <- paste(c(val), sep = '', collapse = '\t')
      readr::write_lines(x = val_,
                         file = file.conn,
                         append = TRUE)
    }
    return(val)
  }
  
  # doParallel::stopImplicitCluster()
  stopCluster(cl)
  # closeAllConnections()
  x <- bind_rows(x) %>% as_tibble()
  return(x)
}

# Function to run PheWAS on a quantitative variable (used for mtscore_sum_resid) 
# Note: USED IN HET PAPER PHEWAS (Main + PathogenicVariant)
process.glm3 <- function(phecodes,
                         both_sexes = NULL,
                         data_,
                         ncores = 1,
                         var = 'var',
                         file.conn = NULL,
                         run.firth = TRUE,
                         back_correct_FirthSE = TRUE,
                         var.is_factor = FALSE,
                         missing_var = NULL,
                         debug_log = F) {
  require(dplyr)
  require(logistf)
  require(doParallel)
  require(parallel)
  require(splines)
  require(logistf)
  require(bigmemory)
  require(biganalytics)
  require(snow)
  require(doSNOW)
  
  # pre-compute some info about dataset
  # message('Precomputing data: n.var in cases, controls, etc.')
  # sex.check <-
  #   apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) {
  #     ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
  #   })
  # var <- data %>% select(var)
  # n.cases <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x)))
  # n.var.cases <- apply(data %>% dplyr::select(var,all_of(phecodes)), 2, function(x) length(which(x==TRUE & var==1)))
  # n.var.cases <- n.var.cases[-1]
  # n.controls <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x==FALSE)))
  # case.ctrl_ratio <- n.cases / n.controls
  
  # print(str(data[,1:10]))
  # print(str(n.var.cases))
  # print(str(n.var.controls))
  
  # filter phecodes to test
  # filter: number of cases  with variant (`var`) > 0
  # phecodes <- phecodes[-which(n.var.cases == 0)]
  # message('Number of phecodes: ', length(phecodes))
  
  # if streaming results to file, write a file with a header
  if (!is.null(file.conn)) {
    if (var.is_factor) {
      header = c(
        'phecode',
        'beta',
        'SE',
        'chisq',
        'OR',
        'p',
        'n.cases',
        'n.var.cases',
        'n.controls',
        'n.var.controls',
        'case.ctrl_ratio',
        'firth_correct'
      )
    } else {
      header = c(
        'phecode',
        'beta',
        'SE',
        'chisq',
        'OR',
        'p',
        'n.cases',
        'mean.var.cases',
        'n.controls',
        'mean.var.controls',
        'case.ctrl_ratio',
        'firth_correct'
      )
    }
    readr::write_lines(
      x = paste(header, sep = '\t', collapse = '\t'),
      file = file.conn,
      append = TRUE
    )
  }
  
  # cl = parallel::makeCluster(ncores)
  # doSNOW::registerDoSNOW(cl)
  
  # x <- foreach(i = 1:length(phecodes), .errorhandling = 'pass',
  #              .noexport = c('sex.check'),
  #              # .export = c('both_sexes','phecodes','data_desc','back_correct_FirthSE','file.conn','var'),
  #              # .combine = 'rbind',
  #              .packages = c('splines','logistf','biganalytics')) %dopar% {
  x <- list()
  # data_ <- bigmemory::attach.big.matrix(data_desc)
  for (i in 1:length(phecodes)) {
    phecode <- phecodes[i]
    
    # # Formula with het_count
    # if(is.null(both_sexes)){
    #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center + het_count')#,
    #            # paste0('PC',1:40,separate='',collapse='+'))
    # } else if(!both_sexes[as.character(phecode)]){
    #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age,df=4) + Center + het_count')#,
    #            # paste0('PC',1:40,separate='',collapse='+'))
    # } else {
    #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center + het_count')#,
    #            # paste0('PC',1:40,separate='',collapse='+'))
    # }
    
    # Formula without het_count
    if (is.null(both_sexes)) {
      f = paste0('`',
                 phecode,
                 '`',
                 ' ~ ',
                 var,
                 ' + ns(age, df=4) + sex + Center + SmokingStatus')#,
      # paste0('PC',1:40,separate='',collapse='+'))
    } else if (!both_sexes[as.character(phecode)]) {
      f = paste0('`',
                 phecode,
                 '`',
                 ' ~ ',
                 var,
                 ' + ns(age,df=4) + Center + SmokingStatus')#,
      # paste0('PC',1:40,separate='',collapse='+'))
    } else {
      f = paste0('`',
                 phecode,
                 '`',
                 ' ~ ',
                 var,
                 ' + ns(age, df=4) + sex + Center + SmokingStatus')#,
      # paste0('PC',1:40,separate='',collapse='+'))
    }
    message(f)
    
    # Run GLM
    dat = as.data.frame(data_[, c('id',
                                  phecode,
                                  'var',
                                  'age',
                                  'sex',
                                  'Center',
                                  'het_count',
                                  'SmokingStatus')])
    dat$sex <- as.factor(dat$sex)
    dat$Center <- as.factor(dat$Center)
    dat$het_count <- as.factor(dat$het_count)
    dat$SmokingStatus <- as.factor(dat$SmokingStatus)
    if (!is.null(missing_var)) {
      # set missing var to NA if IIDs given
      dat$var[which(dat$id %in% missing_var)] <- NA
    }
    if (var.is_factor) {
      dat$var <- as.factor(dat$var)
    }
    
    fit <-
      glm(
        formula = formula(f),
        data = dat,
        family = binomial(link = 'logit')
      )
    
    # print(summary(fit))
    
    # alternate GLM using biganalytics (doesn't properly remove missing var values, so not using)
    # fc <- c('sex','Center')
    # if(!is.null(both_sexes)) {
    #   if (!both_sexes[as.character(phecode)]) {
    #     fc <- fc[!fc %in% c('sex')]
    #   }
    # }
    # fit <-
    #   biganalytics::bigglm.big.matrix(
    #     formula(f),
    #     data_,
    #     family=gaussian(link='identity'),
    #     fc=fc,
    #     # fc=ifelse(!is.null(both_sexes),
    #     #           yes = ifelse(!both_sexes[as.character(phecode)],
    #     #                        yes=c('Center'),  # <--- both_sexes==TRUE, so !both_sexes==FALSE; thus only one sex present, so remove sex covar
    #     #                        no=c('sex','Center')), # <--- both_sexes==FALSE, so !both_sexes==TRUE, thus both present, so keep sex covar
    #     #           no = c('sex','Center')),
    #     quiet=T
    #   )
    
    # Firth Correct GLM p-value if it is <= 0.05
    firth_correct = FALSE
    if (run.firth) {
      # if(summary(fit)$coef[2,4] <= 0.05){
      if (summary(fit)$coef[2, 4] <= 5e-02) {
        # if(summary(fit)$mat[2,'p'] <= 0.05){
        message(phecode, ' Firth Correcting...')
        firth_correct = TRUE
        fit.firth <-
          logistf(
            formula(f),
            # as.data.frame(data_[, c(phecode, 'var','age','sex','Center',paste0('PC',1:40, separate=''))]),
            # as.data.frame(data_[, c(phecode, 'var','age','sex','Center')]),
            dat,
            control = logistf.control(maxit = 100, maxstep = 50),
            plcontrol = logistpl.control(maxit = 400, maxstep = 20),
            plconf = c(2),
            model = F
          )
      }
    }
    
    # extract values from GLM or Firth model -> write to file, and/or return
    val <- rep(NA, 12)
    if (var.is_factor) {
      names(val) <- c(
        'phecode',
        'beta',
        'SE',
        'chisq',
        'OR',
        'p',
        'n.cases',
        'n.var.cases',
        'n.controls',
        'n.var.controls',
        'case.ctrl_ratio',
        'firth_correct'
      )
    } else {
      names(val) <- c(
        'phecode',
        'beta',
        'SE',
        'chisq',
        'OR',
        'p',
        'n.cases',
        'mean.var.cases',
        'n.controls',
        'mean.var.controls',
        'case.ctrl_ratio',
        'firth_correct'
      )
    }
    val['phecode'] <- phecode
    if (firth_correct) {
      val['beta'] <- fit.firth$coefficients[2]
      val['OR'] <- exp(fit.firth$coefficients[2])
      val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
      val['chisq'] <- qchisq(1 - fit.firth$prob, df = 1)[2]
      val['p'] <- fit.firth$prob[2]
      if (back_correct_FirthSE) {
        val['SE'] <- # recompute SE using beta and Chisq value
          abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
      }
    } else {
      # val['beta'] <- coef(fit)[2]
      # val['OR'] <- exp(coef(fit)[2])
      # val['SE'] <- summary(fit)$coef[2, 2]
      # val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
      # val['p'] <- summary(fit)$coef[2, 4]
      # val['beta'] <- summary(fit)$mat[2,'Coef']
      val['beta'] <- summary(fit)$coef[2, 'Estimate']
      # val['OR'] <- exp(summary(fit)$mat[2,'Coef'])
      val['OR'] <- exp(summary(fit)$coef[2, 'Estimate'])
      # val['SE'] <- summary(fit)$mat[2,'SE']
      val['SE'] <- summary(fit)$coef[2, 'Std. Error']
      # val['chisq'] <- qchisq(1 - summary(fit)$mat[2,'p'], df=1)
      val['chisq'] <- qchisq(1 - summary(fit)$coef[2, 'Pr(>|z|)'], df =
                               1)
      # val['p'] <- summary(fit)$mat[2,'p']
      val['p'] <- summary(fit)$coef[2, 'Pr(>|z|)']
    }
    # val['n.cases'] <- n.cases[phecode]
    # data <- data %>% select(var, all_of(phecode))
    # val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
    # val['n.controls'] <- n.controls[phecode]
    # val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
    # val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
    # val['firth_correct'] <- firth_correct
    # message('table var: ', unique(data_[,'var']))
    val['n.cases'] <- length(which(data_[, phecode] == TRUE &
                                     !is.na(data_[, 'var'])))
    val['n.controls'] <- length(which(data_[, phecode] == FALSE &
                                        !is.na(data_[, 'var'])))
    if (var.is_factor) {
      val['n.var.cases'] <- length(which(data_[, phecode] == TRUE &
                                           data_[, 'var'] == 1))
      val['n.var.controls'] <- length(which(data_[, phecode] ==
                                              FALSE &
                                              data_[, 'var'] == 1))
    } else {
      var_ <- dat[, 'var']
      val['mean.var.cases'] <- mean(var_$var[which(!is.na(data_[, 'var']) &
                                                     data_[, phecode] == TRUE)], na.rm = T)
      val['mean.var.controls'] <- mean(var_$var[which(!is.na(data_[, 'var']) &
                                                        data_[, phecode] == FALSE)], na.rm = T)
    }
    
    val['case.ctrl_ratio'] <-
      (length(which(data_[, phecode] == TRUE &
                      !is.na(data_[, 'var'])))) / (length(which(data_[, phecode] == FALSE &
                                                                  !is.na(data_[, 'var']))))
    val['firth_correct'] <- firth_correct
    
    # optionally, stream to file
    if (!is.null(file.conn)) {
      val_ <- paste(c(val), sep = '', collapse = '\t')
      readr::write_lines(x = val_,
                         file = file.conn,
                         append = TRUE)
    }
    rm(dat)
    # return(val)
    x[[i]] <- val
  }
  
  # doParallel::stopImplicitCluster()
  # stopCluster(cl)
  # closeAllConnections()
  # x <- bind_rows(x) %>% as_tibble()
  return(x)
}

# Function to run PheWAS on a quantitative variable
process.glm3.1 <- function(phecodes,
                           both_sexes = NULL,
                           data_desc,
                           ncores = 1,
                           var = 'var',
                           file.conn = NULL,
                           back_correct_FirthSE = TRUE,
                           missing_var = NULL,
                           debug_log = F) {
  require(dplyr)
  require(logistf)
  require(doParallel)
  require(parallel)
  require(splines)
  require(logistf)
  require(bigmemory)
  require(biganalytics)
  require(snow)
  require(doSNOW)
  
  # pre-compute some info about dataset
  # message('Precomputing data: n.var in cases, controls, etc.')
  # sex.check <-
  #   apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) {
  #     ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
  #   })
  # var <- data %>% select(var)
  # n.cases <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x)))
  # n.var.cases <- apply(data %>% dplyr::select(var,all_of(phecodes)), 2, function(x) length(which(x==TRUE & var==1)))
  # n.var.cases <- n.var.cases[-1]
  # n.controls <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x==FALSE)))
  # case.ctrl_ratio <- n.cases / n.controls
  
  # print(str(data[,1:10]))
  # print(str(n.var.cases))
  # print(str(n.var.controls))
  
  # filter phecodes to test
  # filter: number of cases  with variant (`var`) > 0
  # phecodes <- phecodes[-which(n.var.cases == 0)]
  # message('Number of phecodes: ', length(phecodes))
  
  # if streaming results to file, write a file with a header
  if (!is.null(file.conn)) {
    header = c(
      'phecode',
      'beta',
      'SE',
      'chisq',
      'OR',
      'p',
      'n.cases',
      'mean.var.cases',
      'n.controls',
      'mean.var.controls',
      'case.ctrl_ratio',
      'firth_correct'
    )
    readr::write_lines(
      x = paste(header, sep = '\t', collapse = '\t'),
      file = file.conn,
      append = TRUE
    )
  }
  
  # cl = parallel::makeCluster(ncores)
  # doSNOW::registerDoSNOW(cl)
  
  # x <- foreach(i = 1:length(phecodes), .errorhandling = 'pass',
  #              .noexport = c('sex.check'),
  #              # .export = c('both_sexes','phecodes','data_desc','back_correct_FirthSE','file.conn','var'),
  #              # .combine = 'rbind',
  #              .packages = c('splines','logistf','biganalytics')) %dopar% {
  x <- list()
  data_ <- bigmemory::attach.big.matrix(data_desc)
  for (i in 1:length(phecodes)) {
    phecode <- phecodes[i]
    if (is.null(both_sexes)) {
      f = paste0('`', phecode, '`', ' ~ ', var, ' + ns(age, df=4) + sex + Center')#,
      # paste0('PC',1:40,separate='',collapse='+'))
    } else if (!both_sexes[as.character(phecode)]) {
      f = paste0('`', phecode, '`', ' ~ ', var, ' + ns(age,df=4) + Center')#,
      # paste0('PC',1:40,separate='',collapse='+'))
    } else {
      f = paste0('`', phecode, '`', ' ~ ', var, ' + ns(age, df=4) + sex + Center')#,
      # paste0('PC',1:40,separate='',collapse='+'))
    }
    message(f)
    
    # Run GLM
    dat = as.data.frame(data_[, c('id', phecode, 'var', 'age', 'sex', 'Center', 'het_count')])
    dat$sex <- as.factor(dat$sex)
    dat$Center <- as.factor(dat$Center)
    # dat$het_count <- as.factor(dat$het_count)
    if (!is.null(missing_var)) {
      # set missing var to NA if IIDs given
      dat$var[which(dat$id %in% missing_var)] <- NA
    }
    
    fit <-
      glm(
        formula = formula(f),
        data = dat,
        family = binomial(link = 'logit')
      )
    
    # alternate GLM using biganalytics (doesn't properly remove missing var values, so not using)
    # fc <- c('sex','Center')
    # if(!is.null(both_sexes)) {
    #   if (!both_sexes[as.character(phecode)]) {
    #     fc <- fc[!fc %in% c('sex')]
    #   }
    # }
    # fit <-
    #   biganalytics::bigglm.big.matrix(
    #     formula(f),
    #     data_,
    #     family=gaussian(link='identity'),
    #     fc=fc,
    #     # fc=ifelse(!is.null(both_sexes),
    #     #           yes = ifelse(!both_sexes[as.character(phecode)],
    #     #                        yes=c('Center'),  # <--- both_sexes==TRUE, so !both_sexes==FALSE; thus only one sex present, so remove sex covar
    #     #                        no=c('sex','Center')), # <--- both_sexes==FALSE, so !both_sexes==TRUE, thus both present, so keep sex covar
    #     #           no = c('sex','Center')),
    #     quiet=T
    #   )
    
    # Firth Correct GLM p-value if it is <= 0.05
    firth_correct = FALSE
    if (FALSE) {
      # if(summary(fit)$coef[2,4] <= 0.05){
      if (summary(fit)$coef[2, 4] <= 0.05) {
        # if(summary(fit)$mat[2,'p'] <= 0.05){
        message(phecode, ' Firth Correcting...')
        firth_correct = TRUE
        fit.firth <-
          logistf(
            formula(f),
            # as.data.frame(data_[, c(phecode, 'var','age','sex','Center',paste0('PC',1:40, separate=''))]),
            # as.data.frame(data_[, c(phecode, 'var','age','sex','Center')]),
            dat,
            control = logistf.control(maxit = 100, maxstep = 50),
            plcontrol = logistpl.control(maxit = 400, maxstep = 20),
            plconf = c(2),
            model = F
          )
      }
    }
    
    # extract values from GLM or Firth model -> write to file, and/or return
    val <- rep(NA, 12)
    names(val) <- c(
      'phecode',
      'beta',
      'SE',
      'chisq',
      'OR',
      'p',
      'n.cases',
      'mean.var.cases',
      'n.controls',
      'mean.var.controls',
      'case.ctrl_ratio',
      'firth_correct'
    )
    val['phecode'] <- phecode
    if (firth_correct) {
      val['beta'] <- fit.firth$coefficients[2]
      val['OR'] <- exp(fit.firth$coefficients[2])
      val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
      val['chisq'] <- qchisq(1 - fit.firth$prob, df = 1)[2]
      val['p'] <- fit.firth$prob[2]
      if (back_correct_FirthSE) {
        val['SE'] <- # recompute SE using beta and Chisq value
          abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
      }
    } else {
      # val['beta'] <- coef(fit)[2]
      # val['OR'] <- exp(coef(fit)[2])
      # val['SE'] <- summary(fit)$coef[2, 2]
      # val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
      # val['p'] <- summary(fit)$coef[2, 4]
      # val['beta'] <- summary(fit)$mat[2,'Coef']
      val['beta'] <- summary(fit)$coef[2, 'Estimate']
      # val['OR'] <- exp(summary(fit)$mat[2,'Coef'])
      val['OR'] <- exp(summary(fit)$coef[2, 'Estimate'])
      # val['SE'] <- summary(fit)$mat[2,'SE']
      val['SE'] <- summary(fit)$coef[2, 'Std. Error']
      # val['chisq'] <- qchisq(1 - summary(fit)$mat[2,'p'], df=1)
      val['chisq'] <- qchisq(1 - summary(fit)$coef[2, 'Pr(>|z|)'], df =
                               1)
      # val['p'] <- summary(fit)$mat[2,'p']
      val['p'] <- summary(fit)$coef[2, 'Pr(>|z|)']
    }
    # val['n.cases'] <- n.cases[phecode]
    # data <- data %>% select(var, all_of(phecode))
    # val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
    # val['n.controls'] <- n.controls[phecode]
    # val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
    # val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
    # val['firth_correct'] <- firth_correct
    # message('table var: ', unique(data_[,'var']))
    var_ <- data[, 'var']
    val['n.cases'] <- length(which(data_[, phecode] == 1 &
                                     !is.na(data_[, 'var'])))
    val['mean.var.cases'] <- mean(var_$var[which(!is.na(data_[, 'var']) &
                                                   data_[, phecode] == 1)], na.rm = T)
    val['n.controls'] <- length(which(data_[, phecode] == 0 &
                                        !is.na(data_[, 'var'])))
    val['mean.var.controls'] <- mean(var_$var[which(!is.na(data_[, 'var']) &
                                                      data_[, phecode] == 0)], na.rm = T)
    
    val['case.ctrl_ratio'] <-
      (length(which(data_[, phecode] == 1 &
                      !is.na(data_[, 'var'])))) / (length(which(data_[, phecode] == 0 &
                                                                  !is.na(data_[, 'var']))))
    val['firth_correct'] <- firth_correct
    
    # optionally, stream to file
    if (!is.null(file.conn)) {
      val_ <- paste(c(val), sep = '', collapse = '\t')
      readr::write_lines(x = val_,
                         file = file.conn,
                         append = TRUE)
    }
    rm(dat)
    # return(val)
    x[[i]] <- val
  }
  
  # doParallel::stopImplicitCluster()
  # stopCluster(cl)
  # closeAllConnections()
  # x <- bind_rows(x) %>% as_tibble()
  return(x)
}

# Function to run PheWAS on a quantitative variable
process.glm3.2 <- function(phecodes,
                           both_sexes = NULL,
                           data_,
                           ncores = 1,
                           var = 'var',
                           file.conn = NULL,
                           back_correct_FirthSE = TRUE,
                           missing_var = NULL,
                           debug_log = F) {
  require(dplyr)
  require(logistf)
  require(doParallel)
  require(parallel)
  require(splines)
  require(logistf)
  require(bigmemory)
  require(biganalytics)
  require(snow)
  require(doSNOW)
  
  # pre-compute some info about dataset
  # message('Precomputing data: n.var in cases, controls, etc.')
  # sex.check <-
  #   apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) {
  #     ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
  #   })
  # var <- data %>% select(var)
  # n.cases <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x)))
  # n.var.cases <- apply(data %>% dplyr::select(var,all_of(phecodes)), 2, function(x) length(which(x==TRUE & var==1)))
  # n.var.cases <- n.var.cases[-1]
  # n.controls <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x==FALSE)))
  # case.ctrl_ratio <- n.cases / n.controls
  
  # print(str(data[,1:10]))
  # print(str(n.var.cases))
  # print(str(n.var.controls))
  
  # filter phecodes to test
  # filter: number of cases  with variant (`var`) > 0
  # phecodes <- phecodes[-which(n.var.cases == 0)]
  # message('Number of phecodes: ', length(phecodes))
  
  # if streaming results to file, write a file with a header
  if (!is.null(file.conn)) {
    header = c(
      'phecode',
      'beta',
      'SE',
      'chisq',
      'OR',
      'p',
      'n.cases',
      'mean.var.cases',
      'n.controls',
      'mean.var.controls',
      'case.ctrl_ratio',
      'firth_correct'
    )
    readr::write_lines(
      x = paste(header, sep = '\t', collapse = '\t'),
      file = file.conn,
      append = TRUE
    )
  }
  
  # cl = parallel::makeCluster(ncores)
  # doSNOW::registerDoSNOW(cl)
  
  # x <- foreach(i = 1:length(phecodes), .errorhandling = 'pass',
  #              .noexport = c('sex.check'),
  #              # .export = c('both_sexes','phecodes','data_desc','back_correct_FirthSE','file.conn','var'),
  #              # .combine = 'rbind',
  #              .packages = c('splines','logistf','biganalytics')) %dopar% {
  x <- list()
  # data_ <- bigmemory::attach.big.matrix(data_desc)
  for (i in 1:length(phecodes)) {
    phecode <- phecodes[i]
    # if(is.null(both_sexes)){
    #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center')#,
    #   # paste0('PC',1:40,separate='',collapse='+'))
    # } else if(!both_sexes[as.character(phecode)]){
    #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age,df=4) + Center')#,
    #   # paste0('PC',1:40,separate='',collapse='+'))
    # } else {
    #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center')#,
    #   # paste0('PC',1:40,separate='',collapse='+'))
    # }
    if (is.null(both_sexes)) {
      f = paste0(
        '`',
        phecode,
        '`',
        ' ~ ',
        var,
        ' + ns(age, df=2) + sex + Center + ',
        paste0('PC', 1:40, separate = '', collapse = '+')
      )
    } else if (!both_sexes[as.character(phecode)]) {
      f = paste0(
        '`',
        phecode,
        '`',
        ' ~ ',
        var,
        ' + ns(age,df=2) + Center + ',
        paste0('PC', 1:40, separate = '', collapse = '+')
      )
    } else {
      f = paste0(
        '`',
        phecode,
        '`',
        ' ~ ',
        var,
        ' + ns(age, df=2) + sex + Center + ',
        paste0('PC', 1:40, separate = '', collapse = '+')
      )
    }
    message(f)
    
    # Run GLM
    # dat=as.data.frame(data_[, c('id', phecode, 'var','age','sex','Center','het_count')])
    dat = as.data.frame(data_[, c('id',
                                  phecode,
                                  'var',
                                  'age',
                                  'sex',
                                  'Center',
                                  paste0('PC', 1:40))])
    dat$sex <- as.factor(dat$sex)
    dat$Center <- as.factor(dat$Center)
    # dat$het_count <- as.factor(dat$het_count)
    if (!is.null(missing_var)) {
      # set missing var to NA if IIDs given
      dat$var[which(dat$id %in% missing_var)] <- NA
    }
    
    fit <-
      glm(
        formula = formula(f),
        data = dat,
        family = binomial(link = 'logit')
      )
    
    # alternate GLM using biganalytics (doesn't properly remove missing var values, so not using)
    # fc <- c('sex','Center')
    # if(!is.null(both_sexes)) {
    #   if (!both_sexes[as.character(phecode)]) {
    #     fc <- fc[!fc %in% c('sex')]
    #   }
    # }
    # fit <-
    #   biganalytics::bigglm.big.matrix(
    #     formula(f),
    #     data_,
    #     family=gaussian(link='identity'),
    #     fc=fc,
    #     # fc=ifelse(!is.null(both_sexes),
    #     #           yes = ifelse(!both_sexes[as.character(phecode)],
    #     #                        yes=c('Center'),  # <--- both_sexes==TRUE, so !both_sexes==FALSE; thus only one sex present, so remove sex covar
    #     #                        no=c('sex','Center')), # <--- both_sexes==FALSE, so !both_sexes==TRUE, thus both present, so keep sex covar
    #     #           no = c('sex','Center')),
    #     quiet=T
    #   )
    
    # Firth Correct GLM p-value if it is <= 0.05
    firth_correct = FALSE
    if (FALSE) {
      # if(summary(fit)$coef[2,4] <= 0.05){
      if (summary(fit)$coef[2, 4] <= 0.05) {
        # if(summary(fit)$mat[2,'p'] <= 0.05){
        message(phecode, ' Firth Correcting...')
        firth_correct = TRUE
        fit.firth <-
          logistf(
            formula(f),
            # as.data.frame(data_[, c(phecode, 'var','age','sex','Center',paste0('PC',1:40, separate=''))]),
            # as.data.frame(data_[, c(phecode, 'var','age','sex','Center')]),
            dat,
            control = logistf.control(maxit = 100, maxstep = 50),
            plcontrol = logistpl.control(maxit = 400, maxstep = 20),
            plconf = c(2),
            model = F
          )
      }
    }
    
    # extract values from GLM or Firth model -> write to file, and/or return
    val <- rep(NA, 12)
    names(val) <- c(
      'phecode',
      'beta',
      'SE',
      'chisq',
      'OR',
      'p',
      'n.cases',
      'mean.var.cases',
      'n.controls',
      'mean.var.controls',
      'case.ctrl_ratio',
      'firth_correct'
    )
    val['phecode'] <- phecode
    if (firth_correct) {
      val['beta'] <- fit.firth$coefficients[2]
      val['OR'] <- exp(fit.firth$coefficients[2])
      val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
      val['chisq'] <- qchisq(1 - fit.firth$prob, df = 1)[2]
      val['p'] <- fit.firth$prob[2]
      if (back_correct_FirthSE) {
        val['SE'] <- # recompute SE using beta and Chisq value
          abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
      }
    } else {
      # val['beta'] <- coef(fit)[2]
      # val['OR'] <- exp(coef(fit)[2])
      # val['SE'] <- summary(fit)$coef[2, 2]
      # val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
      # val['p'] <- summary(fit)$coef[2, 4]
      # val['beta'] <- summary(fit)$mat[2,'Coef']
      val['beta'] <- summary(fit)$coef[2, 'Estimate']
      # val['OR'] <- exp(summary(fit)$mat[2,'Coef'])
      val['OR'] <- exp(summary(fit)$coef[2, 'Estimate'])
      # val['SE'] <- summary(fit)$mat[2,'SE']
      val['SE'] <- summary(fit)$coef[2, 'Std. Error']
      # val['chisq'] <- qchisq(1 - summary(fit)$mat[2,'p'], df=1)
      val['chisq'] <- qchisq(1 - summary(fit)$coef[2, 'Pr(>|z|)'], df =
                               1)
      # val['p'] <- summary(fit)$mat[2,'p']
      val['p'] <- summary(fit)$coef[2, 'Pr(>|z|)']
    }
    # val['n.cases'] <- n.cases[phecode]
    # data <- data %>% select(var, all_of(phecode))
    # val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
    # val['n.controls'] <- n.controls[phecode]
    # val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
    # val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
    # val['firth_correct'] <- firth_correct
    # message('table var: ', unique(data_[,'var']))
    var_ <- data[, 'var']
    val['n.cases'] <- length(which(data_[, phecode] == 1 &
                                     !is.na(data_[, 'var'])))
    val['mean.var.cases'] <- mean(var_$var[which(!is.na(data_[, 'var']) &
                                                   data_[, phecode] == 1)], na.rm = T)
    val['n.controls'] <- length(which(data_[, phecode] == 0 &
                                        !is.na(data_[, 'var'])))
    val['mean.var.controls'] <- mean(var_$var[which(!is.na(data_[, 'var']) &
                                                      data_[, phecode] == 0)], na.rm = T)
    
    val['case.ctrl_ratio'] <-
      (length(which(data_[, phecode] == 1 &
                      !is.na(data_[, 'var'])))) / (length(which(data_[, phecode] == 0 &
                                                                  !is.na(data_[, 'var']))))
    val['firth_correct'] <- firth_correct
    
    # optionally, stream to file
    if (!is.null(file.conn)) {
      val_ <- paste(c(val), sep = '', collapse = '\t')
      readr::write_lines(x = val_,
                         file = file.conn,
                         append = TRUE)
    }
    rm(dat)
    # return(val)
    x[[i]] <- val
  }
  
  # doParallel::stopImplicitCluster()
  # stopCluster(cl)
  # closeAllConnections()
  # x <- bind_rows(x) %>% as_tibble()
  return(x)
}

# See `process-firth.R` for alternate version
process.firth <- function(s, ncores=1, gene=NULL) {
  require(logistf)
  require(doParallel)
  require(parallel)
  require(doParallel)
  
  doParallel::registerDoParallel(ncores)
  
  x <- foreach(i = 1:nrow(s), .packages='logistf', .errorhandling = 'pass') %dopar% {
    phecode <- s$phenotype[i]
    if(s$note[i] ==" [Note: Column(s) dropped due to lack of variability:  sex ]"){
      formula=paste0('`',phecode,'`',' ~ ', gene, '.carrier + age + wes.batch + ', paste0('PC',1:40,separate='',collapse='+'))
    } else {
      formula=paste0('`',phecode,'`',' ~ ', gene, '.carrier + age + sex + wes.batch + ', paste0('PC',1:40,separate='',collapse='+'))
    }
    message(formula)
    fit.firth <-
      logistf(
        formula,
        data,
        control = logistf.control(maxit = 100, maxstep = 50),
        plcontrol = logistpl.control(maxit = 400, maxstep = 20),
        plconf = c(2),
        model=F
      )
    
    # extract values from firth model -> return
    val <- rep(NA, 3)
    names(val) <- c('beta','SE','p')
    val['beta'] <- fit.firth$coefficients[2]
    val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
    val['p'] <- fit.firth$prob[2]
  }
  
  doParallel::stopImplicitCluster()
  return(x)
}




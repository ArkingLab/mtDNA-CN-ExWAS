

# ---- PHECODE GROUPS -----
n500k.phecodes <- getobj('~/mito_rare-variant/resources/n502485.icd10_phecodes.rds')
phecode_def <- read_csv('~/mito_rare-variant/resources/phecode_definitions1.2.csv')

testx <- function(x){
  tmp <- 0
  if (all(is.na(x))) {
    tmp = NA
  } else {
    if (any(x,na.rm=T)) {
      tmp <- TRUE
    } else {
      tmp <- FALSE
    }
  }
  return(tmp)
}

cases.by.groups = phecode_def %>%
  # filter(category_number == 1) %>%
  group_by(category) %>%
  group_map(
    .f = function(.x, ...) {
      codes = .x$phecode
      tmp = n500k.phecodes %>% dplyr::select(any_of(codes))
      case = apply(tmp, 1, testx)
      return(case)
    }
  ) %>% 
  setNames(unique(sort(phecode_def$category)))

cases.by.groups <- 
  as.matrix(simplify2array(cases.by.groups))
rownames(cases.by.groups) <- n500k.phecodes$id

saveRDS(object = cases.by.groups, file = '~/mito_rare-variant/resources/n502485.phecode_cases_by_phecode_category.rds')
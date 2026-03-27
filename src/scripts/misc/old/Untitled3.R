f0=formula(paste0('`',phecode,'`',' ~ ','ns(age, df=2) + sex + genotyping.array + Center + MCV + PlatCrit + MSCV + creatinine + Aspartate_aminotransferase + MPV + Platelet + alanine_aminotransferase +', paste0('PC',1:40,separate='',collapse='+')))
f1=formula(paste0('`',phecode, '` ~ haplo + ns(age, df=2) + sex + genotyping.array + Center + MCV + PlatCrit + MSCV + creatinine + Aspartate_aminotransferase + MPV + Platelet + alanine_aminotransferase +', paste0('PC',1:40,separate='',collapse="+")))
data <- ukb.data %>% dplyr::select(all_of(phecode), haplo, age, sex, genotyping.array, Center, 
                                   MCV, PlatCrit, MSCV,creatinine,Aspartate_aminotransferase, MPV, Platelet, alanine_aminotransferase,
                                   paste0('PC',1:40,separate='')) %>% filter(complete.cases(.))
f0.glm_alt <- glm(f0, data, family='binomial'); f1.glm_alt <- glm(f1, data, family='binomial')


# MT73A_G
phecode='411'
var='MT7028C_T'
f0=formula(paste0('`',phecode,'`',' ~ ','ns(age, df=2) + sex + genotyping.array + Center +', paste0('PC',1:40,separate='',collapse='+')))
f1=formula(paste0('`',phecode, '` ~ ',var,' + ns(age, df=2) + sex + genotyping.array + Center +', paste0('PC',1:40,separate='',collapse="+")))
data <- ukb.data %>% dplyr::select(all_of(phecode), haplo, age, sex, genotyping.array, Center, 
                                   MCV, PlatCrit, MSCV,creatinine,Aspartate_aminotransferase, MPV, Platelet, alanine_aminotransferase,
                                   haplo, MT73A_G, MT7028C_T, MT10238T_C, MT12612A_G, MT13617T_C, MT15257G_A,
                                   paste0('PC',1:40,separate='')) %>% filter(complete.cases(.))
f0.glm_alt <- glm(f0, data, family='binomial'); f1.glm_alt <- glm(f1, data, family='binomial')



process.mtSNP <- function(data, phecode_qual, sex_check, ncores=1, var) {
  require(logistf)
  require(doParallel)
  require(parallel)
  require(doParallel)
  require(splines)
  
  doParallel::registerDoParallel(ncores)
  
  x <- foreach(i = 1:length(phecode_qual), .errorhandling = 'pass') %dopar% {
    phecode <- phecode_qual[i]
    
    if(!sex_check[as.character(phecode)]){
      f0=formula(paste0('`',phecode,'`',' ~ ','ns(age, df=2) + genotyping.array + Center +', paste0('PC',1:40,separate='',collapse='+')))
      f1=formula(paste0('`',phecode, '` ~ ',var,' + ns(age, df=2) + genotyping.array + Center +', paste0('PC',1:40,separate='',collapse="+")))
    } else {
      f0=formula(paste0('`',phecode,'`',' ~ ','ns(age, df=2) + sex + genotyping.array + Center +', paste0('PC',1:40,separate='',collapse='+')))
      f1=formula(paste0('`',phecode, '` ~ ',var,' + ns(age, df=2) + sex + genotyping.array + Center +', paste0('PC',1:40,separate='',collapse="+")))
    }
    data <- ukb.data %>% dplyr::select(all_of(phecode), all_of(var), age, sex, genotyping.array, Center,
                                       paste0('PC',1:40,separate='')) %>% filter(complete.cases(.))
    data$Center <- as.factor(data$Center)
    message(f1)
    
    f0.glm <- glm(formula = f0, family = "binomial", data = data)
    f1.glm <- glm(formula = f1, family = "binomial", data = data)
    a <- anova(f0.glm, f1.glm, test="LRT")
    
    # extract values from glm model and anova -> return
    val <- rep(NA, 3)
    names(val) <- c('phecode','Pr(>Chi)',var)
    val[1] <- as.numeric(phecode)
    val[2] <- a$`Pr(>Chi)`[2]
    val[3] <- summary(f1.glm)$coef[2,4]
    return(val)
  }
  
  doParallel::stopImplicitCluster()
  return(x)
}










# ---- PHECODE GROUPS -----
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

cases.by.groups=phecodes %>% group_by(group) %>% # filter(group=='infectious diseases') %>% 
  group_map(.f=function(.x,...){
  codes=.x$phecode
  tmp=ukb.data %>% dplyr::select(all_of(codes))
  case = apply(tmp, 1, testx)
  return(case)
  })






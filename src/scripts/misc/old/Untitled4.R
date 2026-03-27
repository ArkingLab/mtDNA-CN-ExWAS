# data for each phecode to pass into process.glm2
t <- data.unrel_whites %>% dplyr::select(all_of(phecodes)) %>% as.list()
x <- list()
for(i in 1:length(t)){
  x[[i]] <- list(phecode=phecodes[i],
               t=t[[i]])
}
names(x) <- phecodes

# covariate data for GLM/logistf
cov <- 
  data.unrel_whites %>%
  dplyr::select(
    id,
    age,
    sex,
    Center,
    var,
    genotyping.array,
    in.white.British.ancestry.subset,
    used.in.pca.calculation,
    all_of(paste0('PC', 1:40, separate = '')
    ))

# If streaming results to file, open conn, write header
# file.conn=NULL
file.conn='./phewas.samhd1_mac_gteq6_in450k.tsv'
if(!is.null(file.conn)){
  header=c('phecode','beta','SE','chisq','OR','p','n.cases','n.var.cases','n.controls','n.var.controls','case.ctrl_ratio','firth_correct')
  readr::write_lines(x = paste(header, sep='\t', collapse='\t'), file = file.conn, append = TRUE)
}

# run phewas using parallel::mclapply()
phewas <-
  mclapply(sample(x,20),
           process.glm2,
           cov = cov,
           sex_check = sex.check,
           gene='SAMHD1', 
           var='var',
           n.cases=n.cases, 
           n.controls=n.controls, 
           file.conn=file.conn,
           mc.cores = 10)
phewas <- simplify2array(lm.res, higher=F)

# ----
tmp_run <- function(phecodes, cov, sex_check, ncores=1, gene=NULL, var=NULL, 
                    n.cases=NULL, n.controls=NULL, file.conn=NULL,
                    back_correct_FirthSE=TRUE){
  phecode <- phecodes[[1]]
  t <- phecodes[[2]]
  
  data <- cbind(as_tibble(t), cov)
  colnames(data) <- c(phecode, colnames(cov))
  print(colnames(data))
  #message(phecode)
  #print(str(t))
  #tmp <- t[[1]]
  #return(tmp)
}

process.glm2 <- function(phecodes, cov, sex_check, ncores=1, gene=NULL, var=NULL, 
                        n.cases=NULL, n.controls=NULL, file.conn=NULL,
                        back_correct_FirthSE=TRUE) {
  require(splines)
  require(logistf)

  phecode <- phecodes[[1]]
  t <- phecodes[[2]]
  
  if(!sex_check[as.character(phecode)]){
    formula=paste0('`',phecode,'`',' ~ ', var, ' + ns(age,df=2) + Center + ', paste0('PC',1:40,separate='',collapse='+'))
  } else {
    formula=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=2) + sex + Center + ', paste0('PC',1:40,separate='',collapse='+'))
  }
  message(formula)
  
  # form data
  data <- cbind(as_tibble(t), cov)
  colnames(data) <- c(phecode, colnames(cov))
  # print(colnames(data))
  
  # GLM
  fit <-
    glm(
      formula,
      data,
      family='binomial'
    )
    # Optionally, Firth Correct P-value if GLM p-value <= 0.05
    firth_correct=FALSE
    if(summary(fit)$coef[2,4] <= 0.05){
      message(phecode, ' Firth Correcting...')
      firth_correct=TRUE
      fit.firth <-
        logistf(
          formula,
          data,
          control = logistf.control(maxit = 100, maxstep = 50),
          plcontrol = logistpl.control(maxit = 400, maxstep = 20),
          plconf = c(2),
          model=F
        )
    }
    
    # extract values from GLM or Firth model -> write to file, and/or return
    val <- rep(NA, 12)
    names(val) <- c('phecode','beta','SE','chisq','OR','p','n.cases','n.var.cases','n.controls','n.var.controls','case.ctrl_ratio','firth_correct')
    val['phecode'] <- phecode
    if(firth_correct) {
      val['beta'] <- fit.firth$coefficients[2]
      val['OR'] <- exp(fit.firth$coefficients[2])
      val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
      val['chisq'] <- qchisq(1-fit.firth$prob,df=1)[2]
      val['p'] <- fit.firth$prob[2]
      if(back_correct_FirthSE){
        val['SE'] <- # recompute SE using beta and Chisq value
          abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
      }
    } else {
      val['beta'] <- coef(fit)[2]
      val['OR'] <- exp(coef(fit)[2])
      val['SE'] <- summary(fit)$coef[2, 2]
      val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
      val['p'] <- summary(fit)$coef[2, 4]
    }
    val['n.cases'] <- n.cases[phecode]
    data <- data %>% select(var, all_of(phecode))
    val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
    val['n.controls'] <- n.controls[phecode]
    val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
    val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
    val['firth_correct'] <- firth_correct
    
    # optionally, stream to file
    if(!is.null(file.conn)){
      val_ <- paste(c(val), sep='', collapse='\t')
      readr::write_lines(x = val_, file = file.conn, append = TRUE)
    }
    return(val)
  }
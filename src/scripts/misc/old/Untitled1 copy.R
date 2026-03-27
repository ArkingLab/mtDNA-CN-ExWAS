
# read in results from results directory
# output from assocTestAggregate
getobj <- function (Rdata)
{
  require(tools)
  if (tolower(file_ext(Rdata)) == "rds") {
    return(readRDS(Rdata))
  }
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata, 
                  "\nReturning only the first object"))
  }
  return(get(objname))
}


res.fnames <- list.files(path='~/mtrv/resources/res/', pattern='*.Rda')

res.list <- list()

for(i in 1:length(res.fnames)){
  message(i)
  res.list[[i]] <- 
    getobj(res.fnames[i])[[1]]
}
#Penalized linear regression
library(gelnet)

#Inputs: design matrix, boolean for using mDNAsi, boolean for mRNAsi, boolean for raw genetic data, l1 penalty, l2 penalty. 
#Output: fitted model

getPenLinReg <- function(dat, mDNAsi=F, mRNAsi=F, raw=T, l1=1, l2=1, max.iter = 100, eps = 1e-05){
  #Extract outcome
  daysToDeath <- dat[,"Y"]
  
  #Extract design matrix
  design <- dat[,c("gender", "age", "cancer.class")]
  
  if(mDNAsi){
    design <- cbind(design, dat[,"mDNAsi"])
    colnames(design)[length(colnames(design))] <- "mDNAsi"
  }
  
  if(mDNAsi){
    design <- cbind(design, dat[,"mRNAsi"])
    colnames(design)[length(colnames(design))] <- "mRNAsi"
  }
  
  if(raw){
    rawColNames <- colnames(dat)[grepl("cg", names(dat))]
    design <- cbind(design, dat[,rawColNames])
  }
  
  
  #sapply(design, class)
  
  design <- as.matrix(design)
  
  #https://cran.r-project.org/web/packages/gelnet/gelnet.pdf
  elasticFit <- gelnet(X=design, y=daysToDeath, l1=l1, l2=l2, max.iter=max.iter, eps=eps)
  
  return(elasticFit)
}

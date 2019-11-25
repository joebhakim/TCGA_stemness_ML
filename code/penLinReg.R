#Penalized linear regression
library(gelnet)
#library(glmnet)

#Inputs: design matrix, boolean for using mDNAsi, boolean for mRNAsi, boolean for raw genetic data, l1 penalty, l2 penalty. 
#Output: fitted model

getPenLinReg <- function(design, mDNAsi=F, mRNAsi=F, raw=T, l1=1, l2=1, max.iter = 100, eps = 1e-05){
  #Remove rows with missing values
  #design <- na.omit(design)
  
  #Extract outcome
  daysToDeath <- design[,"Y"]
  
  #Extract design matrix
  cancerCols <- colnames(design)[grepl("cancer", colnames(design))]
  designReduced <- design[,c("genderMALE", "age", cancerCols)]
  
  if(mDNAsi){
    designReduced  <- cbind(designReduced , design[,"mDNAsi"])
    colnames(designReduced)[length(colnames(designReduced))] <- "mDNAsi"
  }
  
  if(raw){
    rawColNames <- colnames(design)[grepl("cg", colnames(design))]
    designReduced  <- cbind(designReduced, design[,rawColNames])
  }
  
  #https://cran.r-project.org/web/packages/gelnet/gelnet.pdf
  elasticFit <- gelnet(X=designReduced, y=daysToDeath, l1=l1, l2=l2, max.iter=max.iter, eps=eps)
  
  #elasticFit <- glmnet(x=designReduced, y=daysToDeath, alpha=0.5)
  
  return(elasticFit)
}

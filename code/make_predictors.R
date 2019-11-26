#Penalized linear regression
library(gelnet)
library(randomForest)

load('data/design-matrix.rda')
design <- as_tibble(design)
#design <- filter(design, censor.indicator==0)

set.seed(1) # Set Seed so that same sample can be reproduced in future also
spec = c(train = .6, test = .2, validate = .2)
g = sample(cut(
  seq(nrow(design)), 
  nrow(design)*cumsum(c(0,spec)),
  labels = names(spec)
))

design_split <- split(design, g)

design_train <- as.matrix(design_split$train)
design_validate <- as.matrix(design_split$validate)
design_test <- as.matrix(design_split$test)

getPenLinReg <- function(design, mDNAsi=T, raw=F, l1=1, l2=1, max.iter = 100, eps = 1e-05){
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
  
  #sapply(design, class)
  
  #design <- model.matrix( ~ ., design)
  
  #design <- as.matrix(design)
  
  #https://cran.r-project.org/web/packages/gelnet/gelnet.pdf
  elasticFit <- gelnet(X=designReduced, y=daysToDeath, l1=l1, l2=l2, max.iter=max.iter, eps=eps)
  
  #elasticFit <- glmnet(x=designReduced, y=daysToDeath, alpha=0.5)
  
  return(elasticFit)
}


getRFModel <- function(design_in, mDNAsi=T, raw=F, l1=1, l2=1, max.iter = 100, eps = 1e-05){
  #Remove rows with missing values
  #design <- na.omit(design)
  
  #Extract design matrix
  cancerCols <- colnames(design_in)[grepl("cancer", colnames(design_in))]
  designReduced <- design_in[,c("genderMALE", "age",'Y', cancerCols)]
  
  if(mDNAsi){
    designReduced  <- cbind(designReduced , design_in[,"mDNAsi"])
    colnames(designReduced)[length(colnames(designReduced))] <- "mDNAsi"
  }
  
  if(raw){
    rawColNames <- colnames(design_in)[grepl("cg", colnames(design_in))]
    designReduced  <- cbind(designReduced, design_in[,rawColNames])
  }
  
  rfFit <- randomForest(Y~., data=data.frame(designReduced), ntree=15)
  
  return(rfFit)

}

getPredLinReg <- function(model, design){
  #Subset design to only use those variables used in the model
  designReduced <- design[,names(model$w)]
  
  #Return score
  return(designReduced %*% model$w + model$b)
}

getPredRFModel <- function(model, design){
  predict(model, newdata=data.frame(design)) 
}


#Get the RMSE
getRMSE <- function(observed, predicted){
  return(sqrt(mean((observed - predicted)^2)))
}

#USAGE:

lin_model <- getPenLinReg(design_train)
RF_model <- getRFModel(design_train)

lin_model_predictions <- getPredLinReg(lin_model, design_validate)
RF_model_predictions_train <- getPredRFModel(RF_model, design_train)
RF_model_predictions_validate <- getPredRFModel(RF_model, design_validate)

getRMSE(design_validate[,'Y'], lin_model_predictions)
getRMSE(design_train[,'Y'], RF_model_predictions_train)
getRMSE(design_validate[,'Y'], RF_model_predictions_validate)

plot(design_train[,'Y'], RF_model_predictions)
plot(design_validate[,'Y'], RF_model_predictions_validate)






# survival forests
rfsrc_fit <- rfsrc(Surv(Y, censor.indicator) ~ ., data=design_train_df, ntree=100)
y_hat <- predict(rfsrc_fit, data.frame(design_validate))
getRMSE(y_hat, design_validate[,'Y'])
# -- Libraries
library(randomForest)
library(tidyverse)
library(gelnet)
library(randomForestSRC)

##
split_data <- function(dat, proportions, seed)
{
  ### dat         : dataset to be split into train, validate, and test sets
  ### proportions : this should be vector with spliting proportions. 
  ###               Specifically, the first element corresponds to train proportion, the
  ###               second to validate proportion, and the last to test proportion
  set.seed(seed)
  
  # -- Number of observations
  n <- nrow(dat)
  
  # -- Train/Validate/Test indices
  idx_train    <- sample(1:n, round(n*proportions[1]))
  idx_validate <- sample((1:n)[!1:n %in% idx_train], round(nrow(dat) * proportions[2]))
  idx_test     <- sample((1:n)[!1:n %in% c(idx_train, idx_validate)], nrow(dat) - (length(idx_train) + length(idx_validate)))
  
  # -- Train/Validate/Test sets
  train    <- dat[idx_train,]
  validate <- dat[idx_validate,]
  test     <- dat[idx_test,]
  
  return(list("train"=as.matrix(train), "validate"=as.matrix(validate), "test"=as.matrix(test)))
  
  #####################################
  # g = sample(cut(
  #   seq(nrow(dat)),
  #   nrow(dat)*cumsum(c(0,spec)),
  #   labels = names(spec)
  # ))
  # 
  # dat_split <- split(dat, g)
  # 
  # dat_train <- as.matrix(dat_split$train)
  # dat_validate <- as.matrix(dat_split$validate)
  # dat_test <- as.matrix(dat_split$test)
  #####################################
}

##
getPenLinReg <- function(dat, mDNAsi=T, raw=T, mutations=T, l1=1, l2=1, max.iter = 100, eps = 1e-05)
{
  ### dat      : should include the outcome and possible covariates
  ### mDNAsi   : flag which if true includes the mDNAsi biomarker as a covariate in the design matrix
  ### raw      : flag which if true includes RNA expression values as covariates in the design matrix
  ### l1 & l2  : hyperparameters for the elastic net regression
  ### max.iter : maximum number of iterations for the iterative algorithm used to fing mle estimates of elastic net
  ### eps      : tolerance parameter for elastic net, this is related to convergence of parameters of the model
  
  # -- Extract outcome
  daysToDeath <- dat[,"Y"]
  
  # -- Extract dat matrix
  cancerCols <- colnames(dat)[grepl("cancer", colnames(dat))]
  datReduced <- dat[,c("genderMALE", "age", cancerCols)]
  
  # -- If mDNAsi = T then add to design matrix
  if(mDNAsi)
  {
    datReduced  <- cbind(datReduced , dat[,"mDNAsi"])
    colnames(datReduced)[length(colnames(datReduced))] <- "mDNAsi"
  }
  
  # -- If raw = T then add to design matrix
  if(raw)
  {
    rawColNames <- colnames(dat)[grepl("cg", colnames(dat))]
    datReduced  <- cbind(datReduced, dat[,rawColNames])
  }
  
  # -- If mutaions = T, then add to design matrix
  if(mutations)
  {
    geneColNames   <- colnames(dat)[grepl("gene", colnames(dat))]
    datReduced     <- cbind(datReduced, dat[, geneColNames])
  }
  #sapply(dat, class)
  #dat <- model.matrix( ~ ., dat)
  #dat <- as.matrix(dat)
  #elasticFit <- glmnet(x=datReduced, y=daysToDeath, alpha=0.5)
  #https://cran.r-project.org/web/packages/gelnet/gelnet.pdf
  
  # -- Running elastic net
  elasticFit <- gelnet(X=datReduced, y=daysToDeath, l1=l1, l2=l2, max.iter=max.iter, eps=eps, silent=TRUE)
  
  return(elasticFit)
}

##
getPredLinReg <- function(model, dat)
{
  ### model : elastice net model object from getPenLinReg
  ### dat   : dataset to make predictions
  
  # -- Subset design to only use those variables used in the model
  datReduced <- dat[,names(model$w)]
  
  #Return score
  return(datReduced %*% model$w + model$b)
}

##
getRFModel <- function(dat, mDNAsi=T, raw=T, mutations=T, ntree = 30, nodesize = 1, maxnodes = NULL)
{
  ### dat    : should include the outcome and possible covariates
  ### mDNAsi : flag which if true includes the mDNAsi biomarker as a covariate in the design matrix
  ### raw    : flag which if true includes RNA expression values as covariates in the design matrix
  ### ntree  : number of trees to be used in RF
  
  # -- Extract design matrix
  cancerCols <- colnames(dat)[grepl("cancer", colnames(dat))]
  datReduced <- dat[,c("genderMALE", "age",'Y', cancerCols)]
  
  # -- If mDNAsi = T then add to design matrix
  if(mDNAsi)
  {
    datReduced  <- cbind(datReduced , dat[,"mDNAsi"])
    colnames(datReduced)[length(colnames(datReduced))] <- "mDNAsi"
  }
  
  # -- If raw = T then add to design matrix
  if(raw)
  {
    rawColNames <- colnames(dat)[grepl("cg", colnames(dat))]
    datReduced  <- cbind(datReduced, dat[,rawColNames])
  }
  
  # -- If mutaions = T, then add to design matrix
  if(mutations)
  {
    geneColNames   <- colnames(dat)[grepl("gene", colnames(dat))]
    datReduced     <- cbind(datReduced, dat[, geneColNames])
  }
  # -- Fitting RF model
  rfFit <- randomForest(Y~., data=data.frame(datReduced), ntree=ntree, nodesize = nodesize, maxnodes = maxnodes)
  
  return(list('rfFit'=rfFit, 'chosenCols'=colnames(datReduced)))
}

##
getPredRFModel <- function(model, dat)
{
  ### model : Random forest model object
  ### dat   : data to make predictions 
  # predict(model, newdata=data.frame(dat)) 
  predict(model, newdata=dat)
}

##
getRMSE <- function(observed, predicted)
{
  ### observed  : Vector of observed outcomes
  ### predicted : Vector of predicted outcomes
  return(sqrt(mean((observed - predicted)^2)))
}

##
assess_RF <- function(train, validate, number_trees, mDNAsi=T, raw=T, mutations=T, verbose=T)
{
  # -- Matrix to save results
  resMat           <- matrix(NA, nrow=length(number_trees), ncol=3)
  colnames(resMat) <- c("trees", "train.mse", "validate.mse")
  resMat[,1]       <- number_trees
  names            <- colnames(train)

  for(i in 1:length(number_trees))
  {
    ##
    colnames(train)    <- names
    colnames(validate) <- names
    
    ##
    if(verbose){cat("Currently at iter =",i, "out of", nrow(resMat), "\n")}

    #  -- Fiting RF model
    tmp_model          <- getRFModel(dat=train, ntree=resMat[i,1], mDNAsi=mDNAsi, raw=raw, mutations=mutations)
    
    # # -- Predictions in the train/validate set
    tmp_preds_train    <- getPredRFModel(tmp_model$rfFit, dat=train)
    tmp_preds_validate <- getPredRFModel(tmp_model$rfFit, dat=validate)
    
    # -- MSE in the train/validate set
    resMat[i,2] <- getRMSE(train[,"Y"], tmp_preds_train)
    resMat[i,3] <- getRMSE(validate[,"Y"], tmp_preds_validate)
  }
  
  
  p <- resMat %>%
    as_tibble() %>%
    gather(set, mse, -trees) %>%
    ggplot(aes(trees, mse, color=set)) +
    geom_line(size=0.80) + 
    geom_point(size=3, alpha=0.90) +
    geom_point(size=3, pch=1, color="black") +
    ylab("RMSE") +
    xlab("# of trees") +
    scale_color_manual(name="",
                       values = c("black", "red3"),
                       labels = c("Train", "Validate")) +
    theme_classic() +
    theme(axis.title = element_text(face="bold"),
          axis.text  = element_text(face="bold"),
          legend.text = element_text(face="bold"),
          legend.position = "bottom")
  
  return(list("resMat" = resMat, "viz" = p))
}

##
assess_EN <- function(train, validate, l1s, l2s, mDNAsi=T, raw=T, mutations=T, verbose=T)
{
  # -- Matrix to save results
  resMat           <- cbind(as.matrix(expand.grid(l1s,l2s)), NA, NA)
  colnames(resMat) <- c("l1", "l2", "train.mse", "validate.mse")
  colnames(validate) <- colnames(train)
  
  for(i in 1:nrow(resMat))
  {
    if(verbose){cat("Currently at iter =",i, "out of", nrow(resMat), "\n")}
  
    #  -- Fiting EN model
    tmp_model <- getPenLinReg(dat=train, l1=resMat[i,1], l2=resMat[i,2], mDNAsi=mDNAsi, raw=raw, mutations=mutations)
    
    # -- Predictions in the train/validate set
    tmp_preds_train    <- getPredLinReg(tmp_model, dat=train)
    tmp_preds_validate <- getPredLinReg(tmp_model, dat=validate)
    
    # -- MSE in the train/validate set
    resMat[i,3] <- getRMSE(train[,"Y"], tmp_preds_train)
    resMat[i,4] <- getRMSE(validate[,"Y"], tmp_preds_validate)
  }
  
  middle <- (max(c(resMat[,3], resMat[,4])) + min(c(resMat[,3], resMat[,4])))/2
  
  p <- resMat %>%
    as_tibble() %>%
    gather(set, mse, -l1, -l2) %>%
    mutate(set = ifelse(set=="train.mse", "Train", "Validate")) %>%
    group_by(set) %>%
    mutate(mse.scaled = scale(mse)) %>%
    ungroup() %>%
    ggplot(aes(l1, l2, fill=mse)) +
    geom_tile(color="black",size=0.50) +
    scale_fill_gradient2(name="RMSE",
                         low="red", high="black", mid="white", midpoint = middle) +
    facet_wrap(~set) +
    theme_bw() +
    theme(axis.title = element_text(face="bold"),
          axis.text  = element_text(face="bold"),
          legend.text = element_text(face="bold"),
          legend.title = element_text(face="bold"),
          strip.text = element_text(face="bold"))
  
  return(list("resMat" = resMat, "viz" = p))
}

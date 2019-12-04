# -- Libraries
library(randomForest)
library(tidyverse)
library(gelnet)
library(randomForestSRC)

##
split_data <- function(dat, proportions)
{
  ### dat         : dataset to be split into train, validate, and test sets
  ### proportions : this should be vector with spliting proportions. 
  ###               Specifically, the first element corresponds to train proportion, the
  ###               second to validate proportion, and the last to test proportion
  
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

get_RF_survival_model <- function(dat, mDNAsi=T, raw=F, mutations=T, ntree = 150)
{
  ### dat    : should include the outcome and possible covariates
  ### mDNAsi : flag which if true includes the mDNAsi biomarker as a covariate in the design matrix
  ### raw    : flag which if true includes RNA expression values as covariates in the design matrix
  ### ntree  : number of trees to be used in RF
  
  # -- Extract design matrix
  cancerCols <- colnames(dat)[grepl("cancer", colnames(dat))]
  datReduced <- dat[,c("genderMALE", "age",'Y', 'censor.indicator', cancerCols)]
  
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
    first_gene_ind <- match('A2M',colnames(dat_all))
    geneColNames <- colnames(dat)[-0:-(first_gene_ind - 1)] 
    datReduced <- cbind(datReduced, dat[, geneColNames])
  }
  # -- Fitting RF Surival model
  rfsrc_fit <- rfsrc(Surv(Y, censor.indicator) ~ ., data=data.frame(datReduced), ntree=ntree)
  
  return(list('rfsrv_fit'=rfsrc_fit, 'chosenCols'=colnames(datReduced)))
}

##
get_pred_RF_survival_model <- function(model, dat)
{
  ### model : Random forest model object
  ### dat   : data to make predictions 
  predict(model, newdata=data.frame(dat))$predicted
}


##
getRMSE <- function(observed, predicted)
{
  ### observed  : Vector of observed outcomes
  ### predicted : Vector of predicted outcomes
  return(sqrt(mean((observed - predicted)^2)))
}

assess_RF_survival <- function(number_trees){
  # -- Matrix to save results
  resMat           <- matrix(NA, nrow=length(number_tress), ncol=3)
  colnames(resMat) <- c("trees", "train.mse", "validate.mse")
  resMat[,1]       <- number_tress
  
  for(i in 1:length(number_tress))
  {
    cat("Currently at iter =",i, "out of", nrow(resMat), "\n")
    #  -- Fiting RF model
    tmp_model          <- get_RF_survival_model(dat=train, ntree=resMat[i,1])
    
    # -- Predictions in the train/validate set
    tmp_preds_train    <- get_pred_RF_survival_model(tmp_model, dat=train)
    tmp_preds_validate <- get_pred_RF_survival_model(tmp_model, dat=validate)
    
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
    ylab("MSE") +
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




#Get predicted values given a linear regression model from gelnet and design matrix

#Input: model, design matrix
#Output: predicted values

getPredLinReg <- function(model, design){
  #Subset design to only use those variables used in the model
  designReduced <- design[,names(model$w)]
  
  #Return score
  return(designReduced %*% model$w + model$b)
}
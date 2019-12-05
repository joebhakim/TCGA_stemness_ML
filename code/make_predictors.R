# -- Init code. Any global variables/functions/libraries should be initialized in init.R
source("code/init.R")

# -- To ensure reproducible results
set.seed(1)

# -- Loading data
load('data/design-matrix.rda')
design <- as_tibble(design)

# -- Train/Validate/Test proportions (Maybe we can turns this into a function)
set_proportions <- c(train = 0.6, validate = 0.2, test = 0.2)

# -- Splitting the data
design_nocensored <- filter(design, censor.indicator == 0)
dat_list          <- split_data(design_nocensored, set_proportions)
train             <- dat_list$train
validate          <- dat_list$validate
test              <- dat_list$test
names             <- colnames(train)
colnames(validate) <- names

# -- Fitting EN and RF models
lin_model <- getPenLinReg(train, mDNAsi=T, raw=F, mutations=F)
RF_model_AND_chosenCols <- getRFModel(train, raw=F, mutations=F)
RF_model                <- RF_model_AND_chosenCols$rfFit      # -- Why are we doing this?
chosenCols              <- RF_model_AND_chosenCols$chosenCols # -- Why are we doing this?

# -- Computing predictions in the validation set
lin_preds <- getPredLinReg(lin_model, validate)

# which(colnames(validate) == "cancer.classneural crest")
RF_preds  <- getPredRFModel(RF_model, validate)
# RF_preds  <- getPredRFModel(RF_model, validate[,chosenCols])

# -- Computing MSE
cat("MSE for EN:",getRMSE(validate[,'Y'], lin_preds),"\n")
cat("MSE for RF:",getRMSE(validate[,'Y'], RF_preds))

# -- Multiple values for RF hyperparameters
number_trees    <- seq(0, 250, by = 50)
number_trees[1] <- 5
# -- Multiple values for EN hyperparameters
l1s <- seq(0, 2, by=0.1)
l2s <- seq(0, 2, by=0.1)

# -- Creating table for results
mDNAsi.vec <- c(T, F)
raw.vec    <- c(T, F)
mut.vec    <- c(T, F)
resMat <- expand.grid(mDNAsi.vec, raw.vec, mut.vec)
colnames(resMat) <- c("mDNAsi", "raw", "mutations")
resMat[,"MSE.EN"] <- NA
resMat[,"MSE.RF"] <- NA
resMat[,"MSE.NN"] <- NA

for(i in 1:nrow(resMat))
{
  # cat("Iter", i, "of", nrow(resMat), "\n")
  
  # # -- Random Forest part of the matrix
  # cat(">> Random Forest\n")
  # tmp_RF <- assess_RF(train, validate, number_trees, mDNAsi=resMat[i,"mDNAsi"], raw=resMat[i,"raw"], mutations=resMat[i,"mutations"], verbose=FALSE)$resMat
  # best_mse <- tmp_RF %>%
  #                   as_tibble() %>%
  #                   filter(validate.mse == min(validate.mse)) %>%
  #                   .$validate.mse
  # best_hyper <- tmp_RF %>%
  #                   as_tibble() %>%
  #                   filter(validate.mse == min(validate.mse)) %>%
  #                   .$trees
  # resMat[i,"MSE.RF"] <- paste0(round(best_mse,2)," (",best_hyper,")")
  
  # -- Random Forest part of the matrix
  cat(">> Elastic Net\n")
  tmp_EN <- assess_EN(train, validate, l1s, l2s, mDNAsi=resMat[i,"mDNAsi"], raw=resMat[i,"raw"], mutations=resMat[i,"mutations"], verbose=F)$resMat
  best_mse <- tmp_EN %>%
                    as_tibble() %>%
                    filter(validate.mse == min(validate.mse)) %>%
                    .$validate.mse
  best_hyper1 <- tmp_EN %>%
                    as_tibble() %>%
                    filter(validate.mse == min(validate.mse)) %>%
                    .$l1
  best_hyper2 <- tmp_EN %>%
                    as_tibble() %>%
                    filter(validate.mse == min(validate.mse)) %>%
                    .$l2
  resMat[i,"MSE.EN"] <- paste0(round(best_mse,2)," (",best_hyper1,",",best_hyper2,")")
}
resMat





















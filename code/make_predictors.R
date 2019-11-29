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
dat_list <- split_data(design, set_proportions)
train    <- dat_list$train
validate <- dat_list$validate
test     <- dat_list$test

# -- Fitting EN and RF models
lin_model <- getPenLinReg(train)
RF_model  <- getRFModel(train)

# -- Computing predictions in the validation set
lin_preds <- getPredLinReg(lin_model, validate)
RF_preds  <- getPredRFModel(RF_model, validate)

# -- Computing MSE
cat("MSE for EN:",getRMSE(validate[,'Y'], lin_preds),"\n")
cat("MSE for RF:",getRMSE(validate[,'Y'], RF_preds))


# -- Multiple values for RF hyperparameters
number_tress    <- seq(0, 400, by = 20)
number_tress[1] <- 1

# -- Multiple values for EN hyperparameters
l1s <- seq(0, 5, by=0.10)
l2s <- seq(0, 5, by=0.10)

# -- Iterating over mutiple values of ntree
resRF <- assess_RF(number_tress)
resEN <- assess_EN(l1s, l2s)

resRF$viz
resEN$viz




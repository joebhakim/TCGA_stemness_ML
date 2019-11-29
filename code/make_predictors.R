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
dat_list <- split_data(design, set_proportions)
#dat_list <- split_data(design_nocensored, set_proportions)
train    <- dat_list$train
validate <- dat_list$validate
test     <- dat_list$test

# -- Fitting EN and RF models
lin_model <- getPenLinReg(train)
RF_model  <- getRFModel(train)
RF_survival_model_AND_chosenCols <- get_RF_survival_model(train)
RF_survival_model <- RF_survival_model_AND_chosenCols$rfsrv_fit
chosenCols <- RF_survival_model_AND_chosenCols$chosenCols

# -- Computing predictions in the validation set
lin_preds <- getPredLinReg(lin_model, validate)
RF_preds  <- getPredRFModel(RF_model, validate)
RF_survival_preds <- get_pred_RF_survival_model(RF_survival_model, validate[,chosenCols])

# -- Computing MSE
cat("MSE for EN:",getRMSE(validate[,'Y'], lin_preds),"\n")
cat("MSE for RF:",getRMSE(validate[,'Y'], RF_preds))
cat('MSE for RF survival', getRMSE(validate[,'Y'], RF_survival_preds$predicted))

# -- Multiple values for RF hyperparameters
number_tress    <- seq(0, 150, by = 50)
number_tress[1] <- 1

# -- Multiple values for EN hyperparameters
l1s <- seq(0, 5, by=0.10)
l2s <- seq(0, 5, by=0.10)

# -- Iterating over mutiple values of ntree
resRF <- assess_RF(number_tress)
resEN <- assess_EN(l1s, l2s)

resRF$viz
resEN$viz




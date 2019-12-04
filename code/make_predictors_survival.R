# -- Init code. Any global variables/functions/libraries should be initialized in init.R
source("code/init_survival.R")

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
RF_survival_model_AND_chosenCols <- get_RF_survival_model(train)
RF_survival_model <- RF_survival_model_AND_chosenCols$rfsrv_fit
chosenCols <- RF_survival_model_AND_chosenCols$chosenCols

# -- Computing predictions in the validation set
RF_survival_preds <- get_pred_RF_survival_model(RF_survival_model, validate[,chosenCols])

# -- Computing MSE
cat('MSE for RF survival', getRMSE(validate[,'Y'], RF_survival_preds))

# -- Multiple values for RF hyperparameters
number_tress    <- seq(5, 20, by = 5)
#number_tress[1] <- 1

# -- Multiple values for EN hyperparameters
l1s <- seq(0, 5, by=1)
l2s <- seq(0, 0.1, by=.02)

# -- Iterating over mutiple values of ntree
resRFSurvival <- assess_RF(number_tress)

resRF$viz
resEN$viz


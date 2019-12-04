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

# -- Need to do this for some stupid reason (doing it the old-fashion way)
colnames(validate)[8]  <- "cancer.classdevelopmental.GI"
colnames(validate)[9]  <- "cancer.classhead.and.neck"
colnames(validate)[12] <- "cancer.classneural.crest"
colnames(validate)[13] <- "cancer.classsoft.tissue"
colnames(validate)[14] <- "cancer.classsolic..urologic"
colnames(validate)[15] <- "cancer.classsolid..core.GI"
colnames(validate)[16] <- "cancer.classsolid..endocrine"
colnames(validate)[17] <- "cancer.classsolid..gynecologic"
# which(colnames(validate) == "cancer.classneural crest")
RF_preds  <- getPredRFModel(RF_model, validate)
# RF_preds  <- getPredRFModel(RF_model, validate[,chosenCols])

# -- Computing MSE
cat("MSE for EN:",getRMSE(validate[,'Y'], lin_preds),"\n")
cat("MSE for RF:",getRMSE(validate[,'Y'], RF_preds))

# -- Multiple values for RF hyperparameters
number_trees    <- seq(50, 250, by = 50)
#number_tress[1] <- 1

# -- Iterating over mutiple values of ntree
resRF <- assess_RF(number_trees, mDNAsi=T, raw=F, mutations=F) # -- When mutations=T it gives problems
resRF$viz


# -- Multiple values for EN hyperparameters
l1s <- seq(0, 6, by=0.1)
l2s <- seq(0, 6, by=0.1)

# -- Iterating over l1s and l2s
resEN <- assess_EN(l1s, l2s, mDNAsi=T, raw=T, mutations=T)
resEN$viz

resEN$resMat %>%
  as_tibble() %>%
  # filter(l1==5)
  filter(validate.mse == min(validate.mse))
  # filter(train.mse == min(train.mse))
# T,T,T: 4.3   0.9      845.         911
# T,T,F: 5.8   0.7      874.         914.
# T,F,F: 1     0      945.         905.









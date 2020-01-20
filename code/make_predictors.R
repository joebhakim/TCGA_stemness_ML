# -- Init code. Any global variables/functions/libraries should be initialized in init.R
source("code/init.R")

# -- rf TTF results
res <- read_csv("data/results-rf2.csv")
res %>%
  group_by(trees, code) %>%
  summarize(train.sd     = sd(train.mse),
            validate.sd  = sd(validate.mse),
            train.mse    = mean(train.mse),
            validate.mse = mean(validate.mse)) %>%
  ungroup() %>%
  filter(code == "TTF") %>%
  ggplot() +
  geom_errorbar(aes(trees, 
                    ymin=train.mse-train.sd,
                    ymax=train.mse+train.sd, color="Train"), width=10) +
  geom_line(aes(trees, train.mse, color="Train")) +
  geom_point(aes(trees, train.mse, color="Train"), size=2) +
  geom_point(aes(trees, train.mse), size=2, pch=1) +
  geom_errorbar(aes(trees, 
                    ymin=validate.mse-validate.sd,
                    ymax=validate.mse+validate.sd, color="Validate"), width=10) +
  geom_line(aes(trees, validate.mse, color="Validate")) +
  geom_point(aes(trees, validate.mse, color="Validate"), size=2) +
  geom_point(aes(trees, validate.mse), size=2, pch=1) +
  xlab("# of trees") +
  ylab("Avg. RMSE") +
  scale_color_manual(name="Set",
                     values=c("#525252", "#2171b5")) +
  theme_minimal() +
  theme(axis.title   = element_text(face="bold"),
        axis.text    = element_text(face="bold"),
        legend.text  = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        strip.text   = element_text(face="bold", color="black"))

# -- rf TFF results
res <- read_csv("data/results-rf2.csv")
res %>%
  group_by(trees, code) %>%
  summarize(train.sd     = sd(train.mse),
            validate.sd  = sd(validate.mse),
            train.mse    = mean(train.mse),
            validate.mse = mean(validate.mse)) %>%
  ungroup() %>%
  filter(code == "TFF") %>%
  ggplot() +
  geom_errorbar(aes(trees, 
                    ymin=train.mse-train.sd,
                    ymax=train.mse+train.sd, color="Train"), width=10) +
  geom_line(aes(trees, train.mse, color="Train")) +
  geom_point(aes(trees, train.mse, color="Train"), size=2) +
  geom_point(aes(trees, train.mse), size=2, pch=1) +
  geom_errorbar(aes(trees, 
                    ymin=validate.mse-validate.sd,
                    ymax=validate.mse+validate.sd, color="Validate"), width=10) +
  geom_line(aes(trees, validate.mse, color="Validate")) +
  geom_point(aes(trees, validate.mse, color="Validate"), size=2) +
  geom_point(aes(trees, validate.mse), size=2, pch=1) +
  xlab("# of trees") +
  ylab("Avg. RMSE") +
  scale_color_manual(name="Set",
                     values=c("#525252", "#2171b5")) +
  theme_minimal() +
  theme(axis.title   = element_text(face="bold"),
        axis.text    = element_text(face="bold"),
        legend.text  = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        strip.text   = element_text(face="bold", color="black"))

# -- en-grid-TTF
results_en <- read_csv("data/results-en.csv")
tmp <- results_en %>%
  as_tibble() %>%
  filter(code == "TTF") %>%
  group_by(l1, l2) %>%
  summarize(train.mse = mean(train.mse),
            validate.mse = mean(validate.mse)) %>%
  ungroup() %>%
  gather(set, mse, -l1, -l2)
tmp %>%
  filter(set == "validate.mse") %>%
  filter(mse==min(mse))
middle <- (min(tmp$mse) + max(tmp$mse)) / 2
tmp %>%
  mutate(set = ifelse(set=="train.mse", "Train", "Validate")) %>%
  ggplot(aes(l1, l2, fill=mse)) +
  geom_tile(color="black",size=0.50) +
  scale_fill_gradient2(name="Avg. \n RMSE",
                       low="#3690c0", high="black", mid="white", midpoint = middle) +
  xlab(expression(lambda[1])) +
  ylab(expression(lambda[2])) +
  facet_wrap(~set) +
  theme_bw() +
  theme(axis.title = element_text(face="bold", size=15),
        axis.text  = element_text(face="bold", size=10),
        legend.text = element_text(face="bold", size=10),
        legend.title = element_text(face="bold", size=12),
        strip.text = element_text(face="bold", size=11, color="black"))

# -- en-grid-TFF
results_en <- read_csv("data/results-en.csv")
tmp <- results_en %>%
  as_tibble() %>%
  filter(code == "TFF") %>%
  group_by(l1, l2) %>%
  summarize(train.mse = mean(train.mse),
            validate.mse = mean(validate.mse)) %>%
  ungroup() %>%
  gather(set, mse, -l1, -l2)
tmp %>%
  filter(set == "validate.mse") %>%
  filter(mse==min(mse))
middle <- (min(tmp$mse) + max(tmp$mse)) / 2
tmp %>%
  mutate(set = ifelse(set=="train.mse", "Train", "Validate")) %>%
  ggplot(aes(l1, l2, fill=mse)) +
  geom_tile(color="black",size=0.50) +
  scale_fill_gradient2(name="Avg. \n RMSE",
                       low="#3690c0", high="black", mid="white", midpoint = middle) +
  xlab(expression(lambda[1])) +
  ylab(expression(lambda[2])) +
  facet_wrap(~set) +
  theme_bw() +
  theme(axis.title = element_text(face="bold", size=15),
        axis.text  = element_text(face="bold", size=10),
        legend.text = element_text(face="bold", size=10),
        legend.title = element_text(face="bold", size=12),
        strip.text = element_text(face="bold", size=11, color="black"))

# -- Optimal NN grid
nn_train    <- read.csv("data/neural100_meanBestRMSE_train.csv", header=FALSE)
nn_validate <- read.csv("data/neural100_meanBestRMSE_val.csv", header=FALSE)
nn_train    <- nn_train %>%
                setNames(c(1, 10, 100, 1000, 1500)) %>%
                mutate(nodes1 = c(1, 10, 100, 1000, 1500)) %>%
                gather(nodes2, train.mse, -nodes1)
nn_validate <-  nn_validate %>%
                  setNames(c(1, 10, 100, 1000, 1500)) %>%
                  mutate(nodes1 = c(1, 10, 100, 1000, 1500)) %>%
                  gather(nodes2, validate.mse, -nodes1)
tmp         <- nn_train %>%
                left_join(nn_validate, by = c("nodes1", "nodes2")) %>%
                gather(set, mse, -nodes1, -nodes2)
middle      <- (max(tmp[,4]) + min(tmp[,4]))/2
tmp %>%
  mutate(set = ifelse(set=="train.mse", "Train", "Validate")) %>%
  ggplot(aes(as.factor(nodes1), as.factor(nodes2), fill = mse)) +
  geom_tile(color="black",size=0.50) +
  scale_fill_gradient2(name="Avg. RMSE",
                       low="#3690c0", high="black", mid="white", midpoint = middle) +
  xlab("Nodes in hidden layer 1") +
  ylab("Nodes in hidden layer 2") +
  facet_wrap(~set) +
  theme_bw() +
  theme(axis.title = element_text(face="bold", size=12),
        axis.text  = element_text(face="bold", size=10),
        legend.text = element_text(face="bold", size=10),
        legend.title = element_text(face="bold", size=12),
        strip.text = element_text(face="bold", size=11, color="black"))

# -- Loading data
load('data/design-matrix.rda')
design            <- as_tibble(design)
design_nocensored <- filter(design, censor.indicator == 0)

# -- Train/Validate/Test proportions (Maybe we can turns this into a function)
set_proportions <- c(train = 0.6, validate = 0.2, test = 0.2)
R                 <- 200
test_rf           <- rep(NA, R)
design_nocensored <- filter(design, censor.indicator == 0)
for(i in 1:R)
{
  cat(paste0("Iter: ",i,"/",R),"\n")
  
  dat_list          <- split_data(design_nocensored, set_proportions, seed=i)
  train             <- dat_list$train
  validate          <- dat_list$validate
  test              <- dat_list$test
  
  # -- Fitting RF model
  RF_model_AND_chosenCols <- getRFModel(train, mDNAsi=F, raw=F, mutations=F, ntree=50)
  RF_model                <- RF_model_AND_chosenCols$rfFit      # -- Why are we doing this?
  chosenCols              <- RF_model_AND_chosenCols$chosenCols # -- Why are we doing this?
  
  # which(colnames(validate) == "cancer.classneural crest")
  RF_preds  <- getPredRFModel(RF_model, test)
  
  # -- Computing MSE
  test_rf[i] <- getRMSE(test[,'Y'], RF_preds)
}
mean(test_rf)
sd(test_rf)










# -- To ensure reproducible results
seed=2
# -- Splitting the data
design_nocensored <- filter(design, censor.indicator == 0)
dat_list          <- split_data(design_nocensored, set_proportions, seed=3)
train             <- dat_list$train
validate          <- dat_list$validate
test              <- dat_list$test

# -- Fitting EN and RF models
lin_model               <- getPenLinReg(train, mDNAsi=T, raw=T, mutations=T, l1=2,l2=0.60)
RF_model_AND_chosenCols <- getRFModel(train, mDNAsi=T, raw=F, mutations=F, ntree=200)
RF_model                <- RF_model_AND_chosenCols$rfFit      # -- Why are we doing this?
chosenCols              <- RF_model_AND_chosenCols$chosenCols # -- Why are we doing this?

# -- Computing predictions in the validation set
lin_preds <- getPredLinReg(lin_model, test)

# which(colnames(validate) == "cancer.classneural crest")
RF_preds  <- getPredRFModel(RF_model, test)

# -- Computing MSE
cat("MSE for EN:",getRMSE(test[,'Y'], lin_preds),"\n")
cat("MSE for RF:",getRMSE(test[,'Y'], RF_preds))

# -- RF optimal
res_rf <- read.csv("~/Desktop/results-rf.csv")
res_rf %>%
  as_tibble() %>%
  filter(code == "TFF") %>%
  # mutate(lower = mean - 2*se,
  #        upper = mean + 2*se) %>%
  ggplot(aes(trees, mean)) +
  # geom_ribbon(aes(ymin = lower, ymax=upper), fill="#3690c0", color="black", alpha=0.70) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), color="#3690c0", width=10, size=0.80) +
  geom_line(size=0.80, color="#3690c0") +
  geom_point(size=3, color="#3690c0", alpha=0.95) +
  geom_point(size=3, pch=1, color="black") +
  ylab("Avg. RMSE") +
  xlab("# of trees") +
  theme_bw() +
  theme(axis.title = element_text(face="bold", size=12),
        axis.text  = element_text(face="bold", size=10),
        legend.text = element_text(face="bold", size=10),
        legend.title = element_text(face="bold", size=12),
        strip.text = element_text(face="bold", size=11, color="black"))

  
  









res_EN_TTF <- assess_EN(train=train, validate=validate, l1s=l1s,l2s=l2s, mDNAsi=T, raw=F, mutations=F)

res_EN_TTF$viz
res_EN_TTF$resMat %>%
  as_tibble() %>%
  filter(validate.mse == min(validate.mse))



# -- Multiple values for RF hyperparameters
number_trees    <- seq(0, 250, by = 50)
number_trees[1] <- 5

res_RF_TTF <- assess_RF(train, validate, number_trees=number_trees, mDNAsi=T, raw=T, mutations=F)
res_RF_TTF$viz
res_RF_TTF$resMat

res_EN_TTF <- assess_EN(train, validate, l1s=l1s,l2s=l2s, mDNAsi=T, raw=F, mutations=F)
res_EN_TTF$viz
res_EN_TTF$resMat

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
write.csv(resMat, file="data/table-1.csv", row.names = F)








































#########################################################
R <- 100
test_en <- rep(NA, R)
test_rf <- rep(NA, R)
# -- Splitting the data
design_nocensored <- filter(design, censor.indicator == 0)
for(i in 1:R)
{
  cat(paste0("Iter: ",i,"/",R),"\n")
  
  dat_list          <- split_data(design_nocensored, set_proportions, seed=i)
  train             <- dat_list$train
  validate          <- dat_list$validate
  test              <- dat_list$test
  
  # -- Fitting EN and RF models
  lin_model               <- getPenLinReg(train, mDNAsi=T, raw=T, mutations=T, l1=2,l2=0.60)
  RF_model_AND_chosenCols <- getRFModel(train, mDNAsi=F, raw=F, mutations=F, ntree=50)
  RF_model                <- RF_model_AND_chosenCols$rfFit      # -- Why are we doing this?
  chosenCols              <- RF_model_AND_chosenCols$chosenCols # -- Why are we doing this?
  
  # -- Computing predictions in the validation set
  lin_preds <- getPredLinReg(lin_model, test)
  
  # which(colnames(validate) == "cancer.classneural crest")
  RF_preds  <- getPredRFModel(RF_model, test)
  
  # -- Computing MSE
  test_en[i] <- getRMSE(test[,'Y'], lin_preds)
  test_rf[i] <- getRMSE(test[,'Y'], RF_preds)
}
mean(test_en)
sd(test_en)
mean(test_rf)
sd(test_rf)
hist(test_en)
hist(test_rf)
#########################################################

#########################################################
# -- Multiple values for EN hyperparameters
l1s <- seq(0, 20,   by=2)
l2s <- seq(0, 0.60, by=0.10)
R   <- 10
dat         <- design_nocensored
proportions <- set_proportions
res <- lapply(1:R, function(x)
{
  set.seed(x)
  
  cat(paste0("Iteration: ",x,"/",R), "\n")
  # -- Number of observations
  n <- nrow(dat)
  
  # -- Train/Validate/Test indices
  idx_train    <- sample(1:n, round(n*proportions[1]))
  idx_validate <- sample((1:n)[!1:n %in% idx_train], round(nrow(dat) * proportions[2]))
  idx_test     <- sample((1:n)[!1:n %in% c(idx_train, idx_validate)], nrow(dat) - (length(idx_train) + length(idx_validate)))
  
  # -- Train/Validate/Test sets
  tmp.train    <- as.matrix(dat[idx_train,])
  tmp.validate <- as.matrix(dat[idx_validate,])
  
  cat("T:T:T\n") #1 -- T:T:T
  tmp.ttt <- assess_EN(train=tmp.train, validate=tmp.validate, l1s=l1s,l2s=l2s, mDNAsi=T, raw=T, mutations=T,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "TTT")
  
  cat("F:T:T\n") #2 -- F:T:T
  tmp.ftt <- assess_EN(train=tmp.train, validate=tmp.validate, l1s=l1s,l2s=l2s, mDNAsi=F, raw=T, mutations=T,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "FTT")
  
  cat("T:F:T\n") #3 -- T:F:T
  tmp.tft <- assess_EN(train=tmp.train, validate=tmp.validate, l1s=l1s,l2s=l2s, mDNAsi=T, raw=F, mutations=T,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "TFT")
  
  cat("F:F:T\n") #4 -- F:F:T
  tmp.fft <- assess_EN(train=tmp.train, validate=tmp.validate, l1s=l1s,l2s=l2s, mDNAsi=F, raw=F, mutations=T,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "FFT")
  
  cat("T:T:F\n") #5 -- T:T:F
  tmp.ttf <- assess_EN(train=tmp.train, validate=tmp.validate, l1s=l1s,l2s=l2s, mDNAsi=T, raw=T, mutations=F,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "TTF")
  
  cat("F:T:F\n") #6 -- F:T:F
  tmp.ftf <- assess_EN(train=tmp.train, validate=tmp.validate, l1s=l1s,l2s=l2s, mDNAsi=F, raw=T, mutations=F,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "FTF")
  
  cat("T:F:F\n") #7 -- T:F:F
  tmp.tff <- assess_EN(train=tmp.train, validate=tmp.validate, l1s=l1s,l2s=l2s, mDNAsi=T, raw=F, mutations=F,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "TFF")
  
  cat("F:F:F\n") #8 -- F:F:F
  tmp.fff <- assess_EN(train=tmp.train, validate=tmp.validate, l1s=l1s,l2s=l2s, mDNAsi=F, raw=F, mutations=F,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "FFF")
  
  res <- bind_rows(tmp.ttt, tmp.ftt, tmp.tft, tmp.fft,
                   tmp.ttf, tmp.ftf, tmp.tff, tmp.fff)
  
  return(res)
  
})
res <- do.call(rbind, res)
#########################################################

#########################################################
# -- Multiple values for RF hyperparameters
number_trees <- seq(50,250,by=50)
R   <- 10
dat         <- design_nocensored
proportions <- set_proportions
res <- lapply(1:R, function(x)
{
  set.seed(x)
  
  cat(paste0("Iteration: ",x,"/",R), "\n")
  # -- Number of observations
  n <- nrow(dat)
  
  # -- Train/Validate/Test indices
  idx_train    <- sample(1:n, round(n*proportions[1]))
  idx_validate <- sample((1:n)[!1:n %in% idx_train], round(nrow(dat) * proportions[2]))
  idx_test     <- sample((1:n)[!1:n %in% c(idx_train, idx_validate)], nrow(dat) - (length(idx_train) + length(idx_validate)))
  
  # -- Train/Validate/Test sets
  tmp.train    <- as.matrix(dat[idx_train,])
  tmp.validate <- as.matrix(dat[idx_validate,])
  
  cat("T:T:F\n") #5 -- T:T:F
  tmp.ttf <- assess_RF(train=tmp.train, validate=tmp.validate, number_trees=number_trees, mDNAsi=T, raw=T, mutations=F,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "TTF")
  
  cat("T:F:F\n") #7 -- T:F:F
  tmp.tff <- assess_RF(train=tmp.train, validate=tmp.validate, number_trees=number_trees, mDNAsi=T, raw=F, mutations=F,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "TFF")
  
  res <- bind_rows(tmp.ttf, tmp.tff)
  
  return(res)
})
res <- do.call(rbind, res)
#########################################################

#########################################################
# -- Multiple values for RF hyperparameters
number_trees <- seq(50,250,by=50)
R   <- 10
dat         <- design_nocensored
proportions <- set_proportions
res <- lapply(1:R, function(x)
{
  set.seed(x)
  
  cat(paste0("Iteration: ",x,"/",R), "\n")
  # -- Number of observations
  n <- nrow(dat)
  
  # -- Train/Validate/Test indices
  idx_train    <- sample(1:n, round(n*proportions[1]))
  idx_validate <- sample((1:n)[!1:n %in% idx_train], round(nrow(dat) * proportions[2]))
  idx_test     <- sample((1:n)[!1:n %in% c(idx_train, idx_validate)], nrow(dat) - (length(idx_train) + length(idx_validate)))
  
  # -- Train/Validate/Test sets
  tmp.train    <- as.matrix(dat[idx_train,])
  tmp.validate <- as.matrix(dat[idx_validate,])
  
  cat("T:T:T\n") #1 -- T:T:T
  tmp.ttt <- assess_RF(train=tmp.train, validate=tmp.validate, number_trees=number_trees, mDNAsi=T, raw=T, mutations=T,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "TTT")
  
  cat("F:T:T\n") #2 -- F:T:T
  tmp.ftt <- assess_RF(train=tmp.train, validate=tmp.validate, number_trees=number_trees, mDNAsi=F, raw=T, mutations=T,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "FTT")
  
  cat("T:F:T\n") #3 -- T:F:T
  tmp.tft <- assess_RF(train=tmp.train, validate=tmp.validate, number_trees=number_trees, mDNAsi=T, raw=F, mutations=T,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "TFT")
  
  cat("F:F:T\n") #4 -- F:F:T
  tmp.fft <- assess_RF(train=tmp.train, validate=tmp.validate, number_trees=number_trees, mDNAsi=F, raw=F, mutations=T,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "FFT")
  
  cat("T:T:F\n") #5 -- T:T:F
  tmp.ttf <- assess_RF(train=tmp.train, validate=tmp.validate, number_trees=number_trees, mDNAsi=T, raw=T, mutations=F,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "TTF")
  
  cat("F:T:F\n") #6 -- F:T:F
  tmp.ftf <- assess_RF(train=tmp.train, validate=tmp.validate,number_trees=number_trees, mDNAsi=F, raw=T, mutations=F,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "FTF")
  
  cat("T:F:F\n") #7 -- T:F:F
  tmp.tff <- assess_RF(train=tmp.train, validate=tmp.validate, number_trees=number_trees, mDNAsi=T, raw=F, mutations=F,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "TFF")
  
  cat("F:F:F\n") #8 -- F:F:F
  tmp.fff <- assess_RF(train=tmp.train, validate=tmp.validate, number_trees=number_trees, mDNAsi=F, raw=F, mutations=F,verbose=F)$resMat %>%
    as_tibble() %>%
    mutate(code = "FFF")
  
  res <- bind_rows(tmp.ttt, tmp.ftt, tmp.tft, tmp.fft,
                   tmp.ttf, tmp.ftf, tmp.tff, tmp.fff)
  
  return(res)
  
})
res <- do.call(rbind, res)
#########################################################
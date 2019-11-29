library(gelnet)
library(dplyr)
library(gdata)
library(DT)


load("data/data.pan.Rda")



replace.NA <-function(data,type.info,by = "mean"){
  if(!"group" %in% colnames(type.info)) stop("type.info must have group column")
  if(!"sample" %in% colnames(type.info)) stop("type.info must have a sample column")
  
  # Do we have NAs?
  if(is.na(table(is.na(data))["TRUE"])){
    message("No NAs were found")
    return(data)
  }
  # get NAs index 
  idx <- which(is.na(data) == TRUE,arr.ind=TRUE)
  count <- table(rownames(idx))
  message("======= Status Number of NA in probes ========")
  message("--------------------- Summary------------------")
  print(summary(as.numeric(count)))
  message("\n----------- Probes with more nb of NAs -----------")
  print(head(sort(count,decreasing = T)))
  message("===============================================")
  
  idx <- cbind(idx, mean = NA, median = NA)
  
  # For each NA value calculate the mean for the same probe for the samples
  # where it belongs
  for(line in 1:nrow(idx)){
    row <- idx[line,1]
    col <- idx[line,2]
    probe <- rownames(idx)[line]
    sample <- colnames(data)[col]
    group <- type.info[type.info$sample == sample,"group"]
    samples.in.group <- type.info[type.info$group == group,]$sample
    
    # get the probe value for all samples in the same group 
    aux <- data[rownames(data) %in% probe, colnames(data) %in% samples.in.group] 
    
    idx[line,3] <- mean(as.numeric(aux),na.rm = TRUE)
    idx[line,4] <- median(as.numeric(aux),na.rm = TRUE)
  }
  # Step 2 replace
  for(line in 1:nrow(idx)){
    row <- idx[line,1]
    col <- idx[line,2]
    if(by == "mean"){
      data[idx[line,1],idx[line,2]] <- idx[line,3]  
    } else if(by == "median") { 
      data[idx[line,1],idx[line,2]] <- idx[line,4]
    }
  }
  return(data)
}



#load in PCBC data
load("data/pcbc.data.Rda")
#preprocess PCBC data
m <- apply(pcbc.data, 1, mean)
pcbc.data.2 <- pcbc.data - m

#separate out stem vs not stem
load("data/pcbc.pd.f.Rda")
M1_smp <- pcbc.pd.f[pcbc.pd.f$Diffname_short %in% "SC",] #SC
M2_smp <- pcbc.pd.f[!(pcbc.pd.f$Diffname_short %in% "SC"),] #non-SC
# Select PCBC data
X.tr <- pcbc.data.2[, as.character(M1_smp$UID)] # 44 samples
X.bk <- pcbc.data.2[, as.character(M2_smp$UID)] # 55 samples


#load in TCGA data
load("data/data.pan.Rda") 
load("data/type.info.Rda")  #has tumor type info
testset <- replace.NA(data.pan, type.info, by="median") 

#bhle
#PCA <- prcomp(t(X.tr))
#SC <- data.frame(PCA$x)
#nonSC <- data.frame(predict(PCA, t(X.bk)))
#cancer <- data.frame(predict(PCA, t(data.pan))[0:1000,])

PCA <- prcomp(t(testset))
cancer <- data.frame(PCA$x)
nonSC <- data.frame(predict(PCA, t(X.bk)))
SC <- data.frame(predict(PCA, t(X.tr)))

mDNA_PCd = bind_rows(list(cancer=cancer, nonSC=nonSC, SC=SC), .id='source')

#ggplot(b, aes(PC1, PC2)) + geom_point() + geom_point(data=c, aes(PC1, PC2, color='red'))
ggplot(data=mDNA_PCd, aes(PC1, PC2, color=source)) + geom_point()


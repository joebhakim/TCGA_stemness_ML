#Setup
##Change working directory
setwd("E:/Matt Work/Harvard/Machine Learning MIT 6867/Project")

##Load the survival package, survminer package, dplyr package, and wrangled data set
if(!require(survival)){
  install.packages("survival")
}
if(!require(survminer)){
  install.packages("survminer")
}
if(!require(dplyr)){
  install.packages("dplyr")
}
library(survival)
library(survminer)
library(dplyr)
load("data/tcga_wrangled.rda")

#Rename EREG-mRNAsi
colnames(dat)[colnames(dat)=="EREG-mRNAsi"] <- "EREG_mRNAsi"

#Multiply all stemness indices by 100 for regression interpretations
dat[,c("mRNAsi", "mDNAsi", "EREG_mRNAsi")] <- 100*dat[,c("mRNAsi", "mDNAsi", "EREG_mRNAsi")]


#Survival Analysis
#Referred to: https://www.datacamp.com/community/tutorials/survival-analysis-R
#but it has some typos (namely that it describes how to include the censoring indicator incorrectly)

##Kaplan-Meier Curves (examples with gender, dichotomized age, dichotomized EREG_mRNAsi given)
surv.object <- Surv(time = dat$Y, event = 1 - dat$censor.indicator) #event = 1 if died, 0 if censored (alive at end of follow-up)
dat <- dat %>% mutate(age_group = ifelse(age >=50, "50 and Older", "Younger than 50")) #Dichotomize age
dat <- dat %>% mutate(risk_group = ifelse(EREG_mRNAsi >= median(EREG_mRNAsi), "Above-Median Stemness", "Below-Median Stemness")) #Dichotomize EREG_mRNAsi
dat <- dat %>% mutate(mDNAsi_group = ifelse(mDNAsi >= median(mDNAsi), "Above-Median Stemness from Methylation", "Below-Median Stemness from Methylation")) #Dichotomize mDNAsi
dat <- dat %>% mutate(mRNSsi_group = ifelse(mRNAsi >= median(mRNAsi), "Above-Median Stemness from Gene Expression", "Below-Median Stemness from Gene Expression")) #Dichotomize mDNAsi


km.curve <- survfit(surv.object ~ gender, data = dat)
ggsurvplot(km.curve, data = dat, pval = TRUE, legend.labs=c("Female", "Male")) + xlab("Time (days post diagnosis)") #pval should correspond to log-rank test
surv_pvalue(km.curve, data = dat, method = "survdiff")

km.curve <- survfit(surv.object ~ age_group, data = dat)
ggsurvplot(km.curve, data = dat, pval = TRUE, legend.labs=c("50 and Older", "Younger than 50")) + xlab("Time (days post diagnosis)")
surv_pvalue(km.curve, data = dat, method = "survdiff")

km.curve <- survfit(surv.object ~ risk_group, data = dat)
ggsurvplot(km.curve, data = dat, pval = TRUE, legend.labs=c("Above-Median Stemness", "Below-Median Stemness")) + xlab("Time (days post diagnosis)")
surv_pvalue(km.curve, data = dat, method = "survdiff")

km.curve <- survfit(surv.object ~ mDNAsi_group, data = dat)
ggsurvplot(km.curve, data = dat, pval = TRUE, legend.labs = c("Above-Median Stemness from Methylation", "Below-Median Stemness from Methylation")) + xlab("Time (days post diagnosis)")
surv_pvalue(km.curve, data = dat, method = "survdiff")

##Cox Proportional-Hazards Model
##EREG-mRNAsi could be used instead of the two separate RNA and DNA stemness indices
fit.coxph <- coxph(surv.object ~  mRNAsi + mDNAsi + gender + age + cancer.class, data = dat)
summary(fit.coxph)
#png(filename="figs/forest_plot.png", res=300)
ggforest(fit.coxph, data = dat)
#dev.off()

#Some interpretations of coefficients. Keep in mind that stemness indices are multiplied by 100.

#Men have a hazards rate that is 1.0570 times that of females, holding age, cancer type, and stemness indices constant.
  #Or: The hazards ratio comparing men to women is 1.0570, holding age, cancer type, and stemness indices constant.

#A one year increase in age is associated with a multiplicative increase of 1.0073 in the hazards rate, holding gender, cancer type, and stemness indices constant.
  #or: The hazards ratio comparing individuals one year older than referent individuals is 1.0073, holding gender, cancer type, and stemness indices constant.

#A one unit increase in the DNA methylation stemness index (when on a 0-100 scale) is associated with a multiplicative increase of 1.0121 in the hazards rate, holding gender, age, cancer type, and RNA stemness index constant.
  #or: The hazards ratio comparing individuals with a DNA methylation stemness index one unit higher than referent individuals (on a 0-100 scale) is 1.0121, holding gender, age, cancer type, and RNA stemness index constant.
library(dplyr)
library(survival)
library(survminer)

setwd('~/TCGA_stemness_ML/')

#load('data/pd.450.prim_20170207.Rda')

covariate_names = c('mDNAsi', 'mRNAsi', 'gender','cancer.type')

outcome <- select(pd.450.prim, days_to_death, days_to_last_followup, covariate_names)
outcome <- filter(outcome, !is.na(days_to_death) | !is.na(days_to_last_followup))
outcome <- mutate(outcome,
                  censored = is.na(days_to_death) * 1,
                  time = ifelse(!censored, days_to_death, days_to_last_followup))


outcome <- mutate(outcome, cancer.class = case_when(
  cancer.type %in% c('LAML', 'DLBC', 'THYM') ~ 'hematologic_lymphatic',
  cancer.type %in% c('OV', 'UCEC', 'CESC', 'BRCA') ~ 'solid, gynecologic',
  cancer.type %in% c('BLCA', 'PRAD', 'TGCT', 'KIRC', 'KICH', 'KIRP') ~ 'solic, urologic',
  cancer.type %in% c('THCA', 'ACC') ~ 'solid, endocrine',
  cancer.type %in% c('ESCA', 'STAD', 'COAD', 'READ') ~ 'solid, core GI',
  cancer.type %in% c('LIHC', 'PAAD', 'CHOL') ~ 'developmental GI',
  cancer.type %in% c('HNSC') ~ 'head and neck',
  cancer.type %in% c('LUAD', 'LUSC', 'MESO') ~ 'thoracic',
  cancer.type %in% c('GBM', 'LGG') ~ 'Central nervous system',
  cancer.type %in% c('SARC', 'UCS') ~ 'soft tissue',
  cancer.type %in% c('PCPG') ~ 'neural crest',
  cancer.type %in% c('SKCM', 'UVM') ~ 'melanocytic'
)
)


outcome <- select(outcome, time, censored, covariate_names, 'cancer.class')

coxph_res <- coxph(Surv(time = time, event = censored) ~ mDNAsi + mRNAsi + gender + cancer.class,
                   data=outcome)
km_gender <- survfit(Surv(time = time, event = censored) ~ gender, data=outcome)

ggsurvplot(km_gender)
ggforest(coxph_res, data=outcome)










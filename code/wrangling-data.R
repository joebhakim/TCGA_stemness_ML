# -- Libraries & Data
library(tidyverse)
library(hexbin)
load("data/pd.450.prim_20170207.Rda")

# -- Wrangling data
# Y : Time to event (outcome)
# Status: Right censor indicator
dat <- pd.450.prim %>%
  as_tibble() %>%
  select(id.12, cancer.type, mRNAsi, mDNAsi, `EREG-mRNAsi`, gender, age_at_initial_pathologic_diagnosis) %>%
  na.omit() %>% 
  left_join(select(pd.450.prim, c(id.12, days_to_death, days_to_last_followup)), by="id.12") %>%
  rename(age = age_at_initial_pathologic_diagnosis) %>%
  select(-id.12) %>%
  filter(!(is.na(days_to_death) & is.na(days_to_last_followup))) %>%
  mutate(Y      = ifelse(!is.na(days_to_death), days_to_death, days_to_last_followup),
         censor.indicator = ifelse(!is.na(days_to_death), 0, 1)) %>%
  mutate(cancer.class = case_when(
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
  ) %>%
  select(-days_to_death, -days_to_last_followup)
save(dat, file="data/tcga_wrangled.rda", compress="xz")

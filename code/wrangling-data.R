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
         status = ifelse(!is.na(days_to_death), 1, 0)) %>%
  select(-days_to_death, -days_to_last_followup)
save(dat, file="data/tcga_wrangled.rda", compress="xz")

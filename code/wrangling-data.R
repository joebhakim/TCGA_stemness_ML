# -- Libraries
library(tidyverse)

# -- Loading data
load("data/pd.450.prim_20170207.Rda")
load("data/RNA_subset.Rda")
load("data/data.pan.Rda")

# -- Wrangling data
# Y : Time to event (outcome)
# Status: Right censor indicator
dat <- pd.450.prim %>%
  as_tibble() %>%
  select(TCGAlong.id, cancer.type, mRNAsi, mDNAsi, `EREG-mRNAsi`, gender, age_at_initial_pathologic_diagnosis) %>%
  na.omit() %>%
  left_join(select(pd.450.prim, c(TCGAlong.id, days_to_death, days_to_last_followup)), by="TCGAlong.id") %>%
  rename(age = age_at_initial_pathologic_diagnosis) %>%
  filter(!(is.na(days_to_death) & is.na(days_to_last_followup))) %>%
  mutate(Y = ifelse(!is.na(days_to_death), days_to_death, days_to_last_followup),
         censor.indicator = ifelse(!is.na(days_to_death), 0, 1),
         TCGAlong.id = as.character(TCGAlong.id)) %>%
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


# -- Methylation probe ids
probes <- rownames(data.pan)

# -- Wrangling methylation data per cancer patient and joinning it to the dataset
tmp.data.pan <- data.pan %>%
  as_tibble() %>%
  gather(TCGAlong.id, methylation) %>%
  mutate(probe = rep(probes, 9627)) %>%
  spread(probe, methylation)

# -- Joining data
dat <- left_join(dat, tmp.data.pan, by = "TCGAlong.id") %>%
  filter(Y > 0)

# -- Last wrangled
dat <- dat %>%
  select(-TCGAlong.id, -cancer.type, -mRNAsi, -`EREG-mRNAsi`)

# -- Creating the design matrix
design <- model.matrix(~., dat)

# # -- Methylation probe ids
# genes <- RNA_subset[,1]
# 
# # -- Wrangling RNA subset data per cancer patient and joinning it to the dataset
# tmp.rna <- RNA_subset[,-1] %>%
#   as_tibble() %>%
#   gather(TCGAlong.id, expression) %>%
#   mutate(genes = rep(genes, 11069)) %>%
#   mutate(TCGAlong.id = gsub("\\.", "-", TCGAlong.id)) %>%
#   spread(genes, expression)
# 
# # -- Joining data
# dat <- left_join(dat, tmp.rna, by = "TCGAlong.id") %>%
#   filter(Y > 0)

save(dat, file="data/tcga_wrangled.rda", compress="xz")
save(design, file="data/design-matrix.rda", compress="xz")
# -- Libraries
library(tidyverse)

# -- Loading data
load("data/pd.450.prim_20170207.Rda")
load("data/RNA_subset.Rda")
load("data/data.pan.Rda")
load('data/pd.maf.450.Rda')

# -- Wrangling data
# Y : Time to event (outcome)
# Status: Right censor indicator


varnames_selected = c(
  'TCGAlong.id',
  'cancer.type',
  'mRNAsi',
  'mDNAsi',
  'EREG-mRNAsi',
  'gender',
  'age_at_initial_pathologic_diagnosis',
  'purity',
  'ploidy',
  'radiation_therapy'
)

dat <- pd.450.prim %>%
  as_tibble() %>%
  select(varnames_selected) %>%
  #na.omit() %>%
  left_join(select(pd.450.prim, c(TCGAlong.id, days_to_death, days_to_last_followup)), by="TCGAlong.id") %>%
  rename(age = age_at_initial_pathologic_diagnosis) %>%
  filter(!(is.na(days_to_death) & is.na(days_to_last_followup))) %>%
  mutate(Y = ifelse(!is.na(days_to_death), days_to_death, days_to_last_followup),
         censor.indicator = ifelse(!is.na(days_to_death), 0, 1),
         TCGAlong.id = as.character(TCGAlong.id)) %>%
  mutate(cancer.class = case_when(
      cancer.type %in% c('LAML', 'DLBC', 'THYM') ~ 'hematologic_lymphatic',
      cancer.type %in% c('OV', 'UCEC', 'CESC', 'BRCA') ~ 'solid_gynecologic',
      cancer.type %in% c('BLCA', 'PRAD', 'TGCT', 'KIRC', 'KICH', 'KIRP') ~ 'solic_urologic',
      cancer.type %in% c('THCA', 'ACC') ~ 'solid_endocrine',
      cancer.type %in% c('ESCA', 'STAD', 'COAD', 'READ') ~ 'solid_core_GI',
      cancer.type %in% c('LIHC', 'PAAD', 'CHOL') ~ 'developmental_GI',
      cancer.type %in% c('HNSC') ~ 'head_and_neck',
      cancer.type %in% c('LUAD', 'LUSC', 'MESO') ~ 'thoracic',
      cancer.type %in% c('GBM', 'LGG') ~ 'central_nervous_system',
      cancer.type %in% c('SARC', 'UCS') ~ 'soft_tissue',
      cancer.type %in% c('PCPG') ~ 'neural_crest',
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

#### Adding in raw mutation data (!), there are more columns i don't understand that im keeping for now
maf_data <- pd.maf.450 %>%
    select(-id.16, -id.12, -sample.type, -cancer.type, -germlayer, -Non.SC_310_DNAmeth, -Non.SC_62_DNAmeth, -SC_310_DNAmeth, -SC_Enh, -RNAss, )

gene_names = colnames(maf_data)[-0:-10]

#remove genes that are barely ever mutated
countWT <- function(x){
  return(sum(x == 'WT')/length(x))
}

mutated_enough <- summarize_all(maf_data[gene_names], countWT)
thresh <- quantile(mutated_enough, 0.1)
selected_genes <- Filter(function(x) x < thresh$`10%`, mutated_enough)
maf_data_subset <- maf_data[c('TCGAlong.id',colnames(selected_genes))]
colnames(maf_data_subset) <- paste0("gene",colnames(maf_data_subset))


dat_all <- left_join(dat, maf_data_subset, by = c('TCGAlong.id'='geneTCGAlong.id'))
dat_all <- rename(dat_all, 'geneBIVM_ERCC5' = `geneBIVM-ERCC5`)
dat_all <- rename(dat_all, 'geneLY75_CD302' = `geneLY75-CD302`)
dat_all <- rename(dat_all, 'genePALM2_AKAP2' = `genePALM2-AKAP2`)
dat_all <- rename(dat_all, 'geneSTON1_GTF2A1L' = `geneSTON1-GTF2A1L`)
dat_all <- rename(dat_all, 'geneRP11_1055B8_7' = `geneRP11-1055B8.7`)

# -- Last wrangled
dat_all <- dat_all %>%
  select(-TCGAlong.id, -cancer.type, -mRNAsi, -`EREG-mRNAsi`, -purity, -ploidy, ) %>%
  na.omit()

# -- Creating the design matrix
design <- model.matrix(~., dat_all)

# -- Used below
ys <- design[,"Y"]

design <- design %>%
  as_tibble() %>%
  mutate(genderMALE = factor(genderMALE),
         censor.indicator = factor(censor.indicator)) %>%
  mutate_at(vars(matches("cancer")), as.factor) %>%
  mutate_at(vars(matches("gene")), as.factor) %>%
  mutate_if(is.numeric, scale) %>%
  mutate(`(Intercept)`=1) %>%
  mutate_all(as.character) %>%
  mutate_all(as.numeric) %>%
  mutate(Y = ys)


# colnames(design) <- gsub(" ", "", colnames(design))
# colnames(design) <- gsub("\\.", "", colnames(design))
# colnames(design) <- gsub("_", "", colnames(design))
# colnames(design) <- gsub(",", "", colnames(design))
# colnames(design)

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


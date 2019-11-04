# -- Libraries & Data
library(tidyverse)
library(hexbin)
load("data/tcga_wrangled.rda")

# -- Proportion of censored/observed
dat %>%
  group_by(status) %>%
  summarize(n = n()/ nrow(.) * 100) %>%
  ungroup() %>%
  ggplot(aes(status, n)) +
  geom_col(fill="black", color="black", alpha=0.80, size=1) +
  scale_x_continuous(breaks = c(0,1), labels = c("Censored", "Observed")) +
  scale_y_continuous(breaks = seq(0, 75, by=5)) +
  ylab("% of the data") +
  xlab("Censor Indicator") +
  theme_minimal() +
  theme(axis.text  = element_text(face="bold", color="black"),
        axis.title = element_text(face="bold", color="black"))

# -- Data distribution by some variables
dat %>%
  group_by(cancer.type, gender, status) %>%
  summarize(n = n() / nrow(.) * 100) %>%
  ungroup() %>%
  mutate(status = ifelse(status==1, "Observed", "Censored")) %>%
  ggplot(aes(reorder(cancer.type, n, median), n, fill = gender)) +
  geom_col(alpha=0.80, color="black") +
  facet_wrap(~status) +
  xlab("Cancer type") +
  ylab("% of the data") +
  scale_fill_manual(name="",
                    values = c("red3", "black"),
                    labels = c("Female", "Male")) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text  = element_text(face="bold", color="black"),
        axis.title = element_text(face="bold", color="black"), 
        strip.text = element_text(face="bold", color="black"),
        legend.text = element_text(face="bold", color="black"),
        legend.position = "bottom")

# -- mRNAsi vs mDNAsi
cor(dat$mDNAsi, dat$mRNAsi)
dat %>%
  ggplot(aes(mDNAsi, mRNAsi)) +
  geom_hex(color="black", bins=30) +
  geom_smooth(method="lm", color="red") +
  theme_minimal() +
  theme(axis.text  = element_text(face="bold", color="black"),
        axis.title = element_text(face="bold", color="black"), 
        strip.text = element_text(face="bold", color="black"),
        legend.text = element_text(face="bold", color="black"),
        legend.title = element_text(face="bold", color="black"))

# -- mRNAsi vs EREG-mRNAsi
cor(dat$`EREG-mRNAsi`, dat$mRNAsi)
dat %>%
  ggplot(aes(`EREG-mRNAsi`, mRNAsi)) +
  geom_hex(color="black", bins=30) +
  geom_smooth(method="lm", color="red") +
  theme_minimal() +
  theme(axis.text  = element_text(face="bold", color="black"),
        axis.title = element_text(face="bold", color="black"), 
        strip.text = element_text(face="bold", color="black"),
        legend.text = element_text(face="bold", color="black"),
        legend.title = element_text(face="bold", color="black"))

# -- mDNAsi vs EREG-mRNAsi
cor(dat$mDNAsi, dat$`EREG-mRNAsi`)
dat %>%
  ggplot(aes(`EREG-mRNAsi`, mDNAsi)) +
  geom_hex(color="black", bins=30) +
  geom_smooth(method="lm", color="red") +
  theme_minimal() +
  theme(axis.text  = element_text(face="bold", color="black"),
        axis.title = element_text(face="bold", color="black"), 
        strip.text = element_text(face="bold", color="black"),
        legend.text = element_text(face="bold", color="black"),
        legend.title = element_text(face="bold", color="black"))

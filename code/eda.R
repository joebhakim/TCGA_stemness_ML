# -- Libraries & Data
library(tidyverse)
library(hexbin)
load("data/tcga_wrangled.rda")

# -- Proportion of censored/observed
dat %>%
  group_by(censor.indicator) %>%
  summarize(n = n()/ nrow(.) * 100) %>%
  ungroup() %>%
  ggplot(aes(censor.indicator, n)) +
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
  group_by(cancer.type, gender, censor.indicator) %>%
  summarize(n = n() / nrow(.) * 100) %>%
  ungroup() %>%
  mutate(censor.indicator = ifelse(censor.indicator==1, "Observed", "Censored")) %>%
  ggplot(aes(reorder(cancer.type, n, median), n, fill = gender)) +
  geom_col(alpha=0.80, color="black") +
  facet_wrap(~censor.indicator) +
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

dat %>%
  group_by(cancer.type) %>%
  summarize(mRNAsi = mean(mRNAsi),
            mDNAsi = mean(mDNAsi),
            `EREG-mRNAsi` = mean(`EREG-mRNAsi`)) %>%
  ungroup() %>%
  ggplot(aes(reorder(cancer.type, mRNAsi), mRNAsi, fill=cancer.type)) +
  geom_col()

dat %>%
  arrange(desc(mRNAsi))
  ggplot() +
  geom_col(aes(reorder(1:8392, mRNAsi), mRNAsi, fill=cancer.type)) +
  theme_classic() +
  xlab("") + 
  theme(legend.position = "none",
        axis.text = element_text(face="bold"),
        axis.text.x = element_blank())

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
cor(dat$`EREG-mRNAsi`, dat$mRNAsi, method = "spearman")
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
cor(dat$mDNAsi, dat$`EREG-mRNAsi`, method = "spearman")
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

# -- Age distribution by sex
dat %>%
  ggplot(aes(age)) +
  geom_density(fill="black", alpha=0.80, size=1) +
  facet_wrap(~gender) +
  xlab("Age") +
  ylab("Density") +
  theme_minimal() 

# -- Age distribution by cancer type
dat %>%
  ggplot(aes(age)) +
  geom_density(fill="black", alpha=0.80, size=1) +
  facet_wrap(~cancer.type) +
  xlab("Age") +
  ylab("Density") +
  theme_minimal() 

# -- mRNAsi vs mDNAsi by gender
dat %>%
  ggplot(aes(mRNAsi, mDNAsi)) +
  geom_point(alpha=0.50) +
  facet_wrap(~gender) +
  theme_minimal()

# -- mRNAsi vs mDNAsi by gender
dat %>%
  ggplot(aes(mRNAsi, mDNAsi, color=gender)) +
  geom_point(alpha=0.50) +
  facet_wrap(~cancer.type) +
  theme_minimal() +
  theme(legend.position = "bottom")


# -- Correlation between mRNAsi and mDNAsi by gender and cancer type
dat %>%
  group_by(cancer.type, gender) %>%
  summarize(corr = cor(mRNAsi, mDNAsi)) %>%
  ungroup() %>% 
  arrange(desc(abs(corr))) %>%
  View()








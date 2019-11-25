# -- Libraries
library(tidyverse)

# -- Loading training data 
load("data/pcbc.data.Rda")
train_probes <- rownames(pcbc.data)
train        <- pcbc.data %>%
                  as_tibble() %>%
                  mutate(probe = train_probes) %>%
                  gather(id, train_methylation, -probe) %>%
                  group_by(probe) %>%
                  summarize(train_methylation = mean(train_methylation)) %>%
                  ungroup()

# -- Loading test data
load("data/data.pan.Rda")
test_probes <- rownames(data.pan)
test        <- data.pan %>%
                as_tibble() %>%
                mutate(probe = test_probes) %>%
                gather(id, test_methylation, -probe) %>%
                group_by(probe) %>%
                summarize(test_methylation = mean(test_methylation, na.rm=TRUE)) %>%
                ungroup()

# -- Viz of methylation distribution
train %>%
  left_join(test, by = "probe") %>%
  gather(set, methylation, -probe) %>%
  mutate(set = ifelse(set=="train_methylation", "Train", "Test")) %>%
  ggplot(aes(probe, methylation)) +
  geom_col(fill="black", position = "dodge2") +
  facet_wrap(~set, nrow=2) +
  xlab("Probe") +
  ylab("Methylation") +
  theme(axis.text.x = element_blank())


# -- Viz methylation values by germ layer in cancer samples
data.pan %>%
  as_tibble() %>%
  mutate(probe = test_probes) %>%
  gather(id, test_methylation, -probe) %>%
  left_join(type.info, by = c("id"="sample")) %>%
  na.omit() %>%
  group_by(probe, germlayer) %>%
  summarize(meth = mean(test_methylation)) %>%
  ungroup() %>%
  ggplot(aes(reorder(probe, meth), meth, color=germlayer)) +
  geom_point(alpha=0.50) +
  ylab("Methylation") +
  xlab("Probe") +
  scale_color_manual(name = "Germ Layer",
                     values = c("black", "red3", "blue2")) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
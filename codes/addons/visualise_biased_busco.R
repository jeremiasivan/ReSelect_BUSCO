# This script aims to visualise proportion of biased BUSCO for each set of short reads
# Written by: Jeremias Ivan
# Last updated: 07 August 2025

library(tidyverse)

df_metadata <- data.table:::fread("~/BusIER/data/eucs_shortreads.txt")

dfpipe_bias <- dfpipe %>%
  group_by(read) %>%
  mutate(total=n()) %>%
  filter(lm_p < 0.05 | spearman_p < 0.05) %>%
  mutate(note=ifelse(lm_p < 0.05 & spearman_p >= 0.05, "lm",
                     ifelse(lm_p < 0.05 & spearman_p < 0.05, "both", "sp"))) %>%
  group_by(read, note) %>%
  summarise(n=n()/total*100) %>%
  mutate(method="pipe") %>%
  unique()

dfpipe_bias_noout <- dfpipe %>%
  group_by(read) %>%
  mutate(total=n()) %>%
  filter(lm_p_noout < 0.05 | spearman_p_noout < 0.05) %>%
  mutate(note=ifelse(lm_p_noout < 0.05 & spearman_p_noout >= 0.05, "lm",
                     ifelse(lm_p_noout < 0.05 & spearman_p_noout < 0.05, "both", "sp"))) %>%
  group_by(read, note) %>%
  summarise(n=n()/total*100) %>%
  mutate(method="pipe_noout") %>%
  unique()

df_bias <- rbind(dfpipe_bias, dfpipe_bias_noout)

ggplot(df_bias, aes(x=read, y=n)) +
  geom_bar(aes(group=read, fill=note, colour=note), stat="identity") +
  coord_cartesian(ylim = c(0, 100)) +
  scale_fill_manual(values = c("#482677FF","#2D708EFF","#29AF7FFF")) +
  scale_colour_manual(values = c("#482677FF","#2D708EFF","#29AF7FFF")) +
  facet_wrap(.~method) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=15, angle=45, hjust=1, vjust=1, face='italic'),
    legend.title=element_text(size=20),
    legend.text=element_text(size=20),
    legend.key.size=unit(1,"cm"),
    strip.text=element_text(size=20)
  )

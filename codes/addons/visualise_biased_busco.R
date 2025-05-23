# This script aims to visualise proportion of biased BUSCO for each set of short reads
# Written by: Jeremias Ivan
# Last updated: 23 May 2025

library(tidyverse)

df_metadata <- data.table:::fread("~/BusIER/data/eucs_shortreads.txt")
dfpipe <- data.table::fread("prefix.cor_pipeline.sumtable")
dfcoor <- data.table::fread("prefix.cor_coordinate.sumtable")

pipe_ln_busco <- length(unique(dfpipe$busco))
coor_ln_busco <- length(unique(dfcoor$busco))

dfpipe_bias <- dfpipe %>%
  group_by(read) %>%
  filter(lm_p < 0.05 | spearman_p < 0.05) %>%
  mutate(note=ifelse(lm_p < 0.05 & spearman_p >= 0.05, "lm",
                     ifelse(lm_p < 0.05 & spearman_p < 0.05, "both", "sp"))) %>%
  group_by(read, note) %>%
  summarise(n=n()/pipe_ln_busco*100) %>%
  mutate(method="pipe")

dfpipe_bias_noout <- dfpipe %>%
  group_by(read) %>%
  filter(lm_p_noout < 0.05 | spearman_p_noout < 0.05) %>%
  mutate(note=ifelse(lm_p_noout < 0.05 & spearman_p_noout >= 0.05, "lm",
                     ifelse(lm_p_noout < 0.05 & spearman_p_noout < 0.05, "both", "sp"))) %>%
  group_by(read, note) %>%
  summarise(n=n()/pipe_ln_busco*100) %>%
  mutate(method="pipe_noout")

dfcoor_bias <- dfcoor %>%
  group_by(read) %>%
  filter(lm_p < 0.05 | spearman_p < 0.05) %>%
  mutate(note=ifelse(lm_p < 0.05 & spearman_p >= 0.05, "lm",
                     ifelse(lm_p < 0.05 & spearman_p < 0.05, "both", "sp"))) %>%
  group_by(read, note) %>%
  summarise(n=n()/coor_ln_busco*100) %>%
  mutate(method="coor")

dfcoor_bias_noout <- dfcoor %>%
  group_by(read) %>%
  filter(lm_p_noout < 0.05 | spearman_p_noout < 0.05) %>%
  mutate(note=ifelse(lm_p_noout < 0.05 & spearman_p_noout >= 0.05, "lm",
                     ifelse(lm_p_noout < 0.05 & spearman_p_noout < 0.05, "both", "sp"))) %>%
  group_by(read, note) %>%
  summarise(n=n()/coor_ln_busco*100) %>%
  mutate(method="coor_nobias")

df_bias <- rbind(dfpipe_bias, dfcoor_bias, dfpipe_bias_noout, dfcoor_bias_noout)
df_bias$read <- sapply(df_bias$read, function(x) {
  name <- df_metadata$species[df_metadata$id==x]
  gsub("\\.", ". ", name)
}) 

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
    axis.text.x=element_text(size=15, angle=45, hjust=1, vjust=1),
    legend.title=element_text(size=20),
    legend.text=element_text(size=20),
    legend.key.size=unit(1,"cm"),
    strip.text=element_text(size=20)
  )

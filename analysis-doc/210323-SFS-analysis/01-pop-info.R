# pop info table
library(tidyr)
library(dplyr)
library(readr)
library(stringr)

mxb_samples <- read_tsv("../../resources/genomes-metadata/50Genomes_info.txt") %>% 
  select(Sample_ID) %>% 
  mutate(
    Sample = str_replace(Sample_ID, "MXB", "MXB_"),
    Population = "MXB"
  ) %>% 
  select(-Sample_ID)

one_tg <- read_tsv("../../resources/1TGP-samples-meta-data/igsr-1000genomes.tsv") %>% 
  select(`Sample name`, `Population code`) %>% 
  rename(Sample = `Sample name`, Population = `Population code`)

sample_list <- read_lines("data/samples.txt")

bind_rows(mxb_samples, one_tg) %>% 
  filter(Sample %in% sample_list) %>% 
  write_csv("data/popinfo.txt")

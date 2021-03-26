library(readr)
library(dplyr)
library(tidyr)
library(stringr)


# This script defines the population used as query and reference haplotypes
# parameters

oneTGP_pops <-  snakemake@params$oneTGP_pops
path_oneT_meta <- "resources/1TGP-samples-meta-data/integrated_call_samples_v3.20130502.ALL.panel"
path_mxb_meta <- "resources/genomes-metadata/50Genomes_info.txt"
oneTGP_NATs <- "resources/1TGP-samples-meta-data/native-american.txt"
# outfiles
sample_map_file <- snakemake@output$sample_map_file
query_pops_file <- snakemake@output$query_pops


# we include this samples as reference for NAT plus the MXB genomes
# NOTE this set of samples will be also in the query and reference haplotypes
# at the same time
oneTGP_NATs <- read_lines(oneTGP_NATs)
oneTGP_NATs <- tibble(
  Sample = oneTGP_NATs,
  Population = "NAT"
)
# filter populations ------------------------------------------------------


oneT_meta <- read.table(path_oneT_meta, header = T) %>% 
  as_tibble()

oneTGP_pops <- 
  oneT_meta %>%
  select(sample, pop, super_pop) %>% 
  filter(str_detect(pop, oneTGP_pops)) %>% 
  rename(Sample = sample, Population = super_pop, `Population code` = pop) 

mxb_meta <- read_tsv(path_mxb_meta) %>% 
  mutate(
    Sample = str_replace(Sample_ID, "MXB", "MXB_"),
    Population = "NAT",
    `Population code` = "MXB"
  )


populations <- 
  bind_rows(oneTGP_pops, mxb_meta) %>% 
  arrange(Sample, Population, `Population code`) %>%
  select(Sample, Population)


# Reference haplotypes:
# All non AMR samples
# and the AMR samples that are proxy for NAT

sample_map <- filter(populations, Population != "AMR")
sample_map <- bind_rows(sample_map, oneTGP_NATs)



write_tsv(sample_map, sample_map_file)


query_pop <- filter(populations, Population == "AMR")
write_tsv(query_pop, query_pops_file)

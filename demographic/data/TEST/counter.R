library(tidyverse)

counter <- read.csv('~/tmp/tmp/vep-chr2-UNIQ.txt', sep = '\t') %>% 
  as_tibble()


counter <- counter %>% 
  select(-Gene, -Feature_type)


Q <- c('missense', 'synonymous', 'LOF', 'LOF', 'LOF')
Consequence <- c(
  'missense_variant',
  'synonymous_variant',
  'stop_lost', 'stop_gained', 'start_lost'
  )

Q <- tibble(
  Q = Q,
  Consequence = Consequence
)


cpgs <- c("CGA", "CGT", "CGG", "ACG", "TCG", "GCG", "CCG", "CGC")

d <- inner_join(counter, Q)

d <- 
  d %>% 
  mutate(
    context = str_extract(Uploaded_variation, '_[ACGT]{5}_'),
    context = str_sub(context, start = 3, end = 5)
  )


d <- 
  d %>% 
  mutate(
    is_CpG = context %in% cpgs
  )


d %>% 
  group_by(Q) %>% 
  summarise(
    CpG_proportion = mean(is_CpG)
  )

d_coutns <- d %>% 
  count(Q, Allele, context)

mus <- read_tsv('../211128-compute-mL/data/mutation_rate_methylation_bins.txt') %>% 
  filter(methylation_level == 0)


# get rev complement ------------------------------------------------------


rev_com <- function(x) {
  Biostrings::reverseComplement(Biostrings::DNAString(x)) %>% 
    as.character()
  
}

mus_rev <- mus %>% 
  mutate(
    context = map_chr(context, rev_com)
  )


all_mus <- bind_rows(mus, mus_rev) 

d_coutns %>% 
  rename(alt = Allele) %>% 
  inner_join(all_mus, by=c('context', 'alt')) %>% 
  mutate(
    mL = n * mu_snp
  ) %>% 
  group_by(Q) %>% 
  summarise(
    mL = sum(mL)
  )

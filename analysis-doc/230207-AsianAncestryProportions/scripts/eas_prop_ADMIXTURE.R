library(tidyverse)

admix <- read_csv('../210310-AdmixurePCA-merged-data-with-1TGP/results/admixture.csv')


# get the EAS ancestry bassed on K = 4 ------------------------------------
cohorts <- c('CLM', 'MXL', 'PEL', 'PUR')

admix %>% 
  filter(
    K == 'K = 4'
  ) %>% 
  group_by(Population, cluster_grp) %>% 
  summarise(
    mean_anp = mean(p)
  ) %>% 
  filter(Population %in% cohorts) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = Population, names_from = cluster_grp, values_from = mean_anp) %>% 
  write_csv('results/ADMIXTURE-ancestry-props.csv')


library(tidyverse)
library(ggtern)
library(ggthemes)
library(scales)

theme_set(theme_tufte(base_family = 'Helvetica'))



get_pop_adx_data <- function(population) {
  
  adx <- read_csv('../210310-AdmixurePCA-merged-data-with-1TGP/results/admixture.csv') %>% 
    filter(
      K == 'K = 3',
      Population == population
    )
  
  adx %>% 
    select(Sample, p, cluster_grp) %>% 
    pivot_wider(id_cols = Sample, names_from = cluster_grp, values_from = p) %>% 
    mutate(
      Population = population
    )
  
  
}


adx <- map_df(c('MXL', 'PEL', 'CLM', 'PUR'), get_pop_adx_data)
mxl <- filter(adx, Population == 'MXL')


adx %>% 
  filter(Population != 'MXL') %>% 
  ggtern(aes(AFR, MXB, EUR)) +
  geom_point(
    aes(shape = Population),
    color = 'grey'
  ) +
  geom_point(
    data = mxl,
    shape = 21,
    size = 2,
    fill = '#440154'
  )
ggsave('plots/anc-p-ternary.pdf', height = 3, width = 4)

adx %>% 
  ggtern(aes(AFR, MXB, EUR)) +
  geom_point(
    data = select(adx, -Population),
    color = 'grey', size = 0.5
  ) +
  geom_point(
    aes(color = Population),
    color = '#440154', size = 0.5
  ) +
  facet_wrap(~Population)
ggsave('plots/anc-p-ternary-grid.pdf', height = 3, width = 4)
